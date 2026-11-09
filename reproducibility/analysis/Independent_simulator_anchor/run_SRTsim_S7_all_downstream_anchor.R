#!/usr/bin/env Rscript

Sys.setenv(RGL_USE_NULL = "TRUE")

parse_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}

parse_reps <- function(x) {
  if (grepl(":", x, fixed = TRUE)) {
    parts <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    return(seq(parts[1], parts[2]))
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

if (!requireNamespace("aricode", quietly = TRUE)) {
  stop("The aricode package is required for ARI/NMI.")
}

repo_dir <- normalizePath(parse_arg("repo", "."), mustWork = TRUE)
reps <- parse_reps(parse_arg("reps", "1:5"))
seed <- as.integer(parse_arg("seed", "20260902"))
simulate_missing <- tolower(parse_arg("simulate-missing", "true")) %in% c("true", "t", "1", "yes", "y")
srtsim_mode <- parse_arg("srtsim-mode", "ct-aggregate")
default_srtsim_dir <- if (srtsim_mode == "ct-aggregate") {
  "Independent_simulator_anchor/SRTsim_S7_celltype_aggregate"
} else {
  "Independent_simulator_anchor/SRTsim_S7_segmentation_anchor"
}
srtsim_dir <- file.path(repo_dir, parse_arg("srtsim-dir", default_srtsim_dir))
out_dir <- file.path(repo_dir, parse_arg("out", "Independent_simulator_anchor/SRTsim_S7_all_downstream_anchor"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ct_names <- c("Cancer", "Tcell", "Endothelial", "CAF")
independent_label <- if (srtsim_mode == "ct-aggregate") "SRTsim CT-aggregate S7" else "SRTsim on S7"
cell_per_spot <- c(`1` = 9, `2` = 9, `3` = 10, `4` = 10, `5` = 10)
marker_genes <- list(
  Cancer = paste0("gene_", 41:50),
  Tcell = paste0("gene_", 51:60),
  Endothelial = paste0("gene_", 61:70),
  CAF = paste0("gene_", 71:80)
)

read_count_matrix <- function(file) {
  df <- read.csv(file, check.names = FALSE)
  gene_cols <- grep("^gene_", colnames(df), value = TRUE)
  if (length(gene_cols) > 0) {
    df <- df[, gene_cols, drop = FALSE]
  }
  mat <- data.matrix(df)
  colnames(mat) <- colnames(df)
  rownames(mat) <- paste0("spot_", seq_len(nrow(mat)))
  mat
}

normalize_counts <- function(counts) {
  lib_size <- rowSums(counts)
  target <- stats::median(lib_size[lib_size > 0])
  if (!is.finite(target)) {
    target <- 1
  }
  scaled <- counts / pmax(lib_size, 1) * target
  log1p(scaled)
}

normalize_rows <- function(mat) {
  mat[!is.finite(mat)] <- 0
  mat <- pmax(mat, 0)
  rs <- rowSums(mat)
  zero <- rs <= 0
  if (any(zero)) {
    mat[zero, ] <- 1 / ncol(mat)
    rs <- rowSums(mat)
  }
  mat / rs
}

neighbor_indices <- function(coord, k = 6) {
  d <- as.matrix(stats::dist(coord))
  diag(d) <- Inf
  apply(d, 1, order)[seq_len(k), , drop = FALSE]
}

knn_smooth <- function(mat, nbrs) {
  out <- mat
  for (i in seq_len(nrow(mat))) {
    out[i, ] <- colMeans(mat[c(i, nbrs[, i]), , drop = FALSE])
  }
  out
}

boundary_mask <- function(labels, nbrs) {
  vapply(seq_along(labels), function(i) any(labels[nbrs[, i]] != labels[i]), logical(1))
}

pca_features <- function(mat, n_pc = 20) {
  keep <- apply(mat, 2, stats::sd, na.rm = TRUE) > 0
  mat <- mat[, keep, drop = FALSE]
  n_pc <- min(n_pc, nrow(mat) - 1, ncol(mat))
  if (n_pc < 1) {
    return(scale(mat))
  }
  fit <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
  fit$x[, seq_len(n_pc), drop = FALSE]
}

permute <- function(x) {
  if (length(x) <= 1) {
    return(list(x))
  }
  unlist(lapply(seq_along(x), function(i) {
    lapply(permute(x[-i]), function(rest) c(x[i], rest))
  }), recursive = FALSE)
}

best_label_accuracy <- function(truth, pred, subset = rep(TRUE, length(truth))) {
  truth <- as.character(truth[subset])
  pred <- as.character(pred[subset])
  pred_levels <- unique(pred)
  truth_levels <- unique(truth)
  if (length(pred_levels) > length(truth_levels)) {
    tab <- table(pred, truth)
    mapping <- colnames(tab)[max.col(tab, ties.method = "first")]
    names(mapping) <- rownames(tab)
    return(mean(mapping[pred] == truth))
  }
  maps <- permute(truth_levels)
  best <- 0
  for (m in maps) {
    mapping <- m[seq_along(pred_levels)]
    names(mapping) <- pred_levels
    best <- max(best, mean(mapping[pred] == truth))
  }
  best
}

binary_metrics <- function(pred, truth, universe) {
  pred <- unique(intersect(pred, universe))
  truth <- unique(intersect(truth, universe))
  tp <- sum(pred %in% truth)
  fp <- length(setdiff(pred, truth))
  fn <- length(setdiff(truth, pred))
  tn <- length(universe) - tp - fp - fn
  precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  fpr <- if ((fp + tn) == 0) 0 else fp / (fp + tn)
  f1 <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  jaccard <- if ((tp + fp + fn) == 0) 0 else tp / (tp + fp + fn)
  data.frame(
    n_pred = length(pred),
    tp = tp,
    fp = fp,
    fn = fn,
    precision = precision,
    TPR = recall,
    FPR = fpr,
    F1 = f1,
    Jaccard = jaccard,
    stringsAsFactors = FALSE
  )
}

rbind_fill <- function(dfs) {
  dfs <- Filter(Negate(is.null), dfs)
  cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(df) {
    missing <- setdiff(cols, names(df))
    for (col in missing) {
      df[[col]] <- NA
    }
    df[, cols, drop = FALSE]
  })
  do.call(rbind, dfs)
}

metric_row <- function(task, simulator, method, rep_id, primary_metric, direction, extras = data.frame()) {
  cbind(
    data.frame(
      task = task,
      simulator = simulator,
      method = method,
      rep = rep_id,
      primary_metric = primary_metric,
      direction = direction,
      stringsAsFactors = FALSE
    ),
    extras
  )
}

score_clustering <- function(truth, pred, boundary) {
  ari <- aricode::ARI(as.character(truth), as.character(pred))
  nmi <- aricode::NMI(as.character(truth), as.character(pred))
  acc <- best_label_accuracy(truth, pred)
  boundary_acc <- best_label_accuracy(truth, pred, boundary)
  interior_acc <- best_label_accuracy(truth, pred, !boundary)
  data.frame(
    ARI = ari,
    NMI = nmi,
    accuracy = acc,
    boundary_accuracy = boundary_acc,
    interior_accuracy = interior_acc,
    boundary_penalty = interior_acc - boundary_acc,
    stringsAsFactors = FALSE
  )
}

run_segmentation <- function(counts, coord, labels, nbrs, rep_id, simulator, seed) {
  labels <- as.character(labels)
  k <- length(unique(labels))
  boundary <- boundary_mask(labels, nbrs)
  log_counts <- normalize_counts(counts)
  expr_pc <- pca_features(log_counts)
  smooth_pc <- pca_features(knn_smooth(log_counts, nbrs))
  coord_features <- scale(coord)

  set.seed(seed + rep_id * 100 + 1)
  preds <- list(`PCA k-means` = stats::kmeans(expr_pc, centers = k, nstart = 50, iter.max = 100)$cluster)
  set.seed(seed + rep_id * 100 + 2)
  preds[["spatial-smoothed PCA k-means"]] <- stats::kmeans(smooth_pc, centers = k, nstart = 50, iter.max = 100)$cluster
  set.seed(seed + rep_id * 100 + 3)
  preds[["coordinate k-means"]] <- stats::kmeans(coord_features, centers = k, nstart = 50, iter.max = 100)$cluster
  preds[["PCA Ward clustering"]] <- stats::cutree(stats::hclust(stats::dist(expr_pc), method = "ward.D2"), k = k)

  rows <- lapply(names(preds), function(method) {
    met <- score_clustering(labels, preds[[method]], boundary)
    metric_row("segmentation", simulator, method, rep_id, met$ARI, "higher", met)
  })
  do.call(rbind, rows)
}

moran_scores <- function(log_counts, nbrs) {
  n <- nrow(log_counts)
  k <- nrow(nbrs)
  vapply(seq_len(ncol(log_counts)), function(j) {
    x <- log_counts[, j]
    xc <- x - mean(x)
    denom <- sum(xc^2)
    if (!is.finite(denom) || denom <= 0) {
      return(-Inf)
    }
    num <- sum(xc * colSums(matrix(xc[nbrs], nrow = k)))
    (n / (n * k)) * num / denom
  }, numeric(1))
}

lm_f_scores <- function(log_counts, design_full) {
  design_null <- matrix(1, nrow = nrow(log_counts), ncol = 1)
  df1 <- ncol(design_full) - ncol(design_null)
  df2 <- nrow(log_counts) - ncol(design_full)
  vapply(seq_len(ncol(log_counts)), function(j) {
    y <- log_counts[, j]
    fit0 <- stats::lm.fit(design_null, y)
    fit1 <- stats::lm.fit(design_full, y)
    sse0 <- sum(fit0$residuals^2)
    sse1 <- sum(fit1$residuals^2)
    if (!is.finite(sse1) || sse1 <= 0 || df2 <= 0) {
      return(0)
    }
    max(0, ((sse0 - sse1) / df1) / (sse1 / df2))
  }, numeric(1))
}

run_svg <- function(counts, coord, layer, truth_genes, rep_id, simulator, nbrs) {
  log_counts <- normalize_counts(counts)
  genes <- colnames(counts)
  names_layer <- factor(layer)
  design_layer <- stats::model.matrix(~ names_layer)
  coord_df <- data.frame(x = coord[, 1], y = coord[, 2])
  design_spatial <- stats::model.matrix(~ x + y + I(x^2) + I(y^2) + I(x * y), data = coord_df)
  score_list <- list(
    `Moran I top-k` = moran_scores(log_counts, nbrs),
    `layer ANOVA top-k` = lm_f_scores(log_counts, design_layer),
    `spatial polynomial top-k` = lm_f_scores(log_counts, design_spatial)
  )

  k <- length(truth_genes)
  rows <- lapply(names(score_list), function(method) {
    scores <- score_list[[method]]
    names(scores) <- genes
    scores[!is.finite(scores)] <- -Inf
    pred <- names(sort(scores, decreasing = TRUE))[seq_len(k)]
    met <- binary_metrics(pred, truth_genes, genes)
    metric_row("SVG detection", simulator, method, rep_id, met$F1, "higher", met)
  })
  do.call(rbind, rows)
}

run_ctsvg <- function(counts, prop, truth_by_ct, rep_id, simulator) {
  log_counts <- normalize_counts(counts)
  genes <- colnames(counts)
  dominant <- ct_names[max.col(prop[, ct_names, drop = FALSE], ties.method = "first")]

  method_scores <- list()
  for (ct in ct_names) {
    p <- prop[, ct]
    top <- p >= stats::quantile(p, 0.75)
    low <- p <= stats::quantile(p, 0.25)
    dom <- dominant == ct & p >= 0.5
    rest <- !dom
    if (sum(dom) < 5 || sum(rest) < 5) {
      dom <- top
      rest <- low
    }
    method_scores[[ct]] <- list(
      `proportion correlation top-k` = vapply(seq_len(ncol(log_counts)), function(j) {
        suppressWarnings(stats::cor(log_counts[, j], p, method = "spearman", use = "pairwise.complete.obs"))
      }, numeric(1)),
      `high-low proportion top-k` = colMeans(log_counts[top, , drop = FALSE]) - colMeans(log_counts[low, , drop = FALSE]),
      `dominant-region logFC top-k` = colMeans(log_counts[dom, , drop = FALSE]) - colMeans(log_counts[rest, , drop = FALSE])
    )
  }

  cell_rows <- list()
  for (ct in ct_names) {
    truth <- truth_by_ct[[ct]]
    k <- length(truth)
    for (method in names(method_scores[[ct]])) {
      scores <- method_scores[[ct]][[method]]
      names(scores) <- genes
      scores[!is.finite(scores)] <- -Inf
      pred <- names(sort(scores, decreasing = TRUE))[seq_len(k)]
      met <- binary_metrics(pred, truth, genes)
      cell_rows[[length(cell_rows) + 1]] <- cbind(
        data.frame(rep = rep_id, simulator = simulator, method = method, cell_type = ct, stringsAsFactors = FALSE),
        met
      )
    }
  }
  cell_dt <- do.call(rbind, cell_rows)
  rep_rows <- lapply(split(cell_dt, cell_dt$method), function(df) {
    extras <- data.frame(
      macro_precision = mean(df$precision),
      macro_TPR = mean(df$TPR),
      macro_FPR = mean(df$FPR),
      macro_F1 = mean(df$F1),
      macro_Jaccard = mean(df$Jaccard),
      stringsAsFactors = FALSE
    )
    metric_row("CT-SVG detection", simulator, df$method[1], rep_id, extras$macro_F1, "higher", extras)
  })
  list(metrics = do.call(rbind, rep_rows), cell_metrics = cell_dt)
}

reference_profiles <- function(scref, genes) {
  ref_counts <- data.matrix(scref[, genes, drop = FALSE])
  profiles <- sapply(ct_names, function(ct) {
    idx <- scref$Cell_type == ct
    colMeans(ref_counts[idx, , drop = FALSE], na.rm = TRUE)
  })
  rownames(profiles) <- genes
  profiles <- apply(profiles + 1e-8, 2, function(v) v / sum(v))
  profiles
}

estimate_ls_props <- function(counts, profiles, genes) {
  genes <- intersect(genes, intersect(colnames(counts), rownames(profiles)))
  x <- profiles[genes, ct_names, drop = FALSE]
  y_mat <- counts[, genes, drop = FALSE]
  y_mat <- sweep(y_mat + 1e-8, 1, rowSums(y_mat + 1e-8), "/")
  qrx <- qr(x)
  pred <- matrix(0, nrow = nrow(y_mat), ncol = length(ct_names), dimnames = list(rownames(counts), ct_names))
  for (i in seq_len(nrow(y_mat))) {
    coef <- tryCatch(qr.coef(qrx, y_mat[i, ]), error = function(e) rep(0, length(ct_names)))
    coef[!is.finite(coef)] <- 0
    pred[i, ] <- pmax(coef, 0)
  }
  normalize_rows(pred)
}

estimate_marker_score_props <- function(counts) {
  log_counts <- normalize_counts(counts)
  score <- sapply(ct_names, function(ct) {
    genes <- intersect(marker_genes[[ct]], colnames(log_counts))
    rowMeans(log_counts[, genes, drop = FALSE])
  })
  score <- sweep(score, 1, apply(score, 1, max), "-")
  normalize_rows(exp(score))
}

run_deconvolution <- function(counts, prop_truth, scref, rep_id, simulator) {
  genes <- colnames(counts)
  profiles <- reference_profiles(scref, intersect(genes, colnames(scref)))
  marker_only <- unique(unlist(marker_genes))
  preds <- list(
    `all-gene reference LS` = estimate_ls_props(counts, profiles, genes),
    `marker-gene reference LS` = estimate_ls_props(counts, profiles, marker_only),
    `marker-score simplex` = estimate_marker_score_props(counts)
  )

  metric_rows <- list()
  cell_rows <- list()
  truth <- as.matrix(prop_truth[, ct_names, drop = FALSE])
  rownames(truth) <- rownames(counts)
  for (method in names(preds)) {
    pred <- preds[[method]][, ct_names, drop = FALSE]
    err <- pred - truth
    overall_rmse <- sqrt(mean(err^2))
    extras <- data.frame(RMSE = overall_rmse, stringsAsFactors = FALSE)
    metric_rows[[length(metric_rows) + 1]] <- metric_row(
      "deconvolution", simulator, method, rep_id, overall_rmse, "lower", extras
    )
    for (ct in ct_names) {
      cell_rows[[length(cell_rows) + 1]] <- data.frame(
        rep = rep_id,
        simulator = simulator,
        method = method,
        cell_type = ct,
        RMSE = sqrt(mean(err[, ct]^2)),
        correlation = suppressWarnings(stats::cor(pred[, ct], truth[, ct], use = "pairwise.complete.obs")),
        stringsAsFactors = FALSE
      )
    }
  }
  list(metrics = do.call(rbind, metric_rows), cell_metrics = do.call(rbind, cell_rows))
}

edge_id <- function(a, b, gene_rank) {
  swap <- gene_rank[a] > gene_rank[b]
  first <- ifelse(swap, b, a)
  second <- ifelse(swap, a, b)
  paste(first, second, sep = "--")
}

truth_edges_from_modules <- function(truth_dt, gene_universe) {
  gene_rank <- setNames(seq_along(gene_universe), gene_universe)
  modules <- split(truth_dt$gene_id, paste(truth_dt$cell_type, truth_dt$module_id, sep = "::"))
  edges <- unlist(lapply(modules, function(genes) {
    genes <- unique(intersect(genes, gene_universe))
    if (length(genes) < 2) {
      return(character())
    }
    cmb <- utils::combn(genes, 2)
    edge_id(cmb[1, ], cmb[2, ], gene_rank)
  }), use.names = FALSE)
  unique(edges)
}

all_gene_edges <- function(gene_universe) {
  cmb <- utils::combn(gene_universe, 2)
  paste(cmb[1, ], cmb[2, ], sep = "--")
}

top_edges_from_matrix <- function(mat, method, n_edges, gene_universe) {
  cors <- suppressWarnings(stats::cor(mat[, gene_universe, drop = FALSE], method = method, use = "pairwise.complete.obs"))
  cors[!is.finite(cors)] <- 0
  diag(cors) <- 0
  idx <- which(upper.tri(cors), arr.ind = TRUE)
  score <- abs(cors[idx])
  genes <- colnames(cors)
  edges <- paste(genes[idx[, 1]], genes[idx[, 2]], sep = "--")
  edges[order(score, decreasing = TRUE)][seq_len(n_edges)]
}

residualize_expression <- function(log_counts, layer, coord) {
  design <- stats::model.matrix(~ factor(layer) + coord[, 1] + coord[, 2])
  res <- matrix(0, nrow = nrow(log_counts), ncol = ncol(log_counts), dimnames = dimnames(log_counts))
  for (j in seq_len(ncol(log_counts))) {
    res[, j] <- stats::lm.fit(design, log_counts[, j])$residuals
  }
  res
}

run_coexpression <- function(counts, layer, coord, nbrs, truth_dt, rep_id, simulator) {
  log_counts <- normalize_counts(counts)
  genes <- colnames(counts)
  truth_edges <- truth_edges_from_modules(truth_dt, genes)
  universe <- all_gene_edges(genes)
  n_edges <- length(truth_edges)
  mats <- list(
    `Pearson correlation top-k` = list(mat = log_counts, method = "pearson"),
    `Spearman correlation top-k` = list(mat = log_counts, method = "spearman"),
    `spatial-smoothed Pearson top-k` = list(mat = knn_smooth(log_counts, nbrs), method = "pearson"),
    `layer/coordinate residual Pearson top-k` = list(mat = residualize_expression(log_counts, layer, coord), method = "pearson")
  )
  rows <- lapply(names(mats), function(method_name) {
    pred <- top_edges_from_matrix(mats[[method_name]]$mat, mats[[method_name]]$method, n_edges, genes)
    met <- binary_metrics(pred, truth_edges, universe)
    met$n_truth_edges <- n_edges
    metric_row("co-expression", simulator, method_name, rep_id, met$F1, "higher", met)
  })
  do.call(rbind, rows)
}

simulate_spot_srtsim_counts <- function(counts, coord, labels, rep_id) {
  if (!requireNamespace("SRTsim", quietly = TRUE)) {
    stop("SRTsim counts are missing and the SRTsim package could not be loaded.")
  }
  set.seed(seed + rep_id)
  rounded_counts <- round(pmax(counts, 0))
  loc <- data.frame(y = coord[, 1], x = coord[, 2], label = as.character(labels), row.names = rownames(counts))
  sim_srt <- SRTsim::createSRT(count_in = t(rounded_counts), loc_in = loc)
  sim_srt <- SRTsim::srtsim_fit(sim_srt, sim_scheme = "tissue", min_nonzero_num = 2, maxiter = 500)
  sim_srt <- SRTsim::srtsim_count(sim_srt, numCores = 1, verbose = FALSE)
  sim_counts <- t(as.matrix(SRTsim::simCounts(sim_srt)))
  colnames(sim_counts) <- rownames(SRTsim::simCounts(sim_srt))
  rownames(sim_counts) <- rownames(counts)
  sim_counts
}

simulate_ct_aggregate_srtsim_counts <- function(scref, coord, labels, prop_truth, rep_id) {
  if (!requireNamespace("SRTsim", quietly = TRUE)) {
    stop("SRTsim counts are missing and the SRTsim package could not be loaded.")
  }
  genes <- grep("^gene_", colnames(scref), value = TRUE)
  spot_ids <- paste0("spot_", seq_len(nrow(coord)))
  loc <- data.frame(y = coord[, 1], x = coord[, 2], label = as.character(labels), row.names = spot_ids)
  agg <- matrix(0, nrow = nrow(coord), ncol = length(genes), dimnames = list(spot_ids, genes))
  n_cells <- cell_per_spot[as.character(rep_id)]
  if (is.na(n_cells)) {
    n_cells <- 10
  }
  for (ct in ct_names) {
    message(sprintf("  SRTsim cell-type fit: rep %s, %s", rep_id, ct))
    ct_mat <- round(pmax(data.matrix(scref[scref$Cell_type == ct, genes, drop = FALSE]), 0))
    rownames(ct_mat) <- spot_ids
    colnames(ct_mat) <- genes
    set.seed(seed + rep_id * 1000 + match(ct, ct_names))
    sim_srt <- SRTsim::createSRT(count_in = t(ct_mat), loc_in = loc)
    sim_srt <- SRTsim::srtsim_fit(sim_srt, sim_scheme = "tissue", min_nonzero_num = 2, maxiter = 500)
    sim_srt <- SRTsim::srtsim_count(sim_srt, numCores = 1, verbose = FALSE)
    ct_counts <- t(as.matrix(SRTsim::simCounts(sim_srt)))
    colnames(ct_counts) <- rownames(SRTsim::simCounts(sim_srt))
    rownames(ct_counts) <- spot_ids
    agg <- agg + ct_counts[, genes, drop = FALSE] * as.numeric(n_cells) * prop_truth[, ct]
  }
  agg
}

truth_marker <- read.csv(file.path(repo_dir, "Data", "S7_gene", "simulated", "S7_marker_gene_truth.csv"), check.names = FALSE)
truth_coexp <- read.csv(file.path(repo_dir, "Data", "S7_gene", "simulated", "S7_coexpression_truth.csv"), check.names = FALSE)

task_metrics <- list()
deconv_cell_metrics <- list()
ctsvg_cell_metrics <- list()

for (rep_id in reps) {
  message(sprintf("Running all downstream anchors for S7 replicate %s", rep_id))
  sim_dir <- file.path(repo_dir, "Data", "S7_gene", "simulated")
  counts_stcos <- read_count_matrix(file.path(sim_dir, sprintf("S7_gene_expression_%s.csv", rep_id)))
  coord_df <- read.csv(file.path(sim_dir, sprintf("S7_layer_%s.csv", rep_id)), check.names = FALSE)
  region_df <- read.csv(file.path(sim_dir, sprintf("S7_region_%s.csv", rep_id)), check.names = FALSE)
  scref <- read.csv(file.path(sim_dir, sprintf("S7_scref_%s.csv", rep_id)), check.names = FALSE)
  coord <- as.matrix(coord_df[, c("x_coord", "y_coord")])
  layer <- coord_df$layer
  prop_truth <- region_df[, ct_names, drop = FALSE]
  rownames(prop_truth) <- rownames(counts_stcos)
  nbrs <- neighbor_indices(coord, k = 6)

  srtsim_file <- file.path(srtsim_dir, if (srtsim_mode == "ct-aggregate") {
    sprintf("S7_CTaggregate_SRTsim_counts_rep%s.csv", rep_id)
  } else {
    sprintf("S7_SRTsim_counts_rep%s.csv", rep_id)
  })
  if (!file.exists(srtsim_file)) {
    if (!simulate_missing) {
      stop("Missing SRTsim count file: ", srtsim_file)
    }
    if (srtsim_mode == "ct-aggregate") {
      counts_srtsim <- simulate_ct_aggregate_srtsim_counts(scref, coord, layer, prop_truth, rep_id)
    } else if (srtsim_mode == "spot") {
      counts_srtsim <- simulate_spot_srtsim_counts(counts_stcos, coord, layer, rep_id)
    } else {
      stop("Unknown --srtsim-mode: ", srtsim_mode)
    }
    dir.create(srtsim_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(counts_srtsim, srtsim_file, row.names = FALSE)
  } else {
    counts_srtsim <- read_count_matrix(srtsim_file)
  }

  truth_rep <- truth_marker[truth_marker$scenario == rep_id, , drop = FALSE]
  svg_truth <- unique(truth_rep$gene_id)
  ctsvg_truth <- split(truth_rep$gene_id, truth_rep$cell_type)
  coexp_rep <- truth_coexp[truth_coexp$scenario == rep_id, , drop = FALSE]

  datasets <- list(`ST-COS S7` = counts_stcos)
  datasets[[independent_label]] <- counts_srtsim
  for (simulator in names(datasets)) {
    counts <- datasets[[simulator]]
    task_metrics[[length(task_metrics) + 1]] <- run_segmentation(counts, coord, layer, nbrs, rep_id, simulator, seed)
    deconv <- run_deconvolution(counts, prop_truth, scref, rep_id, simulator)
    task_metrics[[length(task_metrics) + 1]] <- deconv$metrics
    deconv_cell_metrics[[length(deconv_cell_metrics) + 1]] <- deconv$cell_metrics
    task_metrics[[length(task_metrics) + 1]] <- run_svg(counts, coord, layer, svg_truth, rep_id, simulator, nbrs)
    ctsvg <- run_ctsvg(counts, prop_truth, ctsvg_truth, rep_id, simulator)
    task_metrics[[length(task_metrics) + 1]] <- ctsvg$metrics
    ctsvg_cell_metrics[[length(ctsvg_cell_metrics) + 1]] <- ctsvg$cell_metrics
    task_metrics[[length(task_metrics) + 1]] <- run_coexpression(counts, layer, coord, nbrs, coexp_rep, rep_id, simulator)
  }
}

task_metrics <- rbind_fill(task_metrics)
deconv_cell_metrics <- rbind_fill(deconv_cell_metrics)
ctsvg_cell_metrics <- rbind_fill(ctsvg_cell_metrics)

write.csv(task_metrics, file.path(out_dir, "SRTsim_S7_all_downstream_metrics_by_rep.csv"), row.names = FALSE)
write.csv(deconv_cell_metrics, file.path(out_dir, "SRTsim_S7_deconvolution_celltype_metrics.csv"), row.names = FALSE)
write.csv(ctsvg_cell_metrics, file.path(out_dir, "SRTsim_S7_ctsvg_celltype_metrics.csv"), row.names = FALSE)

summary_stats <- do.call(rbind, lapply(split(task_metrics, list(task_metrics$task, task_metrics$simulator, task_metrics$method), drop = TRUE), function(df) {
  data.frame(
    task = df$task[1],
    simulator = df$simulator[1],
    method = df$method[1],
    direction = df$direction[1],
    n_reps = nrow(df),
    primary_mean = mean(df$primary_metric, na.rm = TRUE),
    primary_sd = stats::sd(df$primary_metric, na.rm = TRUE),
    primary_se = stats::sd(df$primary_metric, na.rm = TRUE) / sqrt(nrow(df)),
    stringsAsFactors = FALSE
  )
}))
summary_stats <- summary_stats[order(summary_stats$task, summary_stats$simulator, summary_stats$method), ]
write.csv(summary_stats, file.path(out_dir, "SRTsim_S7_all_downstream_summary.csv"), row.names = FALSE)

rank_agreement <- function(summary_dt) {
  rows <- list()
  for (task in unique(summary_dt$task)) {
    sub <- summary_dt[summary_dt$task == task, ]
    st <- sub[sub$simulator == "ST-COS S7", c("method", "primary_mean", "direction")]
    sr <- sub[sub$simulator == independent_label, c("method", "primary_mean", "direction")]
    comp <- merge(st, sr, by = "method", suffixes = c("_STCOS", "_SRTsim"))
    if (nrow(comp) == 0) {
      next
    }
    direction <- comp$direction_STCOS[1]
    st_score <- if (direction == "lower") -comp$primary_mean_STCOS else comp$primary_mean_STCOS
    sr_score <- if (direction == "lower") -comp$primary_mean_SRTsim else comp$primary_mean_SRTsim
    rows[[length(rows) + 1]] <- data.frame(
      task = task,
      comparison = "method primary-metric rank",
      n_items = nrow(comp),
      spearman = suppressWarnings(stats::cor(st_score, sr_score, method = "spearman")),
      best_STCOS = comp$method[which.max(st_score)],
      best_SRTsim = comp$method[which.max(sr_score)],
      worst_STCOS = comp$method[which.min(st_score)],
      worst_SRTsim = comp$method[which.min(sr_score)],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

agreement <- rank_agreement(summary_stats)

deconv_ct_summary <- aggregate(RMSE ~ simulator + method + cell_type, deconv_cell_metrics, mean)
dc <- merge(
  deconv_ct_summary[deconv_ct_summary$simulator == "ST-COS S7", c("method", "cell_type", "RMSE")],
  deconv_ct_summary[deconv_ct_summary$simulator == independent_label, c("method", "cell_type", "RMSE")],
  by = c("method", "cell_type"),
  suffixes = c("_STCOS", "_SRTsim")
)
agreement <- rbind(agreement, data.frame(
  task = "deconvolution",
  comparison = "method-cell-type RMSE rank",
  n_items = nrow(dc),
  spearman = suppressWarnings(stats::cor(-dc$RMSE_STCOS, -dc$RMSE_SRTsim, method = "spearman")),
  best_STCOS = paste(dc$method[which.min(dc$RMSE_STCOS)], dc$cell_type[which.min(dc$RMSE_STCOS)], sep = " / "),
  best_SRTsim = paste(dc$method[which.min(dc$RMSE_SRTsim)], dc$cell_type[which.min(dc$RMSE_SRTsim)], sep = " / "),
  worst_STCOS = paste(dc$method[which.max(dc$RMSE_STCOS)], dc$cell_type[which.max(dc$RMSE_STCOS)], sep = " / "),
  worst_SRTsim = paste(dc$method[which.max(dc$RMSE_SRTsim)], dc$cell_type[which.max(dc$RMSE_SRTsim)], sep = " / "),
  stringsAsFactors = FALSE
))

ctsvg_summary <- aggregate(F1 ~ simulator + method + cell_type, ctsvg_cell_metrics, mean)
cc <- merge(
  ctsvg_summary[ctsvg_summary$simulator == "ST-COS S7", c("method", "cell_type", "F1")],
  ctsvg_summary[ctsvg_summary$simulator == independent_label, c("method", "cell_type", "F1")],
  by = c("method", "cell_type"),
  suffixes = c("_STCOS", "_SRTsim")
)
agreement <- rbind(agreement, data.frame(
  task = "CT-SVG detection",
  comparison = "method-cell-type F1 rank",
  n_items = nrow(cc),
  spearman = suppressWarnings(stats::cor(cc$F1_STCOS, cc$F1_SRTsim, method = "spearman")),
  best_STCOS = paste(cc$method[which.max(cc$F1_STCOS)], cc$cell_type[which.max(cc$F1_STCOS)], sep = " / "),
  best_SRTsim = paste(cc$method[which.max(cc$F1_SRTsim)], cc$cell_type[which.max(cc$F1_SRTsim)], sep = " / "),
  worst_STCOS = paste(cc$method[which.min(cc$F1_STCOS)], cc$cell_type[which.min(cc$F1_STCOS)], sep = " / "),
  worst_SRTsim = paste(cc$method[which.min(cc$F1_SRTsim)], cc$cell_type[which.min(cc$F1_SRTsim)], sep = " / "),
  stringsAsFactors = FALSE
))

write.csv(agreement, file.path(out_dir, "SRTsim_S7_all_downstream_failure_mode_agreement.csv"), row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  p <- ggplot(task_metrics, aes(x = method, y = primary_metric, fill = simulator)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
    geom_point(
      aes(color = simulator),
      position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65),
      size = 1.2,
      alpha = 0.85
    ) +
    coord_flip() +
    facet_wrap(~ task, scales = "free", ncol = 2) +
    labs(x = NULL, y = "Primary metric", fill = NULL, color = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "SRTsim_S7_all_downstream_primary_metrics.pdf"), p, width = 11, height = 8)
}

message("Finished. Results written to: ", out_dir)
