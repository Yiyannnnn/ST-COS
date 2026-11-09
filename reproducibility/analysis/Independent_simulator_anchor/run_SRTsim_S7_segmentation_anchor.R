#!/usr/bin/env Rscript

Sys.setenv(RGL_USE_NULL = "TRUE")

suppressPackageStartupMessages({
  library(SRTsim)
  library(aricode)
})

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

repo_dir <- normalizePath(parse_arg("repo", "."), mustWork = TRUE)
out_dir <- file.path(repo_dir, parse_arg("out", "Independent_simulator_anchor/SRTsim_S7_segmentation_anchor"))
reps <- parse_reps(parse_arg("reps", "1:5"))
maxiter <- as.integer(parse_arg("maxiter", "500"))
seed <- as.integer(parse_arg("seed", "20260902"))
save_counts <- tolower(parse_arg("save-counts", "true")) %in% c("true", "t", "1", "yes", "y")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normalize_counts <- function(counts) {
  lib_size <- rowSums(counts)
  target <- stats::median(lib_size[lib_size > 0])
  if (!is.finite(target)) {
    target <- 1
  }
  scaled <- counts / pmax(lib_size, 1) * target
  log1p(scaled)
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

knn_smooth <- function(mat, coord, k = 6) {
  d <- as.matrix(stats::dist(coord))
  diag(d) <- Inf
  nbrs <- apply(d, 1, order)[seq_len(k), , drop = FALSE]
  out <- mat
  for (i in seq_len(nrow(mat))) {
    out[i, ] <- colMeans(mat[c(i, nbrs[, i]), , drop = FALSE])
  }
  out
}

boundary_mask <- function(coord, labels, k = 6) {
  d <- as.matrix(stats::dist(coord))
  diag(d) <- Inf
  nbrs <- apply(d, 1, order)[seq_len(k), , drop = FALSE]
  vapply(seq_along(labels), function(i) any(labels[nbrs[, i]] != labels[i]), logical(1))
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
  padded_truth <- truth_levels
  maps <- permute(padded_truth)
  best <- 0
  for (m in maps) {
    mapping <- m[seq_along(pred_levels)]
    names(mapping) <- pred_levels
    best <- max(best, mean(mapping[pred] == truth))
  }
  best
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

run_simple_segmentation <- function(counts, coord, labels, rep_id, simulator, seed) {
  labels <- as.character(labels)
  k <- length(unique(labels))
  boundary <- boundary_mask(coord, labels)
  log_counts <- normalize_counts(counts)
  expr_pc <- pca_features(log_counts)
  smooth_pc <- pca_features(knn_smooth(log_counts, coord, k = 6))
  coord_features <- scale(coord)

  set.seed(seed + rep_id * 100 + 1)
  km_expr <- stats::kmeans(expr_pc, centers = k, nstart = 50, iter.max = 100)$cluster

  set.seed(seed + rep_id * 100 + 2)
  km_smooth <- stats::kmeans(smooth_pc, centers = k, nstart = 50, iter.max = 100)$cluster

  set.seed(seed + rep_id * 100 + 3)
  km_coord <- stats::kmeans(coord_features, centers = k, nstart = 50, iter.max = 100)$cluster

  hc_expr <- stats::cutree(stats::hclust(stats::dist(expr_pc), method = "ward.D2"), k = k)

  preds <- list(
    `PCA k-means` = km_expr,
    `spatial-smoothed PCA k-means` = km_smooth,
    `coordinate k-means` = km_coord,
    `PCA Ward clustering` = hc_expr
  )

  rows <- lapply(names(preds), function(method) {
    cbind(
      data.frame(rep = rep_id, simulator = simulator, method = method, stringsAsFactors = FALSE),
      score_clustering(labels, preds[[method]], boundary),
      data.frame(
        n_spots = length(labels),
        n_boundary = sum(boundary),
        n_layers = k,
        stringsAsFactors = FALSE
      )
    )
  })
  do.call(rbind, rows)
}

simulate_srtsim <- function(counts, coord, labels, rep_id, seed, maxiter) {
  set.seed(seed + rep_id)
  spot_ids <- paste0("spot_", seq_len(nrow(counts)))
  gene_ids <- colnames(counts)
  rounded_counts <- round(pmax(counts, 0))
  rownames(rounded_counts) <- spot_ids
  colnames(rounded_counts) <- gene_ids

  loc <- data.frame(
    y = coord[, 1],
    x = coord[, 2],
    label = as.character(labels),
    row.names = spot_ids
  )

  sim_srt <- createSRT(count_in = t(rounded_counts), loc_in = loc)
  sim_srt <- srtsim_fit(
    sim_srt,
    sim_scheme = "tissue",
    min_nonzero_num = 2,
    maxiter = maxiter
  )
  sim_srt <- srtsim_count(sim_srt, numCores = 1, verbose = FALSE)
  sim_counts <- as.matrix(simCounts(sim_srt))
  sim_counts <- t(sim_counts)
  rownames(sim_counts) <- spot_ids
  colnames(sim_counts) <- rownames(simCounts(sim_srt))
  sim_counts
}

all_metrics <- list()
timing <- list()

for (rep_id in reps) {
  message(sprintf("Running S7 replicate %s", rep_id))
  expr_file <- file.path(repo_dir, "Data", "S7_gene", "simulated", sprintf("S7_gene_expression_%s.csv", rep_id))
  layer_file <- file.path(repo_dir, "Data", "S7_gene", "simulated", sprintf("S7_layer_%s.csv", rep_id))

  counts <- as.matrix(read.csv(expr_file, check.names = FALSE))
  coord_df <- read.csv(layer_file, check.names = FALSE)
  coord <- as.matrix(coord_df[, c("x_coord", "y_coord")])
  labels <- coord_df$layer

  spot_ids <- paste0("spot_", seq_len(nrow(counts)))
  rownames(counts) <- spot_ids
  rownames(coord) <- spot_ids

  all_metrics[[length(all_metrics) + 1]] <- run_simple_segmentation(
    counts = counts,
    coord = coord,
    labels = labels,
    rep_id = rep_id,
    simulator = "ST-COS S7",
    seed = seed
  )

  t0 <- proc.time()[["elapsed"]]
  srt_counts <- simulate_srtsim(counts, coord, labels, rep_id, seed, maxiter)
  elapsed <- proc.time()[["elapsed"]] - t0
  timing[[length(timing) + 1]] <- data.frame(rep = rep_id, simulator = "SRTsim on S7", seconds = elapsed)

  if (save_counts) {
    write.csv(
      srt_counts,
      file = file.path(out_dir, sprintf("S7_SRTsim_counts_rep%s.csv", rep_id)),
      row.names = FALSE
    )
  }

  all_metrics[[length(all_metrics) + 1]] <- run_simple_segmentation(
    counts = srt_counts,
    coord = coord,
    labels = labels,
    rep_id = rep_id,
    simulator = "SRTsim on S7",
    seed = seed
  )
}

metrics <- do.call(rbind, all_metrics)
write.csv(metrics, file.path(out_dir, "SRTsim_S7_segmentation_anchor_metrics.csv"), row.names = FALSE)

timing <- do.call(rbind, timing)
write.csv(timing, file.path(out_dir, "SRTsim_S7_runtime.csv"), row.names = FALSE)

summary_stats <- do.call(rbind, lapply(split(metrics, list(metrics$simulator, metrics$method), drop = TRUE), function(df) {
  numeric_cols <- c("ARI", "NMI", "accuracy", "boundary_accuracy", "interior_accuracy", "boundary_penalty")
  out <- data.frame(
    simulator = df$simulator[1],
    method = df$method[1],
    n_reps = nrow(df),
    stringsAsFactors = FALSE
  )
  for (col in numeric_cols) {
    out[[paste0(col, "_mean")]] <- mean(df[[col]], na.rm = TRUE)
    out[[paste0(col, "_sd")]] <- stats::sd(df[[col]], na.rm = TRUE)
    out[[paste0(col, "_se")]] <- stats::sd(df[[col]], na.rm = TRUE) / sqrt(nrow(df))
  }
  out
}))
summary_stats <- summary_stats[order(summary_stats$simulator, -summary_stats$ARI_mean), ]
write.csv(summary_stats, file.path(out_dir, "SRTsim_S7_segmentation_anchor_summary.csv"), row.names = FALSE)

stcos_summary <- summary_stats[summary_stats$simulator == "ST-COS S7", ]
srtsim_summary <- summary_stats[summary_stats$simulator == "SRTsim on S7", ]
comparison <- merge(
  stcos_summary,
  srtsim_summary,
  by = "method",
  suffixes = c("_STCOS", "_SRTsim")
)
comparison$ARI_mean_delta_SRTsim_minus_STCOS <- comparison$ARI_mean_SRTsim - comparison$ARI_mean_STCOS
comparison$boundary_penalty_delta_SRTsim_minus_STCOS <- comparison$boundary_penalty_mean_SRTsim - comparison$boundary_penalty_mean_STCOS
comparison <- comparison[order(-comparison$ARI_mean_STCOS), ]
write.csv(comparison, file.path(out_dir, "SRTsim_S7_failure_mode_comparison.csv"), row.names = FALSE)

expression_comparison <- comparison[comparison$method != "coordinate k-means", ]
spearman <- function(df, lhs, rhs) {
  as.character(suppressWarnings(stats::cor(df[[lhs]], df[[rhs]], method = "spearman")))
}

agreement <- data.frame(
  statistic = c(
    "ARI rank Spearman across simple methods",
    "ARI rank Spearman across expression-based methods",
    "boundary penalty Spearman across simple methods",
    "boundary penalty Spearman across expression-based methods",
    "best simple method on ST-COS S7",
    "best simple method on SRTsim-on-S7",
    "lowest simple method on ST-COS S7",
    "lowest simple method on SRTsim-on-S7",
    "mean SRTsim runtime seconds"
  ),
  value = c(
    spearman(comparison, "ARI_mean_STCOS", "ARI_mean_SRTsim"),
    spearman(expression_comparison, "ARI_mean_STCOS", "ARI_mean_SRTsim"),
    spearman(comparison, "boundary_penalty_mean_STCOS", "boundary_penalty_mean_SRTsim"),
    spearman(expression_comparison, "boundary_penalty_mean_STCOS", "boundary_penalty_mean_SRTsim"),
    stcos_summary$method[which.max(stcos_summary$ARI_mean)],
    srtsim_summary$method[which.max(srtsim_summary$ARI_mean)],
    stcos_summary$method[which.min(stcos_summary$ARI_mean)],
    srtsim_summary$method[which.min(srtsim_summary$ARI_mean)],
    as.character(mean(timing$seconds))
  ),
  stringsAsFactors = FALSE
)
write.csv(agreement, file.path(out_dir, "SRTsim_S7_anchor_agreement.csv"), row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  metrics$method <- factor(metrics$method, levels = unique(summary_stats$method[order(-summary_stats$ARI_mean)]))
  p <- ggplot(metrics, aes(x = method, y = ARI, fill = simulator)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
    geom_point(
      aes(color = simulator),
      position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65),
      size = 1.5,
      alpha = 0.85
    ) +
    coord_flip() +
    labs(x = NULL, y = "Adjusted Rand index", fill = NULL, color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "SRTsim_S7_segmentation_anchor_ARI.pdf"), p, width = 8, height = 4.8)

  p2 <- ggplot(comparison, aes(x = ARI_mean_STCOS, y = ARI_mean_SRTsim, label = method)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(size = 2.5) +
    geom_text(nudge_y = 0.03, size = 3) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Mean ARI on ST-COS S7", y = "Mean ARI on SRTsim-on-S7") +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, "SRTsim_S7_failure_mode_comparison.pdf"), p2, width = 5.2, height = 5.2)
}

message("Finished. Results written to: ", out_dir)
