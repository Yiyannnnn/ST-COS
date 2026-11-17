# Shared evaluation rules for the simulated benchmarks (Supplementary Methods S2).
#
# The S1-S7 evaluators source this file, so truth definitions, decision rules,
# alignment, and metrics are implemented once.
# Paths are relative to the analysis root that contains `Data/` and the task
# output folders.

suppressPackageStartupMessages(library(data.table))

EVAL_FDR <- 0.05
COEXP_EDGE_THRESHOLD <- 0.1
COEXP_CONTINUOUS_METHODS <- c("spMOCA", "SpaceX", "CSSN", "Smoothie", "binSpect")

# Scenarios whose truth files describe replicate-invariant design constants:
# the S7 generator fixes svg_gene_sets, marker_gene_sets, and
# background_gene_sets in code.
REPLICATE_INVARIANT_TRUTH <- "S7"

## Identifiers -----------------------------------------------------------------

gene_index_from_names <- function(x) {
  x <- as.character(x)
  hit <- !is.na(x) & grepl("^gene[_.-]?[0-9]+$", x, ignore.case = TRUE)
  out <- rep(NA_integer_, length(x))
  out[hit] <- as.integer(sub("^gene[_.-]?", "", x[hit], ignore.case = TRUE))
  out
}

sim_dir_for <- function(base_dir, scenario) {
  file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated")
}

## Truth -------------------------------------------------------------------------

# Reads `<scenario>_<kind>_<rep>.csv` when present, otherwise the combined
# `<scenario>_<kind>.csv` restricted to the replicate. Separate replicate jobs
# overwrite the combined S7 file, so for replicate-invariant designs the rows
# of one available replicate are used and the source is recorded.
read_replicate_truth <- function(sim_dir, scenario, kind, rep_id) {
  per_rep <- file.path(sim_dir, sprintf("%s_%s_%s.csv", scenario, kind, rep_id))
  if (file.exists(per_rep)) {
    dt <- fread(per_rep)
    attr(dt, "truth_source") <- basename(per_rep)
    return(dt)
  }
  combined <- file.path(sim_dir, sprintf("%s_%s.csv", scenario, kind))
  if (!file.exists(combined)) return(NULL)
  dt <- fread(combined)
  attr(dt, "truth_source") <- basename(combined)
  if (!"scenario" %in% names(dt)) return(dt)

  hit <- dt[dt$scenario == rep_id]
  if (nrow(hit) > 0L) {
    attr(hit, "truth_source") <- basename(combined)
    return(hit)
  }
  if (!scenario %in% REPLICATE_INVARIANT_TRUTH) return(NULL)
  available <- sort(unique(dt$scenario))
  if (length(available) == 0L) return(NULL)
  hit <- dt[dt$scenario == available[1L]]
  attr(hit, "truth_source") <- sprintf(
    "%s (replicate %s rows; replicate-invariant design truth)",
    basename(combined), available[1L]
  )
  hit
}

svg_truth_genes <- function(sim_dir, scenario, rep_id) {
  dt <- read_replicate_truth(sim_dir, scenario, "marker_gene_truth", rep_id)
  if (!is.null(dt) && "gene_id" %in% names(dt)) {
    return(unique(stats::na.omit(as.character(dt$gene_id))))
  }
  # S1-S6 generators: hotspot genes 1-15 (genes 1-5, 6-10, and 11-15 for
  # ct1, ct2, and ct3). S2 and S5 contain no gene-level spatial pattern.
  if (scenario %in% c("S2", "S5")) return(character(0))
  if (scenario %in% c("S1", "S3", "S4", "S6")) return(paste0("gene_", 1:15))
  NULL
}

ctsvg_truth_sets <- function(sim_dir, scenario, rep_id) {
  dt <- read_replicate_truth(sim_dir, scenario, "marker_gene_truth", rep_id)
  if (!is.null(dt) && all(c("cell_type", "gene_id") %in% names(dt))) {
    return(lapply(split(as.character(dt$gene_id), dt$cell_type), unique))
  }
  if (scenario %in% c("S2", "S5")) {
    return(list(ct1 = character(0), ct2 = character(0), ct3 = character(0)))
  }
  if (scenario %in% c("S1", "S3", "S4", "S6")) {
    return(list(ct1 = paste0("gene_", 1:5), ct2 = paste0("gene_", 6:10),
                ct3 = paste0("gene_", 11:15)))
  }
  NULL
}

# Layer-level composition implied by the S1-S6 generators. S2, S3, S5, and S6
# use layer_weights ct1 = (1/2, 0, 0), ct2 = (1/4, 3/4, 0), ct3 = (1/3, 0, 1),
# which generate_cell_prop() normalizes row-wise; S1 and S4 are uniform.
scenario_layer_proportions <- function(scenario, layer) {
  n <- length(layer)
  if (scenario %in% c("S1", "S4")) {
    return(data.frame(ct1 = rep(1 / 3, n), ct2 = rep(1 / 3, n), ct3 = rep(1 / 3, n)))
  }
  if (!scenario %in% c("S2", "S3", "S5", "S6")) {
    stop("No layer-level composition is defined for ", scenario, ".")
  }
  weights <- rbind(c(1 / 2, 1 / 4, 1 / 3), c(0, 3 / 4, 0), c(0, 0, 1))
  weights <- weights / rowSums(weights)
  idx <- as.integer(factor(layer))
  if (max(idx) > nrow(weights)) stop("Unexpected layer label for ", scenario, ".")
  out <- as.data.frame(weights[idx, , drop = FALSE])
  names(out) <- c("ct1", "ct2", "ct3")
  out
}

# Realized proportions from the region file when available; otherwise the
# generator-defined layer composition.
truth_proportions <- function(sim_dir, scenario, rep_id) {
  region_file <- file.path(sim_dir, sprintf("%s_region_%s.csv", scenario, rep_id))
  layer_file <- file.path(sim_dir, sprintf("%s_layer_%s.csv", scenario, rep_id))
  if (file.exists(region_file)) {
    region <- fread(region_file)
    ct_cols <- setdiff(names(region), c("x_coord", "y_coord", "region"))
    truth <- as.data.frame(region[, ..ct_cols])
  } else if (file.exists(layer_file)) {
    truth <- scenario_layer_proportions(scenario, fread(layer_file)$layer)
  } else {
    stop(sprintf("No truth file found for %s rep %s.", scenario, rep_id))
  }
  rownames(truth) <- paste0("spot_", seq_len(nrow(truth)))
  truth
}

coexp_truth_matrix <- function(sim_dir, scenario, rep_id, n) {
  truth <- matrix(0L, nrow = n, ncol = n)
  dt <- read_replicate_truth(sim_dir, scenario, "coexpression_truth", rep_id)
  if (!is.null(dt)) {
    idx <- if ("gene_index" %in% names(dt)) {
      as.integer(dt$gene_index)
    } else {
      gene_index_from_names(dt$gene_id)
    }
    module <- if ("module_id" %in% names(dt)) dt$module_id else dt$module_label
    keep <- !is.na(idx) & idx >= 1L & idx <= n
    for (m in unique(module[keep])) {
      ids <- unique(idx[keep & module == m])
      if (length(ids) > 1L) truth[ids, ids] <- 1L
    }
    attr(truth, "truth_source") <- attr(dt, "truth_source")
  } else if (scenario %in% c("S4", "S5", "S6")) {
    # generate_covariance_group() blocks: genes 1-5 (ct1), 6-10 (ct2), 11-15 (ct3).
    for (ids in list(1:5, 6:10, 11:15)) {
      ids <- ids[ids <= n]
      if (length(ids) > 1L) truth[ids, ids] <- 1L
    }
  } else if (scenario %in% REPLICATE_INVARIANT_TRUTH) {
    return(NULL)
  }
  diag(truth) <- 0L
  truth
}

## Binary decision rules -----------------------------------------------------------

adjusted_calls <- function(q, alpha = EVAL_FDR) {
  q <- suppressWarnings(as.numeric(q))
  as.integer(!is.na(q) & q <= alpha)
}

bh_calls <- function(p, alpha = EVAL_FDR) {
  adjusted_calls(stats::p.adjust(suppressWarnings(as.numeric(p)), method = "BH"), alpha)
}

# Each method's documented multiple-testing output at FDR 0.05; BH is applied
# when only unadjusted P values are available. Untested genes are not called.
svg_calls <- function(dt, method) {
  has <- function(col) col %in% names(dt)
  switch(
    method,
    moranI = if (has("pval_norm_fdr_bh")) adjusted_calls(dt$pval_norm_fdr_bh),
    HEARTSVG = if (has("p_adj")) adjusted_calls(dt$p_adj),
    SPARK = if (has("adjusted_pvalue")) adjusted_calls(dt$adjusted_pvalue),
    SPARKX = if (has("adjustedPval")) adjusted_calls(dt$adjustedPval),
    SpatialDE = ,
    SpatialDE2 = if (has("qval")) {
      adjusted_calls(dt$qval)
    } else if (has("pval")) {
      bh_calls(dt$pval)
    },
    NULL
  )
}

## Binary gene-set metrics -----------------------------------------------------------

gene_set_metrics <- function(pred, truth, universe) {
  pred <- unique(as.character(pred))
  truth <- unique(as.character(truth))
  universe <- unique(c(as.character(universe), pred, truth))
  tp <- length(intersect(pred, truth))
  fp <- length(setdiff(pred, truth))
  fn <- length(setdiff(truth, pred))
  tn <- length(setdiff(universe, union(pred, truth)))
  binary_metrics(tp, fp, fn, tn)
}

# Metrics requiring positive truth are NA when none exists (S2 and S5 for
# SVG and CT-SVG; S1-S3 for co-expression).
binary_metrics <- function(tp, fp, fn, tn) {
  has_positive <- (tp + fn) > 0L
  precision <- if (!has_positive) NA_real_ else if ((tp + fp) > 0L) tp / (tp + fp) else NA_real_
  recall <- if (has_positive) tp / (tp + fn) else NA_real_
  f1 <- if (!has_positive) NA_real_ else if (tp == 0L) 0 else 2 * precision * recall / (precision + recall)
  list(
    tp = tp, fp = fp, fn = fn, tn = tn,
    has_positive_truth = has_positive,
    precision = precision,
    TPR = recall,
    FPR = if ((fp + tn) > 0L) fp / (fp + tn) else NA_real_,
    F1 = f1,
    Jaccard = if (has_positive) tp / (tp + fp + fn) else NA_real_
  )
}

## Co-expression ----------------------------------------------------------------------

# Reads a gene-by-gene matrix written with or without a row-name column and
# aligns it to genes 1..n by gene name when names are present.
read_gene_matrix <- function(path, n) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  if (nrow(dt) == 0L || ncol(dt) == 0L) return(NULL)
  row_ids <- NULL
  if (ncol(dt) == nrow(dt) + 1L) {
    row_ids <- as.character(dt[[1L]])
    dt <- dt[, -1L, with = FALSE]
  }
  m <- suppressWarnings(as.matrix(dt))
  storage.mode(m) <- "double"
  m[!is.finite(m)] <- 0

  col_idx <- gene_index_from_names(colnames(m))
  if (anyNA(col_idx) || anyDuplicated(col_idx)) {
    if (ncol(m) != n) {
      warning(sprintf("%s has %d unnamed columns for %d genes; aligned by position.",
                      basename(path), ncol(m), n), call. = FALSE)
    }
    col_idx <- seq_len(ncol(m))
  }
  row_idx <- if (is.null(row_ids)) NULL else gene_index_from_names(row_ids)
  if (is.null(row_idx) || anyNA(row_idx) || anyDuplicated(row_idx)) {
    if (nrow(m) != ncol(m)) {
      warning(sprintf("%s is not square and has no gene row names.", basename(path)),
              call. = FALSE)
      return(NULL)
    }
    row_idx <- col_idx
  }
  keep_r <- row_idx >= 1L & row_idx <= n
  keep_c <- col_idx >= 1L & col_idx <= n
  out <- matrix(0, nrow = n, ncol = n)
  out[row_idx[keep_r], col_idx[keep_c]] <- m[keep_r, keep_c, drop = FALSE]
  diag(out) <- 0
  out
}

membership_adjacency <- function(idx, label, n) {
  out <- matrix(0, nrow = n, ncol = n)
  keep <- !is.na(idx) & idx >= 1L & idx <= n & !is.na(label) & label != ""
  for (lab in unique(label[keep])) {
    ids <- unique(idx[keep & label == lab])
    if (length(ids) > 1L) out[ids, ids] <- 1
  }
  diag(out) <- 0
  out
}

load_coexp_prediction <- function(base_dir, scenario, rep_id, method, n) {
  coexp_dir <- file.path(base_dir, "coexpression")
  data_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"))
  stem <- paste0(scenario, "_", rep_id)

  if (method == "spMOCA") {
    return(read_gene_matrix(file.path(coexp_dir, "spMOCA_Agent", paste0(stem, "_correlation_matrix.csv")), n))
  }
  if (method == "SpaceX") {
    return(read_gene_matrix(file.path(data_dir, paste0("SpaceX_", rep_id, ".csv")), n))
  }
  if (method == "binSpect") {
    candidates <- file.path(coexp_dir, c("binSpect_agent", "binSpect_Agent"),
                            paste0(stem, "_correlation_matrix.csv"))
    hit <- candidates[file.exists(candidates)]
    return(if (length(hit)) read_gene_matrix(hit[1L], n) else NULL)
  }
  if (method == "CSSN") {
    f_mean <- file.path(coexp_dir, "CSSN", paste0(stem, "_mean_coexpression_network.csv"))
    if (file.exists(f_mean)) return(read_gene_matrix(f_mean, n))
    f <- file.path(coexp_dir, "CSSN", paste0(stem, "_CSSN_subset_results.rds"))
    if (!file.exists(f)) return(NULL)
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    corr_list <- if (is.list(obj)) obj[["Correlation Matrix"]] else NULL
    if (!is.list(corr_list)) return(NULL)
    agg <- matrix(0, nrow = n, ncol = n)
    n_ok <- 0L
    for (cur in corr_list) {
      dm <- tryCatch(dim(cur), error = function(e) NULL)
      if (is.null(dm) || length(dm) != 2L || any(dm < 1L)) next
      k <- min(n, dm[1L], dm[2L])
      tmp <- matrix(0, nrow = n, ncol = n)
      tmp[seq_len(k), seq_len(k)] <- abs(as.matrix(cur[seq_len(k), seq_len(k)]))
      agg <- agg + tmp
      n_ok <- n_ok + 1L
    }
    if (n_ok == 0L) return(NULL)
    agg <- agg / n_ok
    diag(agg) <- 0
    return(agg)
  }
  if (method == "hdWGCNA") {
    f <- file.path(coexp_dir, "hdWGCNA_Agent", paste0(stem, "_modules.csv"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    if (nrow(dt) == 0L) return(matrix(0, nrow = n, ncol = n))
    gene_col <- if ("gene_name" %in% names(dt)) "gene_name" else names(dt)[1L]
    mod_col <- if ("module" %in% names(dt)) "module" else names(dt)[2L]
    module <- as.character(dt[[mod_col]])
    module[module == "grey"] <- NA_character_  # hdWGCNA's unassigned genes
    return(membership_adjacency(gene_index_from_names(dt[[gene_col]]), module, n))
  }
  if (method == "Smoothie") {
    f <- file.path(coexp_dir, "smoothie_Agent", paste0(stem, "_edge_list.csv"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    out <- matrix(0, nrow = n, ncol = n)
    if (nrow(dt) == 0L) return(out)
    g1 <- gene_index_from_names(dt[[if ("gene1" %in% names(dt)) "gene1" else names(dt)[1L]]])
    g2 <- gene_index_from_names(dt[[if ("gene2" %in% names(dt)) "gene2" else names(dt)[2L]]])
    score_col <- intersect(c("PCC", "Rescaled_PCC"), names(dt))
    score <- if (length(score_col)) abs(as.numeric(dt[[score_col[1L]]])) else rep(1, nrow(dt))
    score[!is.finite(score)] <- 0
    keep <- !is.na(g1) & !is.na(g2) & g1 >= 1L & g1 <= n & g2 >= 1L & g2 <= n
    for (r in which(keep)) {
      v <- max(out[g1[r], g2[r]], score[r])
      out[g1[r], g2[r]] <- v
      out[g2[r], g1[r]] <- v
    }
    diag(out) <- 0
    return(out)
  }
  if (method == "stIHC") {
    f <- file.path(coexp_dir, "stIHC_Agent", paste0(scenario, "_gene_stIHC_", rep_id, ".RDS"))
    if (!file.exists(f)) return(NULL)
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    labels <- obj$label
    if (is.null(labels)) return(NULL)
    if (is.list(labels) && "label" %in% names(labels)) labels <- labels$label
    labels <- as.character(unlist(labels))
    L <- min(length(labels), n)
    if (L < 1L) return(matrix(0, nrow = n, ncol = n))
    return(membership_adjacency(seq_len(L), labels[seq_len(L)], n))
  }
  NULL
}

# Step-wise average precision evaluated at each distinct score threshold.
average_precision <- function(score, label) {
  n_pos <- sum(label == 1L)
  if (n_pos == 0L) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  score <- score[ord]
  label <- label[ord]
  block_end <- c(score[-1L] != score[-length(score)], TRUE)
  tp <- cumsum(label)[block_end]
  k <- seq_along(label)[block_end]
  precision <- tp / k
  recall <- tp / n_pos
  sum(diff(c(0, recall)) * precision)
}

# Pairwise edge metrics over the upper triangle. Continuous outputs are
# thresholded at |s| >= COEXP_EDGE_THRESHOLD; co-membership adjacencies are 0/1.
coexp_metrics <- function(pred, truth, continuous, threshold = COEXP_EDGE_THRESHOLD) {
  tri <- upper.tri(truth, diag = FALSE)
  score <- abs(pred[tri])
  label <- as.integer(truth[tri] > 0)
  call <- as.integer(score >= threshold)
  out <- binary_metrics(
    tp = sum(call == 1L & label == 1L),
    fp = sum(call == 1L & label == 0L),
    fn = sum(call == 0L & label == 1L),
    tn = sum(call == 0L & label == 0L)
  )
  out$n_positive_truth <- sum(label)
  out$AUPRC <- if (continuous && out$has_positive_truth) average_precision(score, label) else NA_real_
  out
}

## Deconvolution --------------------------------------------------------------------------

DECONV_UNNAMED_COMPONENTS <- c("STdeconvolve")

read_deconv_prediction <- function(base_dir, scenario, rep_id, method) {
  data_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"))
  f <- switch(
    method,
    CARD = file.path(data_dir, paste0("CARD_", rep_id, ".csv")),
    RCTD = file.path(data_dir, paste0("RCTD_", rep_id, ".csv")),
    STdeconvolve = file.path(data_dir, paste0("STdeconvolve_", rep_id, ".csv")),
    STdGCN = file.path(data_dir, paste0("STdGCN_input_", rep_id), "output", "predict_result.csv"),
    Tangram = file.path(data_dir, paste0("Tangram_", rep_id, ".csv")),
    POLARIS = file.path(data_dir, paste0("POLARIS_", rep_id, ".csv")),
    stop("Unknown deconvolution method: ", method)
  )
  if (!file.exists(f)) return(NULL)
  # All prediction files carry a header. STdeconvolve's topic names are
  # numeric ("1","2","3"), which fread would otherwise read as a data row.
  dt <- fread(f, header = TRUE)
  spot_ids <- NULL
  for (col in c("spot", "barcodes")) {
    if (col %in% names(dt)) {
      spot_ids <- as.character(dt[[col]])
      dt[, (col) := NULL]
    }
  }
  if ("V1" %in% names(dt)) {
    if (!is.numeric(dt$V1)) spot_ids <- as.character(dt$V1)
    if (!is.numeric(dt$V1) || method == "Tangram") dt[, V1 := NULL]
  }
  pred <- as.data.frame(dt)
  pred <- pred[, vapply(pred, is.numeric, logical(1)), drop = FALSE]
  if (ncol(pred) == 0L) return(NULL)
  attr_has_ids <- !is.null(spot_ids)
  rownames(pred) <- if (attr_has_ids) spot_ids else paste0("spot_", seq_len(nrow(pred)))
  attr(pred, "has_spot_ids") <- attr_has_ids
  pred
}

ordered_subsets <- function(pool, k) {
  if (k == 0L) return(matrix(integer(0), nrow = 1L, ncol = 0L))
  do.call(rbind, lapply(seq_along(pool), function(i) {
    cbind(pool[i], ordered_subsets(pool[-i], k - 1L))
  }))
}

# One-to-one assignment of unnamed components to truth cell types that
# minimizes the summed squared error. Truth types without a component receive
# zero; unmatched components are excluded.
align_unnamed_components <- function(truth_m, pred_m) {
  K <- ncol(pred_m)
  C <- ncol(truth_m)
  if (max(K, C) > 8L) stop("Too many components for exhaustive assignment.")
  cost <- matrix(0, nrow = K, ncol = C)
  for (k in seq_len(K)) for (c in seq_len(C)) cost[k, c] <- sum((pred_m[, k] - truth_m[, c])^2)
  zero_cost <- colSums(truth_m^2)

  aligned <- matrix(0, nrow = nrow(truth_m), ncol = C, dimnames = dimnames(truth_m))
  if (K >= C) {
    choices <- ordered_subsets(seq_len(K), C)  # topic assigned to each truth type
    total <- apply(choices, 1L, function(r) sum(cost[cbind(r, seq_len(C))]))
    best <- choices[which.min(total), ]
    aligned[, seq_len(C)] <- pred_m[, best, drop = FALSE]
  } else {
    choices <- ordered_subsets(seq_len(C), K)  # truth type assigned to each topic
    total <- apply(choices, 1L, function(r) {
      sum(cost[cbind(seq_len(K), r)]) + sum(zero_cost[setdiff(seq_len(C), r)])
    })
    best <- choices[which.min(total), ]
    aligned[, best] <- pred_m
  }
  aligned
}

# Aligns a prediction to the truth spots and cell types. Spots are matched by
# identifier, or by position only when the prediction has no identifiers and
# the same number of spots.
align_deconv_prediction <- function(truth, pred, method) {
  fail <- function(status) list(truth = NULL, pred = NULL, n_spots = 0L, status = status)
  if (is.null(pred)) return(fail("missing"))
  if (isTRUE(attr(pred, "has_spot_ids"))) {
    spots <- intersect(rownames(truth), rownames(pred))
  } else if (nrow(pred) == nrow(truth)) {
    spots <- rownames(truth)
    rownames(pred) <- spots
  } else {
    return(fail("spot_count_mismatch_without_ids"))
  }
  if (length(spots) == 0L) return(fail("no_common_spots"))
  truth_m <- as.matrix(truth[spots, , drop = FALSE])
  pred_m <- as.matrix(pred[spots, , drop = FALSE])

  if (method %in% DECONV_UNNAMED_COMPONENTS) {
    aligned <- align_unnamed_components(truth_m, pred_m)
  } else {
    matched <- intersect(colnames(truth_m), colnames(pred_m))
    if (length(matched) == 0L) return(fail("no_matching_cell_types"))
    aligned <- matrix(0, nrow = nrow(truth_m), ncol = ncol(truth_m), dimnames = dimnames(truth_m))
    aligned[, matched] <- pred_m[, matched, drop = FALSE]
  }
  list(
    truth = truth_m,
    pred = aligned,
    n_spots = length(spots),
    status = if (length(spots) == nrow(truth)) "ok" else "partial_spots"
  )
}

# RMSE over all spots and all truth cell types (Supplementary Methods S2).
deconv_rmse <- function(truth, pred, method) {
  a <- align_deconv_prediction(truth, pred, method)
  rmse <- if (is.null(a$truth)) NA_real_ else sqrt(mean((a$truth - a$pred)^2))
  list(RMSE = rmse, n_spots = a$n_spots, status = a$status)
}

## Script helpers -------------------------------------------------------------------------

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# Evaluators source this file with:
#   local({
#     a <- commandArgs(trailingOnly = FALSE)
#     f <- sub("^--file=", "", a[grepl("^--file=", a)])
#     here <- if (length(f) == 1L) dirname(normalizePath(f)) else getwd()
#     rules <- c(Sys.getenv("STCOS_EVALUATION_RULES"),
#                file.path(here, "..", "evaluation_rules.R"),
#                file.path(getwd(), "evaluation_rules.R"))
#     rules <- rules[nzchar(rules) & file.exists(rules)]
#     if (!length(rules)) stop("Cannot find evaluation_rules.R; set STCOS_EVALUATION_RULES.")
#     source(rules[1L], local = globalenv())
#   })
