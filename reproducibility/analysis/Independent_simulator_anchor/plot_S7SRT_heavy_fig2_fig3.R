#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

base_dir <- "."
out_dir <- file.path(base_dir, "Independent_simulator_anchor", "S7SRT_heavy_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- c("S7", "S7SRT")
scenario_labels <- c("S7" = "ST-COS S7", "S7SRT" = "SRTsim-S7 anchor")
reps <- 1:5

label_scenario <- function(x) {
  factor(scenario_labels[as.character(x)], levels = scenario_labels[scenarios])
}

extract_index <- function(x) {
  suppressWarnings(as.integer(sub(".*?(\\d+)$", "\\1", as.character(x))))
}

calc_ari <- function(labels_true, labels_pred) {
  keep <- !(is.na(labels_true) | is.na(labels_pred))
  labels_true <- as.character(labels_true[keep])
  labels_pred <- as.character(labels_pred[keep])
  n <- length(labels_true)
  if (n < 2) return(NA_real_)
  tab <- table(labels_true, labels_pred)
  comb2 <- function(x) x * (x - 1) / 2
  sum_nij <- sum(comb2(tab))
  sum_ai <- sum(comb2(rowSums(tab)))
  sum_bj <- sum(comb2(colSums(tab)))
  sum_n <- comb2(n)
  expected <- (sum_ai * sum_bj) / sum_n
  max_index <- 0.5 * (sum_ai + sum_bj)
  denom <- max_index - expected
  if (denom == 0) return(NA_real_)
  as.numeric((sum_nij - expected) / denom)
}

calc_f1 <- function(tp, fp, fn) {
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
}

align_matrix <- function(mat, n) {
  if (is.null(mat)) return(NULL)
  m <- suppressWarnings(as.matrix(mat))
  if (!is.numeric(m)) storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  if (is.null(nrow(m)) || is.null(ncol(m)) || nrow(m) < 1 || ncol(m) < 1) return(NULL)
  out <- matrix(0, nrow = n, ncol = n)
  k <- min(n, nrow(m), ncol(m))
  out[seq_len(k), seq_len(k)] <- m[seq_len(k), seq_len(k)]
  diag(out) <- 0
  out
}

common_theme <- theme_bw(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.35),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "#d2e0e3", color = "#6f6f6f", linewidth = 0.35),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.margin = margin(6, 8, 6, 8)
  )

save_plot <- function(plot, stem, width, height) {
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  message("Saving ", stem, " PDF...")
  ggsave(pdf_file, plot, width = width, height = height, bg = "white")
  message("Saving ", stem, " PNG...")
  ggsave(png_file, plot, width = width, height = height, dpi = 300, bg = "white")
  c(pdf = pdf_file, png = png_file)
}

mean_ci <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  mean_x <- if (n > 0) mean(x) else NA_real_
  sd_x <- if (n > 1) sd(x) else NA_real_
  se_x <- if (n > 1) sd_x / sqrt(n) else NA_real_
  q <- if (n > 1) qt((1 + conf) / 2, df = n - 1) else NA_real_
  list(
    mean = mean_x,
    median = if (n > 0) median(x) else NA_real_,
    sd = sd_x,
    n = n,
    ci_low = if (n > 1) mean_x - q * se_x else NA_real_,
    ci_high = if (n > 1) mean_x + q * se_x else NA_real_
  )
}

## Fig. 2-style metric collection -------------------------------------

seg_cfg <- data.table(
  method = c("BASS", "BayesSpace", "GraphST", "MAPLE", "SpaGCN", "STAGATE"),
  prefix = c("BASS", "BayesSpace", "GraphST", "MAPLE", "SpaGCN", "STAGATE"),
  truth_col = c("layer", "layer", "Ground Truth", "layer", "Ground Truth", "Ground Truth"),
  pred_col = c("BASS_layer", "BayesSpace_labels", "domain", "cluster", "refined_pred", "mclust")
)

collect_segmentation <- function() {
  rows <- list()
  i <- 1L
  for (scenario in scenarios) {
    for (rep_id in reps) {
      for (k in seq_len(nrow(seg_cfg))) {
        f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), paste0(seg_cfg$prefix[k], "_", rep_id, ".csv"))
        if (!file.exists(f)) next
        dt <- tryCatch(fread(f), error = function(e) NULL)
        if (is.null(dt)) next
        if (!(seg_cfg$truth_col[k] %in% names(dt)) || !(seg_cfg$pred_col[k] %in% names(dt))) next
        rows[[i]] <- data.table(
          scenario = scenario,
          benchmark = label_scenario(scenario),
          rep = rep_id,
          method = seg_cfg$method[k],
          ARI = calc_ari(dt[[seg_cfg$truth_col[k]]], dt[[seg_cfg$pred_col[k]]])
        )
        i <- i + 1L
      }
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

build_truth_props <- function(scenario, rep_id) {
  region_f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_region_", rep_id, ".csv"))
  region <- fread(region_f)
  ct_cols <- intersect(c("Cancer", "Tcell", "Endothelial", "CAF"), names(region))
  truth <- as.data.frame(region[, ..ct_cols])
  rownames(truth) <- paste0("spot_", seq_len(nrow(truth)))
  truth
}

read_deconv <- function(scenario, rep_id, method) {
  data_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"))
  f <- switch(
    method,
    "CARD" = file.path(data_dir, paste0("CARD_", rep_id, ".csv")),
    "RCTD" = file.path(data_dir, paste0("RCTD_", rep_id, ".csv")),
    "STdeconvolve" = file.path(data_dir, paste0("STdeconvolve_", rep_id, ".csv")),
    "STdGCN" = file.path(data_dir, paste0("STdGCN_input_", rep_id), "output", "predict_result.csv"),
    "Tangram" = file.path(data_dir, paste0("Tangram_", rep_id, ".csv")),
    "POLARIS" = file.path(data_dir, paste0("POLARIS_", rep_id, ".csv")),
    NULL
  )
  if (is.null(f) || !file.exists(f)) return(NULL)
  dt <- fread(f)
  spot_ids <- NULL
  if ("spot" %in% names(dt)) {
    spot_ids <- as.character(dt$spot)
    dt[, spot := NULL]
  }
  if ("barcodes" %in% names(dt)) {
    spot_ids <- as.character(dt$barcodes)
    dt[, barcodes := NULL]
  }
  if ("V1" %in% names(dt) && !is.numeric(dt$V1)) {
    spot_ids <- as.character(dt$V1)
    dt[, V1 := NULL]
  }
  if ("V1" %in% names(dt) && method == "Tangram") dt[, V1 := NULL]

  pred <- as.data.frame(dt)
  pred <- pred[, vapply(pred, is.numeric, logical(1)), drop = FALSE]
  if (ncol(pred) == 0) return(NULL)
  if (method == "STdeconvolve" && ncol(pred) == 3) {
    colnames(pred) <- c("Cancer", "Tcell", "Endothelial")
  }
  if (all(c("ct1", "ct2", "ct3") %in% colnames(pred))) {
    colnames(pred)[match(c("ct1", "ct2", "ct3"), colnames(pred))] <- c("Cancer", "Tcell", "Endothelial")
  }
  keep_cols <- intersect(c("Cancer", "Tcell", "Endothelial", "CAF"), colnames(pred))
  pred <- pred[, keep_cols, drop = FALSE]
  if (ncol(pred) == 0) return(NULL)
  rownames(pred) <- if (is.null(spot_ids)) paste0("spot_", seq_len(nrow(pred))) else spot_ids
  pred
}

calc_rmse <- function(truth, pred) {
  if (is.null(pred)) return(NA_real_)
  common_cols <- intersect(colnames(truth), colnames(pred))
  common_spots <- intersect(rownames(truth), rownames(pred))
  if (length(common_cols) == 0 || length(common_spots) == 0) return(NA_real_)
  sqrt(mean((as.matrix(truth[common_spots, common_cols, drop = FALSE]) -
               as.matrix(pred[common_spots, common_cols, drop = FALSE]))^2, na.rm = TRUE))
}

deconv_methods <- c("CARD", "RCTD", "STdeconvolve", "STdGCN", "Tangram", "POLARIS")

collect_deconvolution <- function() {
  rows <- list()
  i <- 1L
  for (scenario in scenarios) {
    for (rep_id in reps) {
      truth <- build_truth_props(scenario, rep_id)
      for (method in deconv_methods) {
        pred <- read_deconv(scenario, rep_id, method)
        rows[[i]] <- data.table(
          scenario = scenario,
          benchmark = label_scenario(scenario),
          rep = rep_id,
          method = method,
          RMSE = calc_rmse(truth, pred),
          available = !is.null(pred)
        )
        i <- i + 1L
      }
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

## Fig. 3-style metric collection -------------------------------------

svg_cfg <- data.table(
  method = c("moranI", "HEARTSVG", "SPARK", "SPARKX", "SpatialDE", "SpatialDE2"),
  prefix = c("moranI", "HEARTSVG", "SPARK", "SPARKX", "spatialDE", "spatialDE2")
)

get_truth_genes <- function(scenario, rep_id) {
  truth_f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_marker_gene_truth.csv"))
  dt <- fread(truth_f)
  if ("scenario" %in% names(dt)) dt <- dt[scenario == rep_id]
  unique(as.character(dt$gene_id))
}

get_gene_col <- function(dt) {
  hit <- intersect(c("gene", "g", "V1", "Unnamed: 0"), names(dt))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

get_svg_pred <- function(dt, method) {
  if (method == "moranI" && "pval_norm_fdr_bh" %in% names(dt)) return(as.integer(dt$pval_norm_fdr_bh <= 0.05))
  if (method == "HEARTSVG" && "p_adj" %in% names(dt)) return(as.integer(dt$p_adj <= 0.05))
  if (method == "SPARK" && "adjusted_pvalue" %in% names(dt)) return(as.integer(dt$adjusted_pvalue <= 0.05))
  if (method == "SPARKX" && "adjustedPval" %in% names(dt)) return(as.integer(dt$adjustedPval <= 0.05))
  if (method == "SpatialDE") {
    if ("qval" %in% names(dt)) return(as.integer(dt$qval <= 0.05))
    if ("pval" %in% names(dt)) return(as.integer(dt$pval <= 0.05))
  }
  if (method == "SpatialDE2" && "FSV" %in% names(dt)) {
    thr <- as.numeric(quantile(dt$FSV, 0.9, na.rm = TRUE))
    return(as.integer(dt$FSV >= thr))
  }
  NULL
}

calc_binary_gene_metrics <- function(pred_genes, truth_genes, universe) {
  pred_genes <- unique(as.character(pred_genes))
  truth_genes <- unique(as.character(truth_genes))
  universe <- unique(c(as.character(universe), pred_genes, truth_genes))
  tp <- length(intersect(pred_genes, truth_genes))
  fp <- length(setdiff(pred_genes, truth_genes))
  fn <- length(setdiff(truth_genes, pred_genes))
  tn <- length(setdiff(universe, union(pred_genes, truth_genes)))
  data.table(
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn,
    TPR = if ((tp + fn) > 0) tp / (tp + fn) else 0,
    FPR = if ((fp + tn) > 0) fp / (fp + tn) else 0,
    F1_score = calc_f1(tp, fp, fn)
  )
}

collect_svg <- function() {
  rows <- list()
  i <- 1L
  for (scenario in scenarios) {
    for (rep_id in reps) {
      truth_genes <- get_truth_genes(scenario, rep_id)
      universe <- colnames(fread(file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_gene_expression_", rep_id, ".csv")), nrows = 0))
      for (k in seq_len(nrow(svg_cfg))) {
        method <- svg_cfg$method[k]
        f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), paste0(svg_cfg$prefix[k], "_", rep_id, ".csv"))
        if (!file.exists(f)) next
        dt <- fread(f)
        gene_col <- get_gene_col(dt)
        pred_vec <- get_svg_pred(dt, method)
        if (is.na(gene_col) || is.null(pred_vec) || length(pred_vec) != nrow(dt)) next
        genes <- as.character(dt[[gene_col]])
        keep <- !is.na(genes) & genes != ""
        met <- calc_binary_gene_metrics(genes[keep][pred_vec[keep] == 1], truth_genes, universe)
        rows[[i]] <- cbind(
          data.table(scenario = scenario, benchmark = label_scenario(scenario), rep = rep_id, method = method),
          met
        )
        i <- i + 1L
      }
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

ctsvg_cfg <- data.table(
  method = c("CELINA", "CELINA_trueprop", "CSIDE", "CSIDE_trueprop", "CTSV_trueprop", "STANCE_trueprop", "spVC_trueprop"),
  prefix = c("CELINA", "CELINA_trueprop", "CSIDE", "CSIDE_trueprop", "CTSV_truth", "STANCE_siggenes", "spVC_trueprop")
)

extract_ctsvg_pred <- function(obj, ct, method) {
  if (is.null(obj)) return(character(0))
  if (method == "STANCE_trueprop") {
    if (is.list(obj) && !is.null(obj[[ct]])) return(unique(as.character(obj[[ct]])))
    if (isS4(obj) && "Test_2" %in% slotNames(obj)) {
      test2 <- obj@Test_2
      if (!is.null(test2[[ct]]) && is.data.frame(test2[[ct]]) && "p_value_adj" %in% names(test2[[ct]])) {
        return(rownames(test2[[ct]])[test2[[ct]]$p_value_adj < 0.05])
      }
    }
    return(character(0))
  }
  if (is.list(obj) && !is.null(obj[[ct]])) return(unique(as.character(obj[[ct]])))
  character(0)
}

collect_ctsvg <- function() {
  cell_rows <- list()
  rep_rows <- list()
  ci <- 1L
  ri <- 1L
  for (scenario in scenarios) {
    for (rep_id in reps) {
      truth_f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_marker_gene_truth.csv"))
      truth_dt <- fread(truth_f)
      if ("scenario" %in% names(truth_dt)) truth_dt <- truth_dt[scenario == rep_id]
      truth_sets <- split(unique(truth_dt$gene_id), truth_dt$cell_type)
      universe <- colnames(fread(file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_gene_expression_", rep_id, ".csv")), nrows = 0))

      for (k in seq_len(nrow(ctsvg_cfg))) {
        method <- ctsvg_cfg$method[k]
        f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), paste0(ctsvg_cfg$prefix[k], "_", rep_id, ".RDS"))
        if (!file.exists(f)) next
        obj <- tryCatch(readRDS(f), error = function(e) NULL)
        if (is.null(obj)) next
        f1_vals <- tpr_vals <- fpr_vals <- numeric(0)
        for (ct in names(truth_sets)) {
          pred <- extract_ctsvg_pred(obj, ct, method)
          met <- calc_binary_gene_metrics(pred, truth_sets[[ct]], universe)
          cell_rows[[ci]] <- cbind(
            data.table(scenario = scenario, benchmark = label_scenario(scenario), rep = rep_id, method = method, cell_type = ct),
            met
          )
          ci <- ci + 1L
          f1_vals <- c(f1_vals, met$F1_score)
          tpr_vals <- c(tpr_vals, met$TPR)
          fpr_vals <- c(fpr_vals, met$FPR)
        }
        rep_rows[[ri]] <- data.table(
          scenario = scenario,
          benchmark = label_scenario(scenario),
          rep = rep_id,
          method = method,
          F1_score = mean(f1_vals, na.rm = TRUE),
          TPR = mean(tpr_vals, na.rm = TRUE),
          FPR = mean(fpr_vals, na.rm = TRUE)
        )
        ri <- ri + 1L
      }
    }
  }
  list(cell = rbindlist(cell_rows, use.names = TRUE, fill = TRUE),
       rep = rbindlist(rep_rows, use.names = TRUE, fill = TRUE))
}

create_truth_matrix <- function(scenario, rep_id, n) {
  truth <- matrix(0, nrow = n, ncol = n)
  f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_coexpression_truth.csv"))
  dt <- fread(f)
  if ("scenario" %in% names(dt)) dt <- dt[scenario == rep_id]
  dt[, idx := if ("gene_index" %in% names(dt)) as.integer(gene_index) else extract_index(gene_id)]
  mod_col <- if ("module_id" %in% names(dt)) "module_id" else "module_label"
  dt <- dt[!is.na(idx) & idx >= 1 & idx <= n]
  for (m in unique(dt[[mod_col]])) {
    ids <- unique(dt[get(mod_col) == m, idx])
    if (length(ids) > 1) truth[ids, ids] <- 1
  }
  diag(truth) <- 0
  truth
}

coexp_methods <- c("spMOCA", "SpaceX", "CSSN", "hdWGCNA", "Smoothie", "binSpect", "stIHC")
coexp_dir <- file.path(base_dir, "coexpression")

load_coexp_matrix <- function(scenario, rep_id, method, n) {
  data_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"))
  if (method == "spMOCA") {
    f <- file.path(coexp_dir, "spMOCA_Agent", paste0(scenario, "_", rep_id, "_correlation_matrix.csv"))
    if (file.exists(f)) return(align_matrix(as.matrix(fread(f)), n))
  }
  if (method == "SpaceX") {
    f <- file.path(data_dir, paste0("SpaceX_", rep_id, ".csv"))
    if (file.exists(f)) return(align_matrix(as.matrix(fread(f)), n))
  }
  if (method == "CSSN") {
    f_mean <- file.path(coexp_dir, "CSSN", paste0(scenario, "_", rep_id, "_mean_coexpression_network.csv"))
    if (file.exists(f_mean)) return(align_matrix(as.matrix(fread(f_mean)), n))
    f <- file.path(coexp_dir, "CSSN", paste0(scenario, "_", rep_id, "_CSSN_subset_results.rds"))
    if (!file.exists(f)) return(NULL)
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    corr_list <- obj[["Correlation Matrix"]]
    if (!is.list(corr_list)) return(align_matrix(obj, n))
    agg <- matrix(0, nrow = n, ncol = n)
    n_ok <- 0L
    for (cur in corr_list) {
      dm <- tryCatch(dim(cur), error = function(e) NULL)
      if (is.null(dm) || length(dm) != 2 || any(dm < 1)) next
      tmp <- matrix(0, nrow = n, ncol = n)
      k <- min(n, dm[1], dm[2])
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
    f <- file.path(coexp_dir, "hdWGCNA_Agent", paste0(scenario, "_", rep_id, "_modules.csv"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    out <- matrix(0, nrow = n, ncol = n)
    if (nrow(dt) == 0) return(out)
    gene_col <- if ("gene_name" %in% names(dt)) "gene_name" else names(dt)[1]
    mod_col <- if ("module" %in% names(dt)) "module" else names(dt)[2]
    dt[, idx := extract_index(get(gene_col))]
    dt <- dt[!is.na(idx) & idx >= 1 & idx <= n & !is.na(get(mod_col)) & get(mod_col) != "" & get(mod_col) != "grey"]
    for (m in unique(dt[[mod_col]])) {
      ids <- unique(dt[get(mod_col) == m, idx])
      if (length(ids) > 1) out[ids, ids] <- 1
    }
    diag(out) <- 0
    return(out)
  }
  if (method == "Smoothie") {
    f <- file.path(coexp_dir, "smoothie_Agent", paste0(scenario, "_", rep_id, "_edge_list.csv"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    out <- matrix(0, nrow = n, ncol = n)
    if (nrow(dt) == 0) return(out)
    g1_col <- if ("gene1" %in% names(dt)) "gene1" else names(dt)[1]
    g2_col <- if ("gene2" %in% names(dt)) "gene2" else names(dt)[2]
    score_col <- if ("PCC" %in% names(dt)) "PCC" else if ("Rescaled_PCC" %in% names(dt)) "Rescaled_PCC" else NULL
    dt[, i := extract_index(get(g1_col))]
    dt[, j := extract_index(get(g2_col))]
    dt[, score := if (!is.null(score_col)) abs(as.numeric(get(score_col))) else 1]
    dt <- dt[!is.na(i) & !is.na(j) & i >= 1 & i <= n & j >= 1 & j <= n]
    for (r in seq_len(nrow(dt))) {
      out[dt$i[r], dt$j[r]] <- max(out[dt$i[r], dt$j[r]], dt$score[r], na.rm = TRUE)
      out[dt$j[r], dt$i[r]] <- out[dt$i[r], dt$j[r]]
    }
    diag(out) <- 0
    return(out)
  }
  if (method == "binSpect") {
    f1 <- file.path(coexp_dir, "binSpect_agent", paste0(scenario, "_", rep_id, "_correlation_matrix.csv"))
    f2 <- file.path(coexp_dir, "binSpect_Agent", paste0(scenario, "_", rep_id, "_correlation_matrix.csv"))
    f <- if (file.exists(f1)) f1 else if (file.exists(f2)) f2 else NA_character_
    if (!is.na(f)) return(align_matrix(as.matrix(fread(f)), n))
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
    out <- matrix(0, nrow = n, ncol = n)
    if (L < 1) return(out)
    labels <- labels[seq_len(L)]
    for (lab in unique(labels)) {
      ids <- which(labels == lab)
      if (length(ids) > 1) out[ids, ids] <- 1
    }
    diag(out) <- 0
    return(out)
  }
  NULL
}

compute_coexp_metrics <- function(pred, truth, threshold = 0.25) {
  pred_bin <- ifelse(abs(pred) >= threshold, 1L, 0L)
  true_bin <- ifelse(truth > 0, 1L, 0L)
  tri <- upper.tri(true_bin)
  p <- pred_bin[tri]
  t <- true_bin[tri]
  tp <- sum(p == 1L & t == 1L)
  fp <- sum(p == 1L & t == 0L)
  tn <- sum(p == 0L & t == 0L)
  fn <- sum(p == 0L & t == 1L)
  data.table(
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn,
    TPR = if ((tp + fn) > 0) tp / (tp + fn) else 0,
    FPR = if ((fp + tn) > 0) fp / (fp + tn) else 0,
    F1_score = calc_f1(tp, fp, fn)
  )
}

collect_coexpression <- function() {
  rows <- list()
  i <- 1L
  for (scenario in scenarios) {
    for (rep_id in reps) {
      expr_f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_gene_expression_", rep_id, ".csv"))
      n <- ncol(fread(expr_f, nrows = 0))
      truth <- create_truth_matrix(scenario, rep_id, n)
      for (method in coexp_methods) {
        pred <- tryCatch(load_coexp_matrix(scenario, rep_id, method, n), error = function(e) NULL)
        if (is.null(pred)) next
        met <- compute_coexp_metrics(pred, truth)
        rows[[i]] <- cbind(
          data.table(scenario = scenario, benchmark = label_scenario(scenario), rep = rep_id, method = method),
          met
        )
        i <- i + 1L
      }
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

## Representative maps -------------------------------------------------

plot_cat_map <- function(dt, value_col, title, colors) {
  ggplot(dt, aes(x = x_coord, y = y_coord, color = .data[[value_col]])) +
    geom_point(size = 0.75) +
    coord_equal() +
    scale_color_manual(values = colors, drop = FALSE) +
    labs(title = title, color = NULL) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(size = 10.5, face = "bold", hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(3, 3, 3, 3)
    )
}

plot_cont_map <- function(dt, value_col, title) {
  ggplot(dt, aes(x = x_coord, y = y_coord, color = .data[[value_col]])) +
    geom_point(size = 0.75) +
    coord_equal() +
    scale_color_viridis_c(option = "magma") +
    labs(title = title, color = NULL) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(size = 10.5, face = "bold", hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(3, 3, 3, 3)
    )
}

make_deconv_maps <- function(rep_id = 1) {
  scenario <- "S7SRT"
  region <- fread(file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_region_", rep_id, ".csv")))
  truth_props <- as.data.frame(region[, .(Cancer, Tcell, Endothelial, CAF)])
  rownames(truth_props) <- paste0("spot_", seq_len(nrow(truth_props)))
  truth_dom <- colnames(truth_props)[max.col(as.matrix(truth_props), ties.method = "first")]
  base_dt <- data.table(x_coord = region$x_coord, y_coord = region$y_coord, value = truth_dom)
  ct_cols <- c("Cancer" = "#d7301f", "Tcell" = "#3182bd", "Endothelial" = "#2ca25f", "CAF" = "#756bb1")
  plots <- list(plot_cat_map(base_dt, "value", "Truth", ct_cols))
  for (method in deconv_methods) {
    pred <- read_deconv(scenario, rep_id, method)
    if (is.null(pred)) next
    common_spots <- intersect(rownames(truth_props), rownames(pred))
    idx <- extract_index(common_spots)
    pred_dom <- colnames(pred)[max.col(as.matrix(pred[common_spots, , drop = FALSE]), ties.method = "first")]
    dt <- data.table(x_coord = region$x_coord[idx], y_coord = region$y_coord[idx], value = pred_dom)
    plots[[length(plots) + 1L]] <- plot_cat_map(dt, "value", method, ct_cols)
  }
  plots
}

make_seg_maps <- function(rep_id = 1) {
  scenario <- "S7SRT"
  colors <- hue_pal()(8)
  names(colors) <- as.character(seq_len(8) - 1L)
  plots <- list()
  for (k in seq_len(nrow(seg_cfg))) {
    f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), paste0(seg_cfg$prefix[k], "_", rep_id, ".csv"))
    dt <- fread(f)
    x_col <- if ("x_coord" %in% names(dt)) "x_coord" else "x"
    y_col <- if ("y_coord" %in% names(dt)) "y_coord" else "y"
    if (k == 1L) {
      truth_dt <- data.table(x_coord = dt[[x_col]], y_coord = dt[[y_col]], value = as.character(dt[[seg_cfg$truth_col[k]]]))
      plots[[1L]] <- plot_cat_map(truth_dt, "value", "Truth", colors)
    }
    pred_dt <- data.table(x_coord = dt[[x_col]], y_coord = dt[[y_col]], value = as.character(dt[[seg_cfg$pred_col[k]]]))
    plots[[length(plots) + 1L]] <- plot_cat_map(pred_dt, "value", seg_cfg$method[k], colors)
  }
  plots
}

plot_adj_matrix <- function(mat, title, threshold = 0.25) {
  bin <- ifelse(abs(mat) >= threshold, 1L, 0L)
  diag(bin) <- 0
  n <- nrow(bin)
  dt <- CJ(i = seq_len(n), j = seq_len(n))
  dt[, edge := factor(ifelse(bin[cbind(i, j)] > 0, "1", "0"), levels = c("0", "1"))]
  ggplot(dt, aes(i, j, fill = edge)) +
    geom_raster() +
    coord_fixed() +
    scale_fill_manual(values = c("0" = "white", "1" = "#222222"), drop = FALSE) +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(size = 10.5, face = "bold", hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(3, 3, 3, 3)
    )
}

make_svg_example_maps <- function(rep_id = 1) {
  scenario <- "S7SRT"
  expr <- fread(file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_gene_expression_", rep_id, ".csv")))
  coord <- fread(file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_coord_", rep_id, ".csv")))
  genes <- c("gene_1", "gene_6", "gene_11", "gene_100")
  titles <- c("Cancer SVG", "T cell SVG", "Endothelial SVG", "Null gene")
  Map(function(g, tt) {
    dt <- data.table(x_coord = coord$x_coord, y_coord = coord$y_coord, expr = expr[[g]])
    plot_cont_map(dt, "expr", paste0(tt, ": ", g))
  }, genes, titles)
}

make_coexp_adj_maps <- function(rep_id = 1, metrics_dt) {
  scenario <- "S7SRT"
  expr_f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_gene_expression_", rep_id, ".csv"))
  n <- ncol(fread(expr_f, nrows = 0))
  truth <- create_truth_matrix(scenario, rep_id, n)
  best <- metrics_dt[scenario == "S7SRT" & rep == rep_id][order(-F1_score, FPR)][1:3, method]
  plots <- list(plot_adj_matrix(truth, "Truth", threshold = 0.5))
  for (method in best) {
    mat <- load_coexp_matrix(scenario, rep_id, method, n)
    if (!is.null(mat)) plots[[length(plots) + 1L]] <- plot_adj_matrix(mat, method)
  }
  plots
}

## Plot builders -------------------------------------------------------

plot_metric_box <- function(dt, y_col, method_levels, method_colors, y_label, title, ylim = NULL) {
  dt <- copy(dt)
  dt[, method := factor(method, levels = method_levels)]
  p <- ggplot(dt, aes(x = method, y = .data[[y_col]], fill = method)) +
    geom_boxplot(width = 0.62, linewidth = 0.35, outlier.shape = NA, na.rm = TRUE) +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.3, alpha = 0.75, na.rm = TRUE) +
    facet_wrap(~benchmark, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(title = title, y = y_label, x = NULL) +
    common_theme
  if (!is.null(ylim)) {
    p <- p + scale_y_continuous(limits = ylim, expand = c(0.01, 0.01))
  }
  p
}

make_fig2 <- function(deconv_dt, seg_dt) {
  deconv_colors <- c("CARD" = "#e19d92", "RCTD" = "#eac47c", "STdeconvolve" = "#7e8cad",
                     "STdGCN" = "#a5c2cd", "Tangram" = "#dadae4", "POLARIS" = "#9ccf9b")
  seg_colors <- c("BASS" = "#e19d92", "BayesSpace" = "#e4c078", "GraphST" = "#edc400",
                  "MAPLE" = "#b39a89", "SpaGCN" = "#7db8bf", "STAGATE" = "#cfcfd8")

  p_rmse <- plot_metric_box(deconv_dt, "RMSE", deconv_methods, deconv_colors,
                            "RMSE (lower is better)", "A. Spot deconvolution")
  p_ari <- plot_metric_box(seg_dt, "ARI", seg_cfg$method, seg_colors,
                           "ARI (higher is better)", "B. Spatial clustering", ylim = c(0, 1.02))
  metrics_row <- plot_grid(p_rmse, p_ari, ncol = 2, rel_widths = c(1, 1))

  deconv_maps <- make_deconv_maps(1)
  seg_maps <- make_seg_maps(1)
  deconv_row <- plot_grid(plotlist = deconv_maps, nrow = 1)
  seg_row <- plot_grid(plotlist = seg_maps, nrow = 1)

  title <- ggdraw() + draw_label("Fig. 2-like heavy benchmark: original S7 vs SRTsim-S7 anchor",
                                 fontface = "bold", size = 16, x = 0, hjust = 0)
  deconv_label <- ggdraw() + draw_label("C. SRTsim-S7 representative deconvolution maps, replicate 1",
                                        fontface = "bold", size = 13, x = 0, hjust = 0)
  seg_label <- ggdraw() + draw_label("D. SRTsim-S7 representative clustering maps, replicate 1",
                                     fontface = "bold", size = 13, x = 0, hjust = 0)

  plot_grid(title, metrics_row, deconv_label, deconv_row, seg_label, seg_row,
            ncol = 1, rel_heights = c(0.06, 0.42, 0.04, 0.25, 0.04, 0.25))
}

make_fig3 <- function(svg_dt, ctsvg_dt, coexp_dt) {
  svg_colors <- c("moranI" = "#d98b5f", "HEARTSVG" = "#a7c5bd", "SPARK" = "#d9b44a",
                  "SPARKX" = "#8ca6db", "SpatialDE" = "#c7b8a3", "SpatialDE2" = "#b8d8b8")
  ctsvg_colors <- c("CELINA" = "#E69F00", "CELINA_trueprop" = "#F6C667", "CSIDE" = "#56B4E9",
                    "CSIDE_trueprop" = "#0072B2", "CTSV_trueprop" = "#CC79A7",
                    "STANCE_trueprop" = "#999999", "spVC_trueprop" = "#009E73")
  coexp_colors <- c("spMOCA" = "#4B9CD3", "SpaceX" = "#4c9390", "CSSN" = "#d65d48",
                    "hdWGCNA" = "#8b5a96", "Smoothie" = "#f4a261",
                    "binSpect" = "#2a9d8f", "stIHC" = "#264653")

  p_svg <- plot_metric_box(svg_dt, "F1_score", svg_cfg$method, svg_colors,
                           "F1 score", "A. SVG detection", ylim = c(0, 1.02))
  p_ctsvg <- plot_metric_box(ctsvg_dt, "F1_score", ctsvg_cfg$method, ctsvg_colors,
                             "F1 score", "B. CT-SVG detection", ylim = c(0, 1.02))
  p_coexp <- plot_metric_box(coexp_dt, "F1_score", coexp_methods, coexp_colors,
                             "F1 score", "C. Co-expression", ylim = c(0, 1.02))
  metric_row <- plot_grid(p_svg, p_ctsvg, p_coexp, nrow = 1, rel_widths = c(1, 1, 1))

  svg_examples <- make_svg_example_maps(1)
  svg_row <- plot_grid(plotlist = svg_examples, nrow = 1)
  coexp_maps <- make_coexp_adj_maps(1, coexp_dt)
  coexp_row <- plot_grid(plotlist = coexp_maps, nrow = 1)

  title <- ggdraw() + draw_label("Fig. 3-like heavy benchmark: original S7 vs SRTsim-S7 anchor",
                                 fontface = "bold", size = 16, x = 0, hjust = 0)
  svg_label <- ggdraw() + draw_label("D. SRTsim-S7 representative expression patterns, replicate 1",
                                     fontface = "bold", size = 13, x = 0, hjust = 0)
  coexp_label <- ggdraw() + draw_label("E. SRTsim-S7 co-expression adjacency, truth and top F1 methods, replicate 1",
                                       fontface = "bold", size = 13, x = 0, hjust = 0)

  plot_grid(title, metric_row, svg_label, svg_row, coexp_label, coexp_row,
            ncol = 1, rel_heights = c(0.06, 0.48, 0.04, 0.20, 0.04, 0.22))
}

message("Collecting metrics...")
seg_dt <- collect_segmentation()
deconv_dt <- collect_deconvolution()
svg_dt <- collect_svg()
ctsvg <- collect_ctsvg()
ctsvg_rep_dt <- ctsvg$rep
coexp_dt <- collect_coexpression()

fwrite(seg_dt, file.path(out_dir, "fig2_segmentation_ARI_by_rep.csv"))
fwrite(deconv_dt, file.path(out_dir, "fig2_deconvolution_RMSE_by_rep.csv"))
fwrite(svg_dt, file.path(out_dir, "fig3_svg_metrics_by_rep.csv"))
fwrite(ctsvg_rep_dt, file.path(out_dir, "fig3_ctsvg_metrics_by_rep.csv"))
fwrite(ctsvg$cell, file.path(out_dir, "fig3_ctsvg_celltype_metrics.csv"))
fwrite(coexp_dt, file.path(out_dir, "fig3_coexpression_metrics_by_rep.csv"))

summarize_metric <- function(dt, value_col, task, higher_better) {
  out <- dt[, as.data.table(mean_ci(get(value_col))), by = .(benchmark, method)]
  out[, `:=`(task = task, metric = value_col, higher_better = higher_better)]
  out
}

summary_parts <- list(
  summarize_metric(deconv_dt, "RMSE", "Deconvolution", FALSE),
  summarize_metric(seg_dt, "ARI", "Segmentation", TRUE),
  summarize_metric(svg_dt, "F1_score", "SVG", TRUE),
  summarize_metric(ctsvg_rep_dt, "F1_score", "CT-SVG", TRUE),
  summarize_metric(coexp_dt, "F1_score", "Coexpression", TRUE)
)
summary_dt <- rbindlist(summary_parts, use.names = TRUE, fill = TRUE)
setcolorder(summary_dt, c("benchmark", "task", "metric", "higher_better", "method",
                          "mean", "median", "sd", "n", "ci_low", "ci_high"))
fwrite(summary_dt, file.path(out_dir, "S7_vs_S7SRT_metric_summary.csv"))

rank_rows <- lapply(unique(summary_dt$task), function(task_i) {
  sub <- copy(summary_dt[task == task_i])
  direction <- unique(sub$higher_better)
  sub[, rank := frank(if (direction) -mean else mean, ties.method = "average"), by = benchmark]
  wide <- dcast(sub, task + metric + method + higher_better ~ benchmark,
                value.var = c("mean", "rank"))
  setnames(
    wide,
    old = c("mean_ST-COS S7", "mean_SRTsim-S7 anchor",
            "rank_ST-COS S7", "rank_SRTsim-S7 anchor"),
    new = c("mean_STCOS_S7", "mean_SRTsim_S7",
            "rank_STCOS_S7", "rank_SRTsim_S7"),
    skip_absent = TRUE
  )
  wide
})
rank_dt <- rbindlist(rank_rows, use.names = TRUE, fill = TRUE)
fwrite(rank_dt, file.path(out_dir, "S7_vs_S7SRT_method_ranks.csv"))

rank_agreement <- rank_dt[, {
  ok <- is.finite(rank_STCOS_S7) & is.finite(rank_SRTsim_S7)
  stcos_best <- method[ok][which.min(rank_STCOS_S7[ok])]
  srt_best <- method[ok][which.min(rank_SRTsim_S7[ok])]
  stcos_worst <- method[ok][which.max(rank_STCOS_S7[ok])]
  srt_worst <- method[ok][which.max(rank_SRTsim_S7[ok])]
  list(
    n_methods = sum(ok),
    spearman_rank_r = if (sum(ok) >= 3) cor(rank_STCOS_S7[ok], rank_SRTsim_S7[ok], method = "spearman") else NA_real_,
    kendall_tau = if (sum(ok) >= 3) cor(rank_STCOS_S7[ok], rank_SRTsim_S7[ok], method = "kendall") else NA_real_,
    best_STCOS_S7 = stcos_best,
    best_SRTsim_S7 = srt_best,
    same_best = identical(stcos_best, srt_best),
    worst_STCOS_S7 = stcos_worst,
    worst_SRTsim_S7 = srt_worst,
    same_worst = identical(stcos_worst, srt_worst)
  )
}, by = .(task, metric, higher_better)]
fwrite(rank_agreement, file.path(out_dir, "S7_vs_S7SRT_rank_concordance_by_task.csv"))

message("Building figures...")
fig2 <- make_fig2(deconv_dt, seg_dt)
fig3 <- make_fig3(svg_dt, ctsvg_rep_dt, coexp_dt)
paths2 <- save_plot(fig2, "S7SRT_fig2_like_deconv_clustering", width = 22, height = 15)
paths3 <- save_plot(fig3, "S7SRT_fig3_like_svg_ctsvg_coexpression", width = 24, height = 15)

message("Saved:")
message(paths2["pdf"])
message(paths2["png"])
message(paths3["pdf"])
message(paths3["png"])
