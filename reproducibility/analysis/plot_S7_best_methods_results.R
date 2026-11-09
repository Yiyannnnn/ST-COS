#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

base_dir <- "."
scenario <- "S7"

seg_cfg <- data.table(
  method = c("BASS", "BayesSpace", "MAPLE", "STAGATE", "GraphST", "SpaGCN"),
  prefix = c("BASS", "BayesSpace", "MAPLE", "STAGATE", "GraphST", "SpaGCN"),
  gt_col = c("layer", "layer", "layer", "Ground Truth", "Ground Truth", "Ground Truth"),
  pred_col = c("BASS_layer", "BayesSpace_labels", "cluster", "mclust", "domain", "refined_pred")
)

extract_index <- function(x) {
  x <- as.character(x)
  suppressWarnings(as.integer(sub(".*?(\\d+)$", "\\1", x)))
}

align_matrix <- function(mat, n) {
  if (is.null(mat)) return(NULL)
  m <- suppressWarnings(as.matrix(mat))
  if (!is.numeric(m)) storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  nr <- nrow(m)
  nc <- ncol(m)
  if (is.null(nr) || is.null(nc) || nr < 1 || nc < 1) return(NULL)
  out <- matrix(0, nrow = n, ncol = n)
  k <- min(n, nr, nc)
  out[1:k, 1:k] <- m[1:k, 1:k]
  diag(out) <- 0
  out
}

read_pred_deconv <- function(method, rep_id) {
  data_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"))
  file_path <- switch(
    method,
    "CARD" = file.path(data_dir, paste0("CARD_", rep_id, ".csv")),
    "RCTD" = file.path(data_dir, paste0("RCTD_", rep_id, ".csv")),
    "Redeconve" = file.path(data_dir, paste0("Redeconve_", rep_id, ".csv")),
    "STdeconvolve" = file.path(data_dir, paste0("STdeconvolve_", rep_id, ".csv")),
    "STdGCN" = file.path(data_dir, paste0("STdGCN_input_", rep_id), "output", "predict_result.csv"),
    "Tangram" = file.path(data_dir, paste0("Tangram_", rep_id, ".csv")),
    "POLARIS" = file.path(data_dir, paste0("POLARIS_", rep_id, ".csv")),
    NULL
  )
  if (is.null(file_path) || !file.exists(file_path)) return(NULL)

  dt <- fread(file_path)
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
  pred <- pred[, sapply(pred, is.numeric), drop = FALSE]
  if (ncol(pred) == 0 || nrow(pred) == 0) return(NULL)

  if (method == "STdeconvolve" && ncol(pred) == 3) {
    colnames(pred) <- c("ct1", "ct2", "ct3")
  }
  if (all(c("ct1", "ct2", "ct3") %in% colnames(pred))) {
    colnames(pred)[match(c("ct1", "ct2", "ct3"), colnames(pred))] <- c("Cancer", "Tcell", "Endothelial")
  }

  keep_cols <- intersect(c("Cancer", "Tcell", "Endothelial", "CAF"), colnames(pred))
  pred <- pred[, keep_cols, drop = FALSE]
  if (ncol(pred) == 0) return(NULL)

  if (is.null(spot_ids)) {
    rownames(pred) <- paste0("spot_", seq_len(nrow(pred)))
  } else {
    rownames(pred) <- spot_ids
  }
  pred
}

load_coexp_matrix <- function(method, rep_id, n) {
  coexp_dir <- file.path(base_dir, "coexpression")
  data_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"))

  if (method == "spMOCA") {
    f <- file.path(coexp_dir, "spMOCA_Agent", paste0(scenario, "_", rep_id, "_correlation_matrix.csv"))
    if (!file.exists(f)) return(NULL)
    return(align_matrix(as.matrix(fread(f)), n))
  }

  if (method == "SpaceX") {
    f <- file.path(data_dir, paste0("SpaceX_", rep_id, ".csv"))
    if (!file.exists(f)) return(NULL)
    return(align_matrix(as.matrix(fread(f)), n))
  }

  if (method == "CSSN") {
    f_mean <- file.path(coexp_dir, "CSSN", paste0(scenario, "_", rep_id, "_mean_coexpression_network.csv"))
    if (file.exists(f_mean)) {
      return(align_matrix(as.matrix(fread(f_mean)), n))
    }
    f <- file.path(coexp_dir, "CSSN", paste0(scenario, "_", rep_id, "_CSSN_subset_results.rds"))
    if (!file.exists(f)) return(NULL)
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj) || !is.list(obj) || is.null(obj[["Correlation Matrix"]])) return(NULL)
    agg <- matrix(0, nrow = n, ncol = n)
    n_ok <- 0L
    for (i in seq_along(obj[["Correlation Matrix"]])) {
      cur <- obj[["Correlation Matrix"]][[i]]
      dm <- tryCatch(dim(cur), error = function(e) NULL)
      if (is.null(dm) || length(dm) != 2 || any(dm < 1)) next
      k <- min(n, dm[1], dm[2])
      sub <- tryCatch(as.matrix(cur[1:k, 1:k]), error = function(e) NULL)
      if (is.null(sub)) next
      tmp <- matrix(0, nrow = n, ncol = n)
      tmp[1:k, 1:k] <- abs(sub)
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
    idx <- extract_index(dt[[gene_col]])
    module <- as.character(dt[[mod_col]])
    keep <- !is.na(idx) & idx >= 1 & idx <= n & !is.na(module) & module != "" & module != "grey"
    if (!any(keep)) return(out)
    dt2 <- data.table(idx = idx[keep], module = module[keep])
    for (m in unique(dt2$module)) {
      ids <- unique(dt2[module == m, idx])
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
    idx1 <- extract_index(dt[[g1_col]])
    idx2 <- extract_index(dt[[g2_col]])
    score <- if (!is.null(score_col)) abs(as.numeric(dt[[score_col]])) else rep(1, nrow(dt))
    keep <- !is.na(idx1) & !is.na(idx2) & idx1 >= 1 & idx1 <= n & idx2 >= 1 & idx2 <= n
    idx1 <- idx1[keep]
    idx2 <- idx2[keep]
    score <- score[keep]
    for (i in seq_along(idx1)) {
      a <- idx1[i]
      b <- idx2[i]
      v <- score[i]
      if (!is.finite(v)) v <- 0
      if (abs(v) > abs(out[a, b])) {
        out[a, b] <- v
        out[b, a] <- v
      }
    }
    diag(out) <- 0
    return(out)
  }

  if (method == "binSpect") {
    f1 <- file.path(coexp_dir, "binSpect_agent", paste0(scenario, "_", rep_id, "_correlation_matrix.csv"))
    f2 <- file.path(coexp_dir, "binSpect_Agent", paste0(scenario, "_", rep_id, "_correlation_matrix.csv"))
    f <- if (file.exists(f1)) f1 else if (file.exists(f2)) f2 else NULL
    if (is.null(f)) return(NULL)
    return(align_matrix(as.matrix(fread(f)), n))
  }

  if (method == "stIHC") {
    f <- file.path(coexp_dir, "stIHC_Agent", paste0(scenario, "_gene_stIHC_", rep_id, ".RDS"))
    if (!file.exists(f)) return(NULL)
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj)) return(NULL)
    labels <- obj$label
    if (is.null(labels)) return(NULL)
    if (is.list(labels) && "label" %in% names(labels)) labels <- labels$label
    labels <- as.character(unlist(labels))
    L <- min(length(labels), n)
    if (L < 1) return(NULL)
    out <- matrix(0, nrow = n, ncol = n)
    labels <- labels[1:L]
    for (lab in unique(labels)) {
      ids <- which(labels == lab)
      if (length(ids) > 1) out[ids, ids] <- 1
    }
    diag(out) <- 0
    return(out)
  }

  NULL
}

plot_spatial <- function(dt, value_col, title_txt, color_map) {
  ggplot(dt, aes(x = x_coord, y = y_coord, color = .data[[value_col]])) +
    geom_point(size = 0.85) +
    coord_equal() +
    scale_color_manual(values = color_map, drop = FALSE) +
    labs(title = title_txt, color = NULL) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 8)
    )
}

plot_binary_matrix <- function(mat, title_txt) {
  if (is.null(rownames(mat))) rownames(mat) <- as.character(seq_len(nrow(mat)))
  if (is.null(colnames(mat))) colnames(mat) <- as.character(seq_len(ncol(mat)))
  dt <- as.data.table(as.table(mat))
  setnames(dt, c("i", "j", "edge"))
  dt[, i := as.integer(i)]
  dt[, j := as.integer(j)]
  edge_chr <- as.character(dt$edge)
  edge_num <- suppressWarnings(as.numeric(edge_chr))
  edge_num[is.na(edge_num)] <- as.integer(edge_chr[is.na(edge_num)] %in% c("TRUE", "T"))
  dt[, edge := factor(ifelse(edge_num > 0, "1", "0"), levels = c("0", "1"))]

  ggplot(dt, aes(x = i, y = j, fill = edge)) +
    geom_raster() +
    scale_fill_manual(values = c("0" = "white", "1" = "#222222"), labels = c("0" = "No edge", "1" = "Edge"), drop = FALSE) +
    coord_fixed() +
    labs(title = title_txt, fill = NULL) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 8)
    )
}

# 1) Pick best method per task and best replicate for each method.
seg_mean <- fread(file.path(base_dir, "Segmentation", "ARI_summary_S1_S7.csv"))[
  scenario == "S7"
]
seg_best <- seg_mean[order(-mean_ARI, method)][1]
seg_rep <- fread(file.path(base_dir, "Segmentation", "ARI_by_rep_S1_S7.csv"))[
  method == seg_best$method
][order(-ARI, rep)][1]

deconv_mean <- fread(file.path(base_dir, "Deconvolution", "RMSE_summary_S1_S7.csv"))[
  scenario == "S7" & n_rmse > 0 & is.finite(mean_RMSE)
]
deconv_best <- deconv_mean[order(mean_RMSE, method)][1]
deconv_rep <- fread(file.path(base_dir, "Deconvolution", "RMSE_by_rep_S1_S7.csv"))[
  scenario == "S7" & method == deconv_best$method & available == TRUE & is.finite(RMSE)
][order(RMSE, rep)][1]

coexp_mean <- fread(file.path(base_dir, "coexpression", "coexpression_method_summary.csv"))[
  scenario == "S7" & n > 0
]
coexp_best <- coexp_mean[order(-mean_F1, mean_FPR, method)][1]
coexp_rep <- fread(file.path(base_dir, "coexpression", "coexpression_method_metrics.csv"))[
  scenario == "S7" & method == coexp_best$method
][order(-F1, replicate)][1]

selection_dt <- rbindlist(list(
  data.table(task = "Segmentation", method = seg_best$method, replicate = seg_rep$rep, metric = "ARI", metric_value = seg_rep$ARI),
  data.table(task = "Deconvolution", method = deconv_best$method, replicate = deconv_rep$rep, metric = "RMSE", metric_value = deconv_rep$RMSE),
  data.table(task = "Coexpression", method = coexp_best$method, replicate = coexp_rep$replicate, metric = "F1", metric_value = coexp_rep$F1)
))

# 2) Segmentation: truth vs inferred map.
seg_method <- as.character(seg_best$method)
seg_rep_id <- as.integer(seg_rep$rep)
seg_row <- seg_cfg[method == seg_method]
if (nrow(seg_row) != 1) stop("Unsupported segmentation method: ", seg_method)

seg_file <- file.path(base_dir, "Data", paste0(scenario, "_gene"), paste0(seg_row$prefix, "_", seg_rep_id, ".csv"))
seg_dt <- fread(seg_file)
if (!(seg_row$gt_col %in% names(seg_dt)) || !(seg_row$pred_col %in% names(seg_dt))) {
  stop("Required segmentation columns missing in ", seg_file)
}
seg_plot_dt <- seg_dt[, .(
  x_coord,
  y_coord,
  truth = as.character(get(seg_row$gt_col)),
  pred = as.character(get(seg_row$pred_col))
)]
seg_levels <- sort(unique(c(seg_plot_dt$truth, seg_plot_dt$pred)))
seg_cols <- setNames(scales::hue_pal()(length(seg_levels)), seg_levels)

p_seg_truth <- plot_spatial(seg_plot_dt, "truth", sprintf("Segmentation GT (S7 rep %d)", seg_rep_id), seg_cols)
p_seg_pred <- plot_spatial(seg_plot_dt, "pred", sprintf("Segmentation %s (ARI=%.3f)", seg_method, seg_rep$ARI), seg_cols)

# 3) Deconvolution: dominant cell type map truth vs inferred.
dec_method <- as.character(deconv_best$method)
dec_rep_id <- as.integer(deconv_rep$rep)

truth_region <- fread(file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0("S7_region_", dec_rep_id, ".csv")))
truth_props <- as.data.frame(truth_region[, .(Cancer, Tcell, Endothelial, CAF)])
rownames(truth_props) <- paste0("spot_", seq_len(nrow(truth_props)))

pred_props <- read_pred_deconv(dec_method, dec_rep_id)
if (is.null(pred_props)) stop("Missing deconvolution predictions for ", dec_method, " rep ", dec_rep_id)

common_spots <- intersect(rownames(truth_props), rownames(pred_props))
if (length(common_spots) == 0) {
  n <- min(nrow(truth_props), nrow(pred_props))
  common_spots <- paste0("spot_", seq_len(n))
  rownames(pred_props) <- paste0("spot_", seq_len(nrow(pred_props)))
}

truth_use <- truth_props[common_spots, , drop = FALSE]
pred_use <- pred_props[common_spots, , drop = FALSE]

dominant_truth <- colnames(truth_use)[max.col(as.matrix(truth_use), ties.method = "first")]
dominant_pred <- colnames(pred_use)[max.col(as.matrix(pred_use), ties.method = "first")]

idx <- suppressWarnings(as.integer(sub("^spot_", "", common_spots)))
coords <- truth_region[idx, .(x_coord, y_coord)]
dec_plot_dt <- data.table(
  x_coord = coords$x_coord,
  y_coord = coords$y_coord,
  truth = dominant_truth,
  pred = dominant_pred
)

ct_levels <- c("Cancer", "Tcell", "Endothelial", "CAF")
ct_cols <- c("Cancer" = "#d7301f", "Tcell" = "#3182bd", "Endothelial" = "#2ca25f", "CAF" = "#756bb1")

p_dec_truth <- plot_spatial(dec_plot_dt, "truth", sprintf("Deconv GT dominant CT (S7 rep %d)", dec_rep_id), ct_cols[ct_levels])
p_dec_pred <- plot_spatial(dec_plot_dt, "pred", sprintf("Deconv %s (RMSE=%.3f)", dec_method, deconv_rep$RMSE), ct_cols[ct_levels])

# 4) Co-expression: truth vs inferred adjacency heatmap.
co_method <- as.character(coexp_best$method)
co_rep_id <- as.integer(coexp_rep$replicate)

truth_coexp <- fread(file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", "S7_coexpression_truth.csv"))
if ("scenario" %in% names(truth_coexp)) {
  truth_coexp <- truth_coexp[scenario == co_rep_id]
}
if ("gene_index" %in% names(truth_coexp)) {
  truth_coexp[, idx := as.integer(gene_index)]
} else {
  truth_coexp[, idx := extract_index(gene_id)]
}
module_col <- if ("module_id" %in% names(truth_coexp)) "module_id" else "module_label"
truth_coexp <- truth_coexp[!is.na(idx)]
truth_coexp <- unique(truth_coexp[, .(idx, module = get(module_col))])
setorder(truth_coexp, module, idx)

order_idx <- unique(truth_coexp$idx)
order_idx <- order_idx[order_idx >= 1]
if (length(order_idx) > 200) order_idx <- order_idx[1:200]
n_gene <- length(order_idx)
if (n_gene < 2) stop("Insufficient co-expression truth genes for plotting.")

truth_mat <- matrix(0, nrow = n_gene, ncol = n_gene)
for (m in unique(truth_coexp$module)) {
  ids <- truth_coexp[module == m, idx]
  pos <- which(order_idx %in% ids)
  if (length(pos) > 1) truth_mat[pos, pos] <- 1
}
diag(truth_mat) <- 0

pred_full <- load_coexp_matrix(co_method, co_rep_id, max(order_idx))
if (is.null(pred_full)) stop("Missing co-expression prediction for ", co_method, " rep ", co_rep_id)

pred_sub <- pred_full[order_idx, order_idx, drop = FALSE]
pred_bin <- ifelse(abs(pred_sub) >= 0.25, 1, 0)
diag(pred_bin) <- 0

p_co_truth <- plot_binary_matrix(truth_mat, sprintf("Coexp GT edges (S7 rep %d)", co_rep_id))
p_co_pred <- plot_binary_matrix(pred_bin, sprintf("Coexp %s inferred (F1=%.3f)", co_method, coexp_rep$F1))

# 5) Compose and save.
row1 <- plot_grid(p_seg_truth, p_seg_pred, ncol = 2, rel_widths = c(1, 1.05))
row2 <- plot_grid(p_dec_truth, p_dec_pred, ncol = 2, rel_widths = c(1, 1.05))
row3 <- plot_grid(p_co_truth, p_co_pred, ncol = 2, rel_widths = c(1, 1.05))
title <- ggdraw() + draw_label("S7 Best Methods: Ground Truth vs Inferred Results", fontface = "bold", size = 15)
final_plot <- plot_grid(title, row1, row2, row3, ncol = 1, rel_heights = c(0.07, 1, 1, 1))

out_pdf <- file.path(base_dir, "S7_best_methods_groundtruth_vs_inferred.pdf")
ggsave(out_pdf, final_plot, width = 14, height = 17, bg = "white")

out_csv <- file.path(base_dir, "S7_best_methods_selection.csv")
fwrite(selection_dt, out_csv)

message("Saved: ", out_pdf)
message("Saved: ", out_csv)
print(selection_dt)
