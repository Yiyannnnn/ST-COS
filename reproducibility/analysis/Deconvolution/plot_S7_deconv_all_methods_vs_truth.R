#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(scatterpie)
})

base_dir <- "."
scenario <- "S7"
data_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"))
sim_dir <- file.path(data_dir, "simulated")
out_dir <- file.path(base_dir, "Deconvolution")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

method_order <- c("CARD", "RCTD", "STdeconvolve", "STdGCN", "Tangram", "POLARIS")
ct_levels <- c("Cancer", "Tcell", "Endothelial", "CAF")
ct_colors <- c("Cancer" = "#ff2e00", "Tcell" = "#29e600", "Endothelial" = "#1149ff", "CAF" = "#8b5cf6")

permute_idx <- function(idx) {
  if (length(idx) <= 1) return(list(idx))
  out <- vector("list", 0)
  for (i in seq_along(idx)) {
    rest <- idx[-i]
    sub <- permute_idx(rest)
    out <- c(out, lapply(sub, function(x) c(idx[i], x)))
  }
  out
}

read_truth <- function(rep_id) {
  region_file <- file.path(sim_dir, paste0(scenario, "_region_", rep_id, ".csv"))
  layer_file <- file.path(sim_dir, paste0(scenario, "_layer_", rep_id, ".csv"))
  region_dt <- fread(region_file)
  layer_dt <- fread(layer_file)
  truth_prop <- as.data.frame(region_dt[, ..ct_levels])
  rownames(truth_prop) <- paste0("spot_", seq_len(nrow(truth_prop)))
  list(
    coord = region_dt[, .(x_coord, y_coord)],
    prop = truth_prop,
    layer = as.character(layer_dt$layer)
  )
}

normalize_rows <- function(mat) {
  rs <- rowSums(mat, na.rm = TRUE)
  rs[!is.finite(rs) | rs <= 0] <- 1
  mat / rs
}

align_stdeconv_topics <- function(pred, truth_prop) {
  if (ncol(pred) != 3) return(pred)
  truth3 <- truth_prop[, c("Cancer", "Tcell", "Endothelial"), drop = FALSE]

  common <- intersect(rownames(pred), rownames(truth3))
  if (length(common) == 0) {
    n <- min(nrow(pred), nrow(truth3))
    if (n < 2) return(pred)
    pred_use <- pred[seq_len(n), , drop = FALSE]
    truth_use <- truth3[seq_len(n), , drop = FALSE]
  } else {
    pred_use <- pred[common, , drop = FALSE]
    truth_use <- truth3[common, , drop = FALSE]
  }
  if (nrow(pred_use) < 2) return(pred)

  cor_mat <- suppressWarnings(cor(as.matrix(pred_use), as.matrix(truth_use), use = "pairwise.complete.obs"))
  if (!all(is.finite(cor_mat))) cor_mat[!is.finite(cor_mat)] <- -1e6

  perms <- permute_idx(1:3)
  best <- perms[[1]]
  best_score <- -Inf
  for (p in perms) {
    score <- sum(cor_mat[cbind(1:3, p)])
    if (is.finite(score) && score > best_score) {
      best_score <- score
      best <- p
    }
  }
  colnames(pred) <- colnames(truth_use)[best]
  pred
}

align_pred_to_truth <- function(pred, n_truth) {
  truth_ids <- paste0("spot_", seq_len(n_truth))
  out <- matrix(0, nrow = n_truth, ncol = length(ct_levels), dimnames = list(truth_ids, ct_levels))
  if (is.null(pred) || nrow(pred) == 0) return(as.data.frame(out))

  for (ct in ct_levels) {
    if (!(ct %in% colnames(pred))) pred[[ct]] <- 0
  }
  pred <- pred[, ct_levels, drop = FALSE]

  common <- intersect(truth_ids, rownames(pred))
  if (length(common) > 0) {
    out[common, ] <- as.matrix(pred[common, ct_levels, drop = FALSE])
  } else {
    n <- min(n_truth, nrow(pred))
    out[seq_len(n), ] <- as.matrix(pred[seq_len(n), ct_levels, drop = FALSE])
  }
  as.data.frame(normalize_rows(out))
}

estimate_pie_radius <- function(coord_dt) {
  xy <- as.matrix(coord_dt[, .(x_coord, y_coord)])
  n <- nrow(xy)
  if (n < 2) return(0.6)

  pick <- seq_len(min(n, 300))
  nn <- rep(NA_real_, length(pick))
  for (k in seq_along(pick)) {
    i <- pick[k]
    d2 <- (xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2
    d2[i] <- Inf
    nn[k] <- sqrt(min(d2, na.rm = TRUE))
  }
  r <- 0.42 * median(nn[is.finite(nn) & nn > 0], na.rm = TRUE)
  if (!is.finite(r) || r <= 0) r <- 0.6
  r
}

read_pred <- function(method, rep_id, truth_prop) {
  file_path <- switch(
    method,
    "CARD" = file.path(data_dir, paste0("CARD_", rep_id, ".csv")),
    "RCTD" = file.path(data_dir, paste0("RCTD_", rep_id, ".csv")),
    "Redeconve" = file.path(data_dir, paste0("Redeconve_", rep_id, ".csv")),
    "STdeconvolve" = file.path(data_dir, paste0("STdeconvolve_", rep_id, ".csv")),
    "STdGCN" = file.path(data_dir, paste0("STdGCN_input_", rep_id), "output", "predict_result.csv"),
    "Tangram" = file.path(data_dir, paste0("Tangram_", rep_id, ".csv")),
    "POLARIS" = file.path(data_dir, paste0("POLARIS_", rep_id, ".csv")),
    stop("Unknown method: ", method)
  )
  if (!file.exists(file_path)) return(NULL)

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
  # Tangram output commonly has an index-like first column V1.
  if ("V1" %in% names(dt) && method == "Tangram") {
    dt[, V1 := NULL]
  }

  pred <- as.data.frame(dt)
  pred <- pred[, sapply(pred, is.numeric), drop = FALSE]
  if (ncol(pred) == 0 || nrow(pred) == 0) return(NULL)

  if (is.null(spot_ids)) {
    rownames(pred) <- paste0("spot_", seq_len(nrow(pred)))
  } else {
    rownames(pred) <- spot_ids
  }

  # Standardize common naming patterns.
  if (all(c("ct1", "ct2", "ct3") %in% colnames(pred))) {
    colnames(pred)[match(c("ct1", "ct2", "ct3"), colnames(pred))] <- c("Cancer", "Tcell", "Endothelial")
  }
  if (method == "STdeconvolve") {
    pred <- align_stdeconv_topics(pred, truth_prop)
  }

  keep <- intersect(colnames(pred), ct_levels)
  pred <- pred[, keep, drop = FALSE]
  if (ncol(pred) == 0) return(NULL)

  # Fill missing cell types with 0 so dominant-cell map is comparable.
  for (ct in ct_levels) {
    if (!(ct %in% colnames(pred))) pred[[ct]] <- 0
  }
  pred <- pred[, ct_levels, drop = FALSE]
  pred[] <- lapply(pred, function(x) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- 0
    x
  })
  pred <- as.data.frame(normalize_rows(as.matrix(pred)))
  rownames(pred) <- if (is.null(spot_ids)) paste0("spot_", seq_len(nrow(pred))) else spot_ids
  pred
}

plot_pie_map <- function(dt, title_txt, pie_r, layer_cols, show_legend = FALSE) {
  ggplot(dt) +
    scatterpie::geom_scatterpie(
      aes(x = x_coord, y = y_coord, r = pie_r, color = layer_group),
      cols = ct_levels, size = 0.28
    ) +
    coord_equal() +
    scale_fill_manual(values = ct_colors, breaks = ct_levels, drop = FALSE, name = "Cell type") +
    scale_color_manual(values = layer_cols, drop = FALSE, name = "Cluster") +
    guides(
      color = guide_legend(override.aes = list(fill = NA, shape = 22, size = 6, stroke = 1.1)),
      fill = guide_legend(override.aes = list(shape = 22, size = 6, color = "black"))
    ) +
    labs(title = title_txt) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(size = 8)
    )
}

plot_missing_panel <- function(title_txt, msg = "Missing output") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = msg, color = "#777777", size = 4.2) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = title_txt) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

shift_second_row_right <- function(plotlist, ncol = 4L) {
  if (length(plotlist) <= ncol) return(plotlist)
  blank <- ggplot() + theme_void()
  c(
    plotlist[seq_len(ncol)],
    list(blank),
    plotlist[seq.int(ncol + 1L, length(plotlist))]
  )
}

reps <- 1:5
for (rep_id in reps) {
  truth_obj <- read_truth(rep_id)
  truth_prop <- truth_obj$prop
  truth_coord <- truth_obj$coord
  layer_group <- factor(truth_obj$layer, levels = sort(unique(truth_obj$layer)))
  layer_cols <- setNames(c("#ff3b1f", "#66d60a", "#00d7de", "#6e3bff")[seq_along(levels(layer_group))], levels(layer_group))
  pie_r <- estimate_pie_radius(truth_coord)

  truth_map <- cbind(
    truth_coord,
    data.table(layer_group = layer_group),
    as.data.table(truth_prop)
  )
  truth_map$pie_r <- pie_r

  legend_plot <- plot_pie_map(truth_map, paste0("S7 rep ", rep_id, " Truth"), pie_r = pie_r, layer_cols = layer_cols, show_legend = TRUE)
  legend_grob <- cowplot::get_legend(legend_plot)

  plot_list <- list(plot_pie_map(truth_map, paste0("S7 rep ", rep_id, " Truth"), pie_r = pie_r, layer_cols = layer_cols, show_legend = FALSE))

  for (method in method_order) {
    pred_prop <- read_pred(method, rep_id, truth_prop)
    if (is.null(pred_prop)) {
      plot_list[[length(plot_list) + 1L]] <- plot_missing_panel(
        paste0("S7 rep ", rep_id, " ", method),
        "Missing output"
      )
      next
    }

    pred_use <- align_pred_to_truth(pred_prop, n_truth = nrow(truth_prop))
    plot_dt <- cbind(
      truth_coord,
      data.table(layer_group = layer_group),
      as.data.table(pred_use)
    )
    plot_dt$pie_r <- pie_r

    plot_list[[length(plot_list) + 1L]] <- plot_pie_map(
      plot_dt, paste0("S7 rep ", rep_id, " ", method), pie_r = pie_r, layer_cols = layer_cols, show_legend = FALSE
    )
  }

  grid <- cowplot::plot_grid(plotlist = plot_list, nrow = 1)
  final <- cowplot::plot_grid(grid, legend_grob, rel_widths = c(1, 0.10))

  out_file <- file.path(out_dir, paste0("S7_rep", rep_id, "_all_deconvolution_methods_vs_truth.pdf"))
  ggsave(out_file, final, width = 26, height = 4.5, bg = "white")
  message("Saved: ", out_file)
}
