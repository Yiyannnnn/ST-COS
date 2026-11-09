#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

strong_palette <- c(
  "#1F77B4", "#D62728", "#2CA02C", "#9467BD", "#FF7F0E",
  "#17BECF", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22"
)

make_strong_cols <- function(levels_chr) {
  lv <- as.character(levels_chr)
  cols <- setNames(rep(strong_palette, length.out = length(lv)), lv)
  # Keep core clusters consistent across all panels; cluster 3 is background.
  if ("1" %in% lv) cols["1"] <- "#1F77B4"
  if ("2" %in% lv) cols["2"] <- "#D62728"
  if ("3" %in% lv) cols["3"] <- "#BDBDBD"
  cols
}

base_dir <- "."
scenario <- "S7"

args <- commandArgs(trailingOnly = TRUE)
rep_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
if (!is.finite(rep_id) || rep_id < 1) rep_id <- 1L

sim_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated")
data_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"))
out_dir <- file.path(base_dir, "S7_example_style_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

coord <- fread(file.path(sim_dir, paste0("S7_coord_", rep_id, ".csv")))
layer <- fread(file.path(sim_dir, paste0("S7_layer_", rep_id, ".csv")))
region <- fread(file.path(sim_dir, paste0("S7_region_", rep_id, ".csv")))
spot_count <- fread(file.path(sim_dir, paste0("S7_gene_expression_", rep_id, ".csv")))
scref <- fread(file.path(sim_dir, paste0("S7_scref_", rep_id, ".csv")))

bass <- fread(file.path(data_dir, paste0("BASS_", rep_id, ".csv")))
bayes <- fread(file.path(data_dir, paste0("BayesSpace_", rep_id, ".csv")))
stdec <- fread(file.path(data_dir, paste0("STdeconvolve_", rep_id, ".csv"))
)
spacex <- as.matrix(fread(file.path(data_dir, paste0("SpaceX_", rep_id, ".csv"))))

ct_truth <- c("Cancer", "Tcell", "Endothelial", "CAF")
ct_truth3 <- c("Cancer", "Tcell", "Endothelial")
gene_cols <- grep("^gene_", names(spot_count), value = TRUE)
if (length(gene_cols) < 50) stop("Need at least 50 genes in S7 spot counts.")
gene50 <- gene_cols[1:50]

plot_corr_heatmap <- function(mat, title_txt) {
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- 0
  if (is.null(rownames(mat))) rownames(mat) <- as.character(seq_len(nrow(mat)))
  if (is.null(colnames(mat))) colnames(mat) <- as.character(seq_len(ncol(mat)))
  dt <- as.data.table(as.table(mat))
  setnames(dt, c("g1", "g2", "rho"))
  dt[, g1 := as.integer(factor(as.character(g1), levels = rownames(mat)))]
  dt[, g2 := as.integer(factor(as.character(g2), levels = colnames(mat)))]
  ggplot(dt, aes(g1, g2, fill = rho)) +
    geom_raster() +
    coord_fixed() +
    scale_fill_gradient2(
      low = "#3b4cc0", mid = "white", high = "#b40426",
      midpoint = 0, limits = c(-1, 1)
    ) +
    labs(title = title_txt, x = NULL, y = NULL, fill = "Spearman") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

plot_prop_facets <- function(dt_long, title_txt) {
  ggplot(dt_long, aes(x = x_coord, y = y_coord, fill = proportion)) +
    geom_point(shape = 21, size = 0.9, stroke = 0) +
    facet_wrap(~cell_type, ncol = 2) +
    coord_equal() +
    scale_fill_gradientn(colors = c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")) +
    labs(title = title_txt, fill = "Prop") +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

plot_label_map <- function(dt, colname, title_txt, cols = NULL, point_size = 0.9, show_legend = TRUE) {
  lv <- sort(unique(as.character(dt[[colname]])))
  if (is.null(cols)) {
    cols <- make_strong_cols(lv)
  }
  legend_breaks <- names(cols)[names(cols) %in% lv]
  ggplot(dt, aes(x = x_coord, y = y_coord, color = as.character(.data[[colname]]))) +
    geom_point(size = point_size, shape = 16) +
    coord_equal() +
    scale_color_manual(values = cols, breaks = legend_breaks, na.value = "#1f1f1f") +
    labs(title = title_txt, color = "Cluster") +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(size = 8)
    )
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

align_labels_to_truth <- function(pred_dt, truth_dt, pred_col, truth_col = "layer") {
  permute_idx <- function(idx) {
    if (length(idx) <= 1) return(list(idx))
    out <- vector("list", 0)
    for (i in seq_along(idx)) {
      rest <- idx[-i]
      subp <- permute_idx(rest)
      out <- c(out, lapply(subp, function(x) c(idx[i], x)))
    }
    out
  }

  pred_use <- copy(pred_dt)
  truth_use <- copy(truth_dt)

  pred_use[, pred_value := as.character(get(pred_col))]
  truth_use[, truth_value := as.character(get(truth_col))]

  ov <- merge(
    pred_use[, .(x_coord, y_coord, pred_value)],
    truth_use[, .(x_coord, y_coord, truth_value)],
    by = c("x_coord", "y_coord"),
    all = FALSE
  )

  truth_levels <- sort(unique(truth_use$truth_value))
  if (nrow(ov) == 0) {
    pred_use[, label_aligned := factor(pred_value, levels = truth_levels)]
    pred_use[, pred_value := NULL]
    return(pred_use)
  }

  pred_levels <- sort(unique(pred_use$pred_value))
  n_pred <- length(pred_levels)
  n_truth <- length(truth_levels)

  conf_mat <- table(
    factor(ov$pred_value, levels = pred_levels),
    factor(ov$truth_value, levels = truth_levels)
  )
  conf_mat <- as.matrix(conf_mat)

  map_vec <- setNames(rep(NA_character_, n_pred), pred_levels)

  # Hungarian-style one-to-one assignment objective (maximize overlap).
  # For n=3 (current S7 use), exhaustive permutation is exact and stable.
  if (n_pred == n_truth && n_pred <= 8) {
    perms <- permute_idx(seq_len(n_truth))
    best_score <- -Inf
    best_perm <- perms[[1]]
    for (p in perms) {
      score <- sum(conf_mat[cbind(seq_len(n_pred), p)])
      if (score > best_score) {
        best_score <- score
        best_perm <- p
      }
    }
    map_vec[pred_levels] <- truth_levels[best_perm]
  } else {
    # Fallback when dimensions differ: greedy unique assignment by overlap.
    conf_dt <- as.data.table(as.table(conf_mat))
    setnames(conf_dt, c("pred_value", "truth_value", "n"))
    setorder(conf_dt, -n)
    used_pred <- character(0)
    used_truth <- character(0)
    for (i in seq_len(nrow(conf_dt))) {
      p <- as.character(conf_dt$pred_value[i])
      t <- as.character(conf_dt$truth_value[i])
      if (!(p %in% used_pred) && !(t %in% used_truth)) {
        map_vec[p] <- t
        used_pred <- c(used_pred, p)
        used_truth <- c(used_truth, t)
      }
    }
    # Any unmapped predicted label falls back to its majority overlap.
    for (p in pred_levels[is.na(map_vec)]) {
      row_counts <- conf_mat[p, , drop = TRUE]
      map_vec[p] <- truth_levels[which.max(row_counts)]
    }
  }

  pred_use[, label_aligned := map_vec[pred_value]]
  pred_use[is.na(label_aligned), label_aligned := pred_value]
  pred_use[, label_aligned := factor(label_aligned, levels = truth_levels)]
  pred_use[, pred_value := NULL]
  pred_use
}

# 1) Cell-type count corr structure from sc reference.
scref_gene_cols <- grep("^gene_", names(scref), value = TRUE)
scref_gene50 <- scref_gene_cols[1:50]
cts <- unique(as.character(scref$Cell_type))
pdf(file.path(out_dir, paste0("S7_rep", rep_id, "_scref_corr_by_celltype_first50.pdf")), width = 6, height = 5)
for (ct in cts) {
  dt_ct <- scref[Cell_type == ct]
  if (nrow(dt_ct) < 5) next
  cor_mat <- suppressWarnings(cor(as.matrix(dt_ct[, ..scref_gene50]), method = "spearman"))
  print(plot_corr_heatmap(cor_mat, paste0("S7 rep ", rep_id, " ", ct, " count corr (first 50 genes)")))
}
dev.off()

# 2) Spot count corr structure.
spot_cor <- suppressWarnings(cor(as.matrix(spot_count[, ..gene50]), method = "spearman"))
p_spot_cor <- plot_corr_heatmap(spot_cor, paste0("S7 rep ", rep_id, " spot count corr (first 50 genes)"))
ggsave(
  filename = file.path(out_dir, paste0("S7_rep", rep_id, "_spot_count_corr_first50.pdf")),
  plot = p_spot_cor,
  width = 6, height = 5, bg = "white"
)

# 3) Cell-type composition and layer segmentation (truth).
truth_prop_long <- melt(
  cbind(region[, .(x_coord, y_coord)], region[, ..ct_truth]),
  id.vars = c("x_coord", "y_coord"),
  variable.name = "cell_type",
  value.name = "proportion"
)
p_truth_prop <- plot_prop_facets(truth_prop_long, paste0("S7 rep ", rep_id, " true cell-type proportion"))
ggsave(
  filename = file.path(out_dir, paste0("S7_rep", rep_id, "_truth_celltype_proportion.pdf")),
  plot = p_truth_prop,
  width = 8, height = 6, bg = "white"
)

layer_dt <- cbind(coord[, .(x_coord, y_coord)], layer = as.factor(layer$layer))
p_layer <- plot_label_map(layer_dt, "layer", paste0("S7 rep ", rep_id, " true layer"))
ggsave(
  filename = file.path(out_dir, paste0("S7_rep", rep_id, "_truth_layer_map.pdf")),
  plot = p_layer,
  width = 6, height = 5, bg = "white"
)

# 4) Domain segmentation: truth + all method predictions.
bass_dt <- bass[, .(x_coord, y_coord, layer = as.factor(layer), BASS_layer = as.factor(BASS_layer))]
bayes_dt <- bayes[, .(x_coord, y_coord, layer = as.factor(layer), BayesSpace_labels = as.factor(BayesSpace_labels))]

truth_levels <- sort(unique(as.character(layer_dt$layer)))
shared_cols <- make_strong_cols(truth_levels)

bass_pred_aligned <- align_labels_to_truth(
  bass_dt[, .(x_coord, y_coord, label = as.character(BASS_layer))],
  layer_dt[, .(x_coord, y_coord, layer = as.character(layer))],
  pred_col = "label",
  truth_col = "layer"
)
bayes_pred_aligned <- align_labels_to_truth(
  bayes_dt[, .(x_coord, y_coord, label = as.character(BayesSpace_labels))],
  layer_dt[, .(x_coord, y_coord, layer = as.character(layer))],
  pred_col = "label",
  truth_col = "layer"
)

p_bass_truth <- plot_label_map(
  layer_dt, "layer", paste0("S7 rep ", rep_id, " truth layer"),
  cols = shared_cols, point_size = 2.8
)
p_bass_pred <- plot_label_map(
  bass_pred_aligned, "label_aligned", paste0("S7 rep ", rep_id, " BASS prediction"),
  cols = shared_cols, point_size = 2.8
)
ggsave(
  filename = file.path(out_dir, paste0("S7_rep", rep_id, "_BASS_vs_truth.pdf")),
  plot = cowplot::plot_grid(p_bass_truth, p_bass_pred, ncol = 2),
  width = 12, height = 5, bg = "white"
)

p_bayes_truth <- plot_label_map(
  layer_dt, "layer", paste0("S7 rep ", rep_id, " truth layer"),
  cols = shared_cols, point_size = 2.8
)
p_bayes_pred <- plot_label_map(
  bayes_pred_aligned, "label_aligned", paste0("S7 rep ", rep_id, " BayesSpace prediction"),
  cols = shared_cols, point_size = 2.8
)
ggsave(
  filename = file.path(out_dir, paste0("S7_rep", rep_id, "_BayesSpace_vs_truth.pdf")),
  plot = cowplot::plot_grid(p_bayes_truth, p_bayes_pred, ncol = 2),
  width = 12, height = 5, bg = "white"
)

graphst_file <- file.path(data_dir, paste0("GraphST_", rep_id, ".csv"))
maple_file <- file.path(data_dir, paste0("MAPLE_", rep_id, ".csv"))
spagcn_file <- file.path(data_dir, paste0("SpaGCN_", rep_id, ".csv"))
stagate_file <- file.path(data_dir, paste0("STAGATE_", rep_id, ".csv"))

pred_panels <- list(
  list(name = "Truth", dt = layer_dt[, .(x_coord, y_coord, label = as.character(layer))]),
  list(name = "BASS", dt = bass_pred_aligned[, .(x_coord, y_coord, label = as.character(label_aligned))]),
  list(name = "BayesSpace", dt = bayes_pred_aligned[, .(x_coord, y_coord, label = as.character(label_aligned))])
)

if (file.exists(graphst_file)) {
  graphst_raw <- fread(graphst_file)[, .(
    x_coord = as.numeric(x),
    y_coord = as.numeric(y),
    label = as.character(domain)
  )]
  graphst_dt <- align_labels_to_truth(
    graphst_raw, layer_dt[, .(x_coord, y_coord, layer = as.character(layer))],
    pred_col = "label", truth_col = "layer"
  )[, .(x_coord, y_coord, label = as.character(label_aligned))]
  pred_panels[[length(pred_panels) + 1L]] <- list(name = "GraphST", dt = graphst_dt)
}

if (file.exists(maple_file)) {
  maple_raw <- fread(maple_file)[, .(
    x_coord = as.numeric(x_coord),
    y_coord = as.numeric(y_coord),
    label = as.character(cluster)
  )]
  maple_dt <- align_labels_to_truth(
    maple_raw, layer_dt[, .(x_coord, y_coord, layer = as.character(layer))],
    pred_col = "label", truth_col = "layer"
  )[, .(x_coord, y_coord, label = as.character(label_aligned))]
  pred_panels[[length(pred_panels) + 1L]] <- list(name = "MAPLE", dt = maple_dt)
}

if (file.exists(spagcn_file)) {
  spagcn_raw <- fread(spagcn_file)[, .(
    x_coord = as.numeric(x),
    y_coord = as.numeric(y),
    label = as.character(refined_pred)
  )]
  spagcn_dt <- align_labels_to_truth(
    spagcn_raw, layer_dt[, .(x_coord, y_coord, layer = as.character(layer))],
    pred_col = "label", truth_col = "layer"
  )[, .(x_coord, y_coord, label = as.character(label_aligned))]
  pred_panels[[length(pred_panels) + 1L]] <- list(name = "SpaGCN", dt = spagcn_dt)
}

if (file.exists(stagate_file)) {
  stagate_raw <- fread(stagate_file)[, .(
    x_coord = as.numeric(x),
    y_coord = as.numeric(y),
    label = as.character(mclust)
  )]
  stagate_dt <- align_labels_to_truth(
    stagate_raw, layer_dt[, .(x_coord, y_coord, layer = as.character(layer))],
    pred_col = "label", truth_col = "layer"
  )[, .(x_coord, y_coord, label = as.character(label_aligned))]
  pred_panels[[length(pred_panels) + 1L]] <- list(name = "STAGATE", dt = stagate_dt)
}

seg_legend_plot <- plot_label_map(
  pred_panels[[1]]$dt,
  "label",
  paste0("S7 rep ", rep_id, " ", pred_panels[[1]]$name),
  cols = shared_cols,
  point_size = 2.8,
  show_legend = TRUE
)
seg_legend <- cowplot::get_legend(seg_legend_plot)

seg_plot_list <- lapply(pred_panels, function(x) {
  plot_label_map(
    x$dt,
    "label",
    paste0("S7 rep ", rep_id, " ", x$name),
    cols = shared_cols,
    point_size = 2.8,
    show_legend = FALSE
  )
})

ggsave(
  filename = file.path(out_dir, paste0("S7_rep", rep_id, "_all_segmentation_methods_vs_truth.pdf")),
  plot = cowplot::plot_grid(
    cowplot::plot_grid(plotlist = seg_plot_list, nrow = 1),
    seg_legend,
    rel_widths = c(1, 0.10)
  ),
  width = 26, height = 4.5, bg = "white"
)

# 5) STdeconvolve composition (3 topics) and truth (3 cell types).
if (ncol(stdec) >= 3) {
  n_use <- min(nrow(coord), nrow(stdec), nrow(region))
  coord_use <- coord[seq_len(n_use), .(x_coord, y_coord)]
  pred_raw <- as.data.frame(stdec[seq_len(n_use), 1:3])
  colnames(pred_raw) <- c("Topic1", "Topic2", "Topic3")
  truth3 <- as.data.frame(region[seq_len(n_use), ..ct_truth3])

  cor_map <- suppressWarnings(cor(pred_raw, truth3))
  perms <- rbind(
    c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
    c(2, 3, 1), c(3, 1, 2), c(3, 2, 1)
  )
  best_idx <- 1L
  best_score <- -Inf
  for (k in seq_len(nrow(perms))) {
    score <- sum(cor_map[cbind(1:3, perms[k, ])], na.rm = TRUE)
    if (is.finite(score) && score > best_score) {
      best_score <- score
      best_idx <- k
    }
  }
  perm_best <- perms[best_idx, ]
  pred_aligned <- pred_raw
  colnames(pred_aligned) <- colnames(truth3)[perm_best]
  pred_aligned <- pred_aligned[, ct_truth3, drop = FALSE]
  rmse <- sqrt(mean((as.matrix(pred_aligned) - as.matrix(truth3))^2, na.rm = TRUE))

  truth3_long <- melt(
    cbind(coord_use, truth3),
    id.vars = c("x_coord", "y_coord"),
    variable.name = "cell_type", value.name = "proportion"
  )
  pred3_long <- melt(
    cbind(coord_use, pred_aligned),
    id.vars = c("x_coord", "y_coord"),
    variable.name = "cell_type", value.name = "proportion"
  )
  p_truth3 <- plot_prop_facets(truth3_long, paste0("S7 rep ", rep_id, " truth (Cancer/Tcell/Endothelial)"))
  p_pred3 <- plot_prop_facets(pred3_long, paste0("S7 rep ", rep_id, " STdeconvolve aligned (RMSE=", sprintf("%.4f", rmse), ")"))

  ggsave(
    filename = file.path(out_dir, paste0("S7_rep", rep_id, "_STdeconvolve_vs_truth.pdf")),
    plot = cowplot::plot_grid(p_truth3, p_pred3, ncol = 2),
    width = 14, height = 6, bg = "white"
  )
}

# 6) SpaceX co-expression heatmap (first 50 genes).
if (nrow(spacex) >= 50 && ncol(spacex) >= 50) {
  spacex50 <- spacex[1:50, 1:50, drop = FALSE]
  p_spacex <- plot_corr_heatmap(spacex50, paste0("S7 rep ", rep_id, " SpaceX co-expression (first 50 genes)"))
  ggsave(
    filename = file.path(out_dir, paste0("S7_rep", rep_id, "_SpaceX_coexpression_first50.pdf")),
    plot = p_spacex,
    width = 6, height = 5, bg = "white"
  )
}

# 7) Optional UMAP, if dependencies are available.
has_umap_pkgs <- all(sapply(
  c("SingleCellExperiment", "scater", "scuttle"),
  requireNamespace, quietly = TRUE
))
if (has_umap_pkgs) {
  tryCatch({
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = t(as.matrix(spot_count)))
    )
    sce$layer <- as.factor(layer$layer)
    sce <- scuttle::logNormCounts(sce)
    sce <- scater::runPCA(sce, exprs_values = "logcounts", scale = TRUE)
    sce <- scater::runUMAP(sce, dimred = "PCA", name = "UMAP")
    um <- as.data.table(SingleCellExperiment::reducedDim(sce, "UMAP"))
    setnames(um, c("UMAP1", "UMAP2"))
    um[, layer := as.factor(layer$layer)]

    p_umap <- ggplot(um, aes(UMAP1, UMAP2, color = layer)) +
      geom_point(size = 0.8) +
      labs(title = paste0("S7 rep ", rep_id, " UMAP by true layer"), color = "Layer") +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    ggsave(
      filename = file.path(out_dir, paste0("S7_rep", rep_id, "_UMAP_by_layer.pdf")),
      plot = p_umap,
      width = 6, height = 5, bg = "white"
    )
  }, error = function(e) {
    message("UMAP plot skipped: ", e$message)
  })
}

message("Saved S7 example-style plots to: ", out_dir)
