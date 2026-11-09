#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

source("./Data/simulator_YZ.R")

set.seed(20260905)

repo_dir <- "."
slide_id <- "151673"
dlpfc_dir <- file.path(repo_dir, "Reference_based", "DLPFC")
out_dir <- file.path(repo_dir, "Independent_simulator_anchor", "DLPFC_layer_SVG_controlled_demo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

coords_file <- file.path(dlpfc_dir, sprintf("%s_simulated_coordinates.csv", slide_id))
stcos_file <- file.path(dlpfc_dir, sprintf("%s_STCOS_simulated_counts.csv", slide_id))
srtsim_file <- file.path(dlpfc_dir, sprintf("%s_SRTsim_counts.csv", slide_id))

read_spot_gene <- function(path) {
  dt <- fread(path, check.names = FALSE)
  ids <- as.character(dt[[1]])
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  colnames(mat) <- names(dt)[-1]
  mat
}

read_gene_spot <- function(path) {
  dt <- fread(path, check.names = FALSE)
  genes <- as.character(dt[[1]])
  spots <- names(dt)[-1]
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes
  colnames(mat) <- spots
  mat <- t(mat)
  colnames(mat) <- genes
  mat
}

fit_nb_fast <- function(y) {
  mu <- mean(y)
  s2 <- var(y)
  theta <- if (is.finite(s2) && s2 > mu && mu > 0) mu^2 / (s2 - mu) else 1e6
  min(max(theta, 1e-4), 1e6)
}

safe_spearman <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

eta_by_group <- function(x, group) {
  x <- log1p(x)
  mu <- mean(x)
  ss_between <- sum(vapply(split(x, group), function(v) length(v) * (mean(v) - mu)^2, numeric(1)))
  ss_total <- sum((x - mu)^2)
  if (ss_total == 0) return(NA_real_)
  ss_between / ss_total
}

dominant_layer <- function(x, group) {
  means <- tapply(log1p(x), group, mean)
  names(means)[which.max(means)]
}

message("Reading DLPFC reference-based outputs ...")
coords <- fread(coords_file, check.names = FALSE)
setnames(coords, names(coords)[1], "spot_id")
stcos_counts <- read_spot_gene(stcos_file)
srtsim_counts <- read_gene_spot(srtsim_file)

common_spots <- base::Reduce(base::intersect, list(coords$spot_id, rownames(stcos_counts), rownames(srtsim_counts)))
common_genes <- base::intersect(colnames(stcos_counts), colnames(srtsim_counts))
if (length(common_spots) == 0 || length(common_genes) == 0) {
  stop("No common spots or genes across coordinates, ST-COS, and SRTsim outputs.")
}

coords <- coords[match(common_spots, spot_id)]
stcos_counts <- stcos_counts[common_spots, common_genes, drop = FALSE]
srtsim_counts <- srtsim_counts[common_spots, common_genes, drop = FALSE]

layer_order <- c("L1", "L2", "L3", "L4", "L5", "L6", "WM")
coords[, label := factor(as.character(label), levels = layer_order)]
coords[, motif := label %in% c("L4", "L5")]
message(sprintf("Layer-like user band covers %d of %d spots (%.1f%%).",
                sum(coords$motif), nrow(coords), 100 * mean(coords$motif)))

message("Selecting target genes from known DLPFC SVG/layer-marker candidates ...")
known_svg <- c(
  "MOBP", "GFAP", "NEFL", "PCP4", "PLP1", "MBP", "SNAP25", "NEFM",
  "RORB", "CUX2", "RELN", "TLE4", "BCL11B", "FEZF2", "AQP4", "PVALB",
  "GAD1", "GAD2", "SST", "VIP", "CALB2"
)
known_svg <- base::intersect(known_svg, common_genes)

gene_stats <- data.table(
  gene = common_genes,
  mean_expr = colMeans(stcos_counts),
  zero_prop = colMeans(stcos_counts == 0)
)
gene_stats[, band_rho := vapply(
  gene,
  function(g) safe_spearman(log1p(stcos_counts[, g]), as.numeric(coords$motif)),
  numeric(1)
)]

known_stats <- gene_stats[gene %in% known_svg]
known_stats[, layer_eta := vapply(
  gene,
  function(g) eta_by_group(stcos_counts[, g], coords$label),
  numeric(1)
)]
known_stats[, dominant_layer := vapply(
  gene,
  function(g) dominant_layer(stcos_counts[, g], coords$label),
  character(1)
)]

target_candidates <- known_stats[
  mean_expr >= 0.10 &
    mean_expr <= 5.00 &
    zero_prop <= 0.80 &
    layer_eta >= 0.15
][order(-layer_eta, zero_prop)]

if (nrow(target_candidates) < 4) {
  target_candidates <- known_stats[
    mean_expr >= 0.05 &
      zero_prop <= 0.90
  ][order(-layer_eta, zero_prop)]
}
preferred_targets <- base::intersect(c("PCP4", "RORB", "NEFL", "NEFM"), common_genes)
fill_targets <- target_candidates[!gene %in% preferred_targets]$gene
target_genes <- head(c(preferred_targets, fill_targets), 4)
if (length(target_genes) < 3) stop("Too few known SVG targets available for the layer-pattern demo.")

target_mean <- mean(gene_stats[gene %in% target_genes, mean_expr])
background_pool <- gene_stats[
  !gene %in% target_genes &
    mean_expr >= 0.05 &
    zero_prop <= 0.90 &
    abs(band_rho) <= 0.08
]
background_pool[, mean_distance := abs(mean_expr - target_mean)]
background_genes <- head(background_pool[order(mean_distance, abs(band_rho))]$gene, min(96, nrow(background_pool)))
selected_genes <- c(target_genes, background_genes)

message(sprintf("Target SVGs: %s", paste(target_genes, collapse = ", ")))
message(sprintf("Background genes: %d matched genes", length(background_genes)))

base_subset <- stcos_counts[, selected_genes, drop = FALSE]

message("Building a layer-aware ST-COS mean surface and adding the user band to target SVGs ...")
mu_mat <- matrix(NA_real_, nrow = nrow(base_subset), ncol = ncol(base_subset),
                 dimnames = list(common_spots, selected_genes))
for (g in selected_genes) {
  layer_means <- tapply(base_subset[, g], coords$label, mean)
  layer_means[is.na(layer_means)] <- mean(base_subset[, g])
  mu_mat[, g] <- as.numeric(layer_means[as.character(coords$label)])
}
mu_mat[!is.finite(mu_mat)] <- 0
mu_mat <- pmax(mu_mat, 1e-6)

fold_change <- 5
mu_mat[coords$motif, target_genes] <- mu_mat[coords$motif, target_genes, drop = FALSE] * fold_change

theta_g <- vapply(seq_len(ncol(base_subset)), function(i) fit_nb_fast(base_subset[, i]), numeric(1))
names(theta_g) <- selected_genes

sigma <- cor(log1p(base_subset), use = "pairwise.complete.obs")
sigma[is.na(sigma)] <- 0
diag(sigma) <- 1
sigma[abs(sigma) < 0.05] <- 0
sigma <- cov2cor(makespd(sigma))
dimnames(sigma) <- list(selected_genes, selected_genes)

gene_rel <- as.data.table(matrix(1, nrow = length(selected_genes), ncol = 1))
setnames(gene_rel, "CT")

cov_gene <- array(
  NA_real_,
  dim = c(length(selected_genes), length(selected_genes), 1),
  dimnames = list(selected_genes, selected_genes, "CT")
)
cov_gene[, , 1] <- sigma

message("Generating ST-COS controlled layer-SVG counts for the selected gene panel ...")
controlled_list <- generate_true_count(
  coord = as.matrix(coords[, .(row, col)]),
  mean_param = list(CT = mu_mat),
  gene_relationship = gene_rel,
  cov_gene = cov_gene,
  lib_size = list(CT = rep(1, nrow(coords))),
  theta_g = list(CT = theta_g)
)

stcos_controlled <- controlled_list[["CT"]]
rownames(stcos_controlled) <- common_spots
colnames(stcos_controlled) <- selected_genes

score_matrix <- function(mat, method_name, genes, gene_class) {
  rbindlist(lapply(genes, function(g) {
    vals <- log1p(mat[, g])
    data.table(
      method = method_name,
      gene_class = gene_class,
      gene = g,
      band_contrast = mean(vals[coords$motif]) - mean(vals[!coords$motif]),
      band_spearman = safe_spearman(vals, as.numeric(coords$motif)),
      mean_inside = mean(vals[coords$motif]),
      mean_outside = mean(vals[!coords$motif])
    )
  }))
}

metric_dt <- rbindlist(list(
  score_matrix(stcos_counts[, selected_genes, drop = FALSE], "ST-COS replay", target_genes, "target SVG"),
  score_matrix(stcos_counts[, selected_genes, drop = FALSE], "ST-COS replay", background_genes, "background"),
  score_matrix(srtsim_counts[, selected_genes, drop = FALSE], "SRTsim replay", target_genes, "target SVG"),
  score_matrix(srtsim_counts[, selected_genes, drop = FALSE], "SRTsim replay", background_genes, "background"),
  score_matrix(stcos_controlled, "ST-COS controlled", target_genes, "target SVG"),
  score_matrix(stcos_controlled, "ST-COS controlled", background_genes, "background")
))

method_levels <- c("ST-COS replay", "SRTsim replay", "ST-COS controlled")
metric_dt[, method := factor(method, levels = method_levels)]
metric_dt[, gene_class := factor(gene_class, levels = c("background", "target SVG"))]

plot_genes <- target_genes
map_dt <- rbindlist(lapply(plot_genes, function(g) {
  rbindlist(list(
    data.table(coords[, .(spot_id, row, col, label, motif)], gene = g, method = "ST-COS replay",
               expr = log1p(stcos_counts[, g])),
    data.table(coords[, .(spot_id, row, col, label, motif)], gene = g, method = "SRTsim replay",
               expr = log1p(srtsim_counts[, g])),
    data.table(coords[, .(spot_id, row, col, label, motif)], gene = g, method = "ST-COS controlled",
               expr = log1p(stcos_controlled[, g]))
  ))
}))
map_dt[, method := factor(method, levels = method_levels)]
map_dt[, gene := factor(gene, levels = plot_genes)]
expr_cap <- quantile(map_dt$expr, 0.99, na.rm = TRUE)
map_dt[, expr_plot := pmin(expr, expr_cap)]

metric_summary <- metric_dt[, .(
  n_genes = .N,
  median_contrast = median(band_contrast, na.rm = TRUE),
  mean_contrast = mean(band_contrast, na.rm = TRUE),
  median_spearman = median(band_spearman, na.rm = TRUE),
  mean_spearman = mean(band_spearman, na.rm = TRUE)
), by = .(method, gene_class)]

target_out <- merge(
  data.table(target_gene = target_genes),
  known_stats,
  by.x = "target_gene",
  by.y = "gene",
  all.x = TRUE
)
target_out[, `:=`(slide_id = slide_id, fold_change = fold_change, user_band = "L4/L5")]

fwrite(target_out, file.path(out_dir, "DLPFC_layer_SVG_target_genes.csv"))
fwrite(metric_dt, file.path(out_dir, "DLPFC_layer_SVG_gene_metrics.csv"))
fwrite(metric_summary, file.path(out_dir, "DLPFC_layer_SVG_metric_summary.csv"))
fwrite(map_dt, file.path(out_dir, "DLPFC_layer_SVG_map_values.csv"))
fwrite(
  data.table(spot_id = rownames(stcos_controlled), as.data.table(stcos_controlled)),
  file.path(out_dir, "DLPFC_layer_SVG_STCOS_selected_counts.csv")
)

caption_text <- sprintf("%s; L4/L5; %d spots; x%d", slide_id, sum(coords$motif), fold_change)

p_band <- ggplot(coords, aes(col, row)) +
  geom_point(aes(color = motif), size = 0.65, alpha = 0.9) +
  scale_color_manual(values = c("FALSE" = "grey83", "TRUE" = "#c43c39"), labels = c("Other layers", "User band")) +
  coord_equal() +
  labs(title = "A. Layer-parallel user band", subtitle = caption_text, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  )

p_maps <- ggplot(map_dt, aes(col, row, color = expr_plot)) +
  geom_point(size = 0.38) +
  coord_equal() +
  facet_grid(gene ~ method) +
  scale_color_gradientn(colors = c("#f2f2f2", "#9fc2d4", "#385991", "#8c2046"), name = "log1p(count)") +
  labs(title = "B. Spatial maps for selected SVG/layer-marker genes") +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    strip.text = element_text(size = 9),
    legend.position = "right"
  )

p_contrast <- ggplot(metric_dt, aes(method, band_contrast, fill = gene_class)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.70, width = 0.65) +
  geom_point(aes(color = gene_class), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65),
             size = 1.0, alpha = 0.85) +
  scale_fill_manual(values = c(background = "grey80", `target SVG` = "#d95f02")) +
  scale_color_manual(values = c(background = "grey35", `target SVG` = "#d95f02")) +
  labs(title = "C. User-band contrast", x = NULL, y = "mean log1p inside - outside", fill = NULL, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom"
  )

p_spearman <- ggplot(metric_dt, aes(method, band_spearman, fill = gene_class)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.70, width = 0.65) +
  geom_point(aes(color = gene_class), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65),
             size = 1.0, alpha = 0.85) +
  scale_fill_manual(values = c(background = "grey80", `target SVG` = "#1b9e77")) +
  scale_color_manual(values = c(background = "grey35", `target SVG` = "#1b9e77")) +
  labs(title = "D. Association with user band", x = NULL, y = "Spearman rho", fill = NULL, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom"
  )

top_row <- plot_grid(p_band, p_contrast, p_spearman, ncol = 3, rel_widths = c(0.9, 1.15, 1.15))
combined <- plot_grid(top_row, p_maps, ncol = 1, rel_heights = c(1.05, 2.35))

pdf_out <- file.path(out_dir, sprintf("DLPFC_%s_STCOS_layer_SVG_controllability_demo.pdf", slide_id))
png_out <- file.path(out_dir, sprintf("DLPFC_%s_STCOS_layer_SVG_controllability_demo.png", slide_id))
ggsave(pdf_out, combined, width = 13.5, height = 11.0, bg = "white")
ggsave(png_out, combined, width = 13.5, height = 11.0, dpi = 300, bg = "white")

print(target_out)
print(metric_summary)
message(sprintf("Saved figure: %s", pdf_out))
message(sprintf("Saved figure: %s", png_out))
