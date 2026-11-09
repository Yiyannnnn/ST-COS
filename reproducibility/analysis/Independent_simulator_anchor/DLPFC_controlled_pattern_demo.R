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
out_dir <- file.path(repo_dir, "Independent_simulator_anchor", "DLPFC_controlled_pattern_demo")
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

message("Reading DLPFC reference-based outputs ...")
coords <- fread(coords_file, check.names = FALSE)
setnames(coords, names(coords)[1], "spot_id")
stcos_counts <- read_spot_gene(stcos_file)
srtsim_counts <- read_gene_spot(srtsim_file)

common_spots <- Reduce(intersect, list(coords$spot_id, rownames(stcos_counts), rownames(srtsim_counts)))
common_genes <- intersect(colnames(stcos_counts), colnames(srtsim_counts))
if (length(common_spots) == 0 || length(common_genes) == 0) {
  stop("No common spots or genes across coordinates, ST-COS, and SRTsim outputs.")
}

coords <- coords[match(common_spots, spot_id)]
stcos_counts <- stcos_counts[common_spots, common_genes, drop = FALSE]
srtsim_counts <- srtsim_counts[common_spots, common_genes, drop = FALSE]

coords[, x01 := (col - min(col)) / (max(col) - min(col))]
coords[, y01 := (row - min(row)) / (max(row) - min(row))]
coords[, motif := ((x01 - 0.66) / 0.18)^2 + ((y01 - 0.36) / 0.16)^2 <= 1]
if (mean(coords$motif) < 0.06 || mean(coords$motif) > 0.30) {
  coords[, motif := ((x01 - 0.58) / 0.22)^2 + ((y01 - 0.42) / 0.18)^2 <= 1]
}
message(sprintf("User-defined hotspot covers %d of %d spots (%.1f%%).",
                sum(coords$motif), nrow(coords), 100 * mean(coords$motif)))

message("Selecting genes with moderate expression and weak pre-existing hotspot association ...")
gene_stats <- data.table(
  gene = common_genes,
  mean_expr = colMeans(stcos_counts),
  zero_prop = colMeans(stcos_counts == 0)
)
gene_stats[, hotspot_rho := vapply(
  gene,
  function(g) safe_spearman(log1p(stcos_counts[, g]), as.numeric(coords$motif)),
  numeric(1)
)]

candidate_genes <- gene_stats[
  mean_expr >= quantile(mean_expr, 0.35, na.rm = TRUE) &
    mean_expr <= quantile(mean_expr, 0.90, na.rm = TRUE) &
    zero_prop <= 0.85 &
    abs(hotspot_rho) <= 0.08
][order(abs(hotspot_rho), -mean_expr)]

if (nrow(candidate_genes) < 40) {
  candidate_genes <- gene_stats[
    mean_expr >= quantile(mean_expr, 0.25, na.rm = TRUE) &
      zero_prop <= 0.90
  ][order(abs(hotspot_rho), -mean_expr)]
}
if (nrow(candidate_genes) < 10) stop("Too few candidate genes for the controlled-pattern demo.")

n_target <- min(4, nrow(candidate_genes))
target_genes <- candidate_genes$gene[seq_len(n_target)]
target_mean <- mean(candidate_genes$mean_expr[seq_len(n_target)])

background_pool <- candidate_genes[!gene %in% target_genes]
background_pool[, mean_distance := abs(mean_expr - target_mean)]
background_genes <- head(background_pool[order(mean_distance, abs(hotspot_rho))]$gene, min(96, nrow(background_pool)))
selected_genes <- c(target_genes, background_genes)

message(sprintf("Target genes: %s", paste(target_genes, collapse = ", ")))
message(sprintf("Background genes: %d matched genes", length(background_genes)))

base_subset <- stcos_counts[, selected_genes, drop = FALSE]
lib_factor <- rowSums(base_subset) / mean(rowSums(base_subset))
lib_factor[!is.finite(lib_factor) | lib_factor <= 0] <- 1

mu_g <- pmax(colMeans(base_subset), 1e-6)
theta_g <- vapply(seq_len(ncol(base_subset)), function(i) fit_nb_fast(base_subset[, i]), numeric(1))
names(theta_g) <- selected_genes

mu_mat <- outer(lib_factor, mu_g)
dimnames(mu_mat) <- list(common_spots, selected_genes)

fold_change <- 8
mu_mat[coords$motif, target_genes] <- mu_mat[coords$motif, target_genes, drop = FALSE] * fold_change

sigma <- cor(log1p(base_subset), use = "pairwise.complete.obs")
sigma[is.na(sigma)] <- 0
diag(sigma) <- 1
sigma[abs(sigma) < 0.05] <- 0
sigma <- makespd(sigma)
sigma <- cov2cor(sigma)
dimnames(sigma) <- list(selected_genes, selected_genes)

gene_rel <- as.data.table(matrix(1, nrow = length(selected_genes), ncol = 1))
setnames(gene_rel, "CT")

cov_gene <- array(
  NA_real_,
  dim = c(length(selected_genes), length(selected_genes), 1),
  dimnames = list(selected_genes, selected_genes, "CT")
)
cov_gene[, , 1] <- sigma

message("Generating ST-COS controlled-pattern counts for the selected gene panel ...")
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
      motif_contrast = mean(vals[coords$motif]) - mean(vals[!coords$motif]),
      motif_spearman = safe_spearman(vals, as.numeric(coords$motif)),
      mean_inside = mean(vals[coords$motif]),
      mean_outside = mean(vals[!coords$motif])
    )
  }))
}

metric_dt <- rbindlist(list(
  score_matrix(stcos_counts[, selected_genes, drop = FALSE], "ST-COS replay", target_genes, "target"),
  score_matrix(stcos_counts[, selected_genes, drop = FALSE], "ST-COS replay", background_genes, "background"),
  score_matrix(srtsim_counts[, selected_genes, drop = FALSE], "SRTsim replay", target_genes, "target"),
  score_matrix(srtsim_counts[, selected_genes, drop = FALSE], "SRTsim replay", background_genes, "background"),
  score_matrix(stcos_controlled, "ST-COS controlled", target_genes, "target"),
  score_matrix(stcos_controlled, "ST-COS controlled", background_genes, "background")
))

method_levels <- c("ST-COS replay", "SRTsim replay", "ST-COS controlled")
metric_dt[, method := factor(method, levels = method_levels)]
metric_dt[, gene_class := factor(gene_class, levels = c("background", "target"))]

plot_genes <- target_genes[seq_len(min(3, length(target_genes)))]
map_dt <- rbindlist(lapply(plot_genes, function(g) {
  rbindlist(list(
    data.table(coords[, .(spot_id, row, col, motif)], gene = g, method = "ST-COS replay",
               expr = log1p(stcos_counts[, g])),
    data.table(coords[, .(spot_id, row, col, motif)], gene = g, method = "SRTsim replay",
               expr = log1p(srtsim_counts[, g])),
    data.table(coords[, .(spot_id, row, col, motif)], gene = g, method = "ST-COS controlled",
               expr = log1p(stcos_controlled[, g]))
  ))
}))
map_dt[, method := factor(method, levels = method_levels)]
map_dt[, gene := factor(gene, levels = plot_genes)]

metric_summary <- metric_dt[, .(
  n_genes = .N,
  median_contrast = median(motif_contrast, na.rm = TRUE),
  mean_contrast = mean(motif_contrast, na.rm = TRUE),
  median_spearman = median(motif_spearman, na.rm = TRUE),
  mean_spearman = mean(motif_spearman, na.rm = TRUE)
), by = .(method, gene_class)]

fwrite(data.table(
  slide_id = slide_id,
  fold_change = fold_change,
  target_gene = target_genes,
  baseline_mean = gene_stats[match(target_genes, gene), mean_expr],
  baseline_zero_prop = gene_stats[match(target_genes, gene), zero_prop],
  baseline_hotspot_rho = gene_stats[match(target_genes, gene), hotspot_rho]
), file.path(out_dir, "DLPFC_controlled_pattern_target_genes.csv"))
fwrite(metric_dt, file.path(out_dir, "DLPFC_controlled_pattern_gene_metrics.csv"))
fwrite(metric_summary, file.path(out_dir, "DLPFC_controlled_pattern_metric_summary.csv"))
fwrite(map_dt, file.path(out_dir, "DLPFC_controlled_pattern_map_values.csv"))
fwrite(
  data.table(spot_id = rownames(stcos_controlled), as.data.table(stcos_controlled)),
  file.path(out_dir, "DLPFC_controlled_pattern_STCOS_selected_counts.csv")
)

caption_text <- sprintf(
  "%s; %d hotspot spots; x%d mean",
  slide_id, sum(coords$motif), fold_change
)

p_hotspot <- ggplot(coords, aes(col, row)) +
  geom_point(aes(color = motif), size = 0.65, alpha = 0.9) +
  scale_color_manual(values = c("FALSE" = "grey83", "TRUE" = "#c43c39"), labels = c("Background", "User hotspot")) +
  coord_equal() +
  labs(title = "A. User-specified hotspot", subtitle = caption_text, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  )

p_maps <- ggplot(map_dt, aes(col, row, color = expr)) +
  geom_point(size = 0.45) +
  coord_equal() +
  facet_grid(gene ~ method) +
  scale_color_gradientn(colors = c("#f2f2f2", "#97b8cc", "#364f8b", "#7c1d3f"), name = "log1p(count)") +
  labs(title = "B. Target-gene spatial maps") +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    strip.text = element_text(size = 9),
    legend.position = "right"
  )

p_contrast <- ggplot(metric_dt, aes(method, motif_contrast, fill = gene_class)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.70, width = 0.65) +
  geom_point(aes(color = gene_class), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65),
             size = 1.0, alpha = 0.85) +
  scale_fill_manual(values = c(background = "grey80", target = "#d95f02")) +
  scale_color_manual(values = c(background = "grey35", target = "#d95f02")) +
  labs(title = "C. Hotspot contrast", x = NULL, y = "mean log1p inside - outside", fill = NULL, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom"
  )

p_spearman <- ggplot(metric_dt, aes(method, motif_spearman, fill = gene_class)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.70, width = 0.65) +
  geom_point(aes(color = gene_class), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65),
             size = 1.0, alpha = 0.85) +
  scale_fill_manual(values = c(background = "grey80", target = "#1b9e77")) +
  scale_color_manual(values = c(background = "grey35", target = "#1b9e77")) +
  labs(title = "D. Association with user hotspot", x = NULL, y = "Spearman rho", fill = NULL, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom"
  )

top_row <- plot_grid(p_hotspot, p_contrast, p_spearman, ncol = 3, rel_widths = c(0.85, 1.15, 1.15))
combined <- plot_grid(top_row, p_maps, ncol = 1, rel_heights = c(1.05, 2.15))

pdf_out <- file.path(out_dir, sprintf("DLPFC_%s_STCOS_controllability_hotspot_demo.pdf", slide_id))
png_out <- file.path(out_dir, sprintf("DLPFC_%s_STCOS_controllability_hotspot_demo.png", slide_id))
ggsave(pdf_out, combined, width = 13.5, height = 9.8, bg = "white")
ggsave(png_out, combined, width = 13.5, height = 9.8, dpi = 300, bg = "white")

print(metric_summary)
message(sprintf("Saved figure: %s", pdf_out))
message(sprintf("Saved figure: %s", png_out))
