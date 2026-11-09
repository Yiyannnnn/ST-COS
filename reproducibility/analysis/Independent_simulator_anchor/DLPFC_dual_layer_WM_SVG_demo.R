#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

source("./Data/simulator_YZ.R")

set.seed(20260906)

repo_dir <- "."
slide_id <- "151673"
dlpfc_dir <- file.path(repo_dir, "Reference_based", "DLPFC")
out_dir <- file.path(repo_dir, "Independent_simulator_anchor", "DLPFC_dual_layer_WM_SVG_demo")
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

safe_cor <- function(x, y, method = "spearman") {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = method))
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

generate_panel <- function(mu_mat, sigma, theta_g, coords, seed) {
  set.seed(seed)
  selected_genes <- colnames(mu_mat)
  gene_rel <- as.data.table(matrix(1, nrow = length(selected_genes), ncol = 1))
  setnames(gene_rel, "CT")

  cov_gene <- array(
    NA_real_,
    dim = c(length(selected_genes), length(selected_genes), 1),
    dimnames = list(selected_genes, selected_genes, "CT")
  )
  cov_gene[, , 1] <- sigma

  sim <- generate_true_count(
    coord = as.matrix(coords[, .(row, col)]),
    mean_param = list(CT = mu_mat),
    gene_relationship = gene_rel,
    cov_gene = cov_gene,
    lib_size = list(CT = rep(1, nrow(coords))),
    theta_g = list(CT = theta_g)
  )[["CT"]]

  rownames(sim) <- coords$spot_id
  colnames(sim) <- selected_genes
  sim
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
coords[, l45_band := label %in% c("L4", "L5")]
coords[, wm_band := label == "WM"]
coords[, program := fifelse(wm_band, "WM program", fifelse(l45_band, "L4/L5 program", "Other layers"))]
coords[, program := factor(program, levels = c("Other layers", "L4/L5 program", "WM program"))]

message(sprintf("L4/L5 band: %d spots; WM band: %d spots.", sum(coords$l45_band), sum(coords$wm_band)))

cortical_targets <- base::intersect(c("PCP4", "RORB", "NEFL", "NEFM"), common_genes)
wm_targets <- base::intersect(c("MOBP", "MBP", "PLP1"), common_genes)

wm_candidates <- base::intersect(
  c(
    "MOG", "MAG", "CNP", "CLDN11", "TF", "SOX10", "OLIG1", "OLIG2",
    "OPALIN", "ERMN", "SLC44A1", "MAL", "FA2H", "UGT8", "CA2",
    "AQP4", "GFAP", "VCAN", "PDGFRA"
  ),
  common_genes
)
if (length(cortical_targets) < 3) stop("Too few cortical SVG targets found.")
if (length(wm_targets) < 2) stop("Too few direct WM SVG targets found.")

gene_stats <- data.table(
  gene = common_genes,
  mean_expr = colMeans(stcos_counts),
  zero_prop = colMeans(stcos_counts == 0)
)
gene_stats[, l45_rho := vapply(
  gene,
  function(g) safe_cor(log1p(stcos_counts[, g]), as.numeric(coords$l45_band), "spearman"),
  numeric(1)
)]
gene_stats[, wm_rho := vapply(
  gene,
  function(g) safe_cor(log1p(stcos_counts[, g]), as.numeric(coords$wm_band), "spearman"),
  numeric(1)
)]

wm_stats <- gene_stats[gene %in% wm_candidates]
wm_stats[, layer_eta := vapply(
  gene,
  function(g) eta_by_group(stcos_counts[, g], coords$label),
  numeric(1)
)]
wm_stats[, dominant_layer := vapply(
  gene,
  function(g) dominant_layer(stcos_counts[, g], coords$label),
  character(1)
)]

wm_check_genes <- wm_stats[
  !gene %in% wm_targets &
    dominant_layer == "WM" &
    layer_eta >= 0.10 &
    mean_expr >= 0.05 &
    zero_prop <= 0.95
][order(-layer_eta)]$gene

if (length(wm_check_genes) < 6) {
  wm_check_genes <- wm_stats[
    !gene %in% wm_targets &
      mean_expr >= 0.05 &
      zero_prop <= 0.95
  ][order(-layer_eta)]$gene
}
wm_check_genes <- head(wm_check_genes, 12)

all_focus_genes <- unique(c(cortical_targets, wm_targets, wm_check_genes))
target_mean <- mean(gene_stats[gene %in% all_focus_genes, mean_expr])
background_pool <- gene_stats[
  !gene %in% all_focus_genes &
    mean_expr >= 0.05 &
    zero_prop <= 0.90 &
    abs(l45_rho) <= 0.08 &
    abs(wm_rho) <= 0.08
]
background_pool[, mean_distance := abs(mean_expr - target_mean)]
background_genes <- head(background_pool[order(mean_distance, abs(l45_rho) + abs(wm_rho))]$gene, min(80, nrow(background_pool)))
selected_genes <- unique(c(cortical_targets, wm_targets, wm_check_genes, background_genes))

message(sprintf("Cortical SVG targets: %s", paste(cortical_targets, collapse = ", ")))
message(sprintf("Direct WM SVG targets: %s", paste(wm_targets, collapse = ", ")))
message(sprintf("Other WM SVGs checked: %s", paste(wm_check_genes, collapse = ", ")))
message(sprintf("Background genes: %d matched genes", length(background_genes)))

base_subset <- stcos_counts[, selected_genes, drop = FALSE]
message("Building layer-aware ST-COS mean surfaces ...")
mu_base <- matrix(NA_real_, nrow = nrow(base_subset), ncol = ncol(base_subset),
                  dimnames = list(common_spots, selected_genes))
for (g in selected_genes) {
  layer_means <- tapply(base_subset[, g], coords$label, mean)
  layer_means[is.na(layer_means)] <- mean(base_subset[, g])
  mu_base[, g] <- as.numeric(layer_means[as.character(coords$label)])
}
mu_base[!is.finite(mu_base)] <- 0
mu_base <- pmax(mu_base, 1e-6)

fold_l45 <- 5
fold_wm <- 5
mu_control <- mu_base
mu_control[coords$l45_band, cortical_targets] <- mu_control[coords$l45_band, cortical_targets, drop = FALSE] * fold_l45
mu_control[coords$wm_band, wm_targets] <- mu_control[coords$wm_band, wm_targets, drop = FALSE] * fold_wm

theta_g <- vapply(seq_len(ncol(base_subset)), function(i) fit_nb_fast(base_subset[, i]), numeric(1))
names(theta_g) <- selected_genes

sigma <- cor(log1p(base_subset), use = "pairwise.complete.obs")
sigma[is.na(sigma)] <- 0
diag(sigma) <- 1
sigma[abs(sigma) < 0.05] <- 0
sigma <- cov2cor(makespd(sigma))
dimnames(sigma) <- list(selected_genes, selected_genes)

message("Generating paired ST-COS baseline and dual controlled panels ...")
paired_seed <- 202609061
stcos_param_baseline <- generate_panel(mu_base, sigma, theta_g, coords, paired_seed)
stcos_dual_controlled <- generate_panel(mu_control, sigma, theta_g, coords, paired_seed)

score_region <- function(mat, method_name, genes, gene_class, region_mask, region_name) {
  rbindlist(lapply(genes, function(g) {
    vals <- log1p(mat[, g])
    data.table(
      method = method_name,
      gene_class = gene_class,
      region = region_name,
      gene = g,
      contrast = mean(vals[region_mask]) - mean(vals[!region_mask]),
      spearman = safe_cor(vals, as.numeric(region_mask), "spearman"),
      mean_inside = mean(vals[region_mask]),
      mean_outside = mean(vals[!region_mask])
    )
  }))
}

method_mats <- list(
  "ST-COS replay" = stcos_counts[, selected_genes, drop = FALSE],
  "SRTsim replay" = srtsim_counts[, selected_genes, drop = FALSE],
  "ST-COS parameter baseline" = stcos_param_baseline,
  "ST-COS dual controlled" = stcos_dual_controlled
)

score_all <- function(region_mask, region_name) {
  rbindlist(lapply(names(method_mats), function(m) {
    rbindlist(list(
      score_region(method_mats[[m]], m, cortical_targets, "L4/L5 targets", region_mask, region_name),
      score_region(method_mats[[m]], m, wm_targets, "direct WM targets", region_mask, region_name),
      score_region(method_mats[[m]], m, wm_check_genes, "other WM SVGs", region_mask, region_name),
      score_region(method_mats[[m]], m, background_genes, "background", region_mask, region_name)
    ))
  }))
}

metric_dt <- rbindlist(list(
  score_all(coords$l45_band, "L4/L5"),
  score_all(coords$wm_band, "WM")
))

method_levels <- names(method_mats)
class_levels <- c("background", "other WM SVGs", "L4/L5 targets", "direct WM targets")
metric_dt[, method := factor(method, levels = method_levels)]
metric_dt[, gene_class := factor(gene_class, levels = class_levels)]

baseline_key <- metric_dt[method == "ST-COS parameter baseline" & region == "WM",
                          .(gene, baseline_contrast = contrast, baseline_spearman = spearman)]
control_key <- metric_dt[method == "ST-COS dual controlled" & region == "WM",
                         .(gene, controlled_contrast = contrast, controlled_spearman = spearman)]
spillover_dt <- merge(baseline_key, control_key, by = "gene")
spillover_dt <- merge(
  spillover_dt,
  unique(metric_dt[region == "WM", .(gene, gene_class)]),
  by = "gene",
  all.x = TRUE
)
spillover_dt[, `:=`(
  delta_contrast = controlled_contrast - baseline_contrast,
  delta_spearman = controlled_spearman - baseline_spearman
)]
spillover_dt[, gene_class := factor(gene_class, levels = class_levels)]

pair_cors <- rbindlist(lapply(names(method_mats), function(m) {
  rbindlist(lapply(wm_targets, function(tg) {
    rbindlist(lapply(wm_check_genes, function(pg) {
      data.table(
        method = m,
        target_gene = tg,
        partner_gene = pg,
        pearson_log1p = safe_cor(log1p(method_mats[[m]][, tg]), log1p(method_mats[[m]][, pg]), "pearson"),
        spearman_log1p = safe_cor(log1p(method_mats[[m]][, tg]), log1p(method_mats[[m]][, pg]), "spearman")
      )
    }))
  }))
}))
pair_cors[, method := factor(method, levels = method_levels)]

metric_summary <- metric_dt[, .(
  n_genes = .N,
  median_contrast = median(contrast, na.rm = TRUE),
  mean_contrast = mean(contrast, na.rm = TRUE),
  median_spearman = median(spearman, na.rm = TRUE),
  mean_spearman = mean(spearman, na.rm = TRUE)
), by = .(region, method, gene_class)]

fwrite(data.table(
  slide_id = slide_id,
  gene = c(cortical_targets, wm_targets, wm_check_genes),
  role = c(
    rep("L4/L5 target", length(cortical_targets)),
    rep("direct WM target", length(wm_targets)),
    rep("other WM SVG checked", length(wm_check_genes))
  )
), file.path(out_dir, "DLPFC_dual_layer_gene_roles.csv"))
fwrite(metric_dt, file.path(out_dir, "DLPFC_dual_layer_region_metrics.csv"))
fwrite(metric_summary, file.path(out_dir, "DLPFC_dual_layer_metric_summary.csv"))
fwrite(spillover_dt, file.path(out_dir, "DLPFC_dual_layer_WM_spillover_check.csv"))
fwrite(pair_cors, file.path(out_dir, "DLPFC_dual_layer_WM_pair_correlations.csv"))
fwrite(
  data.table(spot_id = rownames(stcos_dual_controlled), as.data.table(stcos_dual_controlled)),
  file.path(out_dir, "DLPFC_dual_layer_STCOS_selected_counts.csv")
)

wm_plot_genes <- wm_targets
map_dt <- rbindlist(lapply(wm_plot_genes, function(g) {
  rbindlist(list(
    data.table(coords[, .(spot_id, row, col, label, program)], gene = g, method = "ST-COS replay",
               expr = log1p(stcos_counts[, g])),
    data.table(coords[, .(spot_id, row, col, label, program)], gene = g, method = "SRTsim replay",
               expr = log1p(srtsim_counts[, g])),
    data.table(coords[, .(spot_id, row, col, label, program)], gene = g, method = "ST-COS dual controlled",
               expr = log1p(stcos_dual_controlled[, g]))
  ))
}))
map_dt[, method := factor(method, levels = c("ST-COS replay", "SRTsim replay", "ST-COS dual controlled"))]
map_dt[, gene := factor(gene, levels = wm_plot_genes)]
expr_cap <- quantile(map_dt$expr, 0.99, na.rm = TRUE)
map_dt[, expr_plot := pmin(expr, expr_cap)]
fwrite(map_dt, file.path(out_dir, "DLPFC_dual_layer_WM_map_values.csv"))

plot_metric <- metric_dt[
  region == "WM" &
    method %in% c("ST-COS replay", "SRTsim replay", "ST-COS dual controlled") &
    gene_class %in% c("background", "other WM SVGs", "direct WM targets")
]
plot_metric[, method := factor(as.character(method), levels = c("ST-COS replay", "SRTsim replay", "ST-COS dual controlled"))]

p_programs <- ggplot(coords, aes(col, row)) +
  geom_point(aes(color = program), size = 0.58, alpha = 0.9) +
  scale_color_manual(values = c("Other layers" = "grey83", "L4/L5 program" = "#c43c39", "WM program" = "#167a72")) +
  coord_equal() +
  labs(
    title = "A. Two user-defined programs",
    subtitle = sprintf("%s; L4/L5 x%d, WM x%d", slide_id, fold_l45, fold_wm),
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  )

p_maps <- ggplot(map_dt, aes(col, row, color = expr_plot)) +
  geom_point(size = 0.40) +
  coord_equal() +
  facet_grid(gene ~ method) +
  scale_color_gradientn(colors = c("#f2f2f2", "#9fc2d4", "#385991", "#8c2046"), name = "log1p(count)") +
  labs(title = "B. Direct WM target SVGs") +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    strip.text = element_text(size = 9),
    legend.position = "right"
  )

p_wm_contrast <- ggplot(plot_metric, aes(method, contrast, fill = gene_class)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.72, width = 0.66) +
  geom_point(aes(color = gene_class), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.66),
             size = 0.9, alpha = 0.78) +
  scale_fill_manual(values = c(background = "grey80", `other WM SVGs` = "#6baed6", `direct WM targets` = "#d95f02")) +
  scale_color_manual(values = c(background = "grey35", `other WM SVGs` = "#2b6f9e", `direct WM targets` = "#d95f02")) +
  labs(title = "C. WM contrast", x = NULL, y = "mean log1p in WM - outside", fill = NULL, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom"
  )

spill_plot <- spillover_dt[gene_class %in% c("background", "other WM SVGs", "direct WM targets")]
p_spill <- ggplot(spill_plot, aes(gene_class, delta_contrast, fill = gene_class)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.72, width = 0.62) +
  geom_point(aes(color = gene_class), position = position_jitter(width = 0.08), size = 1.05, alpha = 0.82) +
  scale_fill_manual(values = c(background = "grey80", `other WM SVGs` = "#6baed6", `direct WM targets` = "#d95f02")) +
  scale_color_manual(values = c(background = "grey35", `other WM SVGs` = "#2b6f9e", `direct WM targets` = "#d95f02")) +
  labs(
    title = "D. Extra WM signal from edit",
    subtitle = "dual controlled minus paired baseline",
    x = NULL,
    y = "delta WM contrast",
    fill = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 18, hjust = 1),
    legend.position = "none"
  )

top_row <- plot_grid(p_programs, p_wm_contrast, p_spill, ncol = 3, rel_widths = c(0.9, 1.25, 1.05))
combined <- plot_grid(top_row, p_maps, ncol = 1, rel_heights = c(1.05, 2.05))

pdf_out <- file.path(out_dir, sprintf("DLPFC_%s_dual_layer_WM_SVG_demo.pdf", slide_id))
png_out <- file.path(out_dir, sprintf("DLPFC_%s_dual_layer_WM_SVG_demo.png", slide_id))
ggsave(pdf_out, combined, width = 13.5, height = 10.0, bg = "white")
ggsave(png_out, combined, width = 13.5, height = 10.0, dpi = 300, bg = "white")

print(metric_summary[region == "WM"])
print(spillover_dt[, .(
  n_genes = .N,
  median_delta_contrast = median(delta_contrast, na.rm = TRUE),
  mean_delta_contrast = mean(delta_contrast, na.rm = TRUE),
  median_delta_spearman = median(delta_spearman, na.rm = TRUE),
  mean_delta_spearman = mean(delta_spearman, na.rm = TRUE)
), by = gene_class])
print(pair_cors[, .(
  n_pairs = .N,
  median_pearson = median(pearson_log1p, na.rm = TRUE),
  median_spearman = median(spearman_log1p, na.rm = TRUE)
), by = method])
message(sprintf("Saved figure: %s", pdf_out))
message(sprintf("Saved figure: %s", png_out))
