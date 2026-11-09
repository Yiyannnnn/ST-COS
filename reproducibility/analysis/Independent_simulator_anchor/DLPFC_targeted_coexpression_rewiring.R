#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(Matrix)
})

repo_dir <- "."
args <- commandArgs(trailingOnly = TRUE)
slide_id <- if (length(args)) args[1L] else "151673"
dlpfc_dir <- file.path(repo_dir, "Reference_based", "DLPFC")
out_root <- file.path(
  repo_dir, "Independent_simulator_anchor",
  "DLPFC_targeted_coexpression_rewiring"
)
out_dir <- if (length(args)) file.path(out_root, "slides", slide_id) else out_root
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

slide_seed_offset <- as.integer(slide_id) %% 100000L
set.seed(20260908L + slide_seed_offset)
n_reps <- as.integer(Sys.getenv("STCOS_REWIRE_N_REPS", "30"))
edge_thresholds <- c("Raw count" = 0.65, "Layer-adjusted" = 0.55)
requested_latent_floor <- 0.70

read_spot_gene <- function(path) {
  dt <- fread(path, check.names = FALSE)
  ids <- as.character(dt[[1L]])
  mat <- as.matrix(dt[, -1L, with = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- ids
  colnames(mat) <- names(dt)[-1L]
  mat
}

read_gene_spot <- function(path) {
  dt <- fread(path, check.names = FALSE)
  genes <- as.character(dt[[1L]])
  mat <- as.matrix(dt[, -1L, with = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  mat <- t(mat)
  colnames(mat) <- genes
  mat
}

safe_cor <- function(x, y, method) {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = method, use = "pairwise.complete.obs"))
}

fit_theta <- function(x, cap = 1e6) {
  mu <- mean(x)
  v <- var(x)
  if (!is.finite(mu) || !is.finite(v) || mu <= 0 || v <= mu) return(cap)
  min(max(mu^2 / (v - mu), 1e-4), cap)
}

repair_correlation <- function(x, eigen_floor = 1e-8) {
  x <- (x + t(x)) / 2
  diag(x) <- 1
  ee <- eigen(x, symmetric = TRUE)
  ee$values <- ee$values + max(0, eigen_floor - min(ee$values))
  ans <- ee$vectors %*% (ee$values * t(ee$vectors))
  ans <- cov2cor(ans)
  ans <- (ans + t(ans)) / 2
  diag(ans) <- 1
  dimnames(ans) <- dimnames(x)
  ans
}

simulate_copula_counts <- function(mu_mat, theta, rho, seed) {
  set.seed(seed)
  n <- nrow(mu_mat)
  g <- ncol(mu_mat)
  lower_chol <- t(chol(rho))
  latent <- lower_chol %*% matrix(rnorm(g * n), nrow = g)
  uniforms <- pmin(pmax(pnorm(latent), 1e-12), 1 - 1e-12)
  out <- matrix(0, nrow = n, ncol = g, dimnames = dimnames(mu_mat))
  for (idx in seq_len(g)) {
    out[, idx] <- qnbinom(
      uniforms[idx, ], size = theta[idx], mu = mu_mat[, idx]
    )
  }
  out
}

cor_matrix <- function(mat, method, groups = NULL) {
  work <- mat
  if (!is.null(groups)) {
    for (idx in seq_len(ncol(work))) {
      work[, idx] <- work[, idx] - ave(work[, idx], groups, FUN = mean)
    }
  }
  ans <- cor(work, method = method, use = "pairwise.complete.obs")
  ans[!is.finite(ans)] <- 0
  diag(ans) <- 1
  ans
}

average_precision <- function(labels, scores) {
  ok <- is.finite(scores) & !is.na(labels)
  labels <- as.integer(labels[ok])
  scores <- scores[ok]
  n_pos <- sum(labels == 1L)
  if (n_pos == 0L) return(NA_real_)
  ord <- order(scores, decreasing = TRUE)
  labels <- labels[ord]
  tp <- cumsum(labels == 1L)
  fp <- cumsum(labels == 0L)
  precision <- tp / (tp + fp)
  sum(precision[labels == 1L]) / n_pos
}

matrix_to_edges <- function(mat, target_genes, method, replicate,
                            correlation_type, dependence_scale,
                            edge_threshold) {
  idx <- which(upper.tri(mat), arr.ind = TRUE)
  data.table(
    method = method,
    replicate = replicate,
    correlation_type = correlation_type,
    dependence_scale = dependence_scale,
    edge_threshold = edge_threshold,
    gene_1 = rownames(mat)[idx[, 1L]],
    gene_2 = colnames(mat)[idx[, 2L]],
    correlation = mat[idx]
  )[, requested_edge := gene_1 %in% target_genes & gene_2 %in% target_genes]
}

message("Reading DLPFC replay outputs ...")
coords <- fread(
  file.path(dlpfc_dir, sprintf("%s_simulated_coordinates.csv", slide_id)),
  check.names = FALSE
)
setnames(coords, names(coords)[1L], "spot_id")
stcos_replay <- read_spot_gene(file.path(
  dlpfc_dir, sprintf("%s_STCOS_simulated_counts.csv", slide_id)
))
srtsim_replay <- read_gene_spot(file.path(
  dlpfc_dir, sprintf("%s_SRTsim_counts.csv", slide_id)
))

common_spots <- intersect(coords$spot_id, rownames(stcos_replay))
common_genes <- intersect(colnames(stcos_replay), colnames(srtsim_replay))
coords <- coords[match(common_spots, spot_id)]
stcos_replay <- stcos_replay[common_spots, common_genes, drop = FALSE]
srtsim_replay <- srtsim_replay[, common_genes, drop = FALSE]

target_candidates <- c("MOBP", "CNP", "MAG", "MOG", "CLDN11", "ERMN")
target_genes <- intersect(target_candidates, common_genes)
if (length(target_genes) < 4L) {
  stop("At least four prespecified WM target genes are required.")
}

gene_stats <- data.table(
  gene = common_genes,
  mean_expr = colMeans(stcos_replay),
  zero_prop = colMeans(stcos_replay == 0)
)
target_mean <- median(gene_stats[gene %in% target_genes, mean_expr])
background_genes <- gene_stats[
  !gene %in% target_genes & mean_expr >= 0.05 & zero_prop <= 0.92
][order(abs(mean_expr - target_mean))][1:min(.N, 60L), gene]
selected_genes <- c(target_genes, background_genes)

message(sprintf(
  "Target block: %s; background genes: %d",
  paste(target_genes, collapse = ", "), length(background_genes)
))

reference_subset <- stcos_replay[, selected_genes, drop = FALSE]
srtsim_subset <- srtsim_replay[, selected_genes, drop = FALSE]

layer <- factor(coords$label, levels = c("L1", "L2", "L3", "L4", "L5", "L6", "WM"))
mu_mat <- matrix(
  NA_real_, nrow = nrow(reference_subset), ncol = ncol(reference_subset),
  dimnames = dimnames(reference_subset)
)
for (gene in selected_genes) {
  layer_means <- tapply(reference_subset[, gene], layer, mean)
  layer_means[!is.finite(layer_means)] <- mean(reference_subset[, gene])
  mu_mat[, gene] <- layer_means[as.character(layer)]
}
mu_mat[!is.finite(mu_mat)] <- 0
mu_mat <- pmax(mu_mat, 1e-8)

theta <- vapply(selected_genes, function(g) fit_theta(reference_subset[, g]), numeric(1))
rho_base <- cor(log1p(reference_subset), use = "pairwise.complete.obs")
rho_base[!is.finite(rho_base)] <- 0
diag(rho_base) <- 1
rho_base[abs(rho_base) < 0.05] <- 0
rho_base <- repair_correlation(rho_base)

target_block <- rho_base[target_genes, target_genes, drop = FALSE]
minimum_target_rho <- min(target_block[upper.tri(target_block)])
delta <- max(0, (requested_latent_floor - minimum_target_rho) /
  (1 - requested_latent_floor))
target_vector <- as.numeric(selected_genes %in% target_genes)
rho_rewired <- cov2cor(rho_base + delta * tcrossprod(target_vector))
rho_rewired <- repair_correlation(rho_rewired)

non_target <- setdiff(selected_genes, target_genes)
latent_non_target_change <- max(abs(
  rho_rewired[non_target, non_target] - rho_base[non_target, non_target]
))
if (latent_non_target_change > 1e-8) {
  stop("The rewiring changed the non-target latent block.")
}

requested_edges <- as.data.table(which(
  upper.tri(rho_rewired[target_genes, target_genes]), arr.ind = TRUE
))
setnames(requested_edges, c("row", "col"), c("row_index", "col_index"))
requested_edges[, `:=`(
  gene_1 = target_genes[row_index],
  gene_2 = target_genes[col_index]
)]
requested_edges[, `:=`(
  baseline_latent_rho = rho_base[cbind(
    match(gene_1, selected_genes), match(gene_2, selected_genes)
  )],
  requested_latent_rho = rho_rewired[cbind(
    match(gene_1, selected_genes), match(gene_2, selected_genes)
  )]
)]
requested_edges[, c("row_index", "col_index") := NULL]
fwrite(requested_edges, file.path(out_dir, "requested_WM_network_edges.csv"))

message(sprintf(
  "PSD block update delta %.4f; target latent range %.3f--%.3f",
  delta,
  min(requested_edges$requested_latent_rho),
  max(requested_edges$requested_latent_rho)
))

edge_rows <- list()
preservation_rows <- list()
matrix_accumulator <- list(Pearson = list(), Spearman = list())

for (replicate in seq_len(n_reps)) {
  seed <- 202609080L + 100L * slide_seed_offset + replicate
  baseline <- simulate_copula_counts(mu_mat, theta, rho_base, seed)
  rewired <- simulate_copula_counts(mu_mat, theta, rho_rewired, seed)

  for (dependence_scale in names(edge_thresholds)) {
    adjustment <- if (dependence_scale == "Layer-adjusted") layer else NULL
    for (correlation_type in c("Pearson", "Spearman")) {
      method <- tolower(correlation_type)
      c_base <- cor_matrix(baseline, method, adjustment)
      c_rewired <- cor_matrix(rewired, method, adjustment)
      edge_rows[[length(edge_rows) + 1L]] <- matrix_to_edges(
        c_base, target_genes, "ST-COS parameter baseline", replicate,
        correlation_type, dependence_scale, edge_thresholds[[dependence_scale]]
      )
      edge_rows[[length(edge_rows) + 1L]] <- matrix_to_edges(
        c_rewired, target_genes, "ST-COS targeted rewire", replicate,
        correlation_type, dependence_scale, edge_thresholds[[dependence_scale]]
      )
      if (dependence_scale == "Layer-adjusted") {
        matrix_accumulator[[correlation_type]][[paste0("baseline_", replicate)]] <- c_base
        matrix_accumulator[[correlation_type]][[paste0("rewired_", replicate)]] <- c_rewired
      }

      idx_non_target <- which(upper.tri(c_base) &
        outer(selected_genes %in% non_target, selected_genes %in% non_target, `&`))
      idx_cross <- which(upper.tri(c_base) &
        outer(selected_genes %in% target_genes, selected_genes %in% non_target, `|`) &
        !outer(selected_genes %in% target_genes, selected_genes %in% target_genes, `&`))
      preservation_rows[[length(preservation_rows) + 1L]] <- data.table(
        replicate = replicate,
        dependence_scale = dependence_scale,
        correlation_type = correlation_type,
        unedited_edge_rmse = sqrt(mean((c_rewired[idx_non_target] - c_base[idx_non_target])^2)),
        cross_edge_rmse = sqrt(mean((c_rewired[idx_cross] - c_base[idx_cross])^2))
      )
    }
  }

  baseline_mean <- colMeans(baseline)
  rewired_mean <- colMeans(rewired)
  baseline_var <- apply(baseline, 2L, var)
  rewired_var <- apply(rewired, 2L, var)
  baseline_theta <- vapply(seq_along(selected_genes), function(i) {
    fit_theta(baseline[, i])
  }, numeric(1))
  rewired_theta <- vapply(seq_along(selected_genes), function(i) {
    fit_theta(rewired[, i])
  }, numeric(1))

  preservation_rows[[length(preservation_rows) + 1L]] <- data.table(
    replicate = replicate,
    dependence_scale = "Marginals",
    correlation_type = "Marginals",
    mean_abs_log2_ratio = mean(abs(log2((rewired_mean + 0.1) /
      (baseline_mean + 0.1)))),
    variance_abs_log2_ratio = mean(abs(log2((rewired_var + 0.1) /
      (baseline_var + 0.1)))),
    dispersion_abs_log2_ratio = median(abs(log2((rewired_theta + 0.1) /
      (baseline_theta + 0.1))))
  )

  if (replicate == 1L) {
    fwrite(as.data.table(baseline, keep.rownames = "spot_id"),
           file.path(out_dir, "STCOS_parameter_baseline_rep1_counts.csv"))
    fwrite(as.data.table(rewired, keep.rownames = "spot_id"),
           file.path(out_dir, "STCOS_targeted_rewire_rep1_counts.csv"))
  }
}

for (dependence_scale in names(edge_thresholds)) {
  adjustment <- if (dependence_scale == "Layer-adjusted") layer else NULL
  for (correlation_type in c("Pearson", "Spearman")) {
    method <- tolower(correlation_type)
    edge_rows[[length(edge_rows) + 1L]] <- matrix_to_edges(
      cor_matrix(reference_subset, method, adjustment), target_genes,
      "ST-COS replay", 0L, correlation_type, dependence_scale,
      edge_thresholds[[dependence_scale]]
    )
    edge_rows[[length(edge_rows) + 1L]] <- matrix_to_edges(
      cor_matrix(srtsim_subset, method, adjustment), target_genes,
      "SRTsim replay", 0L, correlation_type, dependence_scale,
      edge_thresholds[[dependence_scale]]
    )
  }
}

edge_dt <- rbindlist(edge_rows, fill = TRUE)
preservation_dt <- rbindlist(preservation_rows, fill = TRUE)
if (identical(Sys.getenv("STCOS_REWIRE_CHECKPOINT", "0"), "1")) {
  fwrite(edge_dt, file.path(out_dir, "coexpression_edges_checkpoint.csv"))
  fwrite(preservation_dt, file.path(out_dir, "preservation_checkpoint.csv"))
}

network_metrics <- edge_dt[, {
  labels <- as.integer(requested_edge)
  scores <- correlation
  threshold <- unique(edge_threshold)
  calls <- correlation >= threshold
  tp <- sum(labels == 1L & calls, na.rm = TRUE)
  fp <- sum(labels == 0L & calls, na.rm = TRUE)
  fn <- sum(labels == 1L & !calls, na.rm = TRUE)
  .(
    mean_target_correlation = mean(correlation[requested_edge], na.rm = TRUE),
    median_target_correlation = median(correlation[requested_edge], na.rm = TRUE),
    requested_edge_recall = mean(correlation[requested_edge] >= threshold, na.rm = TRUE),
    requested_edge_precision = ifelse(tp + fp > 0L, tp / (tp + fp), 0),
    requested_edge_F1 = ifelse(2L * tp + fp + fn > 0L,
                               2L * tp / (2L * tp + fp + fn), 0),
    non_target_strong_edge_rate = mean(correlation[!requested_edge] >= threshold, na.rm = TRUE),
    requested_edge_AUPRC = average_precision(labels, scores)
  )
}, by = .(method, replicate, correlation_type, dependence_scale, edge_threshold)]

paired_edges <- merge(
  edge_dt[method == "ST-COS parameter baseline",
          .(replicate, correlation_type, dependence_scale, gene_1, gene_2,
            baseline_correlation = correlation, requested_edge)],
  edge_dt[method == "ST-COS targeted rewire",
          .(replicate, correlation_type, dependence_scale, gene_1, gene_2,
            rewired_correlation = correlation)],
  by = c("replicate", "correlation_type", "dependence_scale", "gene_1", "gene_2")
)
paired_edges[, correlation_change := rewired_correlation - baseline_correlation]

change_metrics <- paired_edges[, .(
  mean_change = mean(correlation_change),
  median_change = median(correlation_change),
  q025 = quantile(correlation_change, 0.025),
  q975 = quantile(correlation_change, 0.975),
  fraction_increased_0.10 = mean(correlation_change >= 0.10)
), by = .(replicate, correlation_type, dependence_scale, requested_edge)]

network_summary <- network_metrics[, .(
  mean_target_correlation = mean(mean_target_correlation),
  ci_low_target_correlation = quantile(mean_target_correlation, 0.025),
  ci_high_target_correlation = quantile(mean_target_correlation, 0.975),
  mean_requested_edge_recall = mean(requested_edge_recall),
  mean_requested_edge_precision = mean(requested_edge_precision),
  mean_requested_edge_F1 = mean(requested_edge_F1),
  mean_requested_edge_AUPRC = mean(requested_edge_AUPRC),
  n_replicates = .N
), by = .(method, correlation_type, dependence_scale, edge_threshold)]

threshold_grid <- seq(0.20, 0.75, by = 0.05)
threshold_curve <- rbindlist(lapply(threshold_grid, function(threshold) {
  edge_dt[dependence_scale == "Layer-adjusted", {
    labels <- as.integer(requested_edge)
    calls <- correlation >= threshold
    tp <- sum(labels == 1L & calls, na.rm = TRUE)
    fp <- sum(labels == 0L & calls, na.rm = TRUE)
    fn <- sum(labels == 1L & !calls, na.rm = TRUE)
    .(
      threshold = threshold,
      precision = ifelse(tp + fp > 0L, tp / (tp + fp), 0),
      recall = ifelse(tp + fn > 0L, tp / (tp + fn), 0),
      F1 = ifelse(2L * tp + fp + fn > 0L, 2L * tp / (2L * tp + fp + fn), 0)
    )
  }, by = .(method, replicate, correlation_type)]
}))

paired_tests <- rbindlist(lapply(names(edge_thresholds), function(scale) {
  rbindlist(lapply(c("Pearson", "Spearman"), function(ct) {
  paired <- dcast(
    network_metrics[
      correlation_type == ct & dependence_scale == scale &
        method %in% c("ST-COS parameter baseline", "ST-COS targeted rewire")
    ],
    replicate ~ method,
    value.var = c("mean_target_correlation", "requested_edge_F1")
  )
  rbindlist(lapply(c("mean_target_correlation", "requested_edge_F1"), function(metric) {
    baseline <- paired[[paste0(metric, "_ST-COS parameter baseline")]]
    rewired <- paired[[paste0(metric, "_ST-COS targeted rewire")]]
    if (length(rewired) >= 2L) {
      tt <- tryCatch(t.test(rewired, baseline, paired = TRUE),
                     error = function(e) NULL)
      wt <- tryCatch(wilcox.test(rewired, baseline, paired = TRUE, exact = FALSE),
                     error = function(e) NULL)
      ci <- if (is.null(tt)) c(mean(rewired - baseline), mean(rewired - baseline)) else
        tt$conf.int
      t_p <- if (is.null(tt)) NA_real_ else tt$p.value
      wilcox_p <- if (is.null(wt)) NA_real_ else wt$p.value
    } else {
      ci <- c(NA_real_, NA_real_)
      t_p <- NA_real_
      wilcox_p <- NA_real_
    }
    data.table(
      correlation_type = ct,
      dependence_scale = scale,
      metric = metric,
      mean_paired_change = mean(rewired - baseline),
      ci_low = ci[1L],
      ci_high = ci[2L],
      paired_t_p_value = t_p,
      paired_wilcoxon_p_value = wilcox_p,
      n_paired_replicates = length(rewired)
    )
  }))
  }))
}))

preservation_summary <- melt(
  preservation_dt,
  id.vars = c("replicate", "dependence_scale", "correlation_type"),
  measure.vars = intersect(
    c("unedited_edge_rmse", "cross_edge_rmse", "mean_abs_log2_ratio",
      "variance_abs_log2_ratio", "dispersion_abs_log2_ratio"),
    names(preservation_dt)
  ),
  variable.name = "metric", value.name = "value", na.rm = TRUE
)[, .(
  mean = mean(value),
  ci_low = quantile(value, 0.025),
  ci_high = quantile(value, 0.975)
), by = .(dependence_scale, correlation_type, metric)]

fwrite(edge_dt, file.path(out_dir, "coexpression_edges_by_method.csv"))
fwrite(network_metrics, file.path(out_dir, "network_metrics_by_replicate.csv"))
fwrite(network_summary, file.path(out_dir, "network_metric_summary.csv"))
fwrite(threshold_curve, file.path(out_dir, "network_threshold_curve.csv"))
fwrite(paired_tests, file.path(out_dir, "paired_significance_tests.csv"))
fwrite(paired_edges, file.path(out_dir, "paired_edge_changes.csv"))
fwrite(change_metrics, file.path(out_dir, "edge_change_metrics_by_replicate.csv"))
fwrite(preservation_dt, file.path(out_dir, "marginal_and_edge_preservation_by_replicate.csv"))
fwrite(preservation_summary, file.path(out_dir, "marginal_and_edge_preservation_summary.csv"))
fwrite(data.table(
  slide_id = slide_id,
  n_spots = nrow(reference_subset),
  n_genes = length(selected_genes),
  n_target_genes = length(target_genes),
  target_genes = paste(target_genes, collapse = ";"),
  requested_latent_floor = requested_latent_floor,
  psd_update_delta = delta,
  maximum_non_target_latent_change = latent_non_target_change,
  n_replicates = n_reps
), file.path(out_dir, "experiment_design.csv"))

mean_matrix <- function(items, prefix) {
  selected <- items[startsWith(names(items), prefix)]
  Reduce(`+`, selected) / length(selected)
}
heatmap_mats <- list(
  "Baseline latent" = rho_base[target_genes, target_genes],
  "Requested latent" = rho_rewired[target_genes, target_genes],
  "ST-COS baseline residual" = mean_matrix(matrix_accumulator$Pearson, "baseline_")[target_genes, target_genes],
  "ST-COS rewired residual" = mean_matrix(matrix_accumulator$Pearson, "rewired_")[target_genes, target_genes],
  "SRTsim replay residual" = cor_matrix(srtsim_subset, "pearson", layer)[target_genes, target_genes]
)
heat_dt <- rbindlist(lapply(names(heatmap_mats), function(label) {
  as.data.table(as.table(heatmap_mats[[label]]))[, panel := label]
}))
setnames(heat_dt, c("gene_1", "gene_2", "correlation", "panel"))
heat_dt[, panel := factor(panel, levels = names(heatmap_mats))]

method_levels <- c(
  "ST-COS replay", "SRTsim replay", "ST-COS parameter baseline",
  "ST-COS targeted rewire"
)
method_colors <- c(
  "ST-COS replay" = "#6C757D",
  "SRTsim replay" = "#D55E00",
  "ST-COS parameter baseline" = "#0072B2",
  "ST-COS targeted rewire" = "#009E73"
)
method_labels <- c(
  "ST-COS replay" = "ST-COS replay",
  "SRTsim replay" = "SRTsim replay",
  "ST-COS parameter baseline" = "ST-COS baseline",
  "ST-COS targeted rewire" = "ST-COS rewired"
)

p_heat <- ggplot(heat_dt, aes(gene_2, gene_1, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  coord_equal() +
  labs(title = "A  Target WM co-expression block (layer-adjusted count correlations)",
       x = NULL, y = NULL, fill = "Correlation") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(), strip.background = element_rect(fill = "#E8EEF0"))

plot_network <- network_metrics[dependence_scale == "Layer-adjusted"]
plot_network[, method := factor(method, levels = method_levels)]
p_target <- ggplot(
  plot_network, aes(method, mean_target_correlation, color = method)
) +
  geom_hline(yintercept = edge_thresholds[["Layer-adjusted"]], linetype = "dotted", color = "#555555") +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.18) +
  geom_jitter(width = 0.08, height = 0, size = 1.2, alpha = 0.6) +
  facet_wrap(~correlation_type) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_x_discrete(labels = method_labels) +
  labs(title = "B  On-target layer-adjusted dependence", x = NULL,
       y = "Mean target-pair correlation") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")

p_recovery <- ggplot(
  plot_network, aes(method, requested_edge_F1, color = method)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.18) +
  geom_jitter(width = 0.08, height = 0, size = 1.2, alpha = 0.6) +
  facet_wrap(~correlation_type) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_x_discrete(labels = method_labels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "C  Requested-edge recovery", x = NULL,
       y = sprintf("Network F1 at residual r >= %.2f", edge_thresholds[["Layer-adjusted"]])) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")

change_plot <- paired_edges[dependence_scale == "Layer-adjusted", .(
  mean_change = mean(correlation_change)
), by = .(replicate, correlation_type, requested_edge)]
change_plot[, edge_class := ifelse(requested_edge, "Edited WM edges", "Unedited edges")]
p_change <- ggplot(
  change_plot, aes(edge_class, mean_change, fill = edge_class)
) +
  geom_hline(yintercept = 0, color = "#555555") +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 1.1, alpha = 0.55) +
  facet_wrap(~correlation_type) +
  scale_fill_manual(values = c("Edited WM edges" = "#B2182B", "Unedited edges" = "#8DA0CB")) +
  labs(title = "D  Correlation change from paired baseline", x = NULL,
       y = "Mean change in count correlation") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none")

pres_plot <- preservation_summary[dependence_scale %in% c("Layer-adjusted", "Marginals")]
pres_plot[, label := fifelse(
  metric == "unedited_edge_rmse", paste0(correlation_type, " unedited-edge RMSE"),
  fifelse(metric == "cross_edge_rmse", paste0(correlation_type, " cross-edge RMSE"),
  fifelse(metric == "mean_abs_log2_ratio", "Mean preservation",
  fifelse(metric == "variance_abs_log2_ratio", "Variance preservation",
          "Dispersion preservation"))))]
p_pres <- ggplot(pres_plot, aes(label, mean)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
  geom_point(size = 2.2, color = "#0072B2") +
  coord_flip() +
  labs(title = "E  Preservation", x = NULL,
       y = "Error relative to paired ST-COS baseline") +
  theme_bw(base_size = 10)

top <- plot_grid(p_heat, ncol = 1)
bottom <- plot_grid(
  p_target, p_recovery, p_change, p_pres,
  nrow = 1, rel_widths = c(1.05, 1.05, 1.0, 1.0)
)
figure <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1.0, 1.15))
ggsave(file.path(out_dir, "DLPFC_targeted_coexpression_rewiring.pdf"),
       figure, width = 19, height = 8.5)
ggsave(file.path(out_dir, "DLPFC_targeted_coexpression_rewiring.png"),
       figure, width = 19, height = 8.5, dpi = 180)

cat(
  "Targeted co-expression rewiring experiment\n",
  sprintf("Slide: %s; spots: %d; genes: %d; replicates: %d\n",
          slide_id, nrow(reference_subset), length(selected_genes), n_reps),
  sprintf("Target genes: %s\n", paste(target_genes, collapse = ", ")),
  sprintf("Requested target latent floor: %.2f; PSD update delta: %.4f\n",
          requested_latent_floor, delta),
  sprintf("Raw-count edge threshold: r >= %.2f\n", edge_thresholds[["Raw count"]]),
  sprintf("Layer-adjusted edge threshold: residual r >= %.2f\n",
          edge_thresholds[["Layer-adjusted"]]),
  sprintf(paste0("Pearson network F1: baseline %.3f, rewired %.3f, ",
                 "SRTsim replay %.3f\n"),
          network_summary[dependence_scale == "Layer-adjusted" &
                            method == "ST-COS parameter baseline" &
                            correlation_type == "Pearson", mean_requested_edge_F1],
          network_summary[dependence_scale == "Layer-adjusted" &
                            method == "ST-COS targeted rewire" &
                            correlation_type == "Pearson", mean_requested_edge_F1],
          network_summary[dependence_scale == "Layer-adjusted" &
                            method == "SRTsim replay" &
                            correlation_type == "Pearson", mean_requested_edge_F1]),
  sprintf(paste0("Spearman network F1: baseline %.3f, rewired %.3f, ",
                 "SRTsim replay %.3f\n"),
          network_summary[dependence_scale == "Layer-adjusted" &
                            method == "ST-COS parameter baseline" &
                            correlation_type == "Spearman", mean_requested_edge_F1],
          network_summary[dependence_scale == "Layer-adjusted" &
                            method == "ST-COS targeted rewire" &
                            correlation_type == "Spearman", mean_requested_edge_F1],
          network_summary[dependence_scale == "Layer-adjusted" &
                            method == "SRTsim replay" &
                            correlation_type == "Spearman", mean_requested_edge_F1]),
  sprintf("Maximum latent change among non-target genes: %.3g\n",
          latent_non_target_change),
  file = file.path(out_dir, "README_results.txt")
)

message("Targeted co-expression rewiring experiment complete: ", out_dir)
