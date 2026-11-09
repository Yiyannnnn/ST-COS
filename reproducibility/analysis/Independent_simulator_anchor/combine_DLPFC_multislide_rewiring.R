#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

repo_dir <- "."
input_root <- file.path(
  repo_dir, "Independent_simulator_anchor",
  "DLPFC_targeted_coexpression_rewiring", "slides"
)
out_dir <- file.path(
  repo_dir, "Independent_simulator_anchor",
  "DLPFC_targeted_coexpression_rewiring", "multislide"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

slide_ids <- as.character(c(151507:151510, 151669:151676))
donor_map <- data.table(
  slide_id = slide_ids,
  donor_group = rep(c("Donor group 1", "Donor group 2", "Donor group 3"), each = 4)
)

read_slide <- function(slide_id, filename) {
  sid <- as.character(slide_id)
  path <- file.path(input_root, slide_id, filename)
  if (!file.exists(path)) stop("Missing slide output: ", path)
  dt <- fread(path)
  if ("slide_id" %in% names(dt)) dt[, slide_id := NULL]
  dt[, slide_id := sid]
  merge(dt, donor_map, by = "slide_id", all.x = TRUE, sort = FALSE)
}

metrics <- rbindlist(lapply(slide_ids, read_slide,
                            filename = "network_metrics_by_replicate.csv"), fill = TRUE)
summaries <- rbindlist(lapply(slide_ids, read_slide,
                              filename = "network_metric_summary.csv"), fill = TRUE)
thresholds <- rbindlist(lapply(slide_ids, read_slide,
                               filename = "network_threshold_curve.csv"), fill = TRUE)
preservation <- rbindlist(lapply(slide_ids, read_slide,
                                 filename = "marginal_and_edge_preservation_summary.csv"), fill = TRUE)

baseline <- metrics[method == "ST-COS parameter baseline", .(
  slide_id, donor_group, replicate, correlation_type, dependence_scale,
  baseline_correlation = mean_target_correlation,
  baseline_F1 = requested_edge_F1
)]
rewired <- metrics[method == "ST-COS targeted rewire", .(
  slide_id, donor_group, replicate, correlation_type, dependence_scale,
  rewired_correlation = mean_target_correlation,
  rewired_F1 = requested_edge_F1
)]
effects <- merge(
  baseline, rewired,
  by = c("slide_id", "donor_group", "replicate", "correlation_type", "dependence_scale")
)
effects[, `:=`(
  correlation_change = rewired_correlation - baseline_correlation,
  F1_change = rewired_F1 - baseline_F1
)]

slide_effects <- effects[, .(
  mean_correlation_change = mean(correlation_change),
  correlation_ci_low = quantile(correlation_change, 0.025),
  correlation_ci_high = quantile(correlation_change, 0.975),
  mean_F1_change = mean(F1_change),
  F1_ci_low = quantile(F1_change, 0.025),
  F1_ci_high = quantile(F1_change, 0.975)
), by = .(slide_id, donor_group, correlation_type, dependence_scale)]

method_by_slide <- summaries[, .(
  mean_target_correlation = mean_target_correlation,
  requested_edge_precision = mean_requested_edge_precision,
  requested_edge_recall = mean_requested_edge_recall,
  requested_edge_F1 = mean_requested_edge_F1,
  requested_edge_AUPRC = mean_requested_edge_AUPRC,
  n_simulation_replicates = n_replicates
), by = .(slide_id, donor_group, method, correlation_type,
          dependence_scale, edge_threshold)]

method_across_slides <- method_by_slide[, .(
  mean_target_correlation = mean(mean_target_correlation),
  correlation_ci_low = quantile(mean_target_correlation, 0.025),
  correlation_ci_high = quantile(mean_target_correlation, 0.975),
  mean_requested_edge_F1 = mean(requested_edge_F1),
  F1_ci_low = quantile(requested_edge_F1, 0.025),
  F1_ci_high = quantile(requested_edge_F1, 0.975),
  mean_requested_edge_AUPRC = mean(requested_edge_AUPRC),
  n_slides = uniqueN(slide_id)
), by = .(method, correlation_type, dependence_scale, edge_threshold)]

cluster_boot_ci <- function(values, donors, n_boot = 10000L, seed = 20260912L) {
  set.seed(seed)
  donor_means <- tapply(values, donors, mean)
  boot <- replicate(n_boot, mean(sample(donor_means, length(donor_means), replace = TRUE)))
  as.numeric(quantile(boot, c(0.025, 0.975)))
}

pooled_tests <- rbindlist(lapply(c("Raw count", "Layer-adjusted"), function(scale) {
  rbindlist(lapply(c("Pearson", "Spearman"), function(ct) {
    current <- slide_effects[dependence_scale == scale & correlation_type == ct]
    rbindlist(lapply(c("mean_correlation_change", "mean_F1_change"), function(metric) {
      x <- current[[metric]]
      tt <- t.test(x, mu = 0)
      wt <- wilcox.test(x, mu = 0, exact = FALSE)
      boot_ci <- cluster_boot_ci(x, current$donor_group)
      data.table(
        dependence_scale = scale,
        correlation_type = ct,
        metric = metric,
        mean_slide_effect = mean(x),
        section_t_ci_low = tt$conf.int[1L],
        section_t_ci_high = tt$conf.int[2L],
        donor_cluster_boot_ci_low = boot_ci[1L],
        donor_cluster_boot_ci_high = boot_ci[2L],
        section_t_p_value = tt$p.value,
        section_wilcoxon_p_value = wt$p.value,
        n_slides = length(x),
        n_donor_groups = uniqueN(current$donor_group)
      )
    }))
  }))
}))

threshold_by_slide <- thresholds[, .(
  precision = mean(precision), recall = mean(recall), F1 = mean(F1)
), by = .(slide_id, donor_group, method, correlation_type, threshold)]
threshold_summary <- threshold_by_slide[, .(
  mean_F1 = mean(F1),
  ci_low = quantile(F1, 0.025),
  ci_high = quantile(F1, 0.975),
  n_slides = uniqueN(slide_id)
), by = .(method, correlation_type, threshold)]

preservation_across_slides <- preservation[, .(
  mean = mean(mean),
  ci_low = quantile(mean, 0.025),
  ci_high = quantile(mean, 0.975),
  n_slides = uniqueN(slide_id)
), by = .(dependence_scale, correlation_type, metric)]

fwrite(effects, file.path(out_dir, "paired_effects_by_replicate_and_slide.csv"))
fwrite(slide_effects, file.path(out_dir, "paired_effects_by_slide.csv"))
fwrite(method_by_slide, file.path(out_dir, "network_metrics_by_slide.csv"))
fwrite(method_across_slides, file.path(out_dir, "network_metrics_across_slides.csv"))
fwrite(pooled_tests, file.path(out_dir, "pooled_section_and_donor_blocked_tests.csv"))
fwrite(threshold_by_slide, file.path(out_dir, "threshold_curves_by_slide.csv"))
fwrite(threshold_summary, file.path(out_dir, "threshold_curve_summary.csv"))
fwrite(preservation_across_slides, file.path(out_dir, "preservation_across_slides.csv"))

plot_effects <- slide_effects[dependence_scale == "Layer-adjusted"]
plot_effects[, slide_id := factor(slide_id, levels = rev(slide_ids))]
donor_colors <- c("Donor group 1" = "#0072B2", "Donor group 2" = "#D55E00",
                  "Donor group 3" = "#009E73")

p_corr <- ggplot(plot_effects, aes(mean_correlation_change, slide_id,
                                   color = donor_group)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_errorbarh(aes(xmin = correlation_ci_low, xmax = correlation_ci_high),
                 height = 0.2) +
  geom_point(size = 2) +
  facet_wrap(~correlation_type, scales = "free_x") +
  scale_color_manual(values = donor_colors) +
  labs(title = "A  Rewiring effect by DLPFC section", x = "Change in target correlation",
       y = "Section", color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

p_f1 <- ggplot(plot_effects, aes(mean_F1_change, slide_id, color = donor_group)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_errorbarh(aes(xmin = F1_ci_low, xmax = F1_ci_high), height = 0.2) +
  geom_point(size = 2) +
  facet_wrap(~correlation_type, scales = "free_x") +
  scale_color_manual(values = donor_colors) +
  labs(title = "B  Requested-edge recovery effect", x = "Change in network F1",
       y = "Section", color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

method_levels <- c("ST-COS replay", "SRTsim replay", "ST-COS parameter baseline",
                   "ST-COS targeted rewire")
method_labels <- c("ST-COS replay", "SRTsim replay", "ST-COS baseline", "ST-COS rewired")
method_colors <- c("#6C757D", "#D55E00", "#0072B2", "#009E73")
plot_methods <- method_by_slide[dependence_scale == "Layer-adjusted"]
plot_methods[, method := factor(method, levels = method_levels, labels = method_labels)]
p_methods <- ggplot(plot_methods, aes(method, mean_target_correlation,
                                      color = donor_group, group = slide_id)) +
  geom_line(alpha = 0.28, color = "#777777") +
  geom_point(size = 1.8) +
  facet_wrap(~correlation_type) +
  scale_color_manual(values = donor_colors) +
  labs(title = "C  Layer-adjusted dependence", x = NULL,
       y = "Mean target-pair correlation", color = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")

curve_methods <- c("ST-COS replay", "SRTsim replay", "ST-COS parameter baseline",
                   "ST-COS targeted rewire")
curve_plot <- threshold_summary[method %in% curve_methods]
curve_plot[, method := factor(method, levels = method_levels, labels = method_labels)]
layer_edge_threshold <- unique(method_by_slide[dependence_scale == "Layer-adjusted", edge_threshold])
stopifnot(length(layer_edge_threshold) == 1L)
p_curve <- ggplot(curve_plot, aes(threshold, mean_F1, color = method, fill = method)) +
  geom_vline(xintercept = layer_edge_threshold,
             linetype = "dashed", color = "#4D4D4D", linewidth = 0.45) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~correlation_type) +
  scale_color_manual(values = setNames(method_colors, method_labels)) +
  scale_fill_manual(values = setNames(method_colors, method_labels)) +
  labs(title = "D  Edge recovery across thresholds",
       subtitle = "Dashed line: prespecified layer-adjusted r = 0.55",
       x = "Layer-adjusted residual correlation threshold",
       y = "Network F1", color = NULL, fill = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 8.5))

pres_plot <- preservation[
  dependence_scale %in% c("Layer-adjusted", "Marginals")
]
pres_plot[, metric_group := fifelse(
  dependence_scale == "Marginals",
  "Marginal absolute log2 ratios",
  "Layer-adjusted correlation RMSE"
)]
pres_plot[, display_metric := fifelse(
  metric == "unedited_edge_rmse", paste(correlation_type, "unedited edges"),
  fifelse(metric == "cross_edge_rmse", paste(correlation_type, "cross edges"),
  fifelse(metric == "mean_abs_log2_ratio", "Gene mean",
  fifelse(metric == "variance_abs_log2_ratio", "Gene variance",
          "NB dispersion"))))]
p_pres <- ggplot(pres_plot, aes(mean, display_metric, color = donor_group)) +
  geom_point(alpha = 0.8) +
  facet_grid(metric_group ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = donor_colors) +
  labs(title = "E  Preservation across sections", x = "Unitless magnitude from paired baseline",
       y = NULL, color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.text.y = element_text(size = 8.2))

top <- plot_grid(p_corr, p_f1, nrow = 1, rel_widths = c(1, 1))
bottom <- plot_grid(p_methods, p_curve, p_pres, nrow = 1,
                    rel_widths = c(1, 1.15, 0.9))
figure <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1.1, 1))
ggsave(file.path(out_dir, "DLPFC_multislide_targeted_rewiring.pdf"),
       figure, width = 18, height = 10)
ggsave(file.path(out_dir, "DLPFC_multislide_targeted_rewiring.png"),
       figure, width = 18, height = 10, dpi = 180)

layer_tests <- pooled_tests[dependence_scale == "Layer-adjusted"]
cat(
  "Multi-slide targeted co-expression rewiring\n",
  sprintf("Sections: %d; donor groups: %d; paired simulations per section: 30.\n",
          length(slide_ids), uniqueN(donor_map$donor_group)),
  paste(apply(layer_tests, 1L, function(row) sprintf(
    "%s %s: mean section effect %s, donor-cluster bootstrap CI [%s, %s].",
    row[["correlation_type"]], row[["metric"]],
    format(as.numeric(row[["mean_slide_effect"]]), digits = 4),
    format(as.numeric(row[["donor_cluster_boot_ci_low"]]), digits = 4),
    format(as.numeric(row[["donor_cluster_boot_ci_high"]]), digits = 4)
  )), collapse = "\n"),
  "\nMonte Carlo replicates quantify simulation variability; section and donor-group summaries quantify cross-section replication.\n",
  file = file.path(out_dir, "README_results.txt")
)

message("Combined multi-slide rewiring outputs written to: ", out_dir)
