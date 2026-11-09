#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

coexp_dir <- "./coexpression"
metrics_file <- file.path(coexp_dir, "coexpression_method_metrics.csv")

if (!file.exists(metrics_file)) {
  stop("Missing metrics file: ", metrics_file)
}

metrics_dt <- fread(metrics_file)
if (nrow(metrics_dt) == 0) {
  stop("coexpression_method_metrics.csv is empty.")
}

scenario_levels <- paste0("S", 1:7)
method_levels <- c("spMOCA", "SpaceX", "CSSN", "hdWGCNA", "Smoothie", "binSpect", "stIHC")
method_colors <- c(
  spMOCA = "#4B9CD3",
  SpaceX = "#4c9390",
  CSSN = "#d65d48",
  hdWGCNA = "#8b5a96",
  Smoothie = "#f4a261",
  binSpect = "#2a9d8f",
  stIHC = "#264653"
)

metrics_dt <- metrics_dt[
  scenario %in% scenario_levels & method %in% method_levels
]
metrics_dt[, scenario := factor(scenario, levels = scenario_levels)]
metrics_dt[, method := factor(method, levels = method_levels)]

if (nrow(metrics_dt) == 0) {
  stop("No rows left after filtering to S1-S7 and known methods.")
}

plot_w <- max(12, 2.7 * length(scenario_levels))
plot_h <- 4.2

base_theme <- theme_bw(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.35),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "#d2e0e3", color = "#1a1a1a", linewidth = 0.9),
    strip.text = element_text(size = 24, face = "bold", color = "#111111"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.14, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.2, color = "black"),
    axis.text.y = element_text(size = 8.6, color = "black"),
    axis.title.y = element_text(size = 10.5, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9.3),
    plot.margin = margin(8, 10, 8, 8)
  )

scenario_annotations <- unique(
  metrics_dt[is.na(F1) | is.na(TPR), .(scenario)]
)[, .(
  scenario,
  x = 0.6,
  xend = length(method_levels) + 0.4,
  y = 0.04,
  yend = 0.96
)]

plot_metric <- function(metric_name, y_lab, out_file, width = plot_w, height = plot_h) {
  dt <- metrics_dt[, .(scenario, method, value = get(metric_name))]
  p <- ggplot(dt, aes(x = method, y = value, fill = method)) +
    geom_boxplot(width = 0.62, linewidth = 0.35, outlier.shape = 95, outlier.size = 1.8) +
    facet_wrap(~scenario, nrow = 1, scales = "free_x", drop = FALSE) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 1.02),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      expand = c(0.01, 0.01)
    ) +
    labs(x = NULL, y = y_lab, fill = "Methods") +
    base_theme

  if (metric_name %in% c("F1", "TPR") && nrow(scenario_annotations) > 0) {
    p <- p + geom_segment(
      data = scenario_annotations,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      color = "#666666",
      linewidth = 0.75,
      lineend = "round"
    )
  }

  ggsave(out_file, p, width = width, height = height, bg = "white")
}

plot_metric(
  metric_name = "F1",
  y_lab = "F1 Score",
  out_file = file.path(coexp_dir, "coexpression_all_methods_F1.pdf")
)

plot_metric(
  metric_name = "TPR",
  y_lab = "True Positive Rate (TPR)",
  out_file = file.path(coexp_dir, "all_methods_TPR_comparison.pdf")
)

plot_metric(
  metric_name = "FPR",
  y_lab = "False Positive Rate (FPR)",
  out_file = file.path(coexp_dir, "all_methods_FPR_comparison.pdf")
)

long_dt <- melt(
  metrics_dt,
  id.vars = c("scenario", "replicate", "method"),
  measure.vars = c("TPR", "FPR"),
  variable.name = "metric",
  value.name = "value"
)
long_dt[, metric := factor(metric, levels = c("TPR", "FPR"))]

p_combined <- ggplot(long_dt, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.62, linewidth = 0.35, outlier.shape = 95, outlier.size = 1.8) +
  facet_grid(metric ~ scenario) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = c(0.01, 0.01)
  ) +
  labs(x = NULL, y = "Rate", fill = "Methods") +
  base_theme

ggsave(
  filename = file.path(coexp_dir, "all_methods_combined_evaluation.pdf"),
  plot = p_combined,
  width = plot_w,
  height = 10,
  bg = "white"
)

message("Updated coexpression F1/TPR/FPR plots with S1-S7 strip style and unified size.")
