#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- "."
old_dir <- file.path(root, "Reviewer3_aggregation_rerun")
source_dir <- file.path(root, "Independent_simulator_anchor", "S7SRT_heavy_figures")
out_dir <- file.path(root, "Independent_simulator_anchor", "S7SRT_spotlevel_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

legacy <- fread(file.path(old_dir, "legacy_vs_corrected_downstream_by_rep.csv"))[
  scenario == "S7" & is.finite(Legacy),
  .(analysis, metric, replicate, method, value = Legacy)
]

task_map <- data.table(
  analysis = c("Deconvolution", "Segmentation", "SVG", "CT-SVG", "Co-expression"),
  task = c("Deconvolution", "Segmentation", "SVG", "CT-SVG", "Coexpression"),
  higher_better = c(FALSE, TRUE, TRUE, TRUE, TRUE)
)
legacy <- merge(legacy, task_map, by = "analysis", all.x = TRUE)
legacy[, benchmark := "ST-COS S7 (spot-level)"]
legacy[, replicate := as.integer(replicate)]

read_srt <- function(filename, task, value_col, metric, higher_better) {
  dt <- fread(file.path(source_dir, filename))[scenario == "S7SRT"]
  dt[, .(
    task = task,
    metric = metric,
    higher_better = higher_better,
    replicate = as.integer(rep),
    method = as.character(method),
    value = as.numeric(get(value_col)),
    benchmark = "SRTsim-S7 anchor"
  )]
}

srt <- rbindlist(list(
  read_srt("fig2_deconvolution_RMSE_by_rep.csv", "Deconvolution", "RMSE", "RMSE", FALSE),
  read_srt("fig2_segmentation_ARI_by_rep.csv", "Segmentation", "ARI", "ARI", TRUE),
  read_srt("fig3_svg_metrics_by_rep.csv", "SVG", "F1_score", "F1", TRUE),
  read_srt("fig3_ctsvg_metrics_by_rep.csv", "CT-SVG", "F1_score", "F1", TRUE),
  read_srt("fig3_coexpression_metrics_by_rep.csv", "Coexpression", "F1_score", "F1", TRUE)
), use.names = TRUE)

legacy <- legacy[, .(task, metric, higher_better, replicate, method, value, benchmark)]
combined <- rbindlist(list(legacy, srt), use.names = TRUE)
setorder(combined, task, benchmark, method, replicate)
fwrite(combined, file.path(out_dir, "S7_vs_S7SRT_spotlevel_metrics_by_rep.csv"))

mean_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  avg <- mean(x)
  se <- if (n > 1L) sd(x) / sqrt(n) else NA_real_
  critical <- if (n > 1L) qt(0.975, n - 1L) else NA_real_
  data.table(mean = avg, median = median(x), sd = if (n > 1L) sd(x) else NA_real_,
             n = n, ci_low = avg - critical * se, ci_high = avg + critical * se)
}

summary_dt <- combined[, mean_ci(value),
                       by = .(benchmark, task, metric, higher_better, method)]
fwrite(summary_dt, file.path(out_dir, "S7_vs_S7SRT_spotlevel_metric_summary.csv"))

summary_dt[, rank := frank(
  if (higher_better[1L]) -mean else mean,
  ties.method = "average"
), by = .(task, benchmark)]

rank_dt <- dcast(
  summary_dt,
  task + metric + higher_better + method ~ benchmark,
  value.var = c("mean", "rank")
)
setnames(
  rank_dt,
  old = c("mean_ST-COS S7 (spot-level)", "mean_SRTsim-S7 anchor",
          "rank_ST-COS S7 (spot-level)", "rank_SRTsim-S7 anchor"),
  new = c("mean_STCOS_S7", "mean_SRTsim_S7",
          "rank_STCOS_S7", "rank_SRTsim_S7")
)
fwrite(rank_dt, file.path(out_dir, "S7_vs_S7SRT_spotlevel_method_ranks.csv"))

agreement <- rank_dt[, {
  ok <- is.finite(rank_STCOS_S7) & is.finite(rank_SRTsim_S7)
  methods <- method[ok]
  old_rank <- rank_STCOS_S7[ok]
  new_rank <- rank_SRTsim_S7[ok]
  list(
    n_methods = sum(ok),
    spearman_rank_r = if (sum(ok) >= 3L) cor(old_rank, new_rank, method = "spearman") else NA_real_,
    kendall_tau = if (sum(ok) >= 3L) cor(old_rank, new_rank, method = "kendall") else NA_real_,
    best_STCOS_S7 = methods[which.min(old_rank)],
    best_SRTsim_S7 = methods[which.min(new_rank)],
    same_best = methods[which.min(old_rank)] == methods[which.min(new_rank)],
    worst_STCOS_S7 = methods[which.max(old_rank)],
    worst_SRTsim_S7 = methods[which.max(new_rank)],
    same_worst = methods[which.max(old_rank)] == methods[which.max(new_rank)]
  )
}, by = .(task, metric, higher_better)]
fwrite(agreement, file.path(out_dir, "S7_vs_S7SRT_spotlevel_rank_concordance_by_task.csv"))

task_labels <- c(
  Deconvolution = "Deconvolution RMSE",
  Segmentation = "Segmentation ARI",
  SVG = "SVG F1",
  `CT-SVG` = "CT-SVG F1",
  Coexpression = "Co-expression F1"
)

plot_dt <- copy(combined)
plot_dt[, task_label := factor(task_labels[task], levels = task_labels)]
plot_dt[, benchmark := factor(benchmark,
                              levels = c("ST-COS S7 (spot-level)", "SRTsim-S7 anchor"))]

p <- ggplot(plot_dt, aes(x = method, y = value, fill = benchmark)) +
  geom_boxplot(position = position_dodge(width = 0.78), width = 0.66,
               linewidth = 0.3, outlier.shape = NA) +
  geom_point(aes(color = benchmark),
             position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.78),
             size = 0.9, alpha = 0.65) +
  facet_wrap(~task_label, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("ST-COS S7 (spot-level)" = "#287271",
                               "SRTsim-S7 anchor" = "#E07A5F")) +
  scale_color_manual(values = c("ST-COS S7 (spot-level)" = "#174A49",
                                "SRTsim-S7 anchor" = "#A8432D")) +
  labs(
    title = "Independent-simulator anchor for S7",
    subtitle = "Five replicates per generator; points show replicates and boxes summarize their distributions",
    x = NULL, y = "Primary metric", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "S7SRT_spotlevel_anchor.png"), p, width = 12, height = 9, dpi = 320)
ggsave(file.path(out_dir, "S7SRT_spotlevel_anchor.pdf"), p, width = 12, height = 9,
       device = cairo_pdf)

print(agreement)
