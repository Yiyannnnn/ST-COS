#!/usr/bin/env Rscript

## Strong gene-pair recovery for the reference-based DLPFC simulations.
## Uses the matched 50,000-pair samples already generated for Figure 4E.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

out_dir <- "./Reference_based/DLPFC"
slides <- c("151507", "151508", "151509", "151510",
            "151669", "151670", "151671", "151672",
            "151673", "151674", "151675", "151676")
method_cols <- c("STCOS", "SRTsim", "scDesign3", "scCube")
method_labels <- c(
  STCOS = "ST-COS",
  SRTsim = "SRTsim",
  scDesign3 = "scDesign3",
  scCube = "scCube"
)
method_levels <- unname(method_labels)
panel_a_thresholds <- c(0.01, 0.05, 0.10, 0.20, 0.50)
thresholds <- sort(unique(c(0.01, 0.05, seq(0.10, 0.35, by = 0.025), 0.40, 0.50)))
primary_threshold <- 0.20

mean_ci <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  s <- if (n > 1) sd(x) else NA_real_
  se <- if (n > 1) s / sqrt(n) else NA_real_
  q <- if (n > 1) qt((1 + conf) / 2, n - 1) else NA_real_
  data.table(
    n = n,
    mean = m,
    median = if (n) median(x) else NA_real_,
    sd = s,
    se = se,
    ci_low = if (n > 1) m - q * se else NA_real_,
    ci_high = if (n > 1) m + q * se else NA_real_
  )
}

average_precision <- function(truth, score) {
  ok <- !is.na(truth) & is.finite(score)
  truth <- as.integer(truth[ok])
  score <- score[ok]
  n_pos <- sum(truth == 1L)
  if (n_pos == 0L) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y <- truth[ord]
  precision <- cumsum(y) / seq_along(y)
  sum(precision[y == 1L]) / n_pos
}

safe_test <- function(x, type = c("t", "wilcox")) {
  type <- match.arg(type)
  x <- x[is.finite(x)]
  if (length(x) < 2L || all(x == x[1])) return(NA_real_)
  if (type == "t") {
    suppressWarnings(t.test(x, mu = 0)$p.value)
  } else {
    suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
  }
}

calc_metrics <- function(ref_r, sim_r, threshold) {
  ok <- is.finite(ref_r) & is.finite(sim_r)
  ref_r <- ref_r[ok]
  sim_r <- sim_r[ok]
  truth <- abs(ref_r) >= threshold
  pred <- abs(sim_r) >= threshold
  tp <- sum(truth & pred)
  fp <- sum(!truth & pred)
  fn <- sum(truth & !pred)
  tn <- sum(!truth & !pred)
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  f1 <- if (!is.finite(recall)) {
    NA_real_
  } else if ((precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    0
  }
  recovered_true <- truth & pred
  data.table(
    n_pairs = length(truth),
    n_reference_strong = sum(truth),
    reference_prevalence = mean(truth),
    n_predicted_strong = sum(pred),
    TP = tp,
    FP = fp,
    FN = fn,
    TN = tn,
    Precision = precision,
    Recall = recall,
    F1 = f1,
    Jaccard = if ((tp + fp + fn) > 0) tp / (tp + fp + fn) else NA_real_,
    Specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_,
    AUPRC = average_precision(truth, abs(sim_r)),
    SignAgreement = if (sum(recovered_true) > 0) mean(sign(ref_r[recovered_true]) == sign(sim_r[recovered_true])) else NA_real_
  )
}

message("Loading matched DLPFC correlation pairs ...")
pairs <- rbindlist(lapply(slides, function(slide_id) {
  f <- file.path(out_dir, paste0(slide_id, "_gene_corr_pairs.csv.gz"))
  if (!file.exists(f)) stop("Missing file: ", f)
  dt <- fread(f)
  dt[, slide := as.character(slide)]
  dt[complete.cases(dt[, c("truth", method_cols), with = FALSE])]
}), use.names = TRUE, fill = TRUE)

rows <- vector("list", length(slides) * length(thresholds) * length(method_cols))
i <- 1L
for (slide_id in slides) {
  slide_dt <- pairs[slide == slide_id]
  for (threshold in thresholds) {
    for (method_col in method_cols) {
      met <- calc_metrics(slide_dt$truth, slide_dt[[method_col]], threshold)
      rows[[i]] <- cbind(
        data.table(
          slide = slide_id,
          method = unname(method_labels[[method_col]]),
          threshold = threshold
        ),
        met
      )
      i <- i + 1L
    }
  }
}
metrics <- rbindlist(rows, use.names = TRUE, fill = TRUE)
metrics[, method := factor(method, levels = method_levels)]

metric_names <- c("Precision", "Recall", "F1", "Jaccard", "Specificity", "AUPRC", "SignAgreement")
metrics_long <- melt(
  metrics,
  id.vars = c("slide", "method", "threshold", "n_pairs", "n_reference_strong",
              "reference_prevalence", "n_predicted_strong", "TP", "FP", "FN", "TN"),
  measure.vars = metric_names,
  variable.name = "metric",
  value.name = "value"
)

summary_dt <- metrics_long[, mean_ci(value), by = .(threshold, metric, method)]
setorder(summary_dt, threshold, metric, method)

primary <- metrics_long[threshold == primary_threshold]
tests <- rbindlist(lapply(setdiff(method_levels, "ST-COS"), function(comparator) {
  rbindlist(lapply(metric_names, function(metric_i) {
    wide <- dcast(primary[metric == metric_i & method %in% c("ST-COS", comparator)],
                  slide ~ method, value.var = "value")
    if (!all(c("ST-COS", comparator) %in% names(wide))) return(NULL)
    difference <- wide[["ST-COS"]] - wide[[comparator]]
    data.table(
      threshold = primary_threshold,
      metric = metric_i,
      reference_method = "ST-COS",
      comparator = comparator,
      n_slides = sum(is.finite(difference)),
      mean_difference = mean(difference, na.rm = TRUE),
      median_difference = median(difference, na.rm = TRUE),
      paired_t_p = safe_test(difference, "t"),
      paired_wilcox_p = safe_test(difference, "wilcox")
    )
  }), use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)
tests[, paired_t_p_BH := p.adjust(paired_t_p, method = "BH"), by = metric]
tests[, paired_wilcox_p_BH := p.adjust(paired_wilcox_p, method = "BH"), by = metric]

fwrite(metrics, file.path(out_dir, "DLPFC_strong_pair_metrics_by_slide.csv"))
fwrite(metrics_long, file.path(out_dir, "DLPFC_strong_pair_metrics_long.csv"))
fwrite(summary_dt, file.path(out_dir, "DLPFC_strong_pair_summary_with_CI.csv"))
fwrite(tests, file.path(out_dir, "DLPFC_strong_pair_paired_tests_vs_STCOS.csv"))

colors <- c("ST-COS" = "#D55E00", "SRTsim" = "#0072B2",
            "scDesign3" = "#009E73", "scCube" = "#7A5195")
theme_clean <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.background = element_rect(fill = "#E8ECEF", color = "#777777"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

panel_a_dt <- metrics_long[
  metric == "F1" & round(threshold, 3) %in% panel_a_thresholds
]
panel_a_dt[, threshold_label := factor(
  sprintf("|r| >= %.2f", round(threshold, 3)),
  levels = sprintf("|r| >= %.2f", panel_a_thresholds)
)]
p_primary <- ggplot(panel_a_dt, aes(method, value, fill = method)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.2, alpha = 0.72) +
  facet_wrap(~threshold_label, nrow = 1) +
  scale_fill_manual(values = colors, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "A. Strong-pair F1 at five selected correlation thresholds",
    subtitle = "The same absolute Pearson threshold defines reference and simulated edges; points are DLPFC sections",
    x = NULL,
    y = "F1 score"
  ) +
  theme_clean

sensitivity_dt <- summary_dt[metric == "F1" & is.finite(mean)]
p_sensitivity <- ggplot(sensitivity_dt,
                        aes(threshold, mean, color = method, fill = method, group = method)) +
  geom_ribbon(aes(ymin = pmax(0, ci_low), ymax = pmin(1, ci_high)), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = colors, drop = FALSE) +
  scale_fill_manual(values = colors, drop = FALSE) +
  scale_x_continuous(
    breaks = c(0.01, 0.05, seq(0.10, 0.40, by = 0.05), 0.50),
    limits = c(0, 0.52)
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "B. F1 sensitivity across an expanded threshold grid",
    subtitle = "Dense grid 0.10-0.35 plus 0.01, 0.05, 0.40, and 0.50; 0.50 retains 1-11 reference edges per section",
    x = "Absolute correlation threshold",
    y = "Mean F1 score",
    color = NULL,
    fill = NULL
  ) +
  theme_clean +
  theme(legend.position = "bottom")

fig <- plot_grid(p_primary, p_sensitivity, ncol = 1,
                 rel_heights = c(1.15, 1),
                 align = "v", axis = "lr")
ggsave(file.path(out_dir, "DLPFC_strong_pair_recovery.pdf"), fig,
       width = 14, height = 8, bg = "white")
ggsave(file.path(out_dir, "DLPFC_strong_pair_recovery.png"), fig,
       width = 14, height = 8, dpi = 300, bg = "white")

message("F1 summary at the five displayed thresholds:")
print(summary_dt[round(threshold, 3) %in% panel_a_thresholds & metric == "F1"])
message("Saved strong-pair recovery tables and figure to: ", out_dir)
