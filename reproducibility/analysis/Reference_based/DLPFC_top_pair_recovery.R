#!/usr/bin/env Rscript

## Rank-based recovery of the strongest observed DLPFC gene pairs.
## Pair identities are represented by their common row index within each
## section's sampled pair table; every simulator is evaluated on that same set.

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
top_specs <- c("Top 100", "Top 250", "Top 500", "Top 1%")
primary_spec <- "Top 1%"

mean_ci <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  s <- if (n > 1L) sd(x) else NA_real_
  se <- if (n > 1L) s / sqrt(n) else NA_real_
  q <- if (n > 1L) qt((1 + conf) / 2, n - 1L) else NA_real_
  data.table(
    n = n,
    mean = m,
    median = if (n) median(x) else NA_real_,
    sd = s,
    se = se,
    ci_low = if (n > 1L) m - q * se else NA_real_,
    ci_high = if (n > 1L) m + q * se else NA_real_
  )
}

average_precision <- function(labels, scores) {
  ok <- is.finite(scores) & !is.na(labels)
  labels <- as.integer(labels[ok])
  scores <- scores[ok]
  n_positive <- sum(labels == 1L)
  if (!n_positive) return(NA_real_)
  ranked <- labels[order(scores, decreasing = TRUE)]
  precision <- cumsum(ranked) / seq_along(ranked)
  sum(precision[ranked == 1L]) / n_positive
}

safe_test <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || all(x == x[1L])) return(NA_real_)
  suppressWarnings(t.test(x, mu = 0)$p.value)
}

top_k_for <- function(spec, n_pairs) {
  switch(spec,
         "Top 100" = min(100L, n_pairs),
         "Top 250" = min(250L, n_pairs),
         "Top 500" = min(500L, n_pairs),
         "Top 1%" = max(1L, ceiling(0.01 * n_pairs)),
         stop("Unknown top-pair specification: ", spec))
}

calc_top_metrics <- function(reference, simulated, spec) {
  ok <- is.finite(reference) & is.finite(simulated)
  reference <- reference[ok]
  simulated <- simulated[ok]
  n_pairs <- length(reference)
  k <- top_k_for(spec, n_pairs)

  reference_top <- order(abs(reference), decreasing = TRUE)[seq_len(k)]
  simulated_top <- order(abs(simulated), decreasing = TRUE)[seq_len(k)]
  common <- intersect(reference_top, simulated_top)
  overlap <- length(common)
  overlap_rate <- overlap / k
  chance_rate <- k / n_pairs

  labels <- seq_len(n_pairs) %in% reference_top
  data.table(
    n_pairs = n_pairs,
    k = k,
    overlap = overlap,
    overlap_rate = overlap_rate,
    jaccard = overlap / (2 * k - overlap),
    enrichment_over_chance = overlap_rate / chance_rate,
    average_precision = average_precision(labels, abs(simulated)),
    sign_agreement = if (overlap) {
      mean(sign(reference[common]) == sign(simulated[common]))
    } else {
      NA_real_
    },
    reference_top_rank_spearman = if (k > 2L) {
      suppressWarnings(cor(abs(reference[reference_top]),
                           abs(simulated[reference_top]), method = "spearman"))
    } else {
      NA_real_
    },
    all_pair_spearman = suppressWarnings(cor(abs(reference), abs(simulated),
                                             method = "spearman"))
  )
}

message("Loading matched DLPFC pair samples ...")
pairs <- rbindlist(lapply(slides, function(slide_id) {
  path <- file.path(out_dir, paste0(slide_id, "_gene_corr_pairs.csv.gz"))
  if (!file.exists(path)) stop("Missing file: ", path)
  result <- fread(path)
  result[, slide := as.character(slide)]
  result[complete.cases(result[, c("truth", method_cols), with = FALSE])]
}), use.names = TRUE, fill = TRUE)

rows <- vector("list", length(slides) * length(method_cols) * length(top_specs))
index <- 1L
for (slide_id in slides) {
  slide_pairs <- pairs[slide == slide_id]
  for (method_col in method_cols) {
    for (spec in top_specs) {
      rows[[index]] <- cbind(
        data.table(
          slide = slide_id,
          method = unname(method_labels[[method_col]]),
          top_set = spec
        ),
        calc_top_metrics(slide_pairs$truth, slide_pairs[[method_col]], spec)
      )
      index <- index + 1L
    }
  }
}

metrics <- rbindlist(rows)
metrics[, method := factor(method, levels = method_levels)]
metrics[, top_set := factor(top_set, levels = top_specs)]

metric_names <- c("overlap_rate", "jaccard", "enrichment_over_chance",
                  "average_precision", "sign_agreement",
                  "reference_top_rank_spearman", "all_pair_spearman")
metrics_long <- melt(
  metrics,
  id.vars = c("slide", "method", "top_set", "n_pairs", "k", "overlap"),
  measure.vars = metric_names,
  variable.name = "metric",
  value.name = "value"
)
summary_dt <- metrics_long[, mean_ci(value), by = .(top_set, metric, method)]
setorder(summary_dt, top_set, metric, method)

primary <- metrics_long[top_set == primary_spec]
tests <- rbindlist(lapply(setdiff(method_levels, "ST-COS"), function(comparator) {
  rbindlist(lapply(metric_names, function(metric_i) {
    wide <- dcast(
      primary[metric == metric_i & method %in% c("ST-COS", comparator)],
      slide ~ method,
      value.var = "value"
    )
    difference <- wide[["ST-COS"]] - wide[[comparator]]
    data.table(
      top_set = primary_spec,
      metric = metric_i,
      reference_method = "ST-COS",
      comparator = comparator,
      n_slides = sum(is.finite(difference)),
      mean_difference = mean(difference, na.rm = TRUE),
      median_difference = median(difference, na.rm = TRUE),
      paired_t_p = safe_test(difference)
    )
  }))
}))
tests[, paired_t_p_BH := p.adjust(paired_t_p, method = "BH"), by = metric]

fwrite(metrics, file.path(out_dir, "DLPFC_top_pair_metrics_by_slide.csv"))
fwrite(metrics_long, file.path(out_dir, "DLPFC_top_pair_metrics_long.csv"))
fwrite(summary_dt, file.path(out_dir, "DLPFC_top_pair_summary_with_CI.csv"))
fwrite(tests, file.path(out_dir, "DLPFC_top_pair_paired_tests_vs_STCOS.csv"))

colors <- c("ST-COS" = "#D55E00", "SRTsim" = "#0072B2",
            "scDesign3" = "#009E73", "scCube" = "#7A5195")
theme_clean <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.background = element_rect(fill = "#E8ECEF", color = "#777777"),
    strip.text = element_text(face = "bold")
  )

primary_plot <- metrics_long[
  top_set == primary_spec &
    metric %in% c("overlap_rate", "jaccard", "average_precision")
]
primary_plot[, metric_label := factor(
  metric,
  levels = c("overlap_rate", "jaccard", "average_precision"),
  labels = c("Top-set overlap", "Jaccard", "Average precision")
)]
p_primary <- ggplot(primary_plot, aes(method, value, fill = method)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.10, height = 0),
             size = 1.2, alpha = 0.72) +
  facet_wrap(~metric_label, nrow = 1) +
  scale_fill_manual(values = colors, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "A. Recovery of the top 1% observed gene pairs",
    subtitle = "Common sampled pair set within each of 12 DLPFC sections",
    x = NULL,
    y = "Recovery metric"
  ) +
  theme_clean +
  theme(legend.position = "none")

sensitivity <- summary_dt[metric %in% c("overlap_rate", "average_precision")]
sensitivity[, metric_label := factor(
  metric,
  levels = c("overlap_rate", "average_precision"),
  labels = c("Top-set overlap", "Average precision")
)]
p_sensitivity <- ggplot(
  sensitivity,
  aes(top_set, mean, color = method, group = method)
) +
  geom_errorbar(aes(ymin = pmax(0, ci_low), ymax = pmin(1, ci_high)),
                width = 0.10, linewidth = 0.35) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2) +
  facet_wrap(~metric_label, nrow = 1) +
  scale_color_manual(values = colors, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "B. Sensitivity to top-set size",
    subtitle = "Means and 95% t confidence intervals across sections",
    x = "Reference top-pair definition",
    y = "Mean recovery metric",
    color = NULL
  ) +
  theme_clean +
  theme(legend.position = "bottom")

figure <- plot_grid(p_primary, p_sensitivity, ncol = 1,
                    rel_heights = c(1.05, 1), align = "v", axis = "lr")
ggsave(file.path(out_dir, "DLPFC_top_pair_recovery.pdf"), figure,
       width = 13.5, height = 8.5, bg = "white")
ggsave(file.path(out_dir, "DLPFC_top_pair_recovery.png"), figure,
       width = 13.5, height = 8.5, dpi = 300, bg = "white")

message("Top-1% summary:")
print(summary_dt[
  top_set == primary_spec &
    metric %in% c("overlap_rate", "jaccard", "average_precision",
                  "sign_agreement", "reference_top_rank_spearman")
])
message("Saved top-pair recovery outputs to: ", out_dir)
