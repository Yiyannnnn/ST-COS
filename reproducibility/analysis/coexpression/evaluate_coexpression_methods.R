#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  here <- if (length(f) == 1L) dirname(normalizePath(f)) else getwd()
  rules <- c(Sys.getenv("STCOS_EVALUATION_RULES"),
             file.path(here, "..", "evaluation_rules.R"),
             file.path(getwd(), "evaluation_rules.R"))
  rules <- rules[nzchar(rules) & file.exists(rules)]
  if (!length(rules)) stop("Cannot find evaluation_rules.R; set STCOS_EVALUATION_RULES.")
  source(rules[1L], local = globalenv())
})

base_dir <- "."
coexp_dir <- file.path(base_dir, "coexpression")

# Edge rule, truth, loaders, and metrics follow Supplementary Methods S2 and are
# defined in evaluation_rules.R.
method_names <- c("spMOCA", "SpaceX", "CSSN", "hdWGCNA", "Smoothie", "binSpect", "stIHC")

scenarios <- paste0("S", 1:7)
replicates <- 1:5

all_rows <- list()
missing_rows <- list()
row_id <- 1L
miss_id <- 1L

for (scenario in scenarios) {
  sim_dir <- sim_dir_for(base_dir, scenario)
  for (rep_id in replicates) {
    coord_file <- file.path(sim_dir, paste0(scenario, "_coord_", rep_id, ".csv"))
    expr_file <- file.path(sim_dir, paste0(scenario, "_gene_expression_", rep_id, ".csv"))
    if (!file.exists(coord_file) || !file.exists(expr_file)) next

    n_eval <- ncol(fread(expr_file, nrows = 0))
    truth <- coexp_truth_matrix(sim_dir, scenario, rep_id, n_eval)
    if (is.null(truth)) {
      for (method_name in method_names) {
        missing_rows[[miss_id]] <- data.table(
          scenario = scenario,
          replicate = rep_id,
          method = method_name,
          status = "missing_truth"
        )
        miss_id <- miss_id + 1L
      }
      next
    }
    truth_source <- attr(truth, "truth_source")
    if (is.null(truth_source)) truth_source <- "generator_definition"

    for (method_name in method_names) {
      pred <- tryCatch(
        load_coexp_prediction(base_dir, scenario, rep_id, method_name, n_eval),
        error = function(e) NULL
      )

      if (is.null(pred)) {
        missing_rows[[miss_id]] <- data.table(
          scenario = scenario,
          replicate = rep_id,
          method = method_name,
          status = "missing_or_failed"
        )
        miss_id <- miss_id + 1L
        next
      }

      met <- coexp_metrics(
        pred, truth,
        continuous = method_name %in% COEXP_CONTINUOUS_METHODS
      )
      all_rows[[row_id]] <- data.table(
        scenario = scenario,
        replicate = rep_id,
        method = method_name,
        F1 = met$F1,
        TPR = met$TPR,
        FPR = met$FPR,
        precision = met$precision,
        Jaccard = met$Jaccard,
        AUPRC = met$AUPRC,
        TP = as.integer(met$tp),
        FP = as.integer(met$fp),
        TN = as.integer(met$tn),
        FN = as.integer(met$fn),
        has_positive_truth = met$has_positive_truth,
        n_positive_truth = as.integer(met$n_positive_truth),
        n_eval = n_eval,
        edge_threshold = COEXP_EDGE_THRESHOLD,
        truth_source = truth_source
      )
      row_id <- row_id + 1L
    }
  }
}

metrics_dt <- if (length(all_rows) > 0) rbindlist(all_rows) else data.table()
missing_dt <- if (length(missing_rows) > 0) rbindlist(missing_rows) else data.table()

fwrite(metrics_dt, file.path(coexp_dir, "coexpression_method_metrics.csv"))
fwrite(missing_dt, file.path(coexp_dir, "coexpression_missing_runs.csv"))

if (nrow(metrics_dt) == 0) {
  stop("No evaluation rows available; nothing to plot.")
}

method_levels <- c("spMOCA", "SpaceX", "CSSN", "hdWGCNA", "Smoothie", "binSpect", "stIHC")
metrics_dt[, method := factor(method, levels = method_levels)]
metrics_dt[, scenario := factor(scenario, levels = scenarios)]

method_colors <- c(
  spMOCA = "#4B9CD3",
  SpaceX = "#4c9390",
  CSSN = "#d65d48",
  hdWGCNA = "#8b5a96",
  Smoothie = "#f4a261",
  binSpect = "#2a9d8f",
  stIHC = "#264653"
)

p_f1 <- ggplot(metrics_dt, aes(x = method, y = F1, fill = method)) +
  geom_boxplot(width = 0.62, linewidth = 0.35, outlier.shape = 95, outlier.size = 1.8) +
  facet_wrap(~ scenario, nrow = 1, scales = "free_x", drop = FALSE) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = c(0.01, 0.01)
  ) +
  labs(x = NULL, y = "F1 Score", fill = "Methods") +
  theme_bw(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.35),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "#d2e0e3", color = "#6f6f6f", linewidth = 0.35),
    strip.text = element_text(size = 11, face = "bold"),
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

ggsave(
  filename = file.path(coexp_dir, "coexpression_all_methods_F1.pdf"),
  plot = p_f1,
  width = 20,
  height = 6,
  bg = "white"
)

long_dt <- melt(
  metrics_dt,
  id.vars = c("scenario", "replicate", "method"),
  measure.vars = c("TPR", "FPR"),
  variable.name = "metric",
  value.name = "value"
)

p_tpr_fpr <- ggplot(long_dt, aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.62, linewidth = 0.35, outlier.shape = 95, outlier.size = 1.8) +
  facet_grid(metric ~ scenario) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = c(0.01, 0.01)
  ) +
  labs(x = NULL, y = "Rate", fill = "Methods") +
  theme_bw(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.35),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "#d2e0e3", color = "#6f6f6f", linewidth = 0.35),
    strip.text = element_text(size = 11, face = "bold"),
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

ggsave(
  filename = file.path(coexp_dir, "coexpression_all_methods_TPR_FPR.pdf"),
  plot = p_tpr_fpr,
  width = 20,
  height = 8,
  bg = "white"
)

summary_metrics <- c("F1", "TPR", "FPR", "precision", "Jaccard", "AUPRC")
summarise_metric <- function(x) {
  n <- sum(!is.na(x))
  m <- mean_or_na(x)
  s <- if (n > 1L) sd(x, na.rm = TRUE) else NA_real_
  half <- if (n > 1L) qt(0.975, n - 1L) * s / sqrt(n) else NA_real_
  list(n = n, mean = m, median = if (n > 0L) median(x, na.rm = TRUE) else NA_real_,
       sd = s, ci_low = m - half, ci_high = m + half)
}
summary_dt <- metrics_dt[, {
  out <- list(n = .N)
  for (metric in summary_metrics) {
    st <- summarise_metric(get(metric))
    names(st) <- paste0(names(st), "_", metric)
    out <- c(out, st)
  }
  out
}, by = .(scenario, method)]
fwrite(summary_dt, file.path(coexp_dir, "coexpression_method_summary.csv"))

cat("Saved:\n")
cat("-", file.path(coexp_dir, "coexpression_method_metrics.csv"), "\n")
cat("-", file.path(coexp_dir, "coexpression_missing_runs.csv"), "\n")
cat("-", file.path(coexp_dir, "coexpression_method_summary.csv"), "\n")
cat("-", file.path(coexp_dir, "coexpression_all_methods_F1.pdf"), "\n")
cat("-", file.path(coexp_dir, "coexpression_all_methods_TPR_FPR.pdf"), "\n")
