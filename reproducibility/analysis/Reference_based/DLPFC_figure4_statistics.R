#!/usr/bin/env Rscript
## ---------------------------------------------------------------
## Formal summary statistics and tests for manuscript Figure 4
## reference-based DLPFC panels.
##
## Figure 4B: gene-level marginal metrics
##   - mean/median distributions by method
##   - per-slide fidelity scores relative to Truth
##   - paired tests across 12 DLPFC slides
##
## Figure 4E: gene-gene correlation RMSE/MAE/Pearson r
##   - mean/median and 95% CI across slides
##   - paired tests across 12 DLPFC slides
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

out_dir <- "./Reference_based/DLPFC"
slides  <- c("151507","151508","151509","151510",
             "151669","151670","151671","151672",
             "151673","151674","151675","151676")

method_levels <- c("Truth", "ST-COS", "SRTsim", "scDesign3", "scCube")
sim_methods <- setdiff(method_levels, "Truth")
corr_methods <- c("ST-COS", "SRTsim", "scDesign3", "scCube")

metric_cols <- c("GeneMean", "GeneVar", "GeneCV", "GeneZeroProp")
metric_labels <- c(
  GeneMean = "Gene Mean",
  GeneVar = "Gene Variance",
  GeneCV = "Gene CV",
  GeneZeroProp = "Zero Proportion"
)
metric_transform <- c(
  GeneMean = "log10(x + 1e-3)",
  GeneVar = "log10(x + 1e-3)",
  GeneCV = "log10(x + 1e-3)",
  GeneZeroProp = "raw"
)

normalize_method <- function(x) {
  x <- as.character(x)
  x[x == "STCOS"] <- "ST-COS"
  x
}

eval_metric <- function(metric_key, value) {
  if (metric_transform[[metric_key]] == "raw") {
    value
  } else {
    log10(value + 1e-3)
  }
}

mean_ci <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- if (n > 1) sd(x) else NA_real_
  se_x <- if (n > 1) sd_x / sqrt(n) else NA_real_
  q <- if (n > 1) qt((1 + conf) / 2, df = n - 1) else NA_real_
  list(
    n = n,
    mean = mean_x,
    median = median(x),
    sd = sd_x,
    se = se_x,
    ci_low = if (n > 1) mean_x - q * se_x else NA_real_,
    ci_high = if (n > 1) mean_x + q * se_x else NA_real_,
    q25 = as.numeric(quantile(x, 0.25, names = FALSE)),
    q75 = as.numeric(quantile(x, 0.75, names = FALSE))
  )
}

safe_ttest <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || sd(x) == 0) return(NA_real_)
  suppressWarnings(t.test(x, mu = 0)$p.value)
}

safe_wilcox <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || all(x == 0)) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}

paired_tests <- function(score_dt, reference_method, comparison_methods,
                         score_col = "score_value", lower_better = TRUE) {
  rbindlist(lapply(unique(score_dt$metric), function(metric_i) {
    rbindlist(lapply(unique(score_dt$score), function(score_i) {
      sub <- score_dt[metric == metric_i & score == score_i]
      wide <- dcast(sub, slide ~ method, value.var = score_col)
      if (!reference_method %in% names(wide)) return(NULL)
      rbindlist(lapply(comparison_methods, function(comp) {
        if (!comp %in% names(wide)) return(NULL)
        diff_vals <- wide[[reference_method]] - wide[[comp]]
        diff_vals <- diff_vals[is.finite(diff_vals)]
        data.table(
          metric = metric_i,
          score = score_i,
          reference_method = reference_method,
          comparator = comp,
          n_slides = length(diff_vals),
          lower_better = lower_better,
          mean_difference = mean(diff_vals),
          median_difference = median(diff_vals),
          reference_better_by_mean = if (lower_better) mean(diff_vals) < 0 else mean(diff_vals) > 0,
          paired_t_p = safe_ttest(diff_vals),
          paired_wilcox_p = safe_wilcox(diff_vals)
        )
      }), use.names = TRUE, fill = TRUE)
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
}

summarize_scores <- function(score_dt, value_col = "score_value") {
  score_dt[, as.data.table(mean_ci(get(value_col))), by = .(metric, score, method)]
}

## =====================================================================
## Figure 4B: gene-level metrics and tests
## =====================================================================
message("Loading Figure 4B gene-level metrics ...")
gene_dt <- rbindlist(lapply(slides, function(sl) {
  f <- file.path(out_dir, paste0(sl, "_gene_metrics.csv"))
  if (!file.exists(f)) stop("Missing file: ", f)
  dt <- fread(f)
  dt[, slide := sl]
  dt
}), use.names = TRUE, fill = TRUE)

gene_dt[, method := normalize_method(method)]
gene_dt <- gene_dt[method %in% method_levels]
gene_dt[, method := factor(method, levels = method_levels)]

fig4b_values <- melt(
  gene_dt[, c("slide", "gene", "method", metric_cols), with = FALSE],
  id.vars = c("slide", "gene", "method"),
  measure.vars = metric_cols,
  variable.name = "metric_key",
  value.name = "raw_value"
)
fig4b_values[, metric := metric_labels[as.character(metric_key)]]
fig4b_values[, transform := metric_transform[as.character(metric_key)]]
fig4b_values[, eval_value := eval_metric(as.character(metric_key[1]), raw_value), by = metric_key]

fig4b_value_summary <- fig4b_values[
  is.finite(raw_value),
  {
    s_raw <- mean_ci(raw_value)
    s_eval <- mean_ci(eval_value)
    data.table(
      n_gene_slide_values = s_raw$n,
      raw_mean = s_raw$mean,
      raw_median = s_raw$median,
      raw_q25 = s_raw$q25,
      raw_q75 = s_raw$q75,
      eval_transform = transform[1],
      eval_mean = s_eval$mean,
      eval_median = s_eval$median,
      eval_q25 = s_eval$q25,
      eval_q75 = s_eval$q75
    )
  },
  by = .(metric_key, metric, method)
]

wide_gene <- dcast(
  fig4b_values,
  slide + gene + metric_key + metric + transform ~ method,
  value.var = "eval_value"
)

fig4b_errors <- rbindlist(lapply(sim_methods, function(m) {
  if (!m %in% names(wide_gene)) return(NULL)
  data.table(
    slide = wide_gene$slide,
    gene = wide_gene$gene,
    metric_key = wide_gene$metric_key,
    metric = wide_gene$metric,
    transform = wide_gene$transform,
    method = m,
    error = wide_gene[[m]] - wide_gene[["Truth"]]
  )
}), use.names = TRUE, fill = TRUE)

fig4b_slide_scores <- fig4b_errors[
  is.finite(error),
  .(
    n_genes = .N,
    RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error)),
    MedianAE = median(abs(error)),
    Bias = mean(error),
    AbsBias = abs(mean(error)),
    KS = {
      truth_vals <- wide_gene[slide == .BY$slide & metric_key == .BY$metric_key][["Truth"]]
      sim_vals <- wide_gene[slide == .BY$slide & metric_key == .BY$metric_key][[as.character(.BY$method)]]
      truth_vals <- truth_vals[is.finite(truth_vals)]
      sim_vals <- sim_vals[is.finite(sim_vals)]
      if (length(truth_vals) > 0 && length(sim_vals) > 0) {
        as.numeric(suppressWarnings(ks.test(sim_vals, truth_vals)$statistic))
      } else {
        NA_real_
      }
    }
  ),
  by = .(slide, metric_key, metric, transform, method)
]
fig4b_slide_scores[, method := factor(method, levels = sim_methods)]

fig4b_score_long <- melt(
  fig4b_slide_scores,
  id.vars = c("slide", "metric_key", "metric", "transform", "method", "n_genes"),
  measure.vars = c("RMSE", "MAE", "MedianAE", "AbsBias", "KS"),
  variable.name = "score",
  value.name = "score_value"
)

fig4b_score_summary <- summarize_scores(fig4b_score_long)
fig4b_score_summary[, method := factor(method, levels = sim_methods)]
setorder(fig4b_score_summary, metric, score, method)

fig4b_tests_vs_stcos <- paired_tests(
  fig4b_score_long,
  reference_method = "ST-COS",
  comparison_methods = setdiff(sim_methods, "ST-COS"),
  score_col = "score_value",
  lower_better = TRUE
)
fig4b_tests_vs_stcos[, paired_t_p_BH := p.adjust(paired_t_p, method = "BH"), by = score]
fig4b_tests_vs_stcos[, paired_wilcox_p_BH := p.adjust(paired_wilcox_p, method = "BH"), by = score]

## =====================================================================
## Figure 4E: gene-gene correlation summaries and tests
## =====================================================================
message("Loading Figure 4E gene-gene-correlation summaries ...")
corr_sum <- rbindlist(lapply(slides, function(sl) {
  f <- file.path(out_dir, paste0(sl, "_gene_corr_summary.csv"))
  if (!file.exists(f)) stop("Missing file: ", f)
  dt <- fread(f)
  if (!"slide" %in% names(dt)) dt[, slide := sl]
  dt
}), use.names = TRUE, fill = TRUE)

corr_sum <- corr_sum[method != "STCOS_Sigma"]
corr_sum[, method := normalize_method(method)]
corr_sum <- corr_sum[method %in% corr_methods]
corr_sum[, method := factor(method, levels = corr_methods)]

fig4e_by_slide <- corr_sum[, .(slide, method, RMSE, MAE, PearsonR)]

fig4e_score_long <- melt(
  fig4e_by_slide,
  id.vars = c("slide", "method"),
  measure.vars = c("RMSE", "MAE", "PearsonR"),
  variable.name = "metric",
  value.name = "score_value"
)
fig4e_score_long[, score := as.character(metric)]
fig4e_score_long[, metric := as.character(metric)]

fig4e_summary <- fig4e_score_long[, as.data.table(mean_ci(score_value)), by = .(metric, score, method)]
fig4e_summary[, direction := fifelse(metric %in% c("RMSE", "MAE"), "lower is better", "higher is better")]
fig4e_summary[, method := factor(method, levels = corr_methods)]
setorder(fig4e_summary, metric, method)

fig4e_tests_lower <- paired_tests(
  fig4e_score_long[metric %in% c("RMSE", "MAE")],
  reference_method = "ST-COS",
  comparison_methods = setdiff(corr_methods, "ST-COS"),
  score_col = "score_value",
  lower_better = TRUE
)
fig4e_tests_higher <- paired_tests(
  fig4e_score_long[metric == "PearsonR"],
  reference_method = "ST-COS",
  comparison_methods = setdiff(corr_methods, "ST-COS"),
  score_col = "score_value",
  lower_better = FALSE
)
fig4e_tests_vs_stcos <- rbindlist(list(fig4e_tests_lower, fig4e_tests_higher), use.names = TRUE, fill = TRUE)
fig4e_tests_vs_stcos[, paired_t_p_BH := p.adjust(paired_t_p, method = "BH"), by = score]
fig4e_tests_vs_stcos[, paired_wilcox_p_BH := p.adjust(paired_wilcox_p, method = "BH"), by = score]

## =====================================================================
## Save outputs
## =====================================================================
fwrite(fig4b_values, file.path(out_dir, "DLPFC_Figure4B_gene_metric_values_long.csv"))
fwrite(fig4b_value_summary, file.path(out_dir, "DLPFC_Figure4B_gene_metric_value_summary.csv"))
fwrite(fig4b_slide_scores, file.path(out_dir, "DLPFC_Figure4B_gene_metric_error_by_slide.csv"))
fwrite(fig4b_score_summary, file.path(out_dir, "DLPFC_Figure4B_gene_metric_error_summary.csv"))
fwrite(fig4b_tests_vs_stcos, file.path(out_dir, "DLPFC_Figure4B_gene_metric_paired_tests_vs_STCOS.csv"))

fwrite(fig4e_by_slide, file.path(out_dir, "DLPFC_Figure4E_gene_corr_by_slide.csv"))
fwrite(fig4e_summary, file.path(out_dir, "DLPFC_Figure4E_gene_corr_summary_with_CI.csv"))
fwrite(fig4e_tests_vs_stcos, file.path(out_dir, "DLPFC_Figure4E_gene_corr_paired_tests_vs_STCOS.csv"))

cat("\nFigure 4B error summary, RMSE rows:\n")
print(fig4b_score_summary[score == "RMSE"][order(metric, method)])

cat("\nFigure 4B paired tests vs ST-COS, RMSE rows:\n")
print(fig4b_tests_vs_stcos[score == "RMSE"][order(metric, comparator)])

cat("\nFigure 4E summary with 95% CI:\n")
print(fig4e_summary[metric %in% c("RMSE", "PearsonR")][order(metric, method)])

cat("\nFigure 4E paired tests vs ST-COS:\n")
print(fig4e_tests_vs_stcos[order(metric, comparator)])

message("Saved Figure 4 statistics under: ", out_dir)
