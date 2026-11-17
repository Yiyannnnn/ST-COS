library(data.table)
library(ggplot2)
library(stringr)

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

method_order <- c("moranI", "HEARTSVG", "SPARK", "SPARKX", "SpatialDE", "SpatialDE2")
method_cfg <- data.table(
  method = method_order,
  prefix = c("moranI", "HEARTSVG", "SPARK", "SPARKX", "spatialDE", "spatialDE2")
)

method_colors <- c(
  "moranI" = "#d98b5f",
  "HEARTSVG" = "#a7c5bd",
  "SPARK" = "#d9b44a",
  "SPARKX" = "#8ca6db",
  "SpatialDE" = "#c7b8a3",
  "SpatialDE2" = "#b8d8b8"
)

get_reps <- function(scenario) {
  sim_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated")
  patt <- paste0("^", scenario, "_gene_expression_([0-9]+)\\.csv$")
  files <- list.files(sim_dir, pattern = patt)
  reps <- suppressWarnings(as.integer(str_match(files, patt)[, 2]))
  reps <- sort(unique(reps[!is.na(reps)]))
  reps
}

# Truth genes and decision rules are defined in evaluation_rules.R
# (Supplementary Methods S2): adjusted P <= 0.05, with BH applied when a
# method reports only unadjusted P values.
get_truth_genes <- function(scenario, rep_id) {
  svg_truth_genes(sim_dir_for(base_dir, scenario), scenario, rep_id)
}

get_gene_col <- function(dt, candidates) {
  hit <- intersect(candidates, colnames(dt))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

get_pred_vec <- function(dt, method) svg_calls(dt, method)

collect_metrics <- function(scenario) {
  reps <- get_reps(scenario)
  if (length(reps) == 0) return(data.table())

  rows <- list()
  idx <- 1L

  for (m in seq_len(nrow(method_cfg))) {
    method_name <- method_cfg$method[m]
    prefix <- method_cfg$prefix[m]

    for (rep_id in reps) {
      f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), paste0(prefix, "_", rep_id, ".csv"))
      if (!file.exists(f)) next

      dt <- tryCatch(fread(f), error = function(e) NULL)
      if (is.null(dt) || nrow(dt) == 0) next

      truth_genes <- get_truth_genes(scenario, rep_id)
      if (is.null(truth_genes)) next
      gene_col <- get_gene_col(dt, c("gene", "g", "V1", "Unnamed: 0"))
      if (is.na(gene_col)) next
      pred <- get_pred_vec(dt, method_name)
      if (is.null(pred) || length(pred) != nrow(dt)) next

      genes <- as.character(dt[[gene_col]])
      keep <- !is.na(genes) & genes != ""
      genes <- genes[keep]
      pred <- pred[keep]
      if (length(genes) == 0 || length(pred) == 0) next

      met <- gene_set_metrics(genes[pred == 1L], truth_genes, genes)

      rows[[idx]] <- data.table(
        scenario = scenario,
        rep = rep_id,
        method = method_name,
        tp = met$tp,
        fp = met$fp,
        fn = met$fn,
        tn = met$tn,
        has_positive_truth = met$has_positive_truth,
        precision = met$precision,
        TPR = met$TPR,
        FPR = met$FPR,
        F1_score = met$F1
      )
      idx <- idx + 1L
    }
  }

  if (length(rows) == 0) return(data.table())
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, method := factor(method, levels = method_order)]
  out
}

plot_single_scenario <- function(metric_dt, scenario, metric_col, metric_label, filename_suffix) {
  out_pdf <- file.path(base_dir, "SVG_detection", paste0(scenario, "_", filename_suffix, ".pdf"))
  pdf(out_pdf, width = 6, height = 6)
  if (nrow(metric_dt) == 0) {
    print(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = paste0(scenario, ": no SVG outputs"), size = 5) +
        xlim(0, 1) + ylim(0, 1) +
        labs(title = paste0(scenario, " ", metric_label, " of SVG Detection Methods")) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5))
    )
  } else {
    print(
      ggplot(metric_dt, aes(x = method, y = .data[[metric_col]], fill = method)) +
        geom_boxplot(width = 0.62, linewidth = 0.35, outlier.shape = 95, outlier.size = 1.8) +
        scale_fill_manual(values = method_colors, drop = FALSE) +
        scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1.0), expand = c(0.01, 0.01)) +
        labs(title = paste0(scenario, " ", metric_label, " of SVG Detection Methods"), x = NULL, y = metric_label, fill = "Methods") +
        theme_bw(base_size = 10) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none"
        )
    )
  }
  dev.off()
}

plot_all_scenarios <- function(all_metrics, scenarios, metric_col, metric_label, out_file) {
  out_pdf <- file.path(base_dir, "SVG_detection", out_file)
  all_metrics[, scenario := factor(scenario, levels = scenarios)]

  n_scen <- length(scenarios)
  plot_w <- max(12, 2.7 * n_scen)

  pdf(out_pdf, width = plot_w, height = 4.2)
  print(
    ggplot(all_metrics, aes(x = method, y = .data[[metric_col]], fill = method)) +
      geom_boxplot(width = 0.62, linewidth = 0.35, outlier.shape = 95, outlier.size = 1.8) +
      facet_wrap(~scenario, nrow = 1, scales = "free_x", drop = FALSE) +
      scale_fill_manual(values = method_colors, drop = FALSE) +
      scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1.0), expand = c(0.01, 0.01)) +
      labs(x = NULL, y = metric_label, fill = "Methods") +
      theme_bw(base_size = 10) +
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
  )
  dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
scenarios <- if (length(args) > 0) args else c("S1", "S2", "S3", "S4", "S5", "S6", "S7")

metric_specs <- list(
  list(col = "F1_score", label = "F1 Score", suffix = "F1_score", all_file = "all_scenarios_F1_score.pdf"),
  list(col = "TPR", label = "TPR", suffix = "TPR_score", all_file = "all_scenarios_TPR_score.pdf"),
  list(col = "FPR", label = "FPR", suffix = "FPR_score", all_file = "all_scenarios_FPR_score.pdf")
)

all_rows <- list()
i <- 1L
for (scenario in scenarios) {
  metrics_dt <- collect_metrics(scenario)
  for (m in metric_specs) {
    plot_single_scenario(metrics_dt, scenario, m$col, m$label, m$suffix)
  }
  if (nrow(metrics_dt) > 0) {
    all_rows[[i]] <- metrics_dt
    i <- i + 1L
  }
}

if (length(all_rows) > 0) {
  all_metrics <- rbindlist(all_rows, use.names = TRUE, fill = TRUE)
  all_metrics[, method := factor(method, levels = method_order)]
  for (m in metric_specs) {
    plot_all_scenarios(all_metrics, scenarios, m$col, m$label, m$all_file)
  }

  scenario_tag <- paste(scenarios, collapse = "_")
  mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  median_or_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
  sd_or_na <- function(x) if (sum(!is.na(x)) > 1L) sd(x, na.rm = TRUE) else NA_real_
  summary_dt <- all_metrics[, .(
    n = .N,
    n_F1 = sum(!is.na(F1_score)),
    n_TPR = sum(!is.na(TPR)),
    n_FPR = sum(!is.na(FPR)),
    mean_F1 = mean_or_na(F1_score),
    median_F1 = median_or_na(F1_score),
    sd_F1 = sd_or_na(F1_score),
    mean_TPR = mean_or_na(TPR),
    median_TPR = median_or_na(TPR),
    sd_TPR = sd_or_na(TPR),
    mean_FPR = mean_or_na(FPR),
    median_FPR = median_or_na(FPR),
    sd_FPR = sd_or_na(FPR)
  ), by = .(scenario, method)]
  summary_dt[, `:=`(
    se_F1 = sd_F1 / sqrt(n_F1),
    se_TPR = sd_TPR / sqrt(n_TPR),
    se_FPR = sd_FPR / sqrt(n_FPR)
  )]
  summary_dt[, `:=`(
    ci_low_F1 = pmax(0, mean_F1 - qt(0.975, pmax(n_F1 - 1L, 1L)) * se_F1),
    ci_high_F1 = pmin(1, mean_F1 + qt(0.975, pmax(n_F1 - 1L, 1L)) * se_F1),
    ci_low_TPR = pmax(0, mean_TPR - qt(0.975, pmax(n_TPR - 1L, 1L)) * se_TPR),
    ci_high_TPR = pmin(1, mean_TPR + qt(0.975, pmax(n_TPR - 1L, 1L)) * se_TPR),
    ci_low_FPR = pmax(0, mean_FPR - qt(0.975, pmax(n_FPR - 1L, 1L)) * se_FPR),
    ci_high_FPR = pmin(1, mean_FPR + qt(0.975, pmax(n_FPR - 1L, 1L)) * se_FPR)
  )]

  fwrite(all_metrics, file.path(base_dir, "SVG_detection", paste0("SVG_metrics_by_rep_", scenario_tag, ".csv")))
  fwrite(summary_dt, file.path(base_dir, "SVG_detection", paste0("SVG_metrics_summary_", scenario_tag, ".csv")))
}
