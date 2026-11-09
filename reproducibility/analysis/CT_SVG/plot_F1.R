library(data.table)
library(ggplot2)
library(stringr)

base_dir <- "."

method_order <- c(
  "CELINA", "CELINA_trueprop", "CSIDE", "CSIDE_trueprop",
  "CTSV_trueprop", "STANCE_trueprop", "spVC_trueprop"
)

method_cfg <- data.table(
  method = method_order,
  prefix = c(
    "CELINA", "CELINA_trueprop", "CSIDE", "CSIDE_trueprop",
    "CTSV_truth", "STANCE_siggenes", "spVC_trueprop"
  )
)

method_colors <- c(
  "CELINA" = "#E69F00",
  "CELINA_trueprop" = "#F6C667",
  "CSIDE" = "#56B4E9",
  "CSIDE_trueprop" = "#0072B2",
  "CTSV_trueprop" = "#CC79A7",
  "STANCE_trueprop" = "#999999",
  "spVC_trueprop" = "#009E73"
)

get_reps <- function(scenario) {
  sim_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated")
  patt <- paste0("^", scenario, "_gene_expression_([0-9]+)\\.csv$")
  files <- list.files(sim_dir, pattern = patt)
  reps <- suppressWarnings(as.integer(str_match(files, patt)[, 2]))
  reps <- sort(unique(reps[!is.na(reps)]))
  reps
}

get_truth_sets <- function(scenario, rep_id) {
  truth_file <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_marker_gene_truth.csv"))
  if (file.exists(truth_file)) {
    truth_dt <- fread(truth_file)
    if (!("cell_type" %in% colnames(truth_dt)) || !("gene_id" %in% colnames(truth_dt))) {
      return(NULL)
    }
    if ("scenario" %in% colnames(truth_dt)) {
      truth_dt <- truth_dt[scenario == rep_id]
    }
    if (nrow(truth_dt) > 0) {
      return(split(unique(truth_dt$gene_id), truth_dt$cell_type))
    }
  }
  if (scenario %in% c("S2", "S5")) {
    return(list(ct1 = character(0), ct2 = character(0), ct3 = character(0)))
  }
  if (scenario %in% c("S1", "S3", "S4", "S6")) {
    return(list(
      ct1 = paste0("gene_", 1:5),
      ct2 = paste0("gene_", 6:10),
      ct3 = paste0("gene_", 11:15)
    ))
  }
  stop("No CT-SVG truth definition found for ", scenario, " replicate ", rep_id, ".")
}

get_gene_universe <- function(scenario, rep_id) {
  expr_file <- file.path(
    base_dir, "Data", paste0(scenario, "_gene"), "simulated",
    paste0(scenario, "_gene_expression_", rep_id, ".csv")
  )
  if (file.exists(expr_file)) {
    genes <- tryCatch(colnames(fread(expr_file, nrows = 0)), error = function(e) character(0))
    genes <- unique(genes[!is.na(genes) & genes != ""])
    if (length(genes) > 0) return(genes)
  }
  paste0("gene_", 1:200)
}

extract_pred_genes <- function(obj, ct, method) {
  if (is.null(obj)) {
    return(character(0))
  }

  if (method == "STANCE_trueprop") {
    if (is.list(obj) && !is.null(obj[[ct]])) {
      return(unique(as.character(obj[[ct]])))
    }
    if (isS4(obj) && "Test_2" %in% slotNames(obj)) {
      test2 <- obj@Test_2
      if (!is.null(test2[[ct]]) && is.data.frame(test2[[ct]]) && "p_value_adj" %in% colnames(test2[[ct]])) {
        return(rownames(test2[[ct]])[test2[[ct]]$p_value_adj < 0.05])
      }
    }
    return(character(0))
  }

  if (is.list(obj) && !is.null(obj[[ct]])) {
    vals <- obj[[ct]]
    if (is.null(vals) || length(vals) == 0) return(character(0))
    return(unique(as.character(vals)))
  }

  character(0)
}

calc_metrics <- function(pred, truth, universe) {
  pred <- unique(as.character(pred))
  truth <- unique(as.character(truth))
  universe <- unique(c(as.character(universe), pred, truth))

  tp <- length(intersect(pred, truth))
  fp <- length(setdiff(pred, truth))
  fn <- length(setdiff(truth, pred))
  tn <- length(setdiff(universe, union(pred, truth)))
  has_positive_truth <- length(truth) > 0L
  precision <- if (has_positive_truth && (tp + fp) > 0) tp / (tp + fp) else if (has_positive_truth) 0 else NA_real_
  recall <- if (has_positive_truth) tp / (tp + fn) else NA_real_
  tpr <- recall
  fpr <- if ((fp + tn) > 0) fp / (fp + tn) else 0
  f1 <- if (!has_positive_truth) NA_real_ else if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  list(tp = tp, fp = fp, fn = fn, tn = tn, precision = precision, recall = recall,
       tpr = tpr, fpr = fpr, f1 = f1, has_positive_truth = has_positive_truth)
}

collect_ctsvg_metrics <- function(scenario) {
  reps <- get_reps(scenario)
  if (length(reps) == 0) {
    return(list(cell = data.table(), rep = data.table()))
  }

  cell_rows <- list()
  rep_rows <- list()
  i <- 1L
  j <- 1L

  for (rep_id in reps) {
    truth_sets <- get_truth_sets(scenario, rep_id)
    if (is.null(truth_sets) || length(truth_sets) == 0) next
    cts <- names(truth_sets)
    universe_genes <- get_gene_universe(scenario, rep_id)

    for (m in seq_len(nrow(method_cfg))) {
      method_name <- method_cfg$method[m]
      prefix <- method_cfg$prefix[m]
      rds_file <- file.path(base_dir, "Data", paste0(scenario, "_gene"), paste0(prefix, "_", rep_id, ".RDS"))
      if (!file.exists(rds_file)) next

      obj <- tryCatch(readRDS(rds_file), error = function(e) NULL)
      if (is.null(obj)) next

      f1_vals <- c()
      tpr_vals <- c()
      fpr_vals <- c()
      for (ct in cts) {
        truth <- unique(as.character(truth_sets[[ct]]))
        pred <- extract_pred_genes(obj, ct, method_name)
        met <- calc_metrics(pred, truth, universe_genes)
        f1_vals <- c(f1_vals, met$f1)
        tpr_vals <- c(tpr_vals, met$tpr)
        fpr_vals <- c(fpr_vals, met$fpr)
        cell_rows[[i]] <- data.table(
          scenario = scenario,
          rep = rep_id,
          method = method_name,
          cell_type = ct,
          n_truth = length(unique(truth)),
          has_positive_truth = met$has_positive_truth,
          n_pred = length(pred),
          tp = met$tp,
          fp = met$fp,
          fn = met$fn,
          tn = met$tn,
          precision = met$precision,
          recall = met$recall,
          TPR = met$tpr,
          FPR = met$fpr,
          F1_score = met$f1
        )
        i <- i + 1L
      }

      rep_rows[[j]] <- data.table(
        scenario = scenario,
        rep = rep_id,
        method = method_name,
        has_positive_truth = any(lengths(truth_sets) > 0L),
        macro_F1 = mean(f1_vals, na.rm = TRUE),
        macro_TPR = mean(tpr_vals, na.rm = TRUE),
        macro_FPR = mean(fpr_vals, na.rm = TRUE)
      )
      j <- j + 1L
    }
  }

  cell_dt <- if (length(cell_rows) > 0) rbindlist(cell_rows, use.names = TRUE, fill = TRUE) else data.table()
  rep_dt <- if (length(rep_rows) > 0) rbindlist(rep_rows, use.names = TRUE, fill = TRUE) else data.table()
  list(cell = cell_dt, rep = rep_dt)
}

plot_all_scenarios <- function(rep_dt, scenarios, metric_col, metric_label, out_file) {
  out_pdf <- file.path(base_dir, "CT_SVG", out_file)
  rep_dt[, scenario := factor(scenario, levels = scenarios)]
  rep_dt[, method := factor(method, levels = method_order)]

  plot_w <- max(12, 2.7 * length(scenarios))
  pdf(out_pdf, width = plot_w, height = 4.2)
  print(
    ggplot(rep_dt, aes(x = method, y = .data[[metric_col]], fill = method)) +
      geom_boxplot(width = 0.62, linewidth = 0.35, outlier.shape = 95, outlier.size = 1.8) +
      facet_wrap(~scenario, nrow = 1, scales = "free_x", drop = FALSE) +
      scale_fill_manual(values = method_colors, drop = FALSE) +
      scale_y_continuous(
        limits = c(0, 1.02),
        breaks = c(0, 0.25, 0.5, 0.75, 1.0),
        expand = c(0.01, 0.01)
      ) +
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

all_cell <- list()
all_rep <- list()
i <- 1L
j <- 1L
for (scenario in scenarios) {
  res <- collect_ctsvg_metrics(scenario)
  if (nrow(res$cell) > 0) {
    all_cell[[i]] <- res$cell
    i <- i + 1L
  }
  if (nrow(res$rep) > 0) {
    all_rep[[j]] <- res$rep
    j <- j + 1L
  }
}

if (length(all_cell) == 0 || length(all_rep) == 0) {
  stop("No CT-SVG outputs found for requested scenarios.")
}

cell_dt <- rbindlist(all_cell, use.names = TRUE, fill = TRUE)
rep_dt <- rbindlist(all_rep, use.names = TRUE, fill = TRUE)
rep_dt[, method := factor(method, levels = method_order)]

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
median_or_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
sd_or_na <- function(x) if (sum(!is.na(x)) > 1L) sd(x, na.rm = TRUE) else NA_real_
summary_dt <- rep_dt[, .(
  n = .N,
  n_macro_F1 = sum(!is.na(macro_F1)),
  n_macro_TPR = sum(!is.na(macro_TPR)),
  n_macro_FPR = sum(!is.na(macro_FPR)),
  mean_macro_F1 = mean_or_na(macro_F1),
  median_macro_F1 = median_or_na(macro_F1),
  sd_macro_F1 = sd_or_na(macro_F1),
  mean_macro_TPR = mean_or_na(macro_TPR),
  median_macro_TPR = median_or_na(macro_TPR),
  sd_macro_TPR = sd_or_na(macro_TPR),
  mean_macro_FPR = mean_or_na(macro_FPR),
  median_macro_FPR = median_or_na(macro_FPR),
  sd_macro_FPR = sd_or_na(macro_FPR)
), by = .(scenario, method)]
summary_dt[, `:=`(
  se_macro_F1 = sd_macro_F1 / sqrt(n_macro_F1),
  se_macro_TPR = sd_macro_TPR / sqrt(n_macro_TPR),
  se_macro_FPR = sd_macro_FPR / sqrt(n_macro_FPR)
)]
summary_dt[, `:=`(
  ci_low_macro_F1 = pmax(0, mean_macro_F1 - qt(0.975, pmax(n_macro_F1 - 1L, 1L)) * se_macro_F1),
  ci_high_macro_F1 = pmin(1, mean_macro_F1 + qt(0.975, pmax(n_macro_F1 - 1L, 1L)) * se_macro_F1),
  ci_low_macro_TPR = pmax(0, mean_macro_TPR - qt(0.975, pmax(n_macro_TPR - 1L, 1L)) * se_macro_TPR),
  ci_high_macro_TPR = pmin(1, mean_macro_TPR + qt(0.975, pmax(n_macro_TPR - 1L, 1L)) * se_macro_TPR),
  ci_low_macro_FPR = pmax(0, mean_macro_FPR - qt(0.975, pmax(n_macro_FPR - 1L, 1L)) * se_macro_FPR),
  ci_high_macro_FPR = pmin(1, mean_macro_FPR + qt(0.975, pmax(n_macro_FPR - 1L, 1L)) * se_macro_FPR)
)]

scenario_tag <- paste(scenarios, collapse = "_")

fwrite(cell_dt, file.path(base_dir, "CT_SVG", paste0("CTSVG_F1_by_celltype_", scenario_tag, ".csv")))
fwrite(rep_dt, file.path(base_dir, "CT_SVG", paste0("CTSVG_F1_by_rep_", scenario_tag, ".csv")))
fwrite(summary_dt, file.path(base_dir, "CT_SVG", paste0("CTSVG_F1_summary_", scenario_tag, ".csv")))
fwrite(cell_dt, file.path(base_dir, "CT_SVG", paste0("CTSVG_metrics_by_celltype_", scenario_tag, ".csv")))
fwrite(rep_dt, file.path(base_dir, "CT_SVG", paste0("CTSVG_metrics_by_rep_", scenario_tag, ".csv")))
fwrite(summary_dt, file.path(base_dir, "CT_SVG", paste0("CTSVG_metrics_summary_", scenario_tag, ".csv")))

plot_all_scenarios(rep_dt, scenarios, "macro_F1", "F1 Score", "all_scenarios_F1_score.pdf")
plot_all_scenarios(rep_dt, scenarios, "macro_TPR", "TPR", "all_scenarios_TPR_score.pdf")
plot_all_scenarios(rep_dt, scenarios, "macro_FPR", "FPR", "all_scenarios_FPR_score.pdf")
