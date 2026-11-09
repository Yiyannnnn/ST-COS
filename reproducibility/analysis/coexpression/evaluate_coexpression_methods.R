#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- "."
coexp_dir <- file.path(base_dir, "coexpression")

a_align_matrix <- function(mat, n) {
  if (is.null(mat)) return(NULL)
  m <- suppressWarnings(as.matrix(mat))
  if (!is.numeric(m)) {
    storage.mode(m) <- "numeric"
  }
  m[!is.finite(m)] <- 0
  nr <- nrow(m)
  nc <- ncol(m)
  if (is.null(nr) || is.null(nc) || nr == 0 || nc == 0) return(NULL)
  k <- min(n, nr, nc)
  out <- matrix(0, nrow = n, ncol = n)
  out[1:k, 1:k] <- m[1:k, 1:k]
  diag(out) <- 0
  out
}

extract_index <- function(x) {
  if (is.null(x)) return(integer(0))
  x <- as.character(x)
  out <- suppressWarnings(as.integer(sub(".*?(\\d+)$", "\\1", x)))
  out
}

create_truth_matrix <- function(scenario, n) {
  truth <- matrix(0, nrow = n, ncol = n)
  truth_file <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_coexpression_truth.csv"))

  if (file.exists(truth_file)) {
    dt <- fread(truth_file)
    if ("gene_index" %in% names(dt)) {
      dt[, idx := as.integer(gene_index)]
    } else if ("gene_id" %in% names(dt)) {
      dt[, idx := extract_index(gene_id)]
    } else {
      dt[, idx := NA_integer_]
    }

    module_col <- if ("module_id" %in% names(dt)) {
      "module_id"
    } else if ("module_label" %in% names(dt)) {
      "module_label"
    } else {
      NULL
    }

    dt <- dt[!is.na(idx) & idx >= 1 & idx <= n]

    if (!is.null(module_col) && nrow(dt) > 0) {
      for (m in unique(dt[[module_col]])) {
        ids <- unique(dt[get(module_col) == m, idx])
        if (length(ids) > 1) {
          truth[ids, ids] <- 1
        }
      }
    }
  } else if (scenario %in% c("S4", "S5", "S6")) {
    blocks <- list(1:5, 6:10, 11:15)
    for (b in blocks) {
      ids <- b[b <= n]
      if (length(ids) > 1) {
        truth[ids, ids] <- 1
      }
    }
  }

  diag(truth) <- 0
  truth
}

compute_metrics <- function(pred, truth, threshold = 0.25) {
  pred_bin <- ifelse(abs(pred) >= threshold, 1L, 0L)
  true_bin <- ifelse(truth > 0, 1L, 0L)

  tri <- upper.tri(true_bin, diag = FALSE)
  p <- pred_bin[tri]
  t <- true_bin[tri]
  n_positive_truth <- sum(t == 1L)
  has_positive_truth <- n_positive_truth > 0L

  tp <- sum(p == 1L & t == 1L)
  fp <- sum(p == 1L & t == 0L)
  tn <- sum(p == 0L & t == 0L)
  fn <- sum(p == 0L & t == 1L)

  precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  recall <- if (has_positive_truth) tp / (tp + fn) else NA_real_
  f1 <- if (has_positive_truth && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else if (has_positive_truth) {
    0
  } else {
    NA_real_
  }
  tpr <- recall
  fpr <- if ((fp + tn) > 0) fp / (fp + tn) else 0

  c(
    F1 = f1,
    TPR = tpr,
    FPR = fpr,
    TP = tp,
    FP = fp,
    TN = tn,
    FN = fn,
    has_positive_truth = as.integer(has_positive_truth),
    n_positive_truth = as.integer(n_positive_truth)
  )
}

load_spacex <- function(scenario, rep_id, n) {
  f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), paste0("SpaceX_", rep_id, ".csv"))
  if (!file.exists(f)) return(NULL)
  mat <- as.matrix(fread(f))
  a_align_matrix(mat, n)
}

load_spmoca <- function(scenario, rep_id, n) {
  f <- file.path(coexp_dir, "spMOCA_Agent", paste0(scenario, "_", rep_id, "_correlation_matrix.csv"))
  if (!file.exists(f)) return(NULL)
  mat <- as.matrix(fread(f))
  a_align_matrix(mat, n)
}

load_cssn <- function(scenario, rep_id, n) {
  f_mean <- file.path(coexp_dir, "CSSN", paste0(scenario, "_", rep_id, "_mean_coexpression_network.csv"))
  if (file.exists(f_mean)) {
    mat <- as.matrix(fread(f_mean))
    return(a_align_matrix(mat, n))
  }

  f <- file.path(coexp_dir, "CSSN", paste0(scenario, "_", rep_id, "_CSSN_subset_results.rds"))
  if (!file.exists(f)) return(NULL)

  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj)) return(NULL)

  corr_list <- obj[["Correlation Matrix"]]
  if (is.list(corr_list)) {
    agg <- matrix(0, nrow = n, ncol = n)
    n_ok <- 0L
    for (i in seq_along(corr_list)) {
      cur <- corr_list[[i]]
      dm <- tryCatch(dim(cur), error = function(e) NULL)
      if (is.null(dm) || length(dm) != 2 || any(dm < 1)) next
      k <- min(n, dm[1], dm[2])
      sub <- tryCatch(as.matrix(cur[1:k, 1:k]), error = function(e) NULL)
      if (is.null(sub)) next
      tmp <- matrix(0, nrow = n, ncol = n)
      tmp[1:k, 1:k] <- sub
      agg <- agg + abs(tmp)
      n_ok <- n_ok + 1L
    }
    if (n_ok == 0L) return(NULL)
    agg <- agg / n_ok
    diag(agg) <- 0
    return(agg)
  }

  a_align_matrix(obj, n)
}

load_hdwgcna <- function(scenario, rep_id, n) {
  f <- file.path(coexp_dir, "hdWGCNA_Agent", paste0(scenario, "_", rep_id, "_modules.csv"))
  if (!file.exists(f)) return(NULL)

  dt <- fread(f)
  out <- matrix(0, nrow = n, ncol = n)
  if (nrow(dt) == 0) return(out)

  gene_col <- if ("gene_name" %in% names(dt)) "gene_name" else names(dt)[1]
  mod_col <- if ("module" %in% names(dt)) "module" else names(dt)[2]

  idx <- extract_index(dt[[gene_col]])
  module <- as.character(dt[[mod_col]])
  keep <- !is.na(idx) & idx >= 1 & idx <= n & !is.na(module) & module != "" & module != "grey"

  if (!any(keep)) return(out)

  dt2 <- data.table(idx = idx[keep], module = module[keep])
  for (m in unique(dt2$module)) {
    ids <- unique(dt2[module == m, idx])
    if (length(ids) > 1) out[ids, ids] <- 1
  }
  diag(out) <- 0
  out
}

load_smoothie <- function(scenario, rep_id, n) {
  f <- file.path(coexp_dir, "smoothie_Agent", paste0(scenario, "_", rep_id, "_edge_list.csv"))
  if (!file.exists(f)) return(NULL)

  dt <- fread(f)
  out <- matrix(0, nrow = n, ncol = n)
  if (nrow(dt) == 0) return(out)

  g1_col <- if ("gene1" %in% names(dt)) "gene1" else names(dt)[1]
  g2_col <- if ("gene2" %in% names(dt)) "gene2" else names(dt)[2]
  score_col <- if ("PCC" %in% names(dt)) "PCC" else if ("Rescaled_PCC" %in% names(dt)) "Rescaled_PCC" else NULL

  idx1 <- extract_index(dt[[g1_col]])
  idx2 <- extract_index(dt[[g2_col]])
  score <- if (!is.null(score_col)) abs(as.numeric(dt[[score_col]])) else rep(1, nrow(dt))

  keep <- !is.na(idx1) & !is.na(idx2) & idx1 >= 1 & idx1 <= n & idx2 >= 1 & idx2 <= n
  if (!any(keep)) return(out)

  idx1 <- idx1[keep]
  idx2 <- idx2[keep]
  score <- score[keep]

  for (i in seq_along(idx1)) {
    a <- idx1[i]
    b <- idx2[i]
    v <- score[i]
    if (!is.finite(v)) v <- 0
    if (abs(v) > abs(out[a, b])) {
      out[a, b] <- v
      out[b, a] <- v
    }
  }

  diag(out) <- 0
  out
}

load_binspect <- function(scenario, rep_id, n) {
  f1 <- file.path(coexp_dir, "binSpect_agent", paste0(scenario, "_", rep_id, "_correlation_matrix.csv"))
  f2 <- file.path(coexp_dir, "binSpect_Agent", paste0(scenario, "_", rep_id, "_correlation_matrix.csv"))
  f <- if (file.exists(f1)) f1 else if (file.exists(f2)) f2 else NULL
  if (is.null(f)) return(NULL)
  mat <- as.matrix(fread(f))
  a_align_matrix(mat, n)
}

load_stihc <- function(scenario, rep_id, n) {
  f <- file.path(coexp_dir, "stIHC_Agent", paste0(scenario, "_gene_stIHC_", rep_id, ".RDS"))
  if (!file.exists(f)) return(NULL)

  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj)) return(NULL)

  labels <- obj$label
  if (is.null(labels)) return(NULL)
  if (is.list(labels) && "label" %in% names(labels)) labels <- labels$label
  if (!is.atomic(labels)) labels <- unlist(labels)

  labels <- as.character(labels)
  L <- min(length(labels), n)
  if (L < 1) return(NULL)

  out <- matrix(0, nrow = n, ncol = n)
  labels <- labels[1:L]
  for (lab in unique(labels)) {
    ids <- which(labels == lab)
    if (length(ids) > 1) out[ids, ids] <- 1
  }
  diag(out) <- 0
  out
}

methods <- list(
  spMOCA = load_spmoca,
  SpaceX = load_spacex,
  CSSN = load_cssn,
  hdWGCNA = load_hdwgcna,
  Smoothie = load_smoothie,
  binSpect = load_binspect,
  stIHC = load_stihc
)

scenarios <- paste0("S", 1:7)
replicates <- 1:5

all_rows <- list()
missing_rows <- list()
row_id <- 1L
miss_id <- 1L

for (scenario in scenarios) {
  for (rep_id in replicates) {
    coord_file <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_coord_", rep_id, ".csv"))
    expr_file <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_gene_expression_", rep_id, ".csv"))
    if (!file.exists(coord_file) || !file.exists(expr_file)) next

    n_eval <- ncol(fread(expr_file))
    truth <- create_truth_matrix(scenario, n_eval)

    for (method_name in names(methods)) {
      pred <- tryCatch(methods[[method_name]](scenario, rep_id, n_eval), error = function(e) NULL)

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

      met <- compute_metrics(pred, truth, threshold = 0.25)
      all_rows[[row_id]] <- data.table(
        scenario = scenario,
        replicate = rep_id,
        method = method_name,
        F1 = as.numeric(met[["F1"]]),
        TPR = as.numeric(met[["TPR"]]),
        FPR = as.numeric(met[["FPR"]]),
        TP = as.integer(met[["TP"]]),
        FP = as.integer(met[["FP"]]),
        TN = as.integer(met[["TN"]]),
        FN = as.integer(met[["FN"]]),
        has_positive_truth = as.logical(met[["has_positive_truth"]]),
        n_positive_truth = as.integer(met[["n_positive_truth"]]),
        n_eval = n_eval
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

summary_dt <- metrics_dt[, .(
  n = .N,
  n_F1 = sum(!is.na(F1)),
  n_TPR = sum(!is.na(TPR)),
  n_FPR = sum(!is.na(FPR)),
  mean_F1 = if (all(is.na(F1))) NA_real_ else mean(F1, na.rm = TRUE),
  median_F1 = if (all(is.na(F1))) NA_real_ else median(F1, na.rm = TRUE),
  sd_F1 = if (sum(!is.na(F1)) > 1L) sd(F1, na.rm = TRUE) else NA_real_,
  mean_TPR = if (all(is.na(TPR))) NA_real_ else mean(TPR, na.rm = TRUE),
  median_TPR = if (all(is.na(TPR))) NA_real_ else median(TPR, na.rm = TRUE),
  sd_TPR = if (sum(!is.na(TPR)) > 1L) sd(TPR, na.rm = TRUE) else NA_real_,
  mean_FPR = if (all(is.na(FPR))) NA_real_ else mean(FPR, na.rm = TRUE),
  median_FPR = if (all(is.na(FPR))) NA_real_ else median(FPR, na.rm = TRUE),
  sd_FPR = if (sum(!is.na(FPR)) > 1L) sd(FPR, na.rm = TRUE) else NA_real_
), by = .(scenario, method)]
summary_dt[, `:=`(
  se_F1 = sd_F1 / sqrt(n_F1),
  se_TPR = sd_TPR / sqrt(n_TPR),
  se_FPR = sd_FPR / sqrt(n_FPR)
)]
summary_dt[, `:=`(
  ci_low_F1 = mean_F1 - qt(0.975, pmax(n_F1 - 1L, 1L)) * se_F1,
  ci_high_F1 = mean_F1 + qt(0.975, pmax(n_F1 - 1L, 1L)) * se_F1,
  ci_low_TPR = mean_TPR - qt(0.975, pmax(n_TPR - 1L, 1L)) * se_TPR,
  ci_high_TPR = mean_TPR + qt(0.975, pmax(n_TPR - 1L, 1L)) * se_TPR,
  ci_low_FPR = mean_FPR - qt(0.975, pmax(n_FPR - 1L, 1L)) * se_FPR,
  ci_high_FPR = mean_FPR + qt(0.975, pmax(n_FPR - 1L, 1L)) * se_FPR
)]
fwrite(summary_dt, file.path(coexp_dir, "coexpression_method_summary.csv"))

cat("Saved:\n")
cat("-", file.path(coexp_dir, "coexpression_method_metrics.csv"), "\n")
cat("-", file.path(coexp_dir, "coexpression_missing_runs.csv"), "\n")
cat("-", file.path(coexp_dir, "coexpression_method_summary.csv"), "\n")
cat("-", file.path(coexp_dir, "coexpression_all_methods_F1.pdf"), "\n")
cat("-", file.path(coexp_dir, "coexpression_all_methods_TPR_FPR.pdf"), "\n")
