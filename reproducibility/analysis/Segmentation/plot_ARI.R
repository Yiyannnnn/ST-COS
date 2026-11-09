library(data.table)
library(ggplot2)
library(stringr)

base_dir <- "."

method_order <- c("BASS", "BayesSpace", "GraphST", "MAPLE", "SpaGCN", "STAGATE")
method_cfg <- data.table(
  method = method_order,
  prefix = method_order,
  gt_col = c("layer", "layer", "Ground Truth", "layer", "Ground Truth", "Ground Truth"),
  pred_col = c("BASS_layer", "BayesSpace_labels", "domain", "cluster", "refined_pred", "mclust")
)

method_colors <- c(
  "BASS" = "#e19d92",
  "BayesSpace" = "#e4c078",
  "GraphST" = "#edc400",
  "MAPLE" = "#b39a89",
  "SpaGCN" = "#7db8bf",
  "STAGATE" = "#cfcfd8"
)

get_reps <- function(scenario) {
  layer_files <- Sys.glob(file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated", paste0(scenario, "_layer_*.csv")))
  reps <- suppressWarnings(as.integer(str_match(basename(layer_files), paste0("^", scenario, "_layer_(\\d+)\\.csv$"))[, 2]))
  reps <- sort(unique(reps[!is.na(reps)]))
  reps
}

calc_ari <- function(labels_true, labels_pred) {
  keep <- !(is.na(labels_true) | is.na(labels_pred))
  labels_true <- as.character(labels_true[keep])
  labels_pred <- as.character(labels_pred[keep])

  n <- length(labels_true)
  if (n < 2) return(NA_real_)

  tab <- table(labels_true, labels_pred)
  comb2 <- function(x) x * (x - 1) / 2

  sum_nij <- sum(comb2(tab))
  sum_ai <- sum(comb2(rowSums(tab)))
  sum_bj <- sum(comb2(colSums(tab)))
  sum_n <- comb2(n)
  if (sum_n == 0) return(NA_real_)

  expected <- (sum_ai * sum_bj) / sum_n
  max_index <- 0.5 * (sum_ai + sum_bj)
  denom <- max_index - expected
  if (denom == 0) return(NA_real_)

  as.numeric((sum_nij - expected) / denom)
}

collect_ari <- function(scenario) {
  reps <- get_reps(scenario)
  if (length(reps) == 0) {
    warning("No replicate layer files found for ", scenario)
    return(data.table())
  }

  rows <- list()
  idx <- 1L
  for (m in seq_len(nrow(method_cfg))) {
    method_name <- method_cfg$method[m]
    prefix <- method_cfg$prefix[m]
    gt_col <- method_cfg$gt_col[m]
    pred_col <- method_cfg$pred_col[m]

    for (rep_id in reps) {
      f <- file.path(base_dir, "Data", paste0(scenario, "_gene"), paste0(prefix, "_", rep_id, ".csv"))
      if (!file.exists(f)) {
        next
      }
      dt <- tryCatch(fread(f), error = function(e) NULL)
      if (is.null(dt)) {
        next
      }
      if (!(gt_col %in% colnames(dt) && pred_col %in% colnames(dt))) {
        next
      }

      ari_value <- tryCatch(calc_ari(dt[[gt_col]], dt[[pred_col]]), error = function(e) NA_real_)
      if (!is.finite(ari_value)) {
        next
      }

      rows[[idx]] <- data.table(
        scenario = scenario,
        rep = rep_id,
        method = method_name,
        ARI = as.numeric(ari_value)
      )
      idx <- idx + 1L
    }
  }

  if (length(rows) == 0) {
    warning("No valid ARI rows found for ", scenario)
    return(data.table())
  }

  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, method := factor(method, levels = method_order)]
  out
}

plot_single_scenario <- function(ari_dt, scenario) {
  if (nrow(ari_dt) == 0) {
    return(invisible(NULL))
  }

  out_pdf <- file.path(base_dir, "Segmentation", paste0(scenario, "_ARI.pdf"))
  pdf(out_pdf, width = 6, height = 6)
  print(
    ggplot(ari_dt, aes(x = method, y = ARI, fill = method)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.18, size = 1.2, alpha = 0.75) +
      scale_fill_manual(values = method_colors, drop = FALSE) +
      labs(title = "ARI of Domain Segmentation Methods", x = "Method", y = "Adjusted Rand Index (ARI)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
}

plot_all_scenarios <- function(all_ari) {
  if (nrow(all_ari) == 0) {
    return(invisible(NULL))
  }

  scenario_levels <- unique(as.character(all_ari$scenario))
  all_ari[, scenario := factor(scenario, levels = scenario_levels)]
  n_scen <- length(scenario_levels)
  # Keep ARI canvas size aligned with RMSE for direct side-by-side comparison.
  plot_w <- max(14, 2.8 * n_scen)

  out_pdf <- file.path(base_dir, "Segmentation", "all_scenarios_ARI.pdf")
  pdf(out_pdf, width = plot_w, height = 4.4)
  print(
    ggplot(all_ari, aes(x = method, y = ARI, fill = method)) +
      geom_boxplot(width = 0.62, linewidth = 0.35, outlier.shape = 95, outlier.size = 1.8) +
      facet_wrap(~scenario, nrow = 1, scales = "free_x") +
      scale_fill_manual(values = method_colors, drop = FALSE) +
      scale_y_continuous(
        limits = c(-0.05, 1.02),
        breaks = c(0, 0.25, 0.5, 0.75, 1.0),
        expand = c(0.01, 0.01)
      ) +
      labs(
        x = NULL,
        y = "Adjusted Rand Index (ARI)",
        fill = "Methods"
      ) +
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
  )
  dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
scenarios <- if (length(args) > 0) args else c("S1", "S2", "S3", "S4", "S5", "S6", "S7")

all_rows <- list()
row_idx <- 1L
for (scenario in scenarios) {
  ari_dt <- collect_ari(scenario)
  plot_single_scenario(ari_dt, scenario)
  if (nrow(ari_dt) > 0) {
    all_rows[[row_idx]] <- ari_dt
    row_idx <- row_idx + 1L
  }
}

if (length(all_rows) > 0) {
  all_ari <- rbindlist(all_rows, use.names = TRUE, fill = TRUE)
  all_ari[, method := factor(method, levels = method_order)]
  ari_summary <- all_ari[, .(
    n = .N,
    mean_ARI = mean(ARI),
    median_ARI = median(ARI),
    sd_ARI = sd(ARI),
    se_ARI = sd(ARI) / sqrt(.N)
  ), by = .(scenario, method)]
  ari_summary[, `:=`(
    ci_low = mean_ARI - qt(0.975, pmax(n - 1L, 1L)) * se_ARI,
    ci_high = mean_ARI + qt(0.975, pmax(n - 1L, 1L)) * se_ARI
  )]
  fwrite(all_ari, file.path(base_dir, "Segmentation", "ARI_by_rep_S1_S7.csv"))
  fwrite(ari_summary, file.path(base_dir, "Segmentation", "ARI_summary_S1_S7.csv"))
  plot_all_scenarios(all_ari)
}
