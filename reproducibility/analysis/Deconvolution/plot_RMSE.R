library(data.table)
library(ggplot2)
library(stringr)

base_dir <- "."
default_scenarios <- paste0("S", 1:7)

method_order <- c("CARD", "RCTD", "STdeconvolve", "STdGCN", "Tangram", "POLARIS")
method_colors <- c(
  "CARD" = "#e19d92",
  "RCTD" = "#eac47c",
  "STdeconvolve" = "#7e8cad",
  "STdGCN" = "#a5c2cd",
  "Tangram" = "#dadae4",
  "POLARIS" = "#9ccf9b"
)

get_reps <- function(scenario) {
  sim_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated")
  patt <- paste0("^", scenario, "_gene_expression_([0-9]+)\\.csv$")
  files <- list.files(sim_dir, pattern = patt)
  reps <- suppressWarnings(as.integer(str_match(files, patt)[, 2]))
  reps <- sort(unique(reps[!is.na(reps)]))
  reps
}

build_truth <- function(scenario, rep_id) {
  sim_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"), "simulated")
  region_file <- file.path(sim_dir, paste0(scenario, "_region_", rep_id, ".csv"))
  layer_file <- file.path(sim_dir, paste0(scenario, "_layer_", rep_id, ".csv"))

  if (file.exists(region_file)) {
    region_dt <- fread(region_file)
    ct_cols <- setdiff(colnames(region_dt), c("x_coord", "y_coord", "region"))
    truth <- as.data.frame(region_dt[, ..ct_cols])
  } else if (file.exists(layer_file)) {
    layer_dt <- fread(layer_file)
    if (scenario %in% c("S1", "S4")) {
      truth <- data.frame(
        ct1 = rep(1 / 3, nrow(layer_dt)),
        ct2 = rep(1 / 3, nrow(layer_dt)),
        ct3 = rep(1 / 3, nrow(layer_dt))
      )
    } else {
      truth <- data.frame(
        ct1 = ifelse(layer_dt$layer == 1, 1 / 2, ifelse(layer_dt$layer == 2, 0, 0)),
        ct2 = ifelse(layer_dt$layer == 1, 1 / 4, ifelse(layer_dt$layer == 2, 3 / 4, 0)),
        ct3 = ifelse(layer_dt$layer == 1, 1 / 4, ifelse(layer_dt$layer == 2, 1 / 4, 1))
      )
    }
  } else {
    stop(sprintf("No truth file found for %s rep %d.", scenario, rep_id))
  }

  rownames(truth) <- paste0("spot_", seq_len(nrow(truth)))
  truth
}

read_pred <- function(scenario, rep_id, method) {
  data_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"))
  file_path <- switch(
    method,
    "CARD" = file.path(data_dir, paste0("CARD_", rep_id, ".csv")),
    "RCTD" = file.path(data_dir, paste0("RCTD_", rep_id, ".csv")),
    "Redeconve" = file.path(data_dir, paste0("Redeconve_", rep_id, ".csv")),
    "STdeconvolve" = file.path(data_dir, paste0("STdeconvolve_", rep_id, ".csv")),
    "STdGCN" = file.path(data_dir, paste0("STdGCN_input_", rep_id), "output", "predict_result.csv"),
    "Tangram" = file.path(data_dir, paste0("Tangram_", rep_id, ".csv")),
    "POLARIS" = file.path(data_dir, paste0("POLARIS_", rep_id, ".csv")),
    stop("Unknown method: ", method)
  )

  if (!file.exists(file_path)) {
    return(NULL)
  }

  dt <- fread(file_path)
  spot_ids <- NULL

  if ("spot" %in% names(dt)) {
    spot_ids <- as.character(dt$spot)
    dt[, spot := NULL]
  }
  if ("barcodes" %in% names(dt)) {
    spot_ids <- as.character(dt$barcodes)
    dt[, barcodes := NULL]
  }
  if ("V1" %in% names(dt) && !is.numeric(dt$V1)) {
    spot_ids <- as.character(dt$V1)
    dt[, V1 := NULL]
  }
  if ("V1" %in% names(dt) && method == "Tangram") {
    dt[, V1 := NULL]
  }

  pred <- as.data.frame(dt)
  pred <- pred[, sapply(pred, is.numeric), drop = FALSE]
  if (ncol(pred) == 0) {
    return(NULL)
  }

  if (method == "STdeconvolve" && ncol(pred) == 3) {
    colnames(pred) <- c("ct1", "ct2", "ct3")
  }

  if (!is.null(spot_ids)) {
    rownames(pred) <- spot_ids
  } else {
    rownames(pred) <- paste0("spot_", seq_len(nrow(pred)))
  }
  pred
}

align_truth_to_pred <- function(truth, pred_cols) {
  # Direct overlap first
  if (all(pred_cols %in% colnames(truth))) {
    return(truth[, pred_cols, drop = FALSE])
  }

  # S7 truth (Cancer/Tcell/Endothelial/CAF) vs 3-component methods (ct1-ct3)
  if (all(c("ct1", "ct2", "ct3") %in% pred_cols) &&
      all(c("Cancer", "Tcell", "Endothelial") %in% colnames(truth))) {
    mapped <- data.frame(
      ct1 = truth$Cancer,
      ct2 = truth$Tcell,
      ct3 = truth$Endothelial,
      row.names = rownames(truth)
    )
    return(mapped[, pred_cols, drop = FALSE])
  }

  # Legacy truth (ct1-ct3) vs 4-component predictions
  if (all(c("Cancer", "Tcell", "Endothelial", "CAF") %in% pred_cols) &&
      all(c("ct1", "ct2", "ct3") %in% colnames(truth))) {
    mapped <- data.frame(
      Cancer = truth$ct1,
      Tcell = truth$ct2,
      Endothelial = truth$ct3,
      CAF = 0,
      row.names = rownames(truth)
    )
    return(mapped[, pred_cols, drop = FALSE])
  }

  common <- intersect(pred_cols, colnames(truth))
  if (length(common) == 0) {
    return(NULL)
  }
  truth[, common, drop = FALSE]
}

calc_rmse <- function(truth, pred) {
  if (is.null(pred) || ncol(pred) == 0) {
    return(NA_real_)
  }

  truth_aligned <- align_truth_to_pred(truth, colnames(pred))
  if (is.null(truth_aligned) || ncol(truth_aligned) == 0) {
    return(NA_real_)
  }

  common_spots <- intersect(rownames(truth_aligned), rownames(pred))
  if (length(common_spots) == 0) {
    return(NA_real_)
  }

  truth_use <- truth_aligned[common_spots, , drop = FALSE]
  pred_use <- pred[common_spots, colnames(truth_use), drop = FALSE]

  # Handle occasional one-row mismatch in methods without spot IDs.
  if (nrow(pred_use) != nrow(truth_use)) {
    n <- min(nrow(pred_use), nrow(truth_use))
    truth_use <- truth_use[seq_len(n), , drop = FALSE]
    pred_use <- pred_use[seq_len(n), , drop = FALSE]
  }

  sqrt(mean((as.matrix(truth_use) - as.matrix(pred_use))^2, na.rm = TRUE))
}

collect_rmse <- function(scenarios) {
  rows <- list()
  idx <- 1L

  for (scenario in scenarios) {
    reps <- get_reps(scenario)
    if (length(reps) == 0) next

    for (rep_id in reps) {
      truth <- build_truth(scenario, rep_id)
      for (method in method_order) {
        pred <- read_pred(scenario, rep_id, method)
        rmse <- calc_rmse(truth, pred)
        rows[[idx]] <- data.table(
          scenario = scenario,
          rep = rep_id,
          method = method,
          RMSE = rmse,
          available = !is.null(pred)
        )
        idx <- idx + 1L
      }
    }
  }

  if (length(rows) == 0) return(data.table())
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, method := factor(method, levels = method_order)]
  out
}

plot_rmse <- function(rmse_dt, scenarios) {
  out_pdf <- file.path(base_dir, "Deconvolution", "all_scenarios_RMSE_S1_S7.pdf")
  rmse_dt[, scenario := factor(scenario, levels = scenarios)]

  plot_w <- max(14, 2.8 * length(scenarios))
  pdf(out_pdf, width = plot_w, height = 4.4)
  print(
    ggplot(rmse_dt, aes(x = method, y = RMSE, fill = method)) +
      geom_boxplot(outlier.shape = NA, width = 0.62, na.rm = TRUE) +
      facet_wrap(~scenario, nrow = 1, scales = "free_x", drop = FALSE) +
      scale_fill_manual(values = method_colors, drop = FALSE) +
      labs(x = NULL, y = "Root Mean Square Error (RMSE)", fill = "Methods") +
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
scenarios <- if (length(args) > 0) args else default_scenarios

rmse_dt <- collect_rmse(scenarios)
if (nrow(rmse_dt) == 0) {
  stop("No deconvolution outputs found for requested scenarios.")
}

summary_dt <- rmse_dt[, .(
  n_total = .N,
  n_available = sum(available, na.rm = TRUE),
  n_rmse = sum(!is.na(RMSE)),
  mean_RMSE = mean(RMSE, na.rm = TRUE),
  median_RMSE = median(RMSE, na.rm = TRUE),
  sd_RMSE = sd(RMSE, na.rm = TRUE),
  se_RMSE = sd(RMSE, na.rm = TRUE) / sqrt(sum(!is.na(RMSE)))
), by = .(scenario, method)]
summary_dt[, `:=`(
  ci_low = mean_RMSE - qt(0.975, pmax(n_rmse - 1L, 1L)) * se_RMSE,
  ci_high = mean_RMSE + qt(0.975, pmax(n_rmse - 1L, 1L)) * se_RMSE
)]

fwrite(rmse_dt, file.path(base_dir, "Deconvolution", "RMSE_by_rep_S1_S7.csv"))
fwrite(summary_dt, file.path(base_dir, "Deconvolution", "RMSE_summary_S1_S7.csv"))

plot_rmse(rmse_dt, scenarios)
