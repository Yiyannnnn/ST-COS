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

# Truth, column alignment, and RMSE follow Supplementary Methods S2 and are
# defined in evaluation_rules.R: named outputs are
# matched to simulated cell types by name, unnamed STdeconvolve topics by an
# optimal one-to-one assignment, and RMSE averages over every truth cell type.
build_truth <- function(scenario, rep_id) {
  truth_proportions(sim_dir_for(base_dir, scenario), scenario, rep_id)
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
        pred <- read_deconv_prediction(base_dir, scenario, rep_id, method)
        res <- deconv_rmse(truth, pred, method)
        rows[[idx]] <- data.table(
          scenario = scenario,
          rep = rep_id,
          method = method,
          RMSE = res$RMSE,
          available = !is.null(pred),
          n_spots = res$n_spots,
          status = res$status
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
  mean_RMSE = mean_or_na(RMSE),
  median_RMSE = if (all(is.na(RMSE))) NA_real_ else median(RMSE, na.rm = TRUE),
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
