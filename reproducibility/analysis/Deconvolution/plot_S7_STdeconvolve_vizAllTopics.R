#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(STdeconvolve)
})

base_dir <- "."
scenario <- "S7"
data_dir <- file.path(base_dir, "Data", paste0(scenario, "_gene"))
sim_dir <- file.path(data_dir, "simulated")
out_dir <- file.path(base_dir, "Deconvolution")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reps <- 1:5
combined_pdf <- file.path(out_dir, paste0(scenario, "_STdeconvolve_vizAllTopics_reps1_5.pdf"))
pdf(combined_pdf, width = 7, height = 7)

for (rep_id in reps) {
  decon_file <- file.path(data_dir, paste0("STdeconvolve_", rep_id, ".csv"))
  coord_file <- file.path(sim_dir, paste0(scenario, "_coord_", rep_id, ".csv"))
  layer_file <- file.path(sim_dir, paste0(scenario, "_layer_", rep_id, ".csv"))

  decon_dt <- fread(decon_file)
  coord_dt <- fread(coord_file)
  layer_dt <- fread(layer_file)

  n_use <- min(nrow(decon_dt), nrow(coord_dt), nrow(layer_dt))
  if (n_use <= 0) next

  deconProp <- as.data.frame(decon_dt[seq_len(n_use)])
  deconProp[] <- lapply(deconProp, as.numeric)

  pos <- as.data.frame(coord_dt[seq_len(n_use), .(x_coord, y_coord)])
  colnames(pos) <- c("x", "y")

  spot_ids <- paste0("spot_", seq_len(n_use))
  rownames(deconProp) <- spot_ids
  rownames(pos) <- spot_ids

  layer_factor <- as.numeric(as.factor(layer_dt$layer[seq_len(n_use)]))

  vizAllTopics(
    deconProp,
    pos,
    groups = layer_factor,
    group_cols = rainbow(length(unique(layer_factor))),
    r = 2,
    lwd = 0.5
  )

  rep_pdf <- file.path(out_dir, paste0(scenario, "_rep", rep_id, "_STdeconvolve_vizAllTopics.pdf"))
  pdf(rep_pdf, width = 7, height = 7)
  vizAllTopics(
    deconProp,
    pos,
    groups = layer_factor,
    group_cols = rainbow(length(unique(layer_factor))),
    r = 2,
    lwd = 0.5
  )
  dev.off()
}

dev.off()

message("Saved: ", combined_pdf)
