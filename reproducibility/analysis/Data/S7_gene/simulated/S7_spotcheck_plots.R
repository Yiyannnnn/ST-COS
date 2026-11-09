#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

ct_names <- c("Cancer", "Tcell", "Endothelial", "CAF")
input_dir <- "."
region_files <- sort(list.files(input_dir, pattern = "^S7_region_\\d+\\.csv$", full.names = TRUE))
if (length(region_files) == 0) {
  stop("No region files found in ", getwd())
}

palette_region <- c(
  background = "#e0e0e0",
  cluster_1 = "#d73027",
  cluster_2 = "#fc8d59",
  cluster_3 = "#fee090",
  cluster_4 = "#91bfdb",
  cluster_5 = "#4575b4",
  stripe_1 = "#1b9e77",
  stripe_2 = "#66a61e",
  ring_1_rim = "#d95f02",
  ring_1_core = "#7570b3",
  ring_2_rim = "#e7298a",
  ring_2_core = "#66a61e",
  ring_3_rim = "#e6ab02",
  ring_3_core = "#a6761d",
  ring_4_rim = "#1b9e77",
  ring_4_core = "#666666"
)

palette_dominant <- c(
  Cancer = "#d73027",
  Tcell = "#1b9e77",
  Endothelial = "#4575b4",
  CAF = "#bdbdbd"
)

plots <- list()
for (file in region_files) {
  data <- fread(file)
  missing_cols <- setdiff(c("x_coord", "y_coord", "region", ct_names), names(data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", file, ": ", paste(missing_cols, collapse = ", "))
  }
  scenario_id <- sub(".*_(\\d+)\\.csv$", "\\1", basename(file))
  scenario_label <- switch(
    scenario_id,
    "1" = "Cold TME",
    "2" = "Mixed TME",
    "3" = "Compartmentalized TME",
    paste0("Scenario ", scenario_id)
  )
  data[, dominant := ct_names[max.col(.SD, ties.method = "first")], .SDcols = ct_names]
  data[, dominant_prop := apply(.SD, 1, max), .SDcols = ct_names]

  p_region <- ggplot(data, aes(x = x_coord, y = y_coord, fill = region)) +
    geom_point(shape = 21, size = 3, colour = "black", stroke = 0.05) +
    scale_fill_manual(values = palette_region, na.value = "#969696") +
    coord_equal() +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    labs(title = paste0("Scenario ", scenario_id, " - ", scenario_label, " (regions)"),
         x = "x", y = "y")

  p_dominant <- ggplot(data, aes(x = x_coord, y = y_coord, fill = dominant)) +
    geom_point(shape = 21, size = 3, colour = "black", stroke = 0.05) +
    scale_fill_manual(values = palette_dominant) +
    coord_equal() +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    labs(title = paste0("Scenario ", scenario_id, " - dominant cell type"),
         x = "x", y = "y")

  p_cancer <- ggplot(data, aes(x = x_coord, y = y_coord, fill = Cancer)) +
    geom_point(shape = 21, size = 3, colour = "black", stroke = 0.05) +
    scale_fill_viridis_c(name = "Cancer prop", option = "magma") +
    coord_equal() +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right") +
    labs(title = paste0("Scenario ", scenario_id, " - Cancer proportion"),
         x = "x", y = "y")

  plot_combined <- plot_grid(p_region, p_dominant, p_cancer, ncol = 3, align = "hv")
  plots[[file]] <- list(plot = plot_combined, scenario_id = scenario_id, label = scenario_label)
}

for (entry in plots) {
  outfile <- file.path(input_dir, sprintf("S7_spotcheck_%s.pdf", entry$scenario_id))
  ggsave(outfile, entry$plot, width = 15, height = 5, units = "in", dpi = 300)
  message("Saved ", outfile)
}
