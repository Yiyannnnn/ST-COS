#!/usr/bin/env Rscript
## ---------------------------------------------------------------
## Plot gene-level simulation metrics across all 12 DLPFC slides,
## style matching all_scenarios_ARI_RMSE_combined.pdf
## Rows: Gene Mean / Gene Variance / Gene CV / Zero Proportion (blue block titles)
## Columns: 12 slides
## x-axis within each panel: method (Truth, ST-COS, SRTsim, scDesign3)
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

out_dir  <- "./Reference_based/DLPFC"
out_file <- file.path(out_dir, "DLPFC_all_slides_metrics_combined.pdf")

slides <- c("151507","151508","151509","151510",
            "151669","151670","151671","151672",
            "151673","151674","151675","151676")

method_levels <- c("Truth", "ST-COS", "SRTsim", "scDesign3", "scCube")
method_colors <- c(
  "Truth"     = "#4d4d4d",
  "ST-COS"    = "#1F77B4",
  "SRTsim"    = "#D62728",
  "scDesign3" = "#2CA02C",
  "scCube"    = "#9467BD"
)

metric_levels <- c("Gene Mean", "Gene Variance", "Gene CV", "Zero Proportion")
metric_map    <- c(
  GeneMean     = "Gene Mean",
  GeneVar      = "Gene Variance",
  GeneCV       = "Gene CV",
  GeneZeroProp = "Zero Proportion"
)

## ---- Load and combine all per-gene metrics ---------------------------
message("Loading gene metrics for all 12 slides ...")
all_dt <- rbindlist(lapply(slides, function(sl) {
  f <- file.path(out_dir, paste0(sl, "_gene_metrics.csv"))
  if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
  dt <- fread(f)
  dt[, slide := sl]
  dt
}), use.names = TRUE)

## Rename method column value "STCOS" -> "ST-COS"
all_dt[method == "STCOS", method := "ST-COS"]
all_dt[method == "scCube", method := "scCube"]
all_dt[, method := factor(method, levels = method_levels)]
all_dt[, slide  := factor(slide,  levels = slides)]

## Melt the 4 metric columns into long form
long_dt <- melt(
  all_dt[, .(gene, method, slide, GeneMean, GeneVar, GeneCV, GeneZeroProp)],
  id.vars      = c("gene", "method", "slide"),
  measure.vars = names(metric_map),
  variable.name = "metric_key",
  value.name    = "value"
)
long_dt[, metric := factor(metric_map[as.character(metric_key)],
                           levels = metric_levels)]
long_dt[, metric_key := NULL]

message(sprintf("Loaded %d rows (%d slides x %d methods x genes x 4 metrics)",
                nrow(long_dt), length(slides), length(method_levels)))

## ---- Build one-row plot per metric -----------------------------------
x_text_size <- 7.5

build_row <- function(metric_name, show_legend = FALSE) {
  plot_dt <- long_dt[metric == metric_name]

  ggplot(plot_dt, aes(x = method, y = value, fill = method)) +
    geom_boxplot(
      width = 0.65, linewidth = 0.3,
      outlier.shape = 95, outlier.size = 1.4,
      na.rm = TRUE
    ) +
    facet_wrap(~slide, nrow = 1, scales = "fixed") +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(
      title = metric_name,
      x = NULL, y = NULL,
      fill = "Method"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.background  = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.3),
      plot.background   = element_rect(fill = "white", color = NA),
      ## Blue block strip for slide labels
      strip.background  = element_rect(fill = "#d2e0e3", color = "#6f6f6f", linewidth = 0.3),
      strip.text        = element_text(size = 9, face = "bold"),
      ## Title text (actual blue block title is added externally via cowplot)
      plot.title        = element_blank(),
      plot.title.position = "plot",
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.spacing.x   = unit(0.12, "lines"),
      axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1,
                                       size = x_text_size, color = "black"),
      axis.text.y       = element_text(size = 8, color = "black"),
      axis.title        = element_blank(),
      legend.position   = if (show_legend) "right" else "none",
      legend.key        = element_rect(fill = "white", color = NA),
      legend.title      = element_text(size = 10, face = "bold"),
      legend.text       = element_text(size = 9),
      plot.margin       = margin(2, 6, 4, 6)
    )
}

## Build 4 metric rows; show legend only on last row
plots <- lapply(seq_along(metric_levels), function(i) {
  build_row(metric_levels[i], show_legend = (i == length(metric_levels)))
})

## ---- Add blue title labels using cowplot draw_label ----------------
## Wrap each plot in a cowplot canvas that adds the blue block title bar
add_blue_title <- function(p, label) {
  title_bar <- ggdraw() +
    draw_label(
      label,
      fontface = "bold", size = 11,
      x = 0.012, hjust = 0,
      color = "white"
    ) +
    theme(plot.background = element_rect(fill = "#2c6e9e", color = NA))

  plot_grid(title_bar, p,
            ncol = 1,
            rel_heights = c(0.055, 1),
            align = "v", axis = "l")
}

titled_plots <- lapply(seq_along(metric_levels), function(i) {
  add_blue_title(plots[[i]], metric_levels[i])
})

## ---- Combine into single figure -------------------------------------
## Extract shared legend from last plot
legend <- get_legend(
  build_row(metric_levels[4], show_legend = TRUE) +
    theme(legend.position = "right",
          legend.margin = margin(0, 4, 0, 8))
)

main_grid <- plot_grid(
  plotlist   = titled_plots,
  ncol       = 1,
  rel_heights = rep(1, length(metric_levels)),
  align      = "v",
  axis       = "l"
)

final_plot <- plot_grid(
  main_grid, legend,
  ncol        = 2,
  rel_widths  = c(1, 0.065)
)

## ---- Save -----------------------------------------------------------
ggsave(out_file, final_plot,
       width = 22, height = 11, bg = "white")
message("Saved: ", out_file)
