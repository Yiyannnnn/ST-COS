#!/usr/bin/env Rscript
## ---------------------------------------------------------------
## Zoomed-in violin + boxplot for Gene Mean and Gene Variance
## across all 12 DLPFC slides, log10 scale to reveal distribution
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

out_dir  <- "./Reference_based/DLPFC"
out_file <- file.path(out_dir, "DLPFC_all_slides_mean_var_zoomed.pdf")

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

## ---- Load all per-gene metrics ----------------------------------------
message("Loading gene metrics ...")
all_dt <- rbindlist(lapply(slides, function(sl) {
  f <- file.path(out_dir, paste0(sl, "_gene_metrics.csv"))
  dt <- fread(f)
  dt[, slide := sl]
  dt
}), use.names = TRUE)

all_dt[method == "STCOS", method := "ST-COS"]
all_dt[, method := factor(method, levels = method_levels)]
all_dt[, slide  := factor(slide,  levels = slides)]

## Log10-transform (shift to avoid log(0))
all_dt[, log10_mean := log10(GeneMean + 1e-3)]
all_dt[, log10_var  := log10(GeneVar  + 1e-3)]

message(sprintf("Loaded %d rows", nrow(all_dt)))

## ---- Build one row plot -----------------------------------------------
build_row <- function(dt, y_col, y_label, title_label,
                      y_limits = NULL, show_legend = FALSE) {

  p <- ggplot(dt, aes(x = method, y = .data[[y_col]], fill = method)) +
    geom_violin(
      trim   = TRUE,
      scale  = "width",
      alpha  = 0.75,
      linewidth = 0.25
    ) +
    geom_boxplot(
      width          = 0.18,
      linewidth      = 0.35,
      outlier.shape  = NA,       # hide outlier dots — violin already shows tails
      fill           = "white",
      alpha          = 0.7,
      color          = "black"
    ) +
    stat_summary(
      fun = median, geom = "point",
      shape = 18, size = 1.6, color = "black"
    ) +
    facet_wrap(~slide, nrow = 1, scales = "fixed") +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    coord_cartesian(ylim = y_limits) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = NULL, y = y_label, fill = "Method") +
    theme_bw(base_size = 10) +
    theme(
      panel.background  = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.3),
      plot.background   = element_rect(fill = "white", color = NA),
      strip.background  = element_rect(fill = "#d2e0e3", color = "#6f6f6f", linewidth = 0.3),
      strip.text        = element_text(size = 9, face = "bold"),
      plot.title        = element_blank(),
      panel.grid.major.y = element_line(color = "#ebebeb", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing.x    = unit(0.12, "lines"),
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = 7.2, color = "black"),
      axis.text.y  = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 9.5, face = "bold"),
      axis.title.x = element_blank(),
      legend.position  = if (show_legend) "right" else "none",
      legend.key       = element_rect(fill = "white", color = NA),
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 9),
      plot.margin      = margin(2, 6, 4, 6)
    )
  p
}

## ---- Determine y-axis limits (1st–99th percentile across all methods/slides)
q_mean <- quantile(all_dt$log10_mean, probs = c(0.005, 0.995), na.rm = TRUE)
q_var  <- quantile(all_dt$log10_var,  probs = c(0.005, 0.995), na.rm = TRUE)

p_mean <- build_row(all_dt, "log10_mean",
                    y_label = expression(log[10](Gene~Mean + 0.001)),
                    title_label = "Gene Mean",
                    y_limits = q_mean,
                    show_legend = FALSE)

p_var  <- build_row(all_dt, "log10_var",
                    y_label = expression(log[10](Gene~Variance + 0.001)),
                    title_label = "Gene Variance",
                    y_limits = q_var,
                    show_legend = TRUE)

## ---- Add blue title bars (same style as combined plot) ----------------
add_blue_title <- function(p, label) {
  title_bar <- ggdraw() +
    draw_label(label, fontface = "bold", size = 11,
               x = 0.012, hjust = 0, color = "white") +
    theme(plot.background = element_rect(fill = "#2c6e9e", color = NA))
  plot_grid(title_bar, p, ncol = 1,
            rel_heights = c(0.07, 1), align = "v", axis = "l")
}

row_mean <- add_blue_title(p_mean, "Gene Mean  [log10 scale]")
row_var  <- add_blue_title(p_var,  "Gene Variance  [log10 scale]")

## ---- Extract shared legend --------------------------------------------
legend <- get_legend(
  p_var + theme(legend.position = "right",
                legend.margin = margin(0, 4, 0, 8))
)

## ---- Assemble final figure --------------------------------------------
main_grid <- plot_grid(
  row_mean, row_var,
  ncol        = 1,
  rel_heights = c(1, 1),
  align = "v", axis = "l"
)

final_plot <- plot_grid(
  main_grid, legend,
  ncol       = 2,
  rel_widths = c(1, 0.065)
)

## ---- Save -------------------------------------------------------------
ggsave(out_file, final_plot, width = 22, height = 8, bg = "white")
message("Saved: ", out_file)
