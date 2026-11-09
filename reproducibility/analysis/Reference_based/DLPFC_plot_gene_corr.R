#!/usr/bin/env Rscript
## ---------------------------------------------------------------
## Gene-gene correlation evaluation plots across all 12 DLPFC slides.
## Panel A: RMSE and MAE bar chart (per method, per slide)
## Panel B: Scatter — truth corr (x) vs simulated corr (y),
##          pooled across all slides, one column per method
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

out_dir  <- "./Reference_based/DLPFC"
out_file <- file.path(out_dir, "DLPFC_gene_corr_eval.pdf")

slides <- c("151507","151508","151509","151510",
            "151669","151670","151671","151672",
            "151673","151674","151675","151676")

method_levels <- c("STCOS", "SRTsim", "scDesign3", "scCube")
method_labels <- c("STCOS" = "ST-COS", "SRTsim" = "SRTsim", "scDesign3" = "scDesign3", "scCube" = "scCube")
method_colors <- c("ST-COS" = "#1F77B4", "SRTsim" = "#D62728", "scDesign3" = "#2CA02C", "scCube" = "#9467BD")

## ---- Load summary metrics --------------------------------------------
message("Loading summary metrics ...")
summary_dt <- rbindlist(lapply(slides, function(sl) {
  f <- file.path(out_dir, paste0(sl, "_gene_corr_summary.csv"))
  if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
  fread(f)
}), use.names = TRUE)

summary_dt[, method_label := method_labels[method]]
summary_dt[, method_label := factor(method_label, levels = names(method_colors))]
summary_dt[, slide := factor(slide, levels = slides)]

## ---- Load sampled pairs for scatter ----------------------------------
message("Loading sampled pairs ...")
pairs_dt <- rbindlist(lapply(slides, function(sl) {
  f <- file.path(out_dir, paste0(sl, "_gene_corr_pairs.csv.gz"))
  if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
  dt <- fread(f)
  dt[, slide := sl]
  dt
}), use.names = TRUE)

## Melt to long for scatter
pairs_long <- melt(
  pairs_dt,
  id.vars       = c("slide", "truth"),
  measure.vars  = c("STCOS", "SRTsim", "scDesign3"),
  variable.name = "method",
  value.name    = "simulated"
)
pairs_long[, method_label := factor(method_labels[as.character(method)],
                                    levels = names(method_colors))]

# Subsample for plotting (max 100K per method across all slides)
set.seed(42)
pairs_plot <- pairs_long[, .SD[sample(.N, min(.N, 100000))], by = method_label]

message(sprintf("Loaded %d slides, %d sampled pairs total", length(slides), nrow(pairs_long)))

## =========================================================================
## PANEL A: RMSE and MAE — faceted bar chart (slides x method)
## =========================================================================
make_bar_row <- function(metric_col, y_label, show_legend = FALSE) {
  ggplot(summary_dt, aes(x = method_label, y = .data[[metric_col]], fill = method_label)) +
    geom_col(width = 0.65) +
    facet_wrap(~slide, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(x = NULL, y = y_label, fill = "Method") +
    theme_bw(base_size = 10) +
    theme(
      panel.background  = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.3),
      plot.background   = element_rect(fill = "white", color = NA),
      strip.background  = element_rect(fill = "#d2e0e3", color = "#6f6f6f", linewidth = 0.3),
      strip.text        = element_text(size = 8.5, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "#ebebeb", linewidth = 0.3),
      panel.spacing.x    = unit(0.12, "lines"),
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 7, color = "black"),
      axis.text.y  = element_text(size = 7.5, color = "black"),
      axis.title.y = element_text(size = 9.5, face = "bold"),
      legend.position  = if (show_legend) "right" else "none",
      legend.key       = element_rect(fill = "white", color = NA),
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 9),
      plot.margin      = margin(2, 6, 4, 6)
    )
}

p_rmse <- make_bar_row("RMSE", "RMSE")
p_mae  <- make_bar_row("MAE",  "MAE", show_legend = TRUE)

add_blue_title <- function(p, label) {
  bar <- ggdraw() +
    draw_label(label, fontface = "bold", size = 11,
               x = 0.012, hjust = 0, color = "white") +
    theme(plot.background = element_rect(fill = "#2c6e9e", color = NA))
  plot_grid(bar, p, ncol = 1, rel_heights = c(0.065, 1), align = "v", axis = "l")
}

row_rmse <- add_blue_title(p_rmse, "RMSE  (gene-gene correlation)")
row_mae  <- add_blue_title(p_mae,  "MAE   (gene-gene correlation)")

legend_bar <- get_legend(
  make_bar_row("RMSE", "", show_legend = TRUE) +
    theme(legend.position = "right", legend.margin = margin(0, 4, 0, 8))
)

bar_grid <- plot_grid(row_rmse, row_mae, ncol = 1, rel_heights = c(1, 1),
                      align = "v", axis = "l")
panel_A <- plot_grid(bar_grid, legend_bar, ncol = 2, rel_widths = c(1, 0.065))

## =========================================================================
## PANEL B: Scatter — truth (x) vs simulated (y), one panel per method
## =========================================================================

# Compute per-method Pearson r for annotation
r_annot <- summary_dt[, .(r_mean = mean(PearsonR), r_sd = sd(PearsonR)),
                       by = method_label]

scatter_panel <- ggplot(pairs_plot,
                        aes(x = truth, y = simulated, color = method_label)) +
  geom_point(alpha = 0.08, size = 0.3, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = "black", linetype = "solid") +
  facet_wrap(~method_label, nrow = 1) +
  scale_color_manual(values = method_colors, guide = "none") +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    x = "Truth gene-gene correlation",
    y = "Simulated gene-gene correlation"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.background  = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.3),
    plot.background   = element_rect(fill = "white", color = NA),
    strip.background  = element_rect(fill = "#d2e0e3", color = "#6f6f6f", linewidth = 0.3),
    strip.text        = element_text(size = 10, face = "bold"),
    panel.grid.major  = element_line(color = "#ebebeb", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(size = 8.5, color = "black"),
    axis.title        = element_text(size = 10, face = "bold"),
    plot.margin       = margin(4, 6, 4, 6)
  )

# Add per-panel Pearson r annotation
r_annot[, label := sprintf("r = %.3f +/- %.3f", r_mean, r_sd)]
scatter_annot <- scatter_panel +
  geom_text(
    data  = r_annot,
    aes(x = -0.95, y = 0.92, label = label),
    inherit.aes = FALSE,
    hjust = 0, size = 3, fontface = "italic", color = "black"
  )

panel_B <- add_blue_title(scatter_annot,
  "Gene-gene correlation scatter  (sampled pairs, pooled across 12 slides)")

## =========================================================================
## Assemble final figure
## =========================================================================
final_plot <- plot_grid(
  panel_A, panel_B,
  ncol        = 1,
  rel_heights = c(1.1, 0.9),
  align = "v", axis = "l"
)

ggsave(out_file, final_plot, width = 22, height = 14, bg = "white")
message("Saved: ", out_file)
