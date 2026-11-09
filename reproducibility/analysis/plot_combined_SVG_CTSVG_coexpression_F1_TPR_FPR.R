#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

base_dir        <- "."
scenario_levels <- paste0("S", 1:7)

## ---- Method levels & colors ------------------------------------------
svg_method_levels <- c("moranI", "HEARTSVG", "SPARK", "SPARKX", "SpatialDE", "SpatialDE2")
svg_method_colors <- c(
  "moranI"     = "#d98b5f",
  "HEARTSVG"   = "#a7c5bd",
  "SPARK"      = "#d9b44a",
  "SPARKX"     = "#8ca6db",
  "SpatialDE"  = "#c7b8a3",
  "SpatialDE2" = "#b8d8b8"
)

ctsvg_method_levels <- c(
  "CELINA", "CELINA_trueprop", "CSIDE", "CSIDE_trueprop",
  "CTSV_trueprop", "STANCE_trueprop", "spVC_trueprop"
)
ctsvg_method_colors <- c(
  "CELINA"          = "#E69F00",
  "CELINA_trueprop" = "#F6C667",
  "CSIDE"           = "#56B4E9",
  "CSIDE_trueprop"  = "#0072B2",
  "CTSV_trueprop"   = "#CC79A7",
  "STANCE_trueprop" = "#999999",
  "spVC_trueprop"   = "#009E73"
)

coexp_method_levels <- c("spMOCA", "SpaceX", "CSSN", "hdWGCNA", "Smoothie", "binSpect", "stIHC")
coexp_method_colors <- c(
  "spMOCA"   = "#4B9CD3",
  "SpaceX"   = "#4c9390",
  "CSSN"     = "#d65d48",
  "hdWGCNA"  = "#8b5a96",
  "Smoothie" = "#f4a261",
  "binSpect" = "#2a9d8f",
  "stIHC"    = "#264653"
)

## ---- Data loading ----------------------------------------------------
load_svg_metrics <- function() {
  f  <- file.path(base_dir, "SVG_detection", "SVG_metrics_by_rep_S1_S2_S3_S4_S5_S6_S7.csv")
  dt <- fread(f)
  long_dt <- melt(
    dt[, .(scenario, rep, method, has_positive_truth, F1_score, TPR, FPR)],
    id.vars      = c("scenario", "rep", "method", "has_positive_truth"),
    variable.name = "metric", value.name = "value"
  )
  long_dt[, source := "SVG"]
  long_dt
}

load_ctsvg_metrics <- function() {
  f  <- file.path(base_dir, "CT_SVG", "CTSVG_metrics_by_rep_S1_S2_S3_S4_S5_S6_S7.csv")
  dt <- fread(f)
  setnames(dt,
    old = c("macro_F1", "macro_TPR", "macro_FPR"),
    new = c("F1_score", "TPR",       "FPR"),
    skip_absent = TRUE
  )
  long_dt <- melt(
    dt[, .(scenario, rep, method, has_positive_truth, F1_score, TPR, FPR)],
    id.vars      = c("scenario", "rep", "method", "has_positive_truth"),
    variable.name = "metric", value.name = "value"
  )
  long_dt[, source := "CT-SVG"]
  long_dt
}

load_coexp_metrics <- function() {
  f  <- file.path(base_dir, "coexpression", "coexpression_method_metrics.csv")
  dt <- fread(f)
  if (!("has_positive_truth" %in% names(dt)))
    dt[, has_positive_truth := !(scenario %in% c("S1", "S2", "S3"))]
  setnames(dt,
    old = c("replicate", "F1"),
    new = c("rep",       "F1_score"),
    skip_absent = TRUE
  )
  long_dt <- melt(
    dt[, .(scenario, rep, method, has_positive_truth, F1_score, TPR, FPR)],
    id.vars      = c("scenario", "rep", "method", "has_positive_truth"),
    variable.name = "metric", value.name = "value"
  )
  long_dt[, source := "Coexpression"]
  long_dt[!has_positive_truth & metric %in% c("F1_score", "TPR"), value := NA_real_]
  long_dt
}

## ---- Cross-out diagonal for missing-truth scenarios ------------------
build_crossout <- function(dt, method_levels) {
  annot <- dt[, .(all_na = all(is.na(value))), by = scenario][all_na == TRUE]
  if (nrow(annot) == 0) return(NULL)
  annot[, `:=`(x = 0.6, xend = length(method_levels) + 0.4, y = 0.04, yend = 0.96)]
  annot
}

## ---- Core plot builder -----------------------------------------------
p1_theme <- theme_bw(base_size = 10) + theme(
  panel.background  = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.35),
  plot.background   = element_rect(fill = "white", color = NA),
  strip.background  = element_rect(fill = "#d2e0e3", color = "#6f6f6f", linewidth = 0.35),
  strip.text        = element_text(size = 10, face = "bold"),
  plot.title        = element_text(size = 10.5, face = "bold", hjust = 0, margin = margin(b = 2)),
  plot.title.position = "plot",
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  panel.spacing.x   = unit(0.14, "lines"),
  axis.text.y       = element_text(size = 8.6, color = "black"),
  axis.title.y      = element_text(size = 10.5, face = "bold"),
  axis.title.x      = element_blank(),
  legend.key        = element_rect(fill = "white", color = NA),
  legend.background = element_rect(fill = "white", color = NA),
  legend.title      = element_text(size = 10.5, face = "bold"),
  legend.text       = element_text(size = 9.3),
  plot.margin       = margin(6, 10, 6, 8)
)

build_panel <- function(dt, method_levels, method_colors,
                        metric_key, y_label, title = NULL,
                        show_legend = FALSE, show_strip = TRUE) {
  pdt <- copy(dt[metric == metric_key])
  pdt[, scenario := factor(scenario, levels = scenario_levels)]
  pdt[, method   := factor(method,   levels = method_levels)]
  annot <- build_crossout(pdt, method_levels)
  x_sz  <- if (length(method_levels) >= 7) 6.6 else 8.2

  p <- ggplot(pdt, aes(x = method, y = value, fill = method)) +
    geom_boxplot(width = 0.62, linewidth = 0.35,
                 outlier.shape = 95, outlier.size = 1.6, na.rm = TRUE) +
    facet_wrap(~scenario, nrow = 1, scales = "free_x", drop = FALSE) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    scale_y_continuous(limits = c(0, 1.02),
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       expand = c(0.01, 0.01)) +
    labs(x = NULL, y = y_label, fill = "Methods",
         title = if (!is.null(title)) title else NULL) +
    p1_theme +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1,
                                     size = x_sz, color = "black"),
      legend.position = if (show_legend) "right" else "none"
    )

  if (!show_strip)
    p <- p + theme(strip.text = element_blank(),
                   strip.background = element_blank())

  if (!is.null(annot))
    p <- p + geom_segment(
      data = annot,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE, color = "#666666",
      linewidth = 0.75, lineend = "round"
    )
  p
}

## ---- Load data -------------------------------------------------------
message("Loading data ...")
svg_dt   <- load_svg_metrics()
ctsvg_dt <- load_ctsvg_metrics()
coexp_dt <- load_coexp_metrics()

## ---- Build 9 panels (3 sources × 3 metrics) -------------------------
## Layout: rows = SVG / CT-SVG / Coexpression   cols = F1 / TPR / FPR
message("Building panels ...")

metrics <- list(
  list(key = "F1_score", label = "F1 Score"),
  list(key = "TPR",      label = "TPR"),
  list(key = "FPR",      label = "FPR")
)

sources <- list(
  list(dt = svg_dt[method   %in% svg_method_levels],   lvl = svg_method_levels,   col = svg_method_colors,   title = "SVG"),
  list(dt = ctsvg_dt[method %in% ctsvg_method_levels], lvl = ctsvg_method_levels, col = ctsvg_method_colors, title = "CT-SVG"),
  list(dt = coexp_dt[method %in% coexp_method_levels], lvl = coexp_method_levels, col = coexp_method_colors, title = "Coexpression")
)

panels <- list()
for (si in seq_along(sources)) {
  for (mi in seq_along(metrics)) {
    show_title  <- (si == 1)                    # column header (metric name) on top row only
    show_strip  <- (si == 1)                    # scenario strip on top row only
    show_legend <- (mi == length(metrics))      # legend on right column only
    y_lbl <- if (mi == 1) sources[[si]]$title else NULL   # row label on left column
    ttl   <- if (show_title) metrics[[mi]]$label else NULL

    panels[[length(panels) + 1]] <- build_panel(
      dt            = sources[[si]]$dt,
      method_levels = sources[[si]]$lvl,
      method_colors = sources[[si]]$col,
      metric_key    = metrics[[mi]]$key,
      y_label       = y_lbl,
      title         = ttl,
      show_legend   = show_legend,
      show_strip    = show_strip
    )
  }
}

## ---- Extract shared legend per source row ----------------------------
leg_svg   <- get_legend(build_panel(sources[[1]]$dt, sources[[1]]$lvl, sources[[1]]$col,
                                    "F1_score", NULL, show_legend = TRUE))
leg_ctsvg <- get_legend(build_panel(sources[[2]]$dt, sources[[2]]$lvl, sources[[2]]$col,
                                    "F1_score", NULL, show_legend = TRUE))
leg_coexp <- get_legend(build_panel(sources[[3]]$dt, sources[[3]]$lvl, sources[[3]]$col,
                                    "F1_score", NULL, show_legend = TRUE))

legend_row <- plot_grid(leg_svg, leg_ctsvg, leg_coexp,
                        ncol = 3, rel_widths = c(1, 1, 1))

## ---- Remove inline legends, arrange 3×3 grid ------------------------
panels_no_leg <- lapply(panels, function(p) p + theme(legend.position = "none"))

grid_3x3 <- plot_grid(
  plotlist   = panels_no_leg,
  nrow       = 3, ncol = 3,
  align      = "hv", axis = "tblr",
  rel_widths = c(1, 1, 1)
)

## ---- Footer note -----------------------------------------------------
footer <- ggdraw() +
  draw_label(
    paste0("Diagonal strike: no positive truth for that scenario/task. ",
           "SVG & CT-SVG in S2/S5; Coexpression in S1/S2/S3."),
    x = 0.01, hjust = 0, size = 8, color = "#444444"
  )

## ---- Assemble --------------------------------------------------------
final_plot <- plot_grid(
  grid_3x3,
  legend_row,
  footer,
  ncol        = 1,
  rel_heights = c(1, 0.13, 0.03),
  align = "v", axis = "l"
)

out_file <- file.path(base_dir, "combined_SVG_CTSVG_coexpression_F1_TPR_FPR.pdf")
ggsave(out_file, final_plot, width = 20, height = 11, bg = "white")
message("Saved: ", out_file)

out_png <- file.path(base_dir, "combined_SVG_CTSVG_coexpression_F1_TPR_FPR.png")
ggsave(out_png, final_plot, width = 20, height = 11, dpi = 300, bg = "white")
message("Saved: ", out_png)

f1_rows <- lapply(seq_along(sources), function(si) {
  build_panel(
    dt = sources[[si]]$dt,
    method_levels = sources[[si]]$lvl,
    method_colors = sources[[si]]$col,
    metric_key = "F1_score",
    y_label = "F1 Score",
    title = sources[[si]]$title,
    show_legend = TRUE,
    show_strip = TRUE
  )
})
f1_footer <- ggdraw() +
  draw_label(
    paste0("Diagonal strike: no positive truth for that scenario/task. ",
           "SVG and CT-SVG in S2/S5; coexpression in S1/S2/S3."),
    x = 0.01, hjust = 0, size = 8, color = "#444444"
  )
f1_vertical <- plot_grid(
  plotlist = c(f1_rows, list(f1_footer)), ncol = 1,
  rel_heights = c(1, 1, 1, 0.06), align = "v", axis = "l"
)

f1_pdf <- file.path(base_dir, "combined_SVG_CTSVG_coexpression_F1_vertical.pdf")
f1_png <- file.path(base_dir, "combined_SVG_CTSVG_coexpression_F1_vertical.png")
ggsave(f1_pdf, f1_vertical, width = 20, height = 12, bg = "white")
ggsave(f1_png, f1_vertical, width = 20, height = 12, dpi = 300, bg = "white")
message("Saved: ", f1_pdf)
message("Saved: ", f1_png)
