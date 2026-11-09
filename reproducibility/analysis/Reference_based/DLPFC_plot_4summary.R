#!/usr/bin/env Rscript
## ---------------------------------------------------------------
## 4 clean summary plots — all slides pooled, no slide facets.
##
## Plot 1  DLPFC_summary_1_metrics_violin.pdf
##   Gene Mean / Var / CV / Zero-Prop  violin+box, genes pooled across 12 slides
## Plot 2  DLPFC_summary_2_metrics_deviation.pdf
##   % deviation from Truth per method, each dot = 1 slide (12 dots),
##   boxplot across slides  — 4 metric panels
## Plot 3  DLPFC_summary_3_gene_corr_scatter.pdf
##   Truth vs simulated gene-gene corr scatter, pooled across 12 slides
##   1 panel per method (4 panels)
## Plot 4  DLPFC_summary_4_gene_corr_rmse.pdf
##   RMSE and Pearson r per method, each dot = 1 slide, 2 panels
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

out_dir <- "./Reference_based/DLPFC"
slides  <- c("151507","151508","151509","151510",
             "151669","151670","151671","151672",
             "151673","151674","151675","151676")

method_levels <- c("Truth","ST-COS","SRTsim","scDesign3","scCube")
method_colors <- c(
  "Truth"     = "#4d4d4d",
  "ST-COS"    = "#1F77B4",
  "SRTsim"    = "#D62728",
  "scDesign3" = "#2CA02C",
  "scCube"    = "#9467BD"
)
sim_methods <- setdiff(method_levels, "Truth")

# Correlation plots evaluate realized simulated counts for every method.
corr_methods <- c("ST-COS","SRTsim","scDesign3","scCube")

base_theme <- theme_bw(base_size = 12) + theme(
  panel.background  = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.35),
  plot.background   = element_rect(fill = "white", color = NA),
  strip.background  = element_rect(fill = "#2c6e9e", color = NA),
  strip.text        = element_text(color = "white", face = "bold", size = 11),
  panel.grid.major.x = element_blank(),
  panel.grid.minor   = element_blank(),
  panel.grid.major.y = element_line(color = "#ebebeb", linewidth = 0.3),
  legend.position   = "none",
  axis.text.x       = element_text(angle = 40, hjust = 1, size = 10, color = "black"),
  axis.text.y       = element_text(size = 9, color = "black"),
  axis.title        = element_text(size = 11, face = "bold"),
  plot.margin       = margin(4, 6, 4, 6)
)

add_legend <- function(p) {
  p + theme(
    legend.position = "right",
    legend.key      = element_rect(fill = "white", color = NA),
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10)
  )
}

signif_p_col <- "paired_t_p_BH"

p_to_stars <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "ns")
      )
    )
  )
}

load_signif_tests <- function(path) {
  if (!file.exists(path)) {
    warning("Missing significance table: ", path)
    return(data.table())
  }
  fread(path)
}

make_stcos_annotations <- function(test_dt, metric_name, score_name, methods,
                                   y_limits = NULL, y_values = NULL,
                                   panel_col = NULL, panel_value = NULL,
                                   p_col = signif_p_col) {
  if (is.null(test_dt) || nrow(test_dt) == 0 || !p_col %in% names(test_dt)) {
    return(data.table())
  }

  ann <- copy(test_dt[
    metric == metric_name & score == score_name &
      comparator %in% setdiff(methods, "ST-COS")
  ])
  if (nrow(ann) == 0) return(data.table())

  if (is.null(y_limits)) {
    vals <- y_values[is.finite(y_values)]
    y_min <- min(vals, na.rm = TRUE)
    y_max <- max(vals, na.rm = TRUE)
    if (!is.finite(y_min)) y_min <- 0
    if (!is.finite(y_max)) y_max <- 1
  } else {
    y_min <- y_limits[1]
    y_max <- y_limits[2]
  }
  rng <- y_max - y_min
  if (!is.finite(rng) || rng <= 0) rng <- max(abs(y_max), 1)

  ann[, p_value := get(p_col)]
  ann[, label := p_to_stars(p_value)]
  ann[, x := match("ST-COS", methods)]
  ann[, xend := match(comparator, methods)]
  ann <- ann[is.finite(x) & is.finite(xend) & label != ""]
  if (nrow(ann) == 0) return(data.table())

  ann[, comp_rank := match(comparator, setdiff(methods, "ST-COS"))]
  ann[, y := y_min + rng * (0.72 + 0.08 * (comp_rank - 1))]
  ann[, y_tip := y - rng * 0.025]
  ann[, y_label := y + rng * 0.018]
  ann[, x_mid := (x + xend) / 2]
  if (!is.null(panel_col) && !is.null(panel_value)) ann[, (panel_col) := panel_value]
  ann[]
}

add_significance_layers <- function(p, ann) {
  if (is.null(ann) || nrow(ann) == 0) return(p)
  p +
    geom_segment(data = ann, aes(x = x, xend = x, y = y_tip, yend = y),
                 inherit.aes = FALSE, linewidth = 0.28) +
    geom_segment(data = ann, aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.28) +
    geom_segment(data = ann, aes(x = xend, xend = xend, y = y_tip, yend = y),
                 inherit.aes = FALSE, linewidth = 0.28) +
    geom_text(data = ann, aes(x = x_mid, y = y_label, label = label),
              inherit.aes = FALSE, size = 3.1, fontface = "bold")
}

## =====================================================================
## DATA LOADING
## =====================================================================
message("Loading gene metrics ...")
gene_dt <- rbindlist(lapply(slides, function(sl) {
  f <- file.path(out_dir, paste0(sl, "_gene_metrics.csv"))
  dt <- fread(f); dt[, slide := sl]; dt
}), use.names = TRUE)
gene_dt[method == "STCOS", method := "ST-COS"]
gene_dt[, method := factor(method, levels = method_levels)]

message("Loading gene summaries ...")
sum_dt <- rbindlist(lapply(slides, function(sl) {
  f <- file.path(out_dir, paste0(sl, "_gene_summary.csv"))
  dt <- fread(f); dt[, slide := sl]; dt
}), use.names = TRUE)
sum_dt[method == "STCOS", method := "ST-COS"]
sum_dt[, method := factor(method, levels = method_levels)]

message("Loading correlation summaries ...")
corr_sum <- rbindlist(lapply(slides, function(sl) {
  f <- file.path(out_dir, paste0(sl, "_gene_corr_summary.csv"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f); dt[, slide := sl]; dt
}), use.names = TRUE)
corr_sum <- corr_sum[method != "STCOS_Sigma"]
corr_sum[method == "STCOS", method := "ST-COS"]
corr_sum[, method := factor(method, levels = corr_methods)]

message("Loading correlation pairs ...")
pairs_dt <- rbindlist(lapply(slides, function(sl) {
  f <- file.path(out_dir, paste0(sl, "_gene_corr_pairs.csv.gz"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f); dt[, slide := sl]; dt
}), use.names = TRUE)
# Use realized output-count correlations for ST-COS.
meas_cols <- intersect(c("STCOS","SRTsim","scDesign3","scCube"), names(pairs_dt))
pairs_long <- melt(pairs_dt, id.vars = c("slide","truth"),
                   measure.vars = meas_cols,
                   variable.name = "method", value.name = "simulated")
pairs_long[method == "STCOS", method := "ST-COS"]
pairs_long[, method := factor(method, levels = corr_methods)]

set.seed(42)
pairs_plot <- pairs_long[!is.na(simulated),
                         .SD[sample(.N, min(.N, 80000))], by = method]

fig4b_tests <- load_signif_tests(file.path(
  out_dir, "DLPFC_Figure4B_gene_metric_paired_tests_vs_STCOS.csv"))
fig4e_tests <- load_signif_tests(file.path(
  out_dir, "DLPFC_Figure4E_gene_corr_paired_tests_vs_STCOS.csv"))

## =====================================================================
## PLOT 1 — Panel A: violin+box 2×2; Panel B: per-gene RMSE 4 metrics
## Two versions: scCube normalized (pseudo-raw counts) vs unnormalized (lognorm)
## =====================================================================
message("Building Plot 1 ...")

metric_meta1 <- list(
  list(col="GeneMean",     label="Gene Mean",       log=TRUE,  ylim=c(0,   0.3)),
  list(col="GeneVar",      label="Gene Variance",   log=TRUE,  ylim=c(0,   0.5)),
  list(col="GeneCV",       label="Gene CV",         log=TRUE,  ylim=c(0,   0.9)),
  list(col="GeneZeroProp", label="Zero Proportion", log=FALSE, ylim=c(0,   0.1))
)

p1_theme <- theme_bw(base_size = 12) + theme(
  panel.background  = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.35),
  plot.background   = element_rect(fill = "white", color = NA),
  strip.background  = element_rect(fill = "#d2e0e3", color = "#6f6f6f", linewidth = 0.35),
  strip.text        = element_text(size = 11, face = "bold", color = "black"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor   = element_blank(),
  panel.grid.major.y = element_line(color = "#ebebeb", linewidth = 0.3),
  legend.position   = "none",
  axis.text.x       = element_text(angle = 40, hjust = 1, size = 10, color = "black"),
  axis.text.y       = element_text(size = 9, color = "black"),
  axis.title        = element_text(size = 11, face = "bold"),
  plot.margin       = margin(4, 6, 4, 6)
)

## ---- Build un-normalized (lognorm) scCube gene_dt --------------------
message("Loading scCube lognorm matrices for un-normalized version ...")
sccube_lognorm_dt <- rbindlist(lapply(slides, function(sl) {
  f <- file.path(out_dir, paste0(sl, "_scCube_lognorm.csv.gz"))
  if (!file.exists(f)) { warning("Missing lognorm: ", f); return(NULL) }
  dt <- fread(f)
  gene_names <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- gene_names
  data.table(
    gene         = gene_names,
    method       = "scCube",
    slide        = sl,
    GeneMean     = rowMeans(mat),
    GeneVar      = apply(mat, 1, var),
    GeneCV       = apply(mat, 1, function(x) { m <- mean(x); if (m == 0) NA_real_ else sd(x) / m }),
    GeneZeroProp = rowMeans(mat == 0)
  )
}), use.names = TRUE)
sccube_lognorm_dt[, method := factor(method, levels = method_levels)]

## gene_dt_norm  = current gene_dt (scCube = pseudo-raw back-transformed counts)
## gene_dt_unnorm = gene_dt with scCube rows replaced by lognorm-based metrics
gene_dt_norm   <- copy(gene_dt)
gene_dt_unnorm <- rbind(
  gene_dt[method != "scCube"],
  sccube_lognorm_dt[gene %in% unique(gene_dt$gene)]
)
gene_dt_unnorm[, method := factor(method, levels = method_levels)]

## ---- Helper: make_violin_panel (uses gdt argument) -------------------
make_violin_panel <- function(gdt, col, label, log_scale = FALSE, show_legend = FALSE) {
  dt <- gdt[, .(method, val = get(col))]
  dt <- dt[is.finite(val)]
  if (log_scale) dt[, val := log10(val + 1e-3)]
  dt[, panel_label := label]
  y_lab <- if (log_scale) bquote(log[10](.(label) + 0.001)) else label

  p <- ggplot(dt, aes(x = method, y = val, fill = method)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.78, linewidth = 0.2) +
    geom_boxplot(width = 0.15, linewidth = 0.35, outlier.shape = NA,
                 fill = "white", alpha = 0.75, color = "black") +
    stat_summary(fun = median, geom = "point", shape = 18,
                 size = 2, color = "black") +
    facet_wrap(~panel_label) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(x = NULL, y = y_lab, fill = "Method") +
    p1_theme
  if (show_legend) p <- p + theme(
    legend.position = "right",
    legend.key      = element_rect(fill = "white", color = NA),
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10)
  )
  p
}

## ---- Helper: save_plot1 — generate all 3 Plot 1 files for one gene_dt
raw_ylims <- list(c(0, 2), c(0, 3.5), c(0, 10), NULL)

save_plot1 <- function(gdt, suffix, sig_tests = NULL) {
  ## violin log10
  vp <- lapply(seq_along(metric_meta1), function(i) {
    m <- metric_meta1[[i]]
    make_violin_panel(gdt, m$col, m$label, m$log, show_legend = (i == length(metric_meta1)))
  })
  leg_v  <- get_legend(vp[[length(vp)]])
  grid_v <- plot_grid(plotlist = lapply(vp, function(p) p + theme(legend.position = "none")),
                      nrow = 2, ncol = 2)
  panel_v <- plot_grid(grid_v, leg_v, ncol = 2, rel_widths = c(1, 0.12))
  violin_file <- file.path(out_dir, paste0("DLPFC_summary_1_metrics_violin_", suffix, ".pdf"))
  ggsave(violin_file, panel_v, width = 14, height = 10, bg = "white")
  if (suffix == "scCube_normalized") {
    ggsave(file.path(out_dir, "DLPFC_summary_1_metrics_violin.pdf"),
           panel_v, width = 14, height = 10, bg = "white")
  }
  message(sprintf("Saved Plot 1 violin [%s]", suffix))

  ## violin raw scale
  vp_raw <- lapply(seq_along(metric_meta1), function(i) {
    m <- metric_meta1[[i]]
    p <- make_violin_panel(gdt, m$col, m$label, log_scale = FALSE,
                           show_legend = (i == length(metric_meta1)))
    if (!is.null(raw_ylims[[i]]))
      p <- p + coord_cartesian(ylim = raw_ylims[[i]], expand = FALSE)
    p
  })
  leg_vr  <- get_legend(vp_raw[[length(vp_raw)]])
  grid_vr <- plot_grid(plotlist = lapply(vp_raw, function(p) p + theme(legend.position = "none")),
                       nrow = 2, ncol = 2)
  panel_vr <- plot_grid(grid_vr, leg_vr, ncol = 2, rel_widths = c(1, 0.12))
  violin_raw_file <- file.path(out_dir, paste0("DLPFC_summary_1_metrics_violin_rawscale_", suffix, ".pdf"))
  ggsave(violin_raw_file, panel_vr, width = 14, height = 10, bg = "white")
  if (suffix == "scCube_normalized") {
    ggsave(file.path(out_dir, "DLPFC_summary_1_metrics_violin_rawscale.pdf"),
           panel_vr, width = 14, height = 10, bg = "white")
  }
  message(sprintf("Saved Plot 1 violin rawscale [%s]", suffix))

  ## per-gene RMSE
  metric_cols   <- c("GeneMean","GeneVar","GeneCV","GeneZeroProp")
  metric_labels <- c("Gene Mean","Gene Variance","Gene CV","Zero Proportion")
  rmse_list <- lapply(seq_along(metric_cols), function(k) {
    col  <- metric_cols[k]
    wide <- dcast(gdt[, .(gene, slide, method, val = get(col))],
                  gene + slide ~ method, value.var = "val")
    sc   <- intersect(sim_methods, names(wide))
    long <- melt(wide, id.vars = c("gene","slide","Truth"),
                 measure.vars = sc, variable.name = "method", value.name = "sim")
    long[, sq := (sim - Truth)^2]
    rmse <- long[!is.na(sq), .(RMSE = sqrt(mean(sq))), by = .(gene, method)]
    rmse[, metric := metric_labels[k]]
    rmse
  })
  rmse_dt <- rbindlist(rmse_list)
  rmse_dt[, method := factor(method, levels = sim_methods)]
  rmse_dt[, metric := factor(metric, levels = metric_labels)]

  make_rmse_panel <- function(met, ylim = NULL, show_legend = FALSE,
                              sig_tests = NULL) {
    dt <- rmse_dt[metric == met]
    dt[, panel_label := met]
    sig_ann <- make_stcos_annotations(
      sig_tests, metric_name = met, score_name = "RMSE", methods = sim_methods,
      y_limits = ylim, y_values = dt$RMSE,
      panel_col = "panel_label", panel_value = met
    )
    p <- ggplot(dt, aes(x = method, y = RMSE, fill = method)) +
      geom_boxplot(width = 0.62, linewidth = 0.35,
                   outlier.shape = 95, outlier.size = 1.2, na.rm = TRUE) +
      facet_wrap(~panel_label) +
      scale_fill_manual(values = method_colors[sim_methods], drop = FALSE) +
      labs(x = NULL, y = "RMSE", fill = "Method") +
      p1_theme +
      theme(strip.background = element_rect(fill = "#2c6e9e", color = NA),
            strip.text = element_text(color = "white", face = "bold", size = 11))
    p <- add_significance_layers(p, sig_ann)
    if (!is.null(ylim))
      p <- p + coord_cartesian(ylim = ylim, expand = FALSE)
    if (show_legend)
      p <- p + theme(legend.position = "right",
                     legend.key   = element_rect(fill = "white", color = NA),
                     legend.title = element_text(size = 11, face = "bold"),
                     legend.text  = element_text(size = 10))
    else
      p <- p + theme(legend.position = "none")
    p
  }

  rp <- lapply(seq_along(metric_labels), function(i)
    make_rmse_panel(metric_labels[i],
                    ylim = metric_meta1[[i]]$ylim,
                    show_legend = (i == length(metric_labels)),
                    sig_tests = sig_tests))
  leg_r  <- get_legend(rp[[length(rp)]])
  grid_r <- plot_grid(plotlist = lapply(rp, function(p) p + theme(legend.position = "none")),
                      nrow = 2, ncol = 2)
  panel_r <- plot_grid(grid_r, leg_r, ncol = 2, rel_widths = c(1, 0.12))
  rmse_file <- file.path(out_dir, paste0("DLPFC_summary_1_metrics_RMSE_", suffix, ".pdf"))
  ggsave(rmse_file, panel_r, width = 14, height = 10, bg = "white")
  if (suffix == "scCube_normalized") {
    ggsave(file.path(out_dir, "DLPFC_summary_1_metrics_RMSE.pdf"),
           panel_r, width = 14, height = 10, bg = "white")
  }
  message(sprintf("Saved Plot 1 RMSE [%s]", suffix))
}

## ---- Generate 6 files ------------------------------------------------
save_plot1(gene_dt_norm,   "scCube_normalized",   sig_tests = fig4b_tests)
save_plot1(gene_dt_unnorm, "scCube_unnormalized", sig_tests = NULL)

## =====================================================================
## PLOT 2 — % deviation from Truth per slide, boxplot across slides
## =====================================================================
message("Building Plot 2 ...")

# For each slide × method compute % diff from truth slide mean
dev_dt <- dcast(sum_dt, slide ~ method, value.var = c("Mean_GeneMean","Mean_GeneVar","Mean_GeneCV","Mean_GeneZeroProp"))
truth_cols <- grep("_Truth$", names(dev_dt), value = TRUE)

make_dev <- function(metric_col, label, is_abs = FALSE) {
  truth_col <- paste0(metric_col, "_Truth")
  sim_cols  <- paste0(metric_col, "_", sim_methods)
  sim_cols  <- intersect(sim_cols, names(dev_dt))

  long <- melt(dev_dt[, c("slide", truth_col, sim_cols), with = FALSE],
               id.vars = c("slide", truth_col),
               measure.vars = sim_cols,
               variable.name = "method_col", value.name = "sim_val")
  long[, method := gsub(paste0("^", metric_col, "_"), "", method_col)]
  long[method == "STCOS", method := "ST-COS"]
  long[, method := factor(method, levels = sim_methods)]
  long[, truth_val := get(truth_col)]

  if (is_abs) {
    long[, dev := sim_val - truth_val]
    y_lab <- paste0(label, "  (sim \u2212 truth)")
    hline <- 0
  } else {
    long[, dev := 100 * (sim_val - truth_val) / truth_val]
    y_lab <- paste0("% deviation from truth\n(", label, ")")
    hline <- 0
  }

  ggplot(long, aes(x = method, y = dev, fill = method)) +
    geom_hline(yintercept = hline, color = "black", linewidth = 0.5, linetype = "dashed") +
    geom_boxplot(width = 0.55, linewidth = 0.35, outlier.shape = NA,
                 alpha = 0.75, color = "black") +
    geom_jitter(aes(color = method), width = 0.15, size = 2.2, alpha = 0.85, shape = 16) +
    scale_fill_manual(values  = method_colors[sim_methods], drop = FALSE) +
    scale_color_manual(values = method_colors[sim_methods], drop = FALSE) +
    labs(x = NULL, y = y_lab, fill = "Method", color = "Method") +
    base_theme
}

dev_panels <- list(
  make_dev("Mean_GeneMean",     "Gene Mean",        is_abs = FALSE),
  make_dev("Mean_GeneVar",      "Gene Variance",    is_abs = FALSE),
  make_dev("Mean_GeneCV",       "Gene CV",          is_abs = FALSE),
  make_dev("Mean_GeneZeroProp", "Zero Proportion",  is_abs = TRUE)
)
dev_panels[[4]] <- dev_panels[[4]] + labs(fill="Method",color="Method") +
  theme(legend.position="right",
        legend.key   = element_rect(fill="white",color=NA),
        legend.title = element_text(size=11,face="bold"),
        legend.text  = element_text(size=10))

legend_d <- get_legend(dev_panels[[4]])
dev_no_leg <- lapply(dev_panels, function(p) p + theme(legend.position="none"))

grid_d <- plot_grid(plotlist = dev_no_leg, nrow = 2, ncol = 2,
                    labels = c("A","B","C","D"), label_size = 13)
plot2 <- plot_grid(grid_d, legend_d, ncol = 2, rel_widths = c(1, 0.14))
ggsave(file.path(out_dir, "DLPFC_summary_2_metrics_deviation.pdf"),
       plot2, width = 14, height = 10, bg = "white")
message("Saved Plot 2")

## =====================================================================
## PLOT 3 — Gene-gene correlation 2×2 density plots, zoomed
## =====================================================================
message("Building Plot 3 ...")

# per-method Pearson r annotated on each panel
r_annot <- corr_sum[, .(r_mean = round(mean(PearsonR, na.rm=TRUE), 3),
                         r_sd   = round(sd(PearsonR,   na.rm=TRUE), 3)), by = method]
r_annot[, label := sprintf("r = %.3f \u00b1 %.3f", r_mean, r_sd)]

xlim3 <- c(-0.2, 0.6)
ylim3 <- c(-0.2, 0.6)

scatter_p <- ggplot(pairs_plot, aes(x = truth, y = simulated, color = method)) +
  geom_point(alpha = 0.06, size = 0.4, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.5, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  geom_text(data = r_annot,
            aes(x = xlim3[1] + 0.02, y = xlim3[2] - 0.04, label = label),
            inherit.aes = FALSE, hjust = 0, size = 3.4,
            fontface = "italic", color = "black") +
  facet_wrap(~method, nrow = 2, ncol = 2) +
  scale_color_manual(values = method_colors, guide = "none") +
  coord_cartesian(xlim = xlim3, ylim = ylim3) +
  labs(x = "Truth gene-gene correlation",
       y = "Simulated gene-gene correlation") +
  base_theme +
  theme(strip.background = element_rect(fill = "#2c6e9e", color = NA),
        strip.text  = element_text(color = "white", face = "bold", size = 10.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
        panel.spacing = unit(0.5, "lines"))

## =====================================================================
## PLOT 4 — Gene-gene correlation RMSE per slide
## =====================================================================
message("Building Plot 4 ...")

fig4e_rmse_ann <- make_stcos_annotations(
  fig4e_tests, metric_name = "RMSE", score_name = "RMSE",
  methods = corr_methods, y_values = corr_sum$RMSE
)

p4_rmse <- ggplot(corr_sum, aes(x = method, y = RMSE, fill = method)) +
  geom_boxplot(width = 0.62, linewidth = 0.35,
               outlier.shape = NA, na.rm = TRUE) +
  geom_jitter(aes(color = method), width = 0.12, size = 2.2,
              alpha = 0.9, shape = 16, show.legend = FALSE) +
  scale_fill_manual(values = method_colors[corr_methods], drop = FALSE) +
  scale_color_manual(values = method_colors[corr_methods], drop = FALSE) +
  labs(x = NULL, y = "RMSE (gene-gene correlation)", fill = "Method") +
  base_theme +
  theme(legend.position = "right",
        legend.key   = element_rect(fill = "white", color = NA),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text  = element_text(size = 10))
p4_rmse <- add_significance_layers(p4_rmse, fig4e_rmse_ann)

legend_c <- get_legend(p4_rmse)
grid_c   <- p4_rmse + theme(legend.position = "none")
plot4_panel <- plot_grid(grid_c, legend_c, ncol = 2, rel_widths = c(1, 0.15))
ggsave(file.path(out_dir, "DLPFC_summary_4_gene_corr_rmse.pdf"),
       plot4_panel, width = 10, height = 4.2, bg = "white")
ggsave(file.path(out_dir, "DLPFC_summary_3_gene_corr_scatter.pdf"),
       scatter_p, width = 10, height = 9, bg = "white")

## =====================================================================
## COMBINED PLOT 3+4 — scatter on top, RMSE/r on bottom
## =====================================================================
plot34 <- plot_grid(
  scatter_p,
  plot4_panel,
  ncol = 1, rel_heights = c(2, 0.8),
  labels = c("A", "B"), label_size = 13
)
ggsave(file.path(out_dir, "DLPFC_summary_3_gene_corr.pdf"),
       plot34, width = 10, height = 13, bg = "white")
message("Saved Plot 3+4")

message("All summary plots complete.")
message("  Plot 1: DLPFC_summary_1_metrics_violin.pdf")
message("  Plot 2: DLPFC_summary_2_metrics_deviation.pdf")
message("  Plot 3+4: DLPFC_summary_3_gene_corr.pdf")
