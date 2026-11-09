#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

repo_dir <- "."
input_root <- file.path(repo_dir, "Reference_based", "DLPFC_parameter_sensitivity", "slides")
out_dir <- file.path(repo_dir, "Reference_based", "DLPFC_parameter_sensitivity", "multislide")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

slide_ids <- as.character(c(151507:151510, 151669:151676))
full_grid_slides <- c("151507", "151669", "151673")
donor_map <- data.table(
  slide_id = slide_ids,
  donor_group = rep(c("Donor group 1", "Donor group 2", "Donor group 3"), each = 4),
  section_range = rep(c("151507-151510", "151669-151672", "151673-151676"), each = 4),
  donor_display = rep(c("Donor 1\n(151507-151510)", "Donor 2\n(151669-151672)",
                        "Donor 3\n(151673-151676)"), each = 4)
)
section_colors <- c("Section 151507" = "#0072B2", "Section 151669" = "#D55E00",
                    "Section 151673" = "#009E73")
donor_colors <- c("Donor 1\n(151507-151510)" = "#0072B2",
                  "Donor 2\n(151669-151672)" = "#D55E00",
                  "Donor 3\n(151673-151676)" = "#009E73")

read_slide <- function(slide_id, filename) {
  sid <- as.character(slide_id)
  path <- file.path(input_root, slide_id, filename)
  if (!file.exists(path)) stop("Missing sensitivity output: ", path)
  dt <- fread(path)
  if ("slide_id" %in% names(dt)) dt[, slide_id := NULL]
  dt[, slide_id := sid]
  merge(dt, donor_map, by = "slide_id", all.x = TRUE, sort = FALSE)
}

defaults <- rbindlist(lapply(slide_ids, read_slide,
                             filename = "sensitivity_default_validation.csv"), fill = TRUE)
matrix_raw <- rbindlist(lapply(full_grid_slides, read_slide,
                               filename = "sensitivity_matrix_repair.csv"), fill = TRUE)
matrix_summary <- matrix_raw[, lapply(.SD, mean, na.rm = TRUE),
  by = .(slide_id, donor_group, tau, epsilon),
  .SDcols = c("correlation_rmse", "mean_log_rmse", "variance_log_rmse",
              "moran_rmse", "median_spatial_spearman", "min_eigen_before",
              "min_eigen_after", "uniform_eigenvalue_shift",
              "n_eigenvalues_below_floor", "n_eigenvalues_modified",
              "relative_frobenius_change", "maximum_entry_change",
              "total_seconds")]
dispersion_summary <- rbindlist(lapply(full_grid_slides, read_slide,
                                       filename = "sensitivity_dispersion_cap_summary.csv"), fill = TRUE)
spatial_summary <- rbindlist(lapply(full_grid_slides, read_slide,
                                    filename = "sensitivity_spatial_summary.csv"), fill = TRUE)
matrix_summary[, section_label := paste("Section", slide_id)]
dispersion_summary[, section_label := paste("Section", slide_id)]
spatial_summary[, section_label := paste("Section", slide_id)]

default_by_slide <- defaults[, .(
  correlation_rmse = mean(correlation_rmse),
  mean_log_rmse = mean(mean_log_rmse),
  variance_log_rmse = mean(variance_log_rmse),
  moran_rmse = mean(moran_rmse),
  median_spatial_spearman = mean(median_spatial_spearman),
  relative_frobenius_change = mean(relative_frobenius_change),
  total_seconds = mean(total_seconds),
  n_dispersion_fallback = unique(n_dispersion_fallback)
), by = .(slide_id, donor_group, section_range, donor_display)]

default_across_slides <- melt(
  default_by_slide,
  id.vars = c("slide_id", "donor_group", "section_range", "donor_display",
              "n_dispersion_fallback"),
  measure.vars = c("correlation_rmse", "mean_log_rmse", "variance_log_rmse",
                   "moran_rmse", "median_spatial_spearman",
                   "relative_frobenius_change", "total_seconds"),
  variable.name = "metric", value.name = "value"
)[, .(
  mean = mean(value),
  ci_low = quantile(value, 0.025),
  ci_high = quantile(value, 0.975),
  n_slides = uniqueN(slide_id),
  n_donor_groups = uniqueN(donor_group)
), by = metric]

epsilon_robustness <- matrix_summary[tau == 0.05, .(
  minimum_correlation_rmse = min(correlation_rmse),
  maximum_correlation_rmse = max(correlation_rmse),
  correlation_rmse_range = diff(range(correlation_rmse)),
  minimum_matrix_perturbation = min(relative_frobenius_change),
  maximum_matrix_perturbation = max(relative_frobenius_change)
), by = .(slide_id, donor_group)]

fallback_summary <- default_by_slide[, .(
  n_genes_using_fallback = n_dispersion_fallback
), by = .(slide_id, donor_group)]

fwrite(defaults, file.path(out_dir, "default_validation_by_replicate.csv"))
fwrite(default_by_slide, file.path(out_dir, "default_validation_by_slide.csv"))
fwrite(default_across_slides, file.path(out_dir, "default_validation_across_slides.csv"))
fwrite(matrix_summary, file.path(out_dir, "full_grid_matrix_three_slides.csv"))
fwrite(dispersion_summary, file.path(out_dir, "full_grid_dispersion_three_slides.csv"))
fwrite(spatial_summary, file.path(out_dir, "full_grid_spatial_three_slides.csv"))
fwrite(epsilon_robustness, file.path(out_dir, "epsilon_robustness_by_slide.csv"))
fwrite(fallback_summary, file.path(out_dir, "dispersion_fallback_by_slide.csv"))
matrix_diagnostics <- copy(matrix_summary)
matrix_diagnostics[, required_minimum_eigenvalue_lift :=
  uniform_eigenvalue_shift]
fwrite(
  matrix_diagnostics[, .(
    slide_id, donor_group, tau, epsilon, correlation_rmse,
    min_eigen_before, required_minimum_eigenvalue_lift, min_eigen_after,
    mean_n_eigenvalues_below_floor = n_eigenvalues_below_floor,
    mean_n_eigenvalues_modified = n_eigenvalues_modified,
    relative_frobenius_change, maximum_entry_change
  )],
  file.path(out_dir, "matrix_repair_diagnostics_three_sections.csv")
)
supplementary_diagnostics <- matrix_diagnostics[tau == 0.05, .(
  section = slide_id,
  epsilon_floor = epsilon,
  count_correlation_rmse = correlation_rmse,
  minimum_eigenvalue_before = min_eigen_before,
  required_minimum_eigenvalue_lift,
  minimum_eigenvalue_after = min_eigen_after,
  mean_number_eigenvalues_below_floor = n_eigenvalues_below_floor,
  mean_number_eigenvalues_modified = n_eigenvalues_modified,
  relative_frobenius_change,
  maximum_entrywise_change = maximum_entry_change
)]
fwrite(
  supplementary_diagnostics,
  file.path(out_dir, "Supplementary_Table_S3_matrix_repair_diagnostics.csv")
)
writeLines(c(
  paste0("Supplementary Table S3. Matrix-repair diagnostics at tau=0.05 in DLPFC ",
         "sections 151507, 151669, and 151673."),
  paste0("Values are means over three Monte Carlo replicates. The required minimum ",
         "eigenvalue lift is max(0, epsilon0 - lambda_min) before unit-diagonal ",
         "renormalization; it is a diagnostic, not a uniform diagonal shift. The ",
         "minimum eigenvalue after repair is reported after renormalization."),
  paste0("Count-correlation RMSE remains essentially unchanged across epsilon0, ",
         "supporting its interpretation as a numerical eigenvalue floor rather than ",
         "a performance-tuning parameter.")
), file.path(out_dir, "Supplementary_Table_S3_legend.txt"))

p_epsilon <- ggplot(
  matrix_summary[tau == 0.05],
  aes(epsilon, correlation_rmse, color = section_label, group = slide_id)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.7) +
  scale_x_log10() +
  scale_color_manual(values = section_colors) +
  labs(title = "A  Epsilon robustness at tau = 0.05",
       subtitle = "Flat responses identify epsilon as a numerical floor",
       x = "Eigenvalue floor",
       y = "Count-correlation RMSE", color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

p_tau <- ggplot(
  matrix_summary[epsilon == 1e-6],
  aes(tau, relative_frobenius_change, color = section_label, group = slide_id)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.7) +
  scale_color_manual(values = section_colors) +
  labs(title = "B  Sparsification perturbation", x = "Correlation threshold (tau)",
       y = "Relative Frobenius change", color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

p_spatial <- ggplot(
  spatial_summary,
  aes(factor(neighbors), moran_rmse, color = factor(spline_basis_dimension),
      group = interaction(slide_id, spline_basis_dimension))
) +
  geom_line(alpha = 0.65) +
  geom_point(size = 1.6) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  facet_wrap(~section_label) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "C  Graph-specific Moran's I recovery",
       subtitle = "Broader graphs lower absolute RMSE by smoothing; dashed 6 = local first-ring default",
       x = "Neighbors",
       y = "Moran's I RMSE", color = "Spline basis") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

stress <- dispersion_summary[evaluation_regime == "Poisson-boundary stress test"]
stress[, cap_order := factor(dispersion_setting,
  levels = c("1e+04", "1e+06", "1e+08", "Poisson limit"))]
p_disp <- ggplot(stress, aes(cap_order, correlation_rmse,
                             color = section_label, group = slide_id)) +
  geom_line() +
  geom_point(size = 1.8) +
  scale_y_continuous(limits = c(0, 1e-4), breaks = c(0, 5e-5, 1e-4),
                     labels = scales::label_scientific()) +
  scale_color_manual(values = section_colors) +
  labs(title = "D  Dispersion boundary: count-correlation error",
       subtitle = "All RMSE < 1e-4; theta = 1e6 is numerically equivalent to Poisson",
       x = "Fallback setting", y = "Count-correlation RMSE", color = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none", plot.margin = margin(5.5, 5.5, 5.5, 14))

default_long <- melt(
  default_by_slide,
  id.vars = c("slide_id", "donor_group", "section_range", "donor_display",
              "n_dispersion_fallback"),
  measure.vars = c("correlation_rmse", "mean_log_rmse", "variance_log_rmse",
                   "moran_rmse", "median_spatial_spearman"),
  variable.name = "metric", value.name = "value"
)
default_long[, metric := factor(metric,
  levels = c("correlation_rmse", "mean_log_rmse", "variance_log_rmse",
             "moran_rmse", "median_spatial_spearman"),
  labels = c("Correlation RMSE", "Mean RMSE", "Variance RMSE",
             "Moran RMSE", "Spatial Spearman"))]
p_default <- ggplot(default_long, aes(donor_display, value, color = donor_display)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.08, size = 1.4) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = donor_colors) +
  labs(title = "E  Default validation across all 12 sections", x = NULL,
       y = "Section-level value", color = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "none")

top <- plot_grid(p_epsilon, p_tau, p_spatial, nrow = 1,
                 rel_widths = c(0.9, 0.9, 1.25), align = "hv", axis = "tblr")
bottom <- plot_grid(p_disp, p_default, nrow = 1, rel_widths = c(0.9, 1.8),
                    align = "hv", axis = "tblr")
figure <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 0.9))
ggsave(file.path(out_dir, "DLPFC_multislide_parameter_sensitivity.pdf"),
       figure, width = 18, height = 9.5)
ggsave(file.path(out_dir, "DLPFC_multislide_parameter_sensitivity.png"),
       figure, width = 18, height = 9.5, dpi = 180)

# Focused multipanel figure for Reviewer 2, Minor Comment 2.1.
m21_a <- p_epsilon +
  geom_vline(xintercept = 1e-6, linetype = "dashed", color = "grey35",
             linewidth = 0.45) +
  labs(
    title = "A  Eigenvalue-floor robustness",
    subtitle = "At tau = 0.05; dashed line marks the default epsilon0 = 1e-6"
  )

m21_b <- p_tau +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey35",
             linewidth = 0.45) +
  labs(
    title = "B  Sparsification changes the requested matrix",
    subtitle = "At epsilon0 = 1e-6; dashed line marks the default tau = 0.05"
  )

m21_c <- p_spatial +
  labs(
    title = "C  Neighborhood and spline sensitivity",
    subtitle = "Graph-specific gene sets; dashed line marks the local k = 6 default"
  )

m21_d <- p_disp +
  geom_vline(xintercept = 2, linetype = "dashed", color = "grey35",
             linewidth = 0.45) +
  labs(
    title = "D  Dispersion-boundary sensitivity",
    subtitle = "All errors < 1e-4; dashed line marks theta = 1e6"
  )

runtime_dt <- rbindlist(list(
  matrix_summary[, .(branch = "tau / epsilon", seconds = total_seconds)],
  spatial_summary[, .(branch = "Moran / spline", seconds = total_seconds)],
  dispersion_summary[, .(branch = "dispersion cap", seconds = total_seconds)]
))
runtime_dt[, branch := factor(
  branch,
  levels = c("dispersion cap", "Moran / spline", "tau / epsilon")
)]
runtime_colors <- c(
  "dispersion cap" = "#E76F51",
  "Moran / spline" = "#2A9D8F",
  "tau / epsilon" = "#457B9D"
)
m21_e <- ggplot(runtime_dt, aes(branch, seconds, fill = branch)) +
  geom_boxplot(width = 0.58, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.45, size = 1) +
  scale_y_log10() +
  scale_fill_manual(values = runtime_colors) +
  labs(
    title = "E  Runtime per parameter-setting replicate",
    subtitle = "Wall-clock time; log scale",
    x = NULL,
    y = "Seconds"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

m21_default <- default_long[
  metric %in% c("Mean RMSE", "Variance RMSE", "Spatial Spearman")
]
m21_default[, metric := droplevels(metric)]
m21_f <- ggplot(
  m21_default,
  aes(donor_display, value, color = donor_display)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.08, size = 1.35) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = donor_colors) +
  labs(
    title = "F  Default-setting validation across all 12 sections",
    subtitle = "Four DLPFC sections per donor group",
    x = NULL,
    y = "Section-level value"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 18, hjust = 1, size = 7.5),
    legend.position = "none"
  )

m21_top <- plot_grid(
  m21_a, m21_b, m21_c,
  nrow = 1, rel_widths = c(0.9, 0.9, 1.25), align = "hv", axis = "tblr"
)
m21_bottom <- plot_grid(
  m21_d, m21_e, m21_f,
  nrow = 1, rel_widths = c(0.9, 0.8, 1.45), align = "hv", axis = "tblr"
)
m21_figure <- plot_grid(
  m21_top, m21_bottom, ncol = 1, rel_heights = c(1, 0.92)
)
ggsave(
  file.path(out_dir, "Reviewer2_Minor2_1_parameter_sensitivity.pdf"),
  m21_figure, width = 18, height = 9.8
)
ggsave(
  file.path(out_dir, "Reviewer2_Minor2_1_parameter_sensitivity.png"),
  m21_figure, width = 18, height = 9.8, dpi = 220
)
writeLines(c(
  paste0("Response Figure R2.M1. Sensitivity of reference-based estimation and ",
         "simulation settings."),
  paste0("(A) Count-correlation RMSE across eigenvalue floors at tau=0.05 in ",
         "DLPFC sections 151507, 151669, and 151673. The flat response identifies ",
         "epsilon0 as a numerical floor; the dashed line marks epsilon0=10^-6."),
  paste0("(B) Relative Frobenius perturbation across sparsification thresholds at ",
         "epsilon0=10^-6. The dashed line marks tau=0.05. Larger tau values yield ",
         "progressively larger departures from the requested matrix."),
  paste0("(C) Graph-specific Moran's I RMSE across 4, 6, 8, and 12 neighbors and ",
         "spline bases 20, 50, and 100. Each graph uses its own selected-gene set, ",
         "so the panel is diagnostic and does not identify a universal optimum. ",
         "The dashed line marks the six-neighbor first-ring default."),
  paste0("(D) Count-correlation RMSE in a 40-gene Poisson-boundary stress test. ",
         "All values are below 10^-4, and theta=10^6 is numerically close to the ",
         "explicit Poisson limit. No gene in the observed DLPFC panels required ",
         "the fallback."),
  paste0("(E) Wall-clock runtime per complete parameter-setting replicate; the ",
         "spatial branch includes spline fitting, simulation, and evaluation."),
  paste0("(F) Marginal and spatial diagnostics at the default settings across all ",
         "12 DLPFC sections, grouped into three donor groups. Points are sections; ",
         "spots within a section are not treated as independent replicates."),
  paste0("Full parameter grids were evaluated in sections 151507, 151669, and ",
         "151673; default settings were validated in all 12 sections. Matrix-repair ",
         "diagnostics are reported separately in Supplementary Table S3.")
), file.path(out_dir, "Reviewer2_Minor2_1_parameter_sensitivity_legend.txt"))

# Focused figure embedded in the response to Reviewer 1, Comment 1.2.
r12_detail <- copy(matrix_summary[slide_id == "151673"])
r12_detail[, epsilon_label := factor(
  format(epsilon, scientific = TRUE),
  levels = format(sort(unique(epsilon)), scientific = TRUE)
)]
r12_a <- ggplot(r12_detail,
                aes(factor(tau), epsilon_label, fill = correlation_rmse)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", correlation_rmse)), size = 3.2) +
  scale_fill_gradient(low = "#f7fbff", high = "#b2182b", limits = c(0, 0.15),
                      oob = scales::squish) +
  labs(title = "A  Realized count-correlation error",
       subtitle = "Section 151673; fixed 0-0.15 scale",
       x = "Sparsification threshold (tau)",
       y = "Eigenvalue floor (epsilon0)", fill = "RMSE") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 8, 8, 14))

r12_b <- ggplot(r12_detail,
                aes(factor(tau), epsilon_label, fill = relative_frobenius_change)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", relative_frobenius_change),
                color = relative_frobenius_change >= 0.06), size = 3.2) +
  scale_fill_viridis_c(option = "B", direction = -1) +
  scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
  labs(title = "B  Perturbation of the requested matrix",
       subtitle = "Section 151673",
       x = "Sparsification threshold (tau)",
       y = "Eigenvalue floor (epsilon0)", fill = "Relative\nFrobenius") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(10, 8, 8, 14))

r12_c <- ggplot(
  matrix_summary[tau == 0.05],
  aes(epsilon, correlation_rmse, color = section_label, group = section_label)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.9) +
  scale_x_log10() +
  scale_color_manual(values = section_colors) +
  labs(title = "C  Epsilon robustness across sections",
       subtitle = "At tau = 0.05; largest within-section RMSE range: 0.000218",
       x = "Eigenvalue floor (epsilon0)", y = "Count-correlation RMSE",
       color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.margin = margin(10, 8, 8, 14))

r12_d <- ggplot(
  matrix_summary[epsilon == 1e-6],
  aes(tau, relative_frobenius_change, color = section_label, group = section_label)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.9) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey35") +
  scale_color_manual(values = section_colors) +
  labs(title = "D  Tau controls matrix perturbation",
       subtitle = "Dashed line: default tau = 0.05",
       x = "Sparsification threshold (tau)", y = "Relative Frobenius change",
       color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.margin = margin(10, 8, 8, 14))

r12_figure <- plot_grid(r12_a, r12_b, r12_c, r12_d, ncol = 2,
                        align = "hv", axis = "tblr")
ggsave(file.path(out_dir, "Reviewer1_Comment1_2_epsilon_sensitivity.pdf"),
       r12_figure, width = 14, height = 9)
ggsave(file.path(out_dir, "Reviewer1_Comment1_2_epsilon_sensitivity.png"),
       r12_figure, width = 14, height = 9, dpi = 220)
writeLines(c(
  paste0("Response Figure R1.2. Sensitivity to the eigenvalue floor epsilon0 and ",
         "sparsification threshold tau."),
  paste0("(A-B) Joint-grid results for DLPFC section 151673. Count-correlation RMSE ",
         "changes negligibly across epsilon0, whereas increasing tau causes progressively ",
         "larger perturbation of the requested latent matrix."),
  paste0("(C-D) Full-grid confirmation in sections 151507, 151669, and 151673. At ",
         "tau=0.05, the largest within-section RMSE range across eleven orders of ",
         "magnitude of epsilon0 is 0.000218. The dashed line in D marks tau=0.05. ",
         "Although tau=0.20 yields slightly lower realized RMSE in section 151673, ",
         "tau=0.05 is retained to limit perturbation of the requested matrix."),
  paste0("Points and heat-map cells summarize three independent copula draws per setting; ",
         "the randomized-PIT estimate is held fixed within each section.")
), file.path(out_dir, "Reviewer1_Comment1_2_epsilon_sensitivity_legend.txt"))

writeLines(c(
  paste0("Supplementary Figure. DLPFC parameter sensitivity. Full grids use sections ",
         paste(full_grid_slides, collapse = ", "),
         "; panel E summarizes all 12 sections, grouped by donor with the four section ",
         "IDs shown beneath each donor label."),
  paste0("(A) Count-correlation RMSE over epsilon0 at tau=0.05. The essentially flat ",
         "responses across eleven orders of magnitude show that epsilon0 is a numerical ",
         "eigenvalue floor rather than a performance-tuning parameter."),
  "(B) Relative Frobenius perturbation over tau at epsilon0=10^-6.",
  paste0("(C) Graph-specific Moran's I RMSE over neighborhood size and spline basis. ",
         "The displayed values are section-specific, not pooled. Twelve neighbors gave ",
         "the lowest absolute RMSE in this grid because broader graphs average over more ",
         "spots, shrink Moran's I magnitudes, and can reward smoothing while obscuring ",
         "local structure. Six remains the default because it represents the immediate ",
         "first-ring neighborhood of a hexagonal Visium array; the dashed line marks ",
         "this default. Users targeting broader spatial scales may select a larger graph."),
  paste0("(D) Count-correlation RMSE in the Poisson-boundary dispersion stress test; ",
         "all values were below 10^-4, and theta=10^6 was numerically indistinguishable ",
         "from the explicit Poisson limit."),
  paste0("(E) Default-setting validation across all 12 sections. Donor labels include ",
         "the four section IDs used for each donor-level group.")
), file.path(out_dir, "DLPFC_multislide_parameter_sensitivity_legend.txt"))

cat(
  "Multi-slide DLPFC parameter sensitivity\n",
  sprintf("Full grids: %s. Default validation: all %d sections.\n",
          paste(full_grid_slides, collapse = ", "), length(slide_ids)),
  sprintf("At tau=0.05, the within-slide epsilon RMSE range was %.3g--%.3g.\n",
          min(epsilon_robustness$correlation_rmse_range),
          max(epsilon_robustness$correlation_rmse_range)),
  sprintf("Dispersion fallback usage across default validations: %d--%d genes per section.\n",
          min(fallback_summary$n_genes_using_fallback),
          max(fallback_summary$n_genes_using_fallback)),
  "Monte Carlo replicates are summarized within section before cross-section summaries.\n",
  file = file.path(out_dir, "README_results.txt")
)

message("Combined multi-slide sensitivity outputs written to: ", out_dir)
