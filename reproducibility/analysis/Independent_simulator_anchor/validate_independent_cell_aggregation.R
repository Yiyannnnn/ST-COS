#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

base_dir <- "."
out_dir <- file.path(
  base_dir, "Independent_simulator_anchor", "spot_aggregation_validation"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(base_dir, "Data", "simulator_YZ.R"))

set.seed(20260905)

cell_numbers <- c(1L, 2L, 5L, 10L, 20L)
latent_rhos <- c(0, 0.2, 0.4, 0.6, 0.8)
n_spots <- 3000L
n_reps <- 30L
mu <- c(5, 5)
theta <- c(5, 5)

sample_nb_pair <- function(n, rho, mu_pair, theta_pair) {
  if (n == 0L) return(matrix(integer(), nrow = 0L, ncol = 2L))
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2L)
  lower_chol <- t(chol(sigma))
  latent <- lower_chol %*% matrix(rnorm(2L * n), nrow = 2L)
  uniforms <- pmin(pmax(pnorm(latent), 1e-12), 1 - 1e-12)
  cbind(
    qnbinom(uniforms[1L, ], size = theta_pair[1L], mu = mu_pair[1L]),
    qnbinom(uniforms[2L, ], size = theta_pair[2L], mu = mu_pair[2L])
  )
}

aggregate_fixed_n <- function(n_spots, n_cells, rho, mu_pair, theta_pair,
                              method = c("Independent cells", "Scaled one draw")) {
  method <- match.arg(method)
  if (method == "Independent cells") {
    cells <- sample_nb_pair(n_spots * n_cells, rho, mu_pair, theta_pair)
    y1 <- colSums(matrix(cells[, 1L], nrow = n_cells, ncol = n_spots))
    y2 <- colSums(matrix(cells[, 2L], nrow = n_cells, ncol = n_spots))
  } else {
    one_draw <- sample_nb_pair(n_spots, rho, mu_pair, theta_pair)
    y1 <- n_cells * one_draw[, 1L]
    y2 <- n_cells * one_draw[, 2L]
  }
  cbind(y1, y2)
}

aggregate_variable_n <- function(cell_counts, rho, mu_pair, theta_pair) {
  n_out <- length(cell_counts)
  ans <- matrix(0, nrow = n_out, ncol = 2L)
  for (n_cells in sort(unique(cell_counts[cell_counts > 0L]))) {
    idx <- which(cell_counts == n_cells)
    cells <- sample_nb_pair(length(idx) * n_cells, rho, mu_pair, theta_pair)
    ans[idx, 1L] <- colSums(matrix(cells[, 1L], nrow = n_cells, ncol = length(idx)))
    ans[idx, 2L] <- colSums(matrix(cells[, 2L], nrow = n_cells, ncol = length(idx)))
  }
  ans
}

message("Simulating homogeneous cell-type aggregation...")
raw_rows <- vector(
  "list", length(cell_numbers) * length(latent_rhos) * n_reps * 2L
)
row_index <- 1L
for (rep_id in seq_len(n_reps)) {
  for (rho in latent_rhos) {
    for (n_cells in cell_numbers) {
      for (method in c("Independent cells", "Scaled one draw")) {
        spot_counts <- aggregate_fixed_n(
          n_spots, n_cells, rho, mu, theta, method
        )
        raw_rows[[row_index]] <- data.table(
          replicate = rep_id,
          latent_rho = rho,
          cell_number = n_cells,
          method = method,
          mean_gene1 = mean(spot_counts[, 1L]),
          variance_gene1 = var(spot_counts[, 1L]),
          pearson = cor(spot_counts[, 1L], spot_counts[, 2L], method = "pearson"),
          spearman = cor(spot_counts[, 1L], spot_counts[, 2L], method = "spearman")
        )
        row_index <- row_index + 1L
      }
    }
  }
}
raw_dt <- rbindlist(raw_rows)

cor_long <- melt(
  raw_dt,
  id.vars = c("replicate", "latent_rho", "cell_number", "method"),
  measure.vars = c("pearson", "spearman"),
  variable.name = "correlation_type",
  value.name = "realized_correlation"
)
cor_long[, correlation_type := fifelse(
  correlation_type == "pearson", "Pearson", "Spearman"
)]
cor_summary <- cor_long[, .(
  mean_correlation = mean(realized_correlation),
  sd_correlation = sd(realized_correlation),
  ci_low = quantile(realized_correlation, 0.025),
  ci_high = quantile(realized_correlation, 0.975),
  n_reps = .N
), by = .(latent_rho, cell_number, method, correlation_type)]

single_cell_variance <- mu[1L] + mu[1L]^2 / theta[1L]
variance_summary <- raw_dt[latent_rho == 0.6, .(
  mean_variance = mean(variance_gene1),
  sd_variance = sd(variance_gene1),
  ci_low = quantile(variance_gene1, 0.025),
  ci_high = quantile(variance_gene1, 0.975),
  theoretical_variance = if (method[1L] == "Independent cells") {
    cell_number[1L] * single_cell_variance
  } else {
    cell_number[1L]^2 * single_cell_variance
  },
  n_reps = .N
), by = .(cell_number, method)]

variance_slopes <- raw_dt[latent_rho == 0.6, {
  fit <- lm(log(variance_gene1) ~ log(cell_number))
  .(log_log_slope = unname(coef(fit)[2L]))
}, by = .(replicate, method)]
variance_slope_summary <- variance_slopes[, .(
  mean_slope = mean(log_log_slope),
  sd_slope = sd(log_log_slope),
  ci_low = quantile(log_log_slope, 0.025),
  ci_high = quantile(log_log_slope, 0.975)
), by = method]

message("Simulating heterogeneous two-cell-type mixtures...")
mixed_rows <- vector("list", length(cell_numbers) * n_reps)
row_index <- 1L
for (rep_id in seq_len(n_reps)) {
  for (n_cells in cell_numbers) {
    # The continuously varying proportions mimic a composition gradient.
    prop_a <- rbeta(n_spots, shape1 = 2, shape2 = 2)
    n_a <- rbinom(n_spots, size = n_cells, prob = prop_a)
    n_b <- n_cells - n_a
    counts_a <- aggregate_variable_n(n_a, 0.6, c(2, 8), theta)
    counts_b <- aggregate_variable_n(n_b, 0.6, c(8, 2), theta)
    mixed_counts <- counts_a + counts_b
    mixed_rows[[row_index]] <- data.table(
      replicate = rep_id,
      latent_rho_within_type = 0.6,
      cell_number = n_cells,
      pearson = cor(mixed_counts[, 1L], mixed_counts[, 2L], method = "pearson"),
      spearman = cor(mixed_counts[, 1L], mixed_counts[, 2L], method = "spearman"),
      realized_mean_prop_a = mean(n_a / n_cells)
    )
    row_index <- row_index + 1L
  }
}
mixed_raw <- rbindlist(mixed_rows)
mixed_long <- melt(
  mixed_raw,
  id.vars = c("replicate", "latent_rho_within_type", "cell_number",
              "realized_mean_prop_a"),
  measure.vars = c("pearson", "spearman"),
  variable.name = "correlation_type",
  value.name = "realized_correlation"
)
mixed_long[, correlation_type := fifelse(
  correlation_type == "pearson", "Pearson", "Spearman"
)]
mixed_summary <- mixed_long[, .(
  mean_correlation = mean(realized_correlation),
  sd_correlation = sd(realized_correlation),
  ci_low = quantile(realized_correlation, 0.025),
  ci_high = quantile(realized_correlation, 0.975),
  n_reps = .N
), by = .(latent_rho_within_type, cell_number, correlation_type)]

# Exercise the production implementation and enforce its key invariants.
smoke_spots <- 400L
smoke_n <- rep(10L, smoke_spots)
smoke_prop <- cbind(
  TypeA = seq(0.1, 0.9, length.out = smoke_spots),
  TypeB = seq(0.9, 0.1, length.out = smoke_spots)
)
smoke_mean <- list(
  TypeA = matrix(rep(c(2, 8), each = smoke_spots), nrow = smoke_spots),
  TypeB = matrix(rep(c(8, 2), each = smoke_spots), nrow = smoke_spots)
)
colnames(smoke_mean$TypeA) <- colnames(smoke_mean$TypeB) <- c("gene_1", "gene_2")
smoke_rel <- data.frame(TypeA = c(1, 1), TypeB = c(1, 1))
smoke_cov <- array(0, dim = c(2, 2, 2),
                   dimnames = list(NULL, NULL, c("TypeA", "TypeB")))
smoke_cov[, , "TypeA"] <- smoke_cov[, , "TypeB"] <- matrix(
  c(1, 0.6, 0.6, 1), nrow = 2
)
smoke <- generate_spot_count_independent(
  mean_param = smoke_mean,
  cell_number = smoke_n,
  celltype_proportion = smoke_prop,
  gene_relationship = smoke_rel,
  cov_gene = smoke_cov,
  theta_g = list(TypeA = theta, TypeB = theta),
  return_details = TRUE,
  seed = 20260905
)
stopifnot(
  all(rowSums(smoke$cell_allocations) == smoke_n),
  all(smoke$cell_allocations >= 0L),
  all(smoke$cell_allocations == round(smoke$cell_allocations)),
  all(smoke$spot_count >= 0L),
  all(smoke$spot_count == round(smoke$spot_count))
)

fwrite(raw_dt, file.path(out_dir, "homogeneous_aggregation_raw.csv"))
fwrite(cor_summary, file.path(out_dir, "homogeneous_correlation_summary.csv"))
fwrite(variance_summary, file.path(out_dir, "variance_scaling_summary.csv"))
fwrite(variance_slopes, file.path(out_dir, "variance_scaling_slopes_by_rep.csv"))
fwrite(variance_slope_summary, file.path(out_dir, "variance_scaling_slope_summary.csv"))
fwrite(mixed_raw, file.path(out_dir, "mixed_composition_aggregation_raw.csv"))
fwrite(mixed_summary, file.path(out_dir, "mixed_composition_correlation_summary.csv"))
fwrite(as.data.table(smoke$cell_allocations),
       file.path(out_dir, "production_smoke_cell_allocations.csv"))

method_colors <- c(
  "Independent cells" = "#0072B2",
  "Scaled one draw" = "#D55E00"
)
correlation_colors <- c("Pearson" = "#0072B2", "Spearman" = "#009E73")
base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e5e5e5", linewidth = 0.3),
    strip.background = element_rect(fill = "#e3eaec", color = "#777777"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "#444444", size = 9.5)
  )

p_variance <- ggplot(
  variance_summary,
  aes(x = cell_number, y = mean_variance, color = method, group = method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.45, linewidth = 0.4) +
  geom_line(aes(y = theoretical_variance), linetype = "dashed", linewidth = 0.65) +
  scale_x_log10(breaks = cell_numbers) +
  scale_y_log10() +
  scale_color_manual(values = method_colors) +
  labs(
    title = "A  Variance scaling under spot aggregation",
    subtitle = "Solid: empirical; dashed: NB expectation",
    x = "Cells per spot (N)", y = "Variance of aggregated gene count",
    color = NULL
  ) +
  base_theme

make_correlation_panel <- function(correlation_name, panel_label) {
  plot_dt <- cor_summary[correlation_type == correlation_name]
  ggplot(
    plot_dt,
    aes(x = latent_rho, y = mean_correlation, color = method, group = method)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#666666") +
    geom_ribbon(
      aes(ymin = ci_low, ymax = ci_high, fill = method),
      alpha = 0.12, color = NA
    ) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.7) +
    facet_wrap(~cell_number, nrow = 1L, labeller = labeller(
      cell_number = function(x) paste0("N = ", x)
    )) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    scale_x_continuous(breaks = latent_rhos) +
    coord_cartesian(ylim = c(-0.05, 0.85)) +
    labs(
      title = paste0(panel_label, "  ", correlation_name, " correlation after aggregation"),
      x = "Specified latent correlation", y = "Realized count correlation",
      color = NULL, fill = NULL
    ) +
    base_theme +
    theme(legend.position = "none", axis.text.x = element_text(size = 7.5))
}

p_pearson <- make_correlation_panel("Pearson", "B")
p_spearman <- make_correlation_panel("Spearman", "C")

p_mixed <- ggplot(
  mixed_summary,
  aes(x = cell_number, y = mean_correlation, color = correlation_type,
      group = correlation_type)
) +
  geom_hline(yintercept = 0, color = "#777777", linewidth = 0.35) +
  geom_ribbon(
    aes(ymin = ci_low, ymax = ci_high, fill = correlation_type),
    alpha = 0.13, color = NA
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.1) +
  scale_x_continuous(breaks = cell_numbers) +
  scale_color_manual(values = correlation_colors) +
  scale_fill_manual(values = correlation_colors) +
  labs(
    title = "D  Correlation in heterogeneous cell-type mixtures",
    subtitle = "Within-type latent correlation is 0.6; cell types have opposing mean profiles",
    x = "Cells per spot (N)", y = "Realized spot-level correlation",
    color = NULL, fill = NULL
  ) +
  base_theme

top_row <- plot_grid(p_variance, p_mixed, nrow = 1L, rel_widths = c(1, 1))
combined <- plot_grid(top_row, p_pearson, p_spearman, ncol = 1L,
                      rel_heights = c(1.05, 1, 1))

ggsave(
  file.path(out_dir, "independent_cell_aggregation_validation.pdf"),
  combined, width = 14, height = 11, device = cairo_pdf
)
ggsave(
  file.path(out_dir, "independent_cell_aggregation_validation.png"),
  combined, width = 14, height = 11, dpi = 300
)

ind_slope <- variance_slope_summary[method == "Independent cells", mean_slope]
legacy_slope <- variance_slope_summary[method == "Scaled one draw", mean_slope]
mixed_p20 <- mixed_summary[
  cell_number == 20L & correlation_type == "Pearson", mean_correlation
]
report_lines <- c(
  "Independent-cell spot aggregation validation",
  sprintf("Simulation: %d spots x %d replicates per condition.", n_spots, n_reps),
  sprintf("Independent aggregation log-log variance slope: %.3f (expected 1).", ind_slope),
  sprintf("Scaled-one-draw log-log variance slope: %.3f (expected 2).", legacy_slope),
  sprintf("Mixed-composition Pearson correlation at N=20: %.3f despite within-type latent rho=0.6.", mixed_p20),
  "The production smoke test confirmed non-negative integer counts and exact per-spot cell totals."
)
writeLines(report_lines, file.path(out_dir, "validation_summary.txt"))

message("Outputs written to: ", out_dir)
