#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

base_dir <- "."
out_dir <- file.path(
  base_dir, "Independent_simulator_anchor", "spot_aggregation_validation"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260920)

n_spots <- 3000L
n_reps <- 30L
cell_numbers <- c(5L, 10L, 20L)
latent_rho <- 0.6
theta <- c(5, 5)

sample_nb_pair <- function(n, rho, mu_pair, theta_pair) {
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2L)
  lower_chol <- t(chol(sigma))
  latent <- lower_chol %*% matrix(rnorm(2L * n), nrow = 2L)
  uniforms <- pmin(pmax(pnorm(latent), 1e-12), 1 - 1e-12)
  cbind(
    qnbinom(uniforms[1L, ], size = theta_pair[1L], mu = mu_pair[1L]),
    qnbinom(uniforms[2L, ], size = theta_pair[2L], mu = mu_pair[2L])
  )
}

rows <- vector("list", length(cell_numbers) * n_reps * 2L)
idx <- 1L

for (rep_id in seq_len(n_reps)) {
  for (n_cells in cell_numbers) {
    homogeneous <- n_cells * sample_nb_pair(
      n_spots, latent_rho, c(5, 5), theta
    )

    rows[[idx]] <- data.table(
      replicate = rep_id,
      cell_number = n_cells,
      setting = "Homogeneous",
      pearson = cor(homogeneous[, 1L], homogeneous[, 2L], method = "pearson"),
      spearman = cor(homogeneous[, 1L], homogeneous[, 2L], method = "spearman")
    )
    idx <- idx + 1L

    prop_a <- rbeta(n_spots, shape1 = 2, shape2 = 2)
    n_a <- n_cells * prop_a
    n_b <- n_cells * (1 - prop_a)
    profile_a <- sample_nb_pair(n_spots, latent_rho, c(2, 8), theta)
    profile_b <- sample_nb_pair(n_spots, latent_rho, c(8, 2), theta)
    mixed <- n_a * profile_a + n_b * profile_b

    rows[[idx]] <- data.table(
      replicate = rep_id,
      cell_number = n_cells,
      setting = "Composition gradient",
      pearson = cor(mixed[, 1L], mixed[, 2L], method = "pearson"),
      spearman = cor(mixed[, 1L], mixed[, 2L], method = "spearman")
    )
    idx <- idx + 1L
  }
}

raw_dt <- rbindlist(rows)
long_dt <- melt(
  raw_dt,
  id.vars = c("replicate", "cell_number", "setting"),
  measure.vars = c("pearson", "spearman"),
  variable.name = "correlation_type",
  value.name = "realized_correlation"
)
long_dt[, correlation_type := fifelse(
  correlation_type == "pearson", "Pearson", "Spearman"
)]

summary_dt <- long_dt[, .(
  mean_correlation = mean(realized_correlation, na.rm = TRUE),
  sd_correlation = sd(realized_correlation, na.rm = TRUE),
  ci_low = quantile(realized_correlation, 0.025, na.rm = TRUE),
  ci_high = quantile(realized_correlation, 0.975, na.rm = TRUE),
  n_reps = sum(is.finite(realized_correlation))
), by = .(cell_number, setting, correlation_type)]

fwrite(
  raw_dt,
  file.path(out_dir, "retained_spotlevel_composition_raw.csv")
)
fwrite(
  summary_dt,
  file.path(out_dir, "retained_spotlevel_composition_summary.csv")
)

print(summary_dt[cell_number == 20L])
