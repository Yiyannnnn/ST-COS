#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- "."
out_dir <- file.path(base_dir, "Independent_simulator_anchor", "count_correlation_attenuation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260902)

rhos <- c(0, 0.2, 0.4, 0.6, 0.8)
n_samples <- 5000L
n_reps <- 50L

mean_grid <- data.table(
  mean_level = factor(c("Low mean", "Medium mean", "High mean"),
                      levels = c("Low mean", "Medium mean", "High mean")),
  mu = c(1, 5, 20)
)

# alpha is the NB overdispersion in Var(Y) = mu + alpha * mu^2.
# The R/ST-COS qnbinom size parameter is theta = 1 / alpha.
disp_grid <- data.table(
  dispersion_level = factor(c("Low dispersion", "Medium dispersion", "High dispersion"),
                            levels = c("Low dispersion", "Medium dispersion", "High dispersion")),
  alpha = c(0.05, 0.20, 1.00)
)
disp_grid[, theta := 1 / alpha]

# Additional independent zero inflation is used to create explicit sparsity levels.
# The realized zero proportion is recorded in the output tables.
zero_grid <- data.table(
  zero_level = factor(c("NB only", "Moderate extra zeros", "High extra zeros"),
                      levels = c("NB only", "Moderate extra zeros", "High extra zeros")),
  zero_inflation = c(0.00, 0.25, 0.50)
)

simulate_pair <- function(rho, mu, theta, zero_inflation, n) {
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  upper_chol <- chol(sigma)
  z <- matrix(rnorm(2L * n), nrow = 2L)
  w <- t(upper_chol) %*% z
  u <- pnorm(w)

  y1 <- qnbinom(u[1, ], size = theta, mu = mu)
  y2 <- qnbinom(u[2, ], size = theta, mu = mu)

  if (zero_inflation > 0) {
    y1[runif(n) < zero_inflation] <- 0L
    y2[runif(n) < zero_inflation] <- 0L
  }

  pearson <- suppressWarnings(cor(y1, y2, method = "pearson"))
  spearman <- suppressWarnings(cor(y1, y2, method = "spearman"))

  data.table(
    pearson = pearson,
    spearman = spearman,
    zero_prop = mean(c(y1, y2) == 0),
    mean_y1 = mean(y1),
    mean_y2 = mean(y2),
    var_y1 = var(y1),
    var_y2 = var(y2)
  )
}

message("Running bivariate NB copula simulations...")
rows <- vector("list", length(rhos) * nrow(mean_grid) * nrow(disp_grid) * nrow(zero_grid) * n_reps)
idx <- 1L

for (rho in rhos) {
  for (m in seq_len(nrow(mean_grid))) {
    for (d in seq_len(nrow(disp_grid))) {
      for (zlev in seq_len(nrow(zero_grid))) {
        for (rep_id in seq_len(n_reps)) {
          sim <- simulate_pair(
            rho = rho,
            mu = mean_grid$mu[m],
            theta = disp_grid$theta[d],
            zero_inflation = zero_grid$zero_inflation[zlev],
            n = n_samples
          )
          rows[[idx]] <- cbind(
            data.table(
              latent_rho = rho,
              rep = rep_id,
              mean_level = mean_grid$mean_level[m],
              mu = mean_grid$mu[m],
              dispersion_level = disp_grid$dispersion_level[d],
              alpha = disp_grid$alpha[d],
              theta = disp_grid$theta[d],
              zero_level = zero_grid$zero_level[zlev],
              zero_inflation = zero_grid$zero_inflation[zlev]
            ),
            sim
          )
          idx <- idx + 1L
        }
      }
    }
  }
}

raw_dt <- rbindlist(rows, use.names = TRUE)
raw_dt[, `:=`(
  mean_level = factor(mean_level, levels = levels(mean_grid$mean_level)),
  dispersion_level = factor(dispersion_level, levels = levels(disp_grid$dispersion_level)),
  zero_level = factor(zero_level, levels = levels(zero_grid$zero_level))
)]

long_dt <- melt(
  raw_dt,
  id.vars = c("latent_rho", "rep", "mean_level", "mu", "dispersion_level",
              "alpha", "theta", "zero_level", "zero_inflation", "zero_prop"),
  measure.vars = c("pearson", "spearman"),
  variable.name = "correlation_type",
  value.name = "realized_r"
)
long_dt[, correlation_type := factor(
  fifelse(correlation_type == "pearson", "Pearson", "Spearman"),
  levels = c("Pearson", "Spearman")
)]

summary_dt <- long_dt[, .(
  mean_r = mean(realized_r, na.rm = TRUE),
  sd_r = sd(realized_r, na.rm = TRUE),
  q025_r = quantile(realized_r, 0.025, na.rm = TRUE),
  q975_r = quantile(realized_r, 0.975, na.rm = TRUE),
  mean_zero_prop = mean(zero_prop, na.rm = TRUE),
  sd_zero_prop = sd(zero_prop, na.rm = TRUE),
  n_reps = .N
), by = .(
  latent_rho, mean_level, mu, dispersion_level, alpha, theta,
  zero_level, zero_inflation, correlation_type
)]

zero_summary <- unique(summary_dt[, .(
  mean_level, mu, dispersion_level, alpha, theta,
  zero_level, zero_inflation, mean_zero_prop, sd_zero_prop
)])

fwrite(raw_dt, file.path(out_dir, "count_correlation_attenuation_raw.csv"))
fwrite(summary_dt, file.path(out_dir, "count_correlation_attenuation_summary.csv"))
fwrite(zero_summary, file.path(out_dir, "count_correlation_zero_proportion_summary.csv"))

plot_theme <- theme_bw(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", color = "#6f6f6f", linewidth = 0.35),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "#d2e0e3", color = "#6f6f6f", linewidth = 0.35),
    strip.text = element_text(size = 10.5, face = "bold"),
    panel.grid.major = element_line(color = "#e5e5e5", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 9.5),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "#444444"),
    plot.margin = margin(8, 8, 8, 8)
  )

zero_colors <- c(
  "NB only" = "#0072B2",
  "Moderate extra zeros" = "#D55E00",
  "High extra zeros" = "#009E73"
)

pearson_dt <- summary_dt[correlation_type == "Pearson"]
p_pearson <- ggplot(
  pearson_dt,
  aes(x = latent_rho, y = mean_r, color = zero_level, group = zero_level)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#555555", linewidth = 0.45) +
  geom_errorbar(aes(ymin = q025_r, ymax = q975_r), width = 0.025, alpha = 0.55, linewidth = 0.35) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.8) +
  facet_grid(mean_level ~ dispersion_level) +
  scale_color_manual(values = zero_colors, drop = FALSE) +
  scale_x_continuous(breaks = rhos, limits = c(0, 0.8)) +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(-0.03, 0.84)) +
  labs(
    title = "Latent Gaussian copula correlation vs realized NB count Pearson correlation",
    subtitle = sprintf(
      "Each point averages %d replicates of %s simulated pairs; intervals are empirical 2.5%%-97.5%% quantiles. Dotted line: y = x.",
      n_reps, format(n_samples, big.mark = ",")
    ),
    x = "Specified latent Gaussian correlation",
    y = "Realized count Pearson r",
    color = "Sparsity setting"
  ) +
  plot_theme

p_both <- ggplot(
  summary_dt,
  aes(x = latent_rho, y = mean_r, color = zero_level,
      linetype = correlation_type, shape = correlation_type,
  group = interaction(zero_level, correlation_type))
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#555555", linewidth = 0.45) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.55) +
  facet_grid(mean_level ~ dispersion_level) +
  scale_color_manual(values = zero_colors, drop = FALSE) +
  scale_linetype_manual(values = c("Pearson" = "solid", "Spearman" = "22"), drop = FALSE) +
  scale_shape_manual(values = c("Pearson" = 16, "Spearman" = 17), drop = FALSE) +
  scale_x_continuous(breaks = rhos, limits = c(0, 0.8)) +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(-0.03, 0.84)) +
  labs(
    title = "Latent Gaussian copula correlation vs realized NB count correlation",
    subtitle = "Attenuation increases with sparsity; additional zeros are most damaging at higher means.",
    x = "Specified latent Gaussian correlation",
    y = "Realized count correlation",
    color = "Sparsity setting",
    linetype = "Correlation",
    shape = "Correlation"
  ) +
  plot_theme

ggsave(file.path(out_dir, "latent_vs_realized_count_pearson.pdf"), p_pearson, width = 11, height = 8.5, bg = "white")
ggsave(file.path(out_dir, "latent_vs_realized_count_pearson.png"), p_pearson, width = 11, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "latent_vs_realized_count_pearson_spearman.pdf"), p_both, width = 11.5, height = 8.8, bg = "white")
ggsave(file.path(out_dir, "latent_vs_realized_count_pearson_spearman.png"), p_both, width = 11.5, height = 8.8, dpi = 300, bg = "white")

message("Saved outputs to: ", out_dir)
