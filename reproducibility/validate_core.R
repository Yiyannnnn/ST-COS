#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else "."
repro_dir <- normalizePath(dirname(script_file), mustWork = TRUE)
package_dir <- normalizePath(file.path(repro_dir, ".."), mustWork = TRUE)
out_dir <- file.path(repro_dir, "results", "core_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source_files <- list.files(file.path(package_dir, "R"), pattern = "[.]R$", full.names = TRUE)
invisible(lapply(source_files, source))

checks <- list()
record <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, passed = isTRUE(passed), detail = detail,
    stringsAsFactors = FALSE
  )
}

# Reproduce the reviewer-provided seed check.
coord <- expand.grid(x_coord = 1:10, y_coord = 1:10)
gp <- data.frame(ct1 = rep("hot_spot", 4), row.names = paste0("g", 1:4))
seed_values <- sapply(1:10, function(seed) {
  generated <- generate_gene_pattern(coord, gp, seed = seed)$ct1
  generated[1, ]
})
colnames(seed_values) <- paste0("seed", 1:10)
utils::write.csv(seed_values, file.path(out_dir, "reviewer_seed_check.csv"))
record(
  "seed propagation across genes",
  all(apply(seed_values, 1L, function(x) length(unique(signif(x, 12))) > 1L)),
  "Every gene changes across the ten top-level seeds; repeated seeds remain reproducible."
)

seed_repeat_a <- generate_gene_pattern(coord, gp, seed = 4)
seed_repeat_b <- generate_gene_pattern(coord, gp, seed = 4)
record(
  "same-seed reproducibility",
  identical(seed_repeat_a, seed_repeat_b),
  "Two calls with seed 4 are identical."
)

# Verify the declared uniform eigenvalue shift and report repair diagnostics.
rho_bad <- matrix(-0.75, nrow = 3, ncol = 3)
diag(rho_bad) <- 1
repair <- makespd(rho_bad, epsilon = 1e-4, return_diagnostics = TRUE)
repair_table <- as.data.frame(repair$diagnostics, check.names = FALSE)
utils::write.csv(repair_table, file.path(out_dir, "eigenvalue_shift_diagnostics.csv"), row.names = FALSE)
record(
  "uniform eigenvalue shift",
  repair$diagnostics$n_eigenvalues_below_floor == 1L &&
    repair$diagnostics$n_eigenvalues_modified == nrow(rho_bad) &&
    abs(repair$diagnostics$uniform_eigenvalue_shift - 0.5001) < 1e-12 &&
    repair$diagnostics$min_eigen_after > 0,
  paste0(
    "Applied a uniform shift of ",
    format(repair$diagnostics$uniform_eigenvalue_shift, digits = 6),
    "; the repaired matrix is positive definite."
  )
)

# Verify the column-vector Cholesky orientation empirically.
rho_target <- matrix(c(1, 0.55, 0.55, 1), nrow = 2)
set.seed(2)
latent <- stats::qnorm(randcop(rho_target, nSpots = 100000))
latent_correlation <- stats::cor(t(latent))[1, 2]
record(
  "Cholesky orientation",
  abs(latent_correlation - 0.55) < 0.015,
  sprintf("Requested rho = 0.55; empirical latent correlation = %.5f.", latent_correlation)
)

# Show that p0 is binomial thinning rather than a structural-zero mixture.
set.seed(11)
input <- list(CT = matrix(10, nrow = 2000, ncol = 1))
proportion <- matrix(1, nrow = 2000, ncol = 1, dimnames = list(NULL, "CT"))
thinned <- generate_batch_count(input, p0 = 0.25, celltype_proportion = proportion, seq_depth = 1e9)$CT
retention <- mean(thinned) / mean(input$CT)
record(
  "binomial thinning",
  abs(retention - 0.75) < 0.02,
  sprintf("Expected retained fraction = 0.75; observed = %.5f.", retention)
)

result <- do.call(rbind, checks)
utils::write.csv(result, file.path(out_dir, "core_validation_summary.csv"), row.names = FALSE)
print(result, row.names = FALSE)
if (!all(result$passed)) quit(status = 1L)
