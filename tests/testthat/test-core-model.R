test_that("spot aggregation is a continuous abundance-weighted mixture", {
  type_a <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE,
                   dimnames = list(c("s1", "s2"), c("g1", "g2")))
  type_b <- matrix(c(5, 6, 7, 8), nrow = 2, byrow = TRUE,
                   dimnames = list(c("s1", "s2"), c("g1", "g2")))
  proportions <- matrix(c(0.25, 0.75, 0.60, 0.40), nrow = 2, byrow = TRUE,
                        dimnames = list(c("s1", "s2"), c("A", "B")))

  observed <- generate_spot_count(
    list(A = type_a, B = type_b), nGenes = 2,
    cell_number = c(10, 5), celltype_proportion = proportions
  )
  expected <- type_a * c(2.5, 3) + type_b * c(7.5, 2)

  expect_equal(observed, expected)
})

test_that("count generators consistently return spots by genes", {
  set.seed(1)
  means <- list(A = matrix(5, nrow = 4, ncol = 3,
                           dimnames = list(paste0("s", 1:4),
                                           paste0("g", 1:3))))
  independent <- generate_true_count(
    matrix(0, nrow = 4, ncol = 2), means,
    lib_size = list(A = rep(20, 4)), theta_g = list(A = rep(5, 3))
  )

  relationship <- data.frame(A = rep(1, 3), row.names = paste0("g", 1:3))
  covariance <- array(
    c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1),
    dim = c(3, 3, 1),
    dimnames = list(paste0("g", 1:3), paste0("g", 1:3), "A")
  )
  correlated <- generate_true_count(
    matrix(0, nrow = 4, ncol = 2), means,
    gene_relationship = relationship, cov_gene = covariance,
    lib_size = list(A = rep(20, 4)), theta_g = list(A = rep(5, 3))
  )

  expect_equal(dim(independent$A), c(4, 3))
  expect_equal(dim(correlated$A), c(4, 3))
})

test_that("Cholesky sampling targets the requested latent correlation", {
  rho <- matrix(c(1, 0.55, 0.55, 1), nrow = 2)
  set.seed(2)
  uniforms <- STCOS:::randcop(rho, nSpots = 10000)
  latent <- stats::qnorm(uniforms)

  expect_equal(stats::cor(t(latent))[1, 2], 0.55, tolerance = 0.03)
})

test_that("reference fitting exposes reviewed defaults and returns dispersion", {
  set.seed(3)
  expression <- matrix(
    stats::rnbinom(180, size = 5, mu = 4), nrow = 6,
    dimnames = list(paste0("g", 1:6), paste0("s", 1:30))
  )
  coordinates <- data.frame(
    row = rep(1:5, 6), col = rep(1:6, each = 5)
  )
  fitted <- get_real_param(
    expression, coordinates, rownames(expression),
    moran_thr = 2, block_size = 2, eigen_floor = 1e-6,
    neighbor_k = 6, spline_k = 50, dispersion_cap = 1e6
  )

  expect_named(fitted, c(
    "mean_param", "theta_g", "gene_relationship", "cov_gene",
    "matrix_diagnostics"
  ))
  expect_length(fitted$theta_g$CT, nrow(expression))
  expect_true(all(c(
    "min_eigen_before", "min_eigen_after", "minimum_eigenvalue_lift",
    "n_eigenvalues_modified", "relative_frobenius_perturbation",
    "maximum_entrywise_change", "epsilon"
  ) %in% names(fitted$matrix_diagnostics)))
})

test_that("gene-pattern generation does not rely on a global gene count", {
  coordinates <- data.frame(x_coord = 1:4, y_coord = 1:4)
  patterns <- data.frame(CT = c("no_pattern", "layer"),
                         row.names = c("g1", "g2"))
  means <- generate_gene_pattern(
    coordinates, patterns, layer_assignments = c("L1", "L1", "L2", "L2"),
    layer_param = c(2, 8), seed = 1
  )

  expect_equal(dim(means$CT), c(4, 2))
  expect_equal(unname(means$CT[, "g2"]), c(2, 2, 8, 8))
})

test_that("image-guided app launcher exposes the hosted application", {
  expect_identical(
    launch_stcos_app(browser = FALSE),
    "https://yiyanz.shinyapps.io/st-cos/"
  )
  expect_error(launch_stcos_app(browser = NA), "TRUE or FALSE")
})

test_that("local image-guided app is a Shiny application", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("png")
  skip_if_not_installed("jpeg")

  expect_s3_class(stcos_shiny_app(), "shiny.appobj")
})
