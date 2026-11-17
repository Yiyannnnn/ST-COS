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

test_that("matrix repair clips only eigenvalues below the floor", {
  rho <- matrix(-0.75, nrow = 3, ncol = 3)
  diag(rho) <- 1
  repaired <- STCOS:::makespd(
    rho, epsilon = 1e-4, return_diagnostics = TRUE
  )

  expect_equal(repaired$diagnostics$n_eigenvalues_below_floor, 1)
  expect_equal(repaired$diagnostics$n_eigenvalues_modified, 1)
  expect_equal(repaired$diagnostics$minimum_eigenvalue_lift, 0.5001)
  expect_equal(repaired$diagnostics$total_eigenvalue_lift, 0.5001)
  expect_identical(
    repaired$diagnostics$repair_method,
    "individual_eigenvalue_clipping"
  )
  expect_gt(min(eigen(repaired$rho, symmetric = TRUE)$values), 0)
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
    "mean_param", "library_size", "theta_g", "gene_relationship", "cov_gene",
    "spatial_genes", "matrix_diagnostics"
  ))
  expect_equal(unname(fitted$library_size), unname(colSums(expression)))
  expect_length(fitted$theta_g$CT, nrow(expression))
  expect_true(all(c(
    "min_eigen_before", "min_eigen_after", "minimum_eigenvalue_lift",
    "total_eigenvalue_lift", "n_eigenvalues_below_floor",
    "n_eigenvalues_modified", "relative_frobenius_perturbation",
    "maximum_entrywise_change", "epsilon", "repair_method"
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

test_that("gene-pattern seeds control every gene without inner resets", {
  coordinates <- expand.grid(x_coord = 1:10, y_coord = 1:10)
  patterns <- data.frame(
    CT = rep("hot_spot", 4), row.names = paste0("g", 1:4)
  )

  first_spot <- sapply(1:10, function(seed) {
    generate_gene_pattern(coordinates, patterns, seed = seed)$CT[1, ]
  })

  expect_true(all(apply(first_spot, 1, function(x) length(unique(x)) > 1L)))
  expect_equal(
    generate_gene_pattern(coordinates, patterns, seed = 4),
    generate_gene_pattern(coordinates, patterns, seed = 4)
  )
})

test_that("pattern parameters can vary by cell type", {
  coordinates <- data.frame(x_coord = 1:10, y_coord = 1)
  patterns <- data.frame(
    A = "streak", B = "streak", row.names = "g1",
    check.names = FALSE
  )
  means <- generate_gene_pattern(
    coordinates,
    patterns,
    streak_x = data.frame(A = 2, B = 8),
    streak_x_param = 10,
    background_param = 1,
    streak_size = 2,
    seed = 1
  )

  expect_equal(unname(which.max(means$A[, "g1"])), 2)
  expect_equal(unname(which.max(means$B[, "g1"])), 8)
})

test_that("pattern parameters can vary by gene and cell type", {
  coordinates <- data.frame(x_coord = 1:10, y_coord = 1)
  patterns <- data.frame(
    A = c("streak", "streak"), B = c("streak", "streak"),
    row.names = c("g1", "g2"), check.names = FALSE
  )
  means <- generate_gene_pattern(
    coordinates,
    patterns,
    streak_x = data.frame(A = c(2, 4), B = c(7, 9)),
    streak_x_param = 10,
    background_param = 1,
    streak_size = 2,
    seed = 1
  )

  expect_equal(
    unname(c(which.max(means$A[, "g1"]), which.max(means$A[, "g2"]))),
    c(2L, 4L)
  )
  expect_equal(
    unname(c(which.max(means$B[, "g1"]), which.max(means$B[, "g2"]))),
    c(7L, 9L)
  )
})

test_that("deterministic pattern construction does not reset the RNG", {
  coordinates <- expand.grid(x_coord = 1:3, y_coord = 1:3)
  set.seed(17)
  state_before <- .Random.seed
  generate_pattern(
    coordinates, "hot_spot",
    center_coord = data.frame(x_coord = 2, y_coord = 2),
    center_param = 10, background_param = 2, hot_spot_size = 2
  )

  expect_identical(.Random.seed, state_before)
})

test_that("randomized PIT uses F(-1) = 0 for zero counts", {
  jitter <- matrix(seq(0.01, 0.99, length.out = 100), nrow = 1)
  counts <- matrix(0, nrow = 1, ncol = ncol(jitter))
  transformed <- STCOS:::counts2z_block(
    1L, counts, mu = 0.2, th = 5,
    rng_block = jitter, eps = 1e-12
  )
  expected_u <- as.numeric(jitter) * stats::pnbinom(0, size = 5, mu = 0.2)

  expect_equal(stats::pnorm(as.numeric(transformed)), expected_u,
               tolerance = 1e-12)
  expect_gt(diff(range(transformed)), 1)
})

test_that("randomized PIT matches analytical intervals for zero-heavy counts", {
  counts_vector <- c(rep(0L, 80), rep(1L, 12), rep(2L, 6), 3L, 4L)
  counts <- matrix(counts_vector, nrow = 1)
  jitter <- matrix(seq(0.005, 0.995, length.out = length(counts_vector)),
                   nrow = 1)
  transformed <- STCOS:::counts2z_block(
    1L, counts, mu = 0.2, th = 5,
    rng_block = jitter, eps = 1e-12
  )
  lower <- stats::pnbinom(counts_vector - 1L, size = 5, mu = 0.2)
  upper <- stats::pnbinom(counts_vector, size = 5, mu = 0.2)
  expected_u <- lower + as.numeric(jitter) * (upper - lower)

  expect_equal(stats::pnorm(as.numeric(transformed)), expected_u,
               tolerance = 1e-12)
  expect_true(all(stats::pnorm(transformed[counts_vector == 0L]) > 0))
  expect_true(all(stats::pnorm(transformed[counts_vector == 0L]) <
                    stats::pnbinom(0, size = 5, mu = 0.2)))
})

test_that("multiple streaks are generated without missing values", {
  coordinates <- data.frame(x_coord = 0:10, y_coord = 0)
  values <- generate_pattern(
    coordinates, "streak",
    background_param = 1,
    streak_x = c(2, 8),
    streak_x_param = c(5, 9),
    streak_size = c(2, 4)
  )

  expect_true(all(is.finite(values)))
  expect_equal(values[coordinates$x_coord == 2], 5)
  expect_equal(values[coordinates$x_coord == 8], 9)
})

test_that("cell-type mean logits are converted by softmax", {
  layer_coordinates <- data.frame(layer = c(1, 1, 2, 2))
  target_weights <- data.frame(A = c(2, 0), B = c(1, 3))
  proportions <- generate_cell_prop(
    layer_coordinates, target_weights, noise_sd = 0
  )

  expect_equal(rowSums(proportions), rep(1, 4))
  expect_equal(
    proportions[1, ],
    c(A = exp(2) / (exp(2) + exp(1)), B = exp(1) / (exp(2) + exp(1)))
  )
  expect_equal(
    proportions[3, ],
    c(A = exp(0) / (exp(0) + exp(3)), B = exp(3) / (exp(0) + exp(3)))
  )
  expect_true(all(proportions >= 0))
})

test_that("log-weight inputs recover intended composition regimes", {
  layer_coordinates <- data.frame(layer = c(1, 2, 3))
  weights <- rbind(
    c(A = 1 / 2, B = 1 / 4, C = 1 / 3),
    c(A = 0, B = 3 / 4, C = 0),
    c(A = 0, B = 0, C = 1)
  )
  logits <- weights
  logits[logits <= 0] <- 1e-8
  logits <- log(logits)
  proportions <- generate_cell_prop(
    layer_coordinates, logits, noise_sd = 0
  )
  expected <- pmax(weights, 1e-8)
  expected <- expected / rowSums(expected)

  expect_equal(proportions, expected, tolerance = 1e-12)
})

test_that("generate_batch_count performs per-molecule binomial thinning", {
  input <- list(CT = matrix(100, nrow = 100, ncol = 10))
  proportions <- matrix(
    1, nrow = 100, ncol = 1,
    dimnames = list(NULL, "CT")
  )

  set.seed(18)
  thinned <- generate_batch_count(input, p0 = 0.25)$CT

  expect_equal(mean(thinned) / mean(input$CT), 0.75, tolerance = 0.01)
  expect_error(generate_batch_count(input, p0 = 1.1), "between zero and one")
  expect_error(generate_batch_count(input, p0 = 0.1, seq_depth = 100),
               "generate_spot_count")
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

test_that("fitted dispersions propagate to the count generator", {
  means <- list(A = matrix(20, nrow = 4000, ncol = 1))
  coordinates <- matrix(0, nrow = 4000, ncol = 2)
  set.seed(21)
  overdispersed <- generate_true_count(coordinates, means, theta_g = list(A = 1))$A
  near_poisson <- generate_true_count(coordinates, means, theta_g = list(A = 1e6))$A

  # NB variance is mu + mu^2 / theta: 420 for theta = 1 and about 20 for 1e6.
  expect_equal(var(overdispersed[, 1]), 420, tolerance = 0.15)
  expect_equal(var(near_poisson[, 1]), 20, tolerance = 0.15)
})

test_that("generate_covariance() output and named lists are accepted", {
  genes <- paste0("gene_", 1:4)
  relationship <- data.frame(A = c(1, 1, NA, NA), B = c(NA, NA, 1, 1),
                             row.names = genes)
  means <- list(A = matrix(5, 6, 4, dimnames = list(NULL, genes)),
                B = matrix(5, 6, 4, dimnames = list(NULL, genes)))
  set.seed(22)
  covariance <- generate_covariance(relationship)
  expect_type(covariance, "list")
  out <- generate_true_count(matrix(0, 6, 2), means,
                             gene_relationship = relationship,
                             cov_gene = covariance)
  expect_named(out, c("A", "B"))
  expect_equal(dim(out$B), c(6L, 4L))

  expect_error(
    generate_true_count(matrix(0, 6, 2), means,
                        gene_relationship = relationship,
                        cov_gene = list(A = diag(3), B = diag(3))),
    "must be a 4 x 4 matrix"
  )
})

test_that("negative or non-finite means are rejected", {
  means <- list(A = matrix(c(1, -1), nrow = 2, ncol = 1))
  expect_error(generate_true_count(matrix(0, 2, 2), means),
               "negative or non-finite")
})

test_that("generate_layer keeps row order and requires full coverage", {
  coordinates <- data.frame(x_coord = c(8, 1, 6, 3), y_coord = c(1, 1, 9, 9))
  layered <- generate_layer(coordinates, xmin = c(0, 5), xmax = c(5, 10),
                            ymin = c(0, 0), ymax = c(10, 10))
  expect_equal(layered$x_coord, coordinates$x_coord)
  expect_equal(layered$layer, c(2L, 1L, 2L, 1L))

  expect_error(generate_layer(coordinates, xmin = 0, xmax = 5, ymin = 0, ymax = 10),
               "2 of 4 spots fall outside")
  expect_error(generate_layer(coordinates, xmin = c(0, 5), xmax = 10,
                              ymin = 0, ymax = 10),
               "same length")
})

test_that("layer patterns require one mean per layer", {
  coordinates <- data.frame(x_coord = 1:4, y_coord = 1:4)
  expect_error(
    generate_pattern(coordinates, "layer", layer = c(1, 1, 2, 3),
                     layer_param = c(5, 6)),
    "one mean per layer"
  )
})

test_that("uniform bounds and layer means can be gene specific", {
  coordinates <- data.frame(x_coord = 1:4, y_coord = 1:4)
  patterns <- data.frame(A = c("no_pattern", "no_pattern", "layer", "layer"),
                         row.names = paste0("g", 1:4))
  means <- generate_gene_pattern(
    coordinates, patterns,
    all_lower_mean = c(1, 100, NA, NA), all_upper_mean = c(2, 101, NA, NA),
    layer_assignments = c("L1", "L1", "L2", "L2"),
    layer_param = list(A = rbind(c(3, 4), c(30, 40))),
    seed = 5
  )
  expect_true(all(means$A[, "g1"] >= 1 & means$A[, "g1"] <= 2))
  expect_true(all(means$A[, "g2"] >= 100 & means$A[, "g2"] <= 101))
  expect_equal(unname(means$A[, "g3"]), c(3, 3, 4, 4))
  expect_equal(unname(means$A[, "g4"]), c(30, 30, 40, 40))

  expect_error(
    generate_gene_pattern(coordinates, patterns[1, , drop = FALSE],
                          all_lower_mean = 5, all_upper_mean = 1),
    "0 <= lower <= upper"
  )
})

test_that("seed arguments reset the generator only when supplied", {
  set.seed(30)
  reference <- runif(1)
  set.seed(30)
  invisible(generate_coordinates(10, 10, "grid", spot_distance = 1))
  expect_identical(runif(1), reference)

  first <- generate_coordinates(10, 10, "random", n_spot = 5)
  second <- generate_coordinates(10, 10, "random", n_spot = 5)
  expect_false(identical(first, second))
  expect_identical(
    generate_coordinates(10, 10, "random", n_spot = 5, seed = 7),
    generate_coordinates(10, 10, "random", n_spot = 5, seed = 7)
  )
})

test_that("spots above the depth threshold are resampled to exactly that depth", {
  profiles <- list(A = matrix(100, nrow = 3, ncol = 10),
                   B = matrix(100, nrow = 3, ncol = 10))
  proportions <- cbind(A = rep(0.5, 3), B = rep(0.5, 3))
  cell_number <- c(10, 0.2, 10)
  set.seed(31)
  mixed <- generate_spot_count(profiles, 10, cell_number = cell_number,
                               celltype_proportion = proportions,
                               seq_depth = 500)

  # Spots 1 and 3 total 10,000 > 500 and are resampled to 500 counts;
  # spot 2 totals 200 <= 500 and keeps its mixed profile.
  expect_equal(unname(rowSums(mixed)), c(500, 200, 500))
  expect_true(all(mixed[c(1, 3), ] == round(mixed[c(1, 3), ])))
  expect_equal(unname(mixed[2, ]), rep(20, 10))
})

test_that("image mapping preserves aspect ratio and exports a y-up frame", {
  layout <- data.frame(x_coord = c(0, 10, 0, 10), y_coord = c(0, 0, 5, 5))
  mapped <- STCOS:::map_coordinates_to_image(layout, width = 400, height = 400)
  expect_equal(diff(range(mapped$x_coord)) / diff(range(mapped$y_coord)), 2)
  expect_equal(range(mapped$x_coord), c(0, 400))
  expect_equal(range(mapped$y_coord), c(100, 300))
  # The bottom-left layout point is lower on the image (larger pixel y).
  expect_gt(mapped$y_coord[1], mapped$y_coord[3])

  exported <- STCOS:::image_to_simulation_coordinates(mapped, height = 400)
  expect_equal(exported$y_coord - min(exported$y_coord),
               (layout$y_coord - min(layout$y_coord)) * 40)
})

test_that("hotspot distances ignore extra coordinate columns", {
  coordinates <- data.frame(x_coord = c(0, 3, 0), y_coord = c(0, 4, 2),
                            cell_type = c("A", "B", "A"))
  expect_equal(STCOS:::cdist(coordinates, data.frame(x_coord = 0, y_coord = 0)),
               c(0, 5, 2))
  patterns <- data.frame(A = "hotspot", row.names = "g1")
  means <- generate_gene_pattern(
    coordinates, patterns,
    center_coord = list(A = data.frame(x_coord = 0, y_coord = 0)),
    center_param = 10, background_param = 2, hot_spot_size = 4
  )
  expect_equal(unname(means$A[, "g1"]), c(10, 2, 6))
})


test_that("composition logits are drawn once per layer", {
  layer_coordinates <- data.frame(layer = rep(1:3, each = 20))
  logits <- data.frame(A = c(0, 1, 2), B = c(1, 0, -1))
  set.seed(40)
  proportions <- generate_cell_prop(layer_coordinates, logits, noise_sd = 0.5)
  for (layer in 1:3) {
    rows <- proportions[layer_coordinates$layer == layer, , drop = FALSE]
    expect_equal(nrow(unique(round(rows, 12))), 1L)
  }
})

test_that("spatial genes are selected by a Moran's I significance test", {
  set.seed(41)
  grid <- expand.grid(row = 1:15, col = 1:15)
  spatial <- rpois(nrow(grid), lambda = ifelse(grid$col > 8, 40, 5))
  noise <- matrix(rpois(4 * nrow(grid), lambda = 10), nrow = 4)
  expression <- rbind(spatial, noise)
  rownames(expression) <- paste0("g", 1:5)
  listw <- spdep::nb2listw(spdep::knn2nb(
    spdep::knearneigh(as.matrix(grid), k = 6)
  ))

  stats <- STCOS:::moran_statistics(expression, listw, rownames(expression))
  expect_lt(stats$p_value[1], 1e-10)
  selected <- suppressMessages(
    STCOS:::detect_spatial_genes(expression, listw, rownames(expression))
  )
  expect_true("g1" %in% selected)
  expect_lte(sum(c("g2", "g3", "g4", "g5") %in% selected), 1L)
})

test_that("observed library sizes are applied by the generator", {
  set.seed(42)
  expression <- matrix(rnbinom(5 * 64, size = 5, mu = 20), nrow = 5,
                       dimnames = list(paste0("g", 1:5), paste0("s", 1:64)))
  expression[, 1:32] <- expression[, 1:32] * 3
  coordinates <- data.frame(row = rep(1:8, 8), col = rep(1:8, each = 8))
  fitted <- suppressMessages(get_real_param(
    expression, coordinates, rownames(expression), spline_k = 20,
    moran_thr = 2
  ))
  # Without spatial genes, the normalized mean of gene g is its sample mean in
  # every spot; the generator restores f_j = L_j / mean(L).
  normalized <- fitted$mean_param$CT
  expect_true(isTRUE(attr(normalized, "library_normalized")))
  expect_equal(unname(normalized[1, ]), unname(rowMeans(expression)))
  expect_lt(max(apply(normalized, 2, sd)), 1e-10)

  means <- list(A = matrix(20, nrow = 2000, ncol = 1))
  libraries <- rep(c(3, 1), each = 1000)
  set.seed(43)
  simulated <- generate_true_count(
    matrix(0, 2000, 2), means, theta_g = list(A = 1e6), library_size = libraries
  )$A
  expect_equal(mean(simulated[1:1000]), 30, tolerance = 0.02)
  expect_equal(mean(simulated[1001:2000]), 10, tolerance = 0.03)

  expect_warning(
    generate_true_count(coordinates, fitted$mean_param, fitted$gene_relationship,
                        fitted$cov_gene, theta_g = fitted$theta_g),
    "library-normalized"
  )
})

test_that("the dispersion cap bounds every moment estimate", {
  y <- c(rep(10, 999), 11, 9)  # variance below the mean: Poisson boundary
  expect_equal(unname(STCOS:::fit_nb_fast(y, dispersion_cap = 1e6)["theta"]), 1e6)
  y <- rep(c(0, 2), 500)       # mean 1, variance ~1.001: theta ~ 1000
  expect_lte(unname(STCOS:::fit_nb_fast(y, dispersion_cap = 100)["theta"]), 100)
  expect_gt(unname(STCOS:::fit_nb_fast(y, dispersion_cap = Inf)["theta"]), 100)
})

test_that("Cholesky orientation reproduces a nonexchangeable 5 x 5 target", {
  rho <- matrix(c(
    1.0, 0.8, 0.5, 0.2, 0.0,
    0.8, 1.0, 0.6, 0.3, 0.1,
    0.5, 0.6, 1.0, 0.4, 0.2,
    0.2, 0.3, 0.4, 1.0, 0.7,
    0.0, 0.1, 0.2, 0.7, 1.0
  ), nrow = 5)
  set.seed(44)
  latent <- stats::qnorm(STCOS:::randcop(rho, nSpots = 1e6))
  correct <- stats::cor(t(latent))
  expect_lt(sqrt(mean((correct - rho)^2)), 0.003)

  upper <- chol(rho)
  wrong <- upper %*% matrix(stats::rnorm(5e6), nrow = 5)
  expect_gt(sqrt(mean((stats::cor(t(wrong)) - rho)^2)), 0.03)
})
