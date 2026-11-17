#' Fit negative binomial parameters
#' 
#' @param y Vector of counts.
#' @param dispersion_cap Upper bound on the negative-binomial size
#'   \eqn{\theta}; also used when the empirical variance is no larger than
#'   the mean (Poisson boundary). Set to `Inf` for no cap and the explicit
#'   Poisson limit.
#' @return Named vector with mu and theta parameters
#' @keywords internal
fit_nb_fast <- function(y, dispersion_cap = 1e6) {
  mu <- mean(y); s2 <- stats::var(y)
  th <- if (s2 > mu) min(mu^2 / (s2 - mu), dispersion_cap) else dispersion_cap
  return(c(mu = mu, theta = th))
}

#' Calculate negative binomial parameters for expression matrix
#' 
#' @param expr Expression matrix
#' @param dispersion_cap Dispersion fallback passed to the moment estimator.
#' @return List with mu and theta vectors
#' @keywords internal
calculate_nb_params <- function(expr, dispersion_cap = 1e6) {
  nb_par <- t(apply(expr, 1, fit_nb_fast,
                    dispersion_cap = dispersion_cap))
  mu_g <- nb_par[, "mu"]
  theta_g <- nb_par[, "theta"]
  return(list(mu = mu_g, theta = theta_g))
}

#' Calculate z-scores for count blocks
#' 
#' @param idx Index of genes
#' @param Y Expression matrix
#' @param mu Mean parameter vector
#' @param th Theta parameter vector  
#' @param rng_block Random number block
#' @param eps Small value to prevent zeros
#' @return Matrix of z-scores
#' @keywords internal
counts2z_block <- function(idx, Y, mu, th, rng_block, eps = 1e-6) {
  out <- matrix(0, nrow = length(idx), ncol = ncol(Y))
  for (ii in seq_along(idx)) {
    g  <- idx[ii]; y <- Y[g, ]
    # For a discrete randomized PIT, the lower endpoint is F(y - 1).
    # In particular, pnbinom(-1, ...) = 0 for an observed zero.
    u0 <- stats::pnbinom(y - 1, size = th[g], mu = mu[g])
    u1 <- stats::pnbinom(y, size = th[g], mu = mu[g])
    u  <- u0 + rng_block[ii, ] * (u1 - u0)               # DT jitter
    u  <- pmin(pmax(u, eps), 1 - eps)                    # keep in (eps,1-eps)
    out[ii, ] <- stats::qnorm(u)
  }
  return(out)
}

#' Generate z-score matrix
#' 
#' @param expr Expression matrix
#' @param mu_g Mean parameters
#' @param theta_g Dispersion parameters
#' @param block Block size for processing
#' @param n_cores Number of cores
#' @param eps Small value to prevent zeros
#' @return Matrix of z-scores
#' @keywords internal
generate_z_matrix <- function(expr, mu_g, theta_g, block = 1024, n_cores = 1, eps = 1e-6) {
  G <- nrow(expr)
  idx_ls <- split(seq_len(G), ceiling(seq_len(G) / block))
  rng <- matrix(runif(G * ncol(expr)), nrow = G)         # pre-draw all jitters
  worker <- function(idx) {
    counts2z_block(
      idx, Y = expr, mu = mu_g, th = theta_g,
      rng_block = rng[idx, , drop = FALSE], eps = eps
    )
  }
  Zparts <- if (n_cores > 1L) {
    parallel::mclapply(idx_ls, worker, mc.cores = n_cores)
  } else {
    lapply(idx_ls, worker)
  }
  return(t(do.call(rbind, Zparts)))
}

#' Calculate correlation matrix
#' 
#' @param Z Z-score matrix
#' @param thr Threshold for correlation
#' @param epsilon Minimum eigenvalue used to repair the thresholded matrix.
#' @return Correlation matrix with repair diagnostics stored in the
#'   `repair_diagnostics` attribute.
#' @keywords internal
calculate_correlation_matrix <- function(Z, thr = 0.05, epsilon = 1e-6) {
  Sigma <- stats::cor(Z, use = "pairwise.complete.obs")
  Sigma[is.na(Sigma)] <- 0                          # undefined → no corr
  diag(Sigma) <- 1                                  # proper 1’s on the diagonal
  Sigma[abs(Sigma) < thr] <- 0                      # optional sparsity
  repaired <- makespd(Sigma, epsilon = epsilon, return_diagnostics = TRUE)
  attr(repaired$rho, "repair_diagnostics") <- repaired$diagnostics
  repaired$rho
}

#' Moran's I statistics
#'
#' Computes Moran's I and the one-sided P value (alternative: positive
#' spatial autocorrelation) for each gene with [spdep::moran.test()].
#'
#' @param expr Genes-by-spots expression matrix.
#' @param listw Spatial weights list.
#' @param gene_names Genes to test.
#' @return Data frame with columns `gene`, `I`, and `p_value`.
#' @keywords internal
moran_statistics <- function(expr, listw, gene_names) {
  stats_list <- lapply(gene_names, function(g) {
    test <- tryCatch(
      spdep::moran.test(as.numeric(expr[g, ]), listw,
                        alternative = "greater", zero.policy = TRUE),
      error = function(e) NULL
    )
    if (is.null(test)) return(c(NA_real_, NA_real_))
    c(unname(test$estimate[["Moran I statistic"]]), test$p.value)
  })
  values <- do.call(rbind, stats_list)
  data.frame(gene = gene_names, I = values[, 1L], p_value = values[, 2L],
             stringsAsFactors = FALSE)
}

#' Detect spatial genes
#'
#' Selects genes with significant positive spatial autocorrelation: one-sided
#' Moran's I tests with Benjamini-Hochberg adjustment across genes.
#'
#' @param expr Expression matrix
#' @param W Spatial weights
#' @param gene_names Gene names
#' @param alpha Significance level for the BH-adjusted P values.
#' @param moran_thr Optional additional minimum Moran's I; `NULL` applies none.
#' @return Vector of spatial gene names
#' @keywords internal
detect_spatial_genes <- function(expr, W, gene_names, alpha = 0.05,
                                 moran_thr = NULL) {
  moran <- moran_statistics(expr, W, gene_names)
  adjusted <- stats::p.adjust(moran$p_value, method = "BH")
  keep <- !is.na(adjusted) & adjusted <= alpha
  if (!is.null(moran_thr)) keep <- keep & !is.na(moran$I) & moran$I > moran_thr
  spatial_genes <- moran$gene[keep]
  message(length(spatial_genes), " spatial genes detected")
  spatial_genes
}

#' Prepare mu matrix
#' 
#' @param expr Expression matrix
#' @param mu_g Mean parameters
#' @param lib_factor Library size factors
#' @return Matrix of mu values
#' @keywords internal
prepare_mu_mat <- function(expr, mu_g, lib_factor) {
  mu_mat <- outer(lib_factor, mu_g)                   # spots × genes
  dimnames(mu_mat) <- list(colnames(expr), rownames(expr))
  return(mu_mat)
}

#' Refit NB-GAM for spatial genes
#' 
#' @param expr Expression matrix
#' @param coord Coordinates data
#' @param lib_factor Library size factors
#' @param spatial_genes Spatial genes
#' @param mu_mat Mu matrix
#' @param theta_g Theta parameters
#' @param spline_k Basis dimension for the thin-plate spline.
#' @param dispersion_cap Upper bound on the refitted negative-binomial size.
#' @return Updated mu_mat and theta_g
#' @keywords internal
refit_nb_gam_for_spatial_genes <- function(expr, coord, lib_factor,
                                           spatial_genes, mu_mat, theta_g,
                                           spline_k = 50, dispersion_cap = Inf) {
  for (g in spatial_genes) {
    dat <- data.frame(
      y       = as.numeric(expr[g, ]),
      row     = coord$row,
      col     = coord$col,
      log_lib = log(lib_factor)
    )
    gam_fit <- mgcv::gam(y ~ offset(log_lib) +
                          s(row, col, bs = "tp", k = spline_k),
                        family = mgcv::nb(link = "log"),
                        data = dat, method = "REML")
    mu_mat[, g] <- stats::fitted(gam_fit)             # spot-wise means
    theta_g[g]  <- min(gam_fit$family$getTheta(TRUE), dispersion_cap)
  }
  return(list(mu_mat = mu_mat, theta_g = theta_g))
}

#' Prepare YZ input objects
#' 
#' @param expr Expression matrix
#' @param Sigma Correlation matrix
#' @param mu_mat Mu matrix
#' @param gene_names Gene names
#' @param theta_g Gene-specific negative-binomial size parameters.
#' @param matrix_diagnostics Diagnostics from correlation-matrix repair.
#' @param library_size Observed library sizes \eqn{L_j}.
#' @param spatial_genes Genes refitted with the NB-GAM.
#' @return List of mean parameters, library sizes, dispersions, dependence
#'   parameters, spatial genes, and matrix-repair diagnostics.
#' @keywords internal
prepare_input <- function(expr, Sigma, mu_mat, gene_names, theta_g,
                          matrix_diagnostics = NULL, library_size = NULL,
                          spatial_genes = character(0)) {
  mean_param <- list(CT = mu_mat)
  
  gene_relationship <- data.table::as.data.table(
    matrix(1, nrow = nrow(expr), ncol = 1,
           dimnames = list(rownames(expr), "CT"))
  )
  
  cov_gene <- array(NA, dim = c(nrow(expr), nrow(expr), 1),
                    dimnames = list(rownames(expr), rownames(expr), "CT"))
  cov_gene[,,"CT"] <- array(Sigma,
                    dim       = c(nrow(expr), nrow(expr),1),
                    dimnames  = list(rownames(expr), rownames(expr), "CT"))
  
  return(list(
    mean_param = mean_param,
    library_size = library_size,
    theta_g = list(CT = theta_g),
    gene_relationship = gene_relationship,
    cov_gene = cov_gene,
    spatial_genes = spatial_genes,
    matrix_diagnostics = matrix_diagnostics
  ))
}

#' Estimate reference-based simulation parameters
#'
#' Fits negative-binomial moments for each gene, estimates one global latent
#' correlation matrix from a randomized probability-integral transform,
#' sparsifies and repairs it, selects genes with significant spatial
#' autocorrelation (one-sided Moran's I test on a k-nearest-neighbour graph of
#' the `row`/`col` grid, BH-adjusted), and refits those genes with a
#' negative-binomial GAM with a library-size offset. The fitted means are
#' returned normalized to the average library size together with the observed
#' library sizes \eqn{L_j}; [generate_true_count()] applies
#' \eqn{f_j = L_j / \bar L} when `library_size` is supplied.
#'
#' @param expr Genes-by-spots count matrix.
#' @param coord Data frame with `row` and `col` columns, one row per spot.
#' @param gene_names Gene names.
#' @param moran_thr Optional minimum Moran's I applied in addition to the
#'   significance test; `NULL` (default) applies none.
#' @param cor_thr Correlations with absolute value below this threshold are set
#'   to zero before repair.
#' @param block_size Number of genes per block in the randomized PIT.
#' @param n_cores Number of cores for the randomized PIT.
#' @param eigen_floor Eigenvalue floor used to repair the latent correlation
#'   matrix by individual eigenvalue clipping.
#' @param neighbor_k Number of nearest neighbours used for Moran's I.
#' @param spline_k Basis dimension of the thin-plate spline.
#' @param dispersion_cap Upper bound on every negative-binomial size
#'   \eqn{\theta}, also used for genes whose sample variance does not exceed
#'   their mean. Set to `Inf` for no cap and the explicit Poisson limit.
#' @param moran_alpha Significance level for the BH-adjusted Moran's I tests.
#' @return List with `mean_param` (library-normalized spot-by-gene means, list
#'   element `CT`), `library_size` (observed \eqn{L_j}), `theta_g`
#'   (gene-specific negative-binomial sizes), `gene_relationship`, `cov_gene`
#'   (repaired latent correlation matrix), `spatial_genes`, and
#'   `matrix_diagnostics`, ready for [generate_true_count()].
#' @export
get_real_param <- function(expr, coord, gene_names, moran_thr = NULL,
                           cor_thr = 0.05, block_size = 1024, n_cores = 1,
                           eigen_floor = 1e-6, neighbor_k = 6,
                           spline_k = 50, dispersion_cap = 1e6,
                           moran_alpha = 0.05) {
  # Step 1: Calculate negative binomial parameters
  nb_params <- calculate_nb_params(expr, dispersion_cap = dispersion_cap)
  mu_g <- nb_params$mu
  theta_g <- nb_params$theta

  # Step 2: Generate z-scores
  Z <- generate_z_matrix(expr, mu_g, theta_g, block = block_size, n_cores = n_cores)

  # Step 3: Calculate the correlation matrix
  Sigma <- calculate_correlation_matrix(
    Z, thr = cor_thr, epsilon = eigen_floor
  )
  matrix_diagnostics <- attr(Sigma, "repair_diagnostics")

  # Step 4: Select genes with significant spatial autocorrelation
  coords_mat <- as.matrix(coord[, c("row", "col")])
  Wnb <- spdep::knearneigh(coords_mat, k = neighbor_k)
  W <- spdep::nb2listw(spdep::knn2nb(Wnb))
  spatial_genes <- detect_spatial_genes(
    expr, W, gene_names, alpha = moran_alpha, moran_thr = moran_thr
  )

  # Step 5: Library-size factors and baseline means
  lib_size <- colSums(expr)
  lib_factor <- lib_size / mean(lib_size)
  mu_mat <- prepare_mu_mat(expr, mu_g, lib_factor)

  # Step 6: Refit NB-GAM for spatial genes
  gam_results <- refit_nb_gam_for_spatial_genes(
    expr, coord, lib_factor, spatial_genes, mu_mat, theta_g,
    spline_k = spline_k, dispersion_cap = dispersion_cap
  )
  mu_mat <- gam_results$mu_mat
  theta_g <- gam_results$theta_g

  # Step 7: Return library-normalized means and the observed library sizes;
  # generate_true_count() restores f_j = L_j / mean(L).
  mu_mat <- sweep(mu_mat, 1L, lib_factor, `/`)
  mu_mat[!is.finite(mu_mat)] <- 0
  attr(mu_mat, "library_normalized") <- TRUE
  names(lib_size) <- colnames(expr)
  prepare_input(
    expr, Sigma, mu_mat, gene_names, theta_g,
    matrix_diagnostics = matrix_diagnostics,
    library_size = lib_size,
    spatial_genes = spatial_genes
  )
}
