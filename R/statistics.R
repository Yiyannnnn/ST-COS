#' Fit negative binomial parameters
#' 
#' @param y Vector of counts.
#' @param dispersion_cap Finite approximation to the Poisson boundary used
#'   when the empirical variance is no larger than the mean. Set to `Inf` to
#'   use the explicit Poisson limit.
#' @return Named vector with mu and theta parameters
#' @keywords internal
fit_nb_fast <- function(y, dispersion_cap = 1e6) {
  mu <- mean(y); s2 <- stats::var(y)
  th <- if (s2 > mu) mu^2 / (s2 - mu) else dispersion_cap
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
    u0 <- pnbinom(pmax(y - 1, 0), size = th[g], mu = mu[g])
    u1 <- pnbinom(y, size = th[g], mu = mu[g])
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

#' Moran's I flag
#' 
#' @param y Response variable
#' @param listw Spatial weights list
#' @param thr Threshold for flagging
#' @return Logical vector indicating flagged genes
#' @keywords internal
moran_flag <- function(y, listw, thr = 0.15) {
  S0  <- spdep::Szero(listw)                          # sum of all weights
  res <- spdep::moran(x  = y,
                      listw = listw,
                      n    = length(y),
                      S0   = S0,
                      zero.policy = TRUE)$I
  return(res > thr)
}

#' Detect spatial genes
#' 
#' @param expr Expression matrix
#' @param W Spatial weights
#' @param gene_names Gene names
#' @param moran_thr Threshold for Moran's I
#' @return Vector of spatial gene names
#' @keywords internal
detect_spatial_genes <- function(expr, W, gene_names, moran_thr = 0.15) {
  spatial_genes <- gene_names[sapply(gene_names, function(g)
    moran_flag(expr[g, ], W, thr = moran_thr))]
  cat(length(spatial_genes), "spatial genes detected\n")
  return(spatial_genes)
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
#' @return Updated mu_mat and theta_g
#' @keywords internal
refit_nb_gam_for_spatial_genes <- function(expr, coord, lib_factor,
                                           spatial_genes, mu_mat, theta_g,
                                           spline_k = 50) {
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
    theta_g[g]  <- gam_fit$family$getTheta(TRUE)      # update dispersion
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
#' @return List of mean parameters, dispersions, dependence parameters, and
#'   matrix-repair diagnostics.
#' @keywords internal
prepare_input <- function(expr, Sigma, mu_mat, gene_names, theta_g,
                          matrix_diagnostics = NULL) {
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
    theta_g = list(CT = theta_g),
    gene_relationship = gene_relationship,
    cov_gene = cov_gene,
    matrix_diagnostics = matrix_diagnostics
  ))
}

#' Main function to estimate real parameters
#' 
#' @param expr Expression matrix
#' @param coord Coordinates data
#' @param gene_names Gene names
#' @param moran_thr Moran's I threshold
#' @param cor_thr Correlation threshold
#' @param block_size Block size for processing
#' @param n_cores Number of cores for parallel processing
#' @param eigen_floor Minimum eigenvalue used to repair the latent correlation
#'   matrix.
#' @param neighbor_k Number of nearest neighbors used for Moran's I.
#' @param spline_k Basis dimension used by the thin-plate spline.
#' @param dispersion_cap Finite negative-binomial size used at the Poisson
#'   boundary. Set to `Inf` for the explicit Poisson limit.
#' @return List of estimated parameters
#' @export
get_real_param <- function(expr, coord, gene_names, moran_thr = 0.15,
                           cor_thr = 0.05, block_size = 1024, n_cores = 1,
                           eigen_floor = 1e-6, neighbor_k = 6,
                           spline_k = 50, dispersion_cap = 1e6) {
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
  
  # Step 4: Create spatial gene detection
  coords_mat <- as.matrix(coord[, c("row", "col")])
  Wnb <- spdep::knearneigh(coords_mat, k = neighbor_k)
  W <- spdep::nb2listw(spdep::knn2nb(Wnb))
  spatial_genes <- detect_spatial_genes(expr, W, gene_names, moran_thr)
  
  # Step 5: Prepare mu_mat
  lib_size <- colSums(expr)
  lib_factor <- lib_size / mean(lib_size)
  mu_mat <- prepare_mu_mat(expr, mu_g, lib_factor)
  
  # Step 6: Refit NB-GAM for spatial genes
  gam_results <- refit_nb_gam_for_spatial_genes(
    expr, coord, lib_factor, spatial_genes, mu_mat, theta_g,
    spline_k = spline_k
  )
  mu_mat <- gam_results$mu_mat
  theta_g <- gam_results$theta_g
  
  # Step 7: Prepare YZ input objects
  real_param <- prepare_input(
    expr, Sigma, mu_mat, gene_names, theta_g,
    matrix_diagnostics = matrix_diagnostics
  )
  
  return(real_param)
}
