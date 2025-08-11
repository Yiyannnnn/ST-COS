#' Fit negative binomial parameters
#' 
#' @param y Vector of counts
#' @return Named vector with mu and theta parameters
#' @keywords internal
fit_nb_fast <- function(y) {
  mu <- mean(y); s2 <- var(y)
  th <- if (s2 > mu) mu^2 / (s2 - mu) else 1e6
  return(c(mu = mu, theta = th))
}

#' Calculate negative binomial parameters for expression matrix
#' 
#' @param expr Expression matrix
#' @return List with mu and theta vectors
#' @keywords internal
calculate_nb_params <- function(expr) {
  nb_par <- t(apply(expr, 1, fit_nb_fast))
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
    out[ii, ] <- qnorm(u)
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
  Zparts <- mclapply(idx_ls, counts2z_block,
                     Y = expr, mu = mu_g, th = theta_g,
                     rng_block = rng[idx_ls[[1]], , drop = FALSE],
                     eps = eps, mc.cores = n_cores)
  return(t(do.call(rbind, Zparts)))
}

#' Calculate correlation matrix
#' 
#' @param Z Z-score matrix
#' @param thr Threshold for correlation
#' @return Correlation matrix
#' @keywords internal
calculate_correlation_matrix <- function(Z, thr = 0.05) {
  Sigma <- cor(Z, use = "pairwise.complete.obs")  # handle constant-variance genes
  Sigma[is.na(Sigma)] <- 0                          # undefined → no corr
  diag(Sigma) <- 1                                  # proper 1’s on the diagonal
  Sigma[abs(Sigma) < thr] <- 0                      # optional sparsity
  return(makespd(Sigma))                         # repair PD if needed
}

#' Moran's I flag
#' 
#' @param y Response variable
#' @param listw Spatial weights list
#' @param thr Threshold for flagging
#' @return Logical vector indicating flagged genes
#' @keywords internal
moran_flag <- function(y, listw, thr = 0.15) {
  S0  <- Szero(listw)                          # sum of all weights
  res <- moran(x  = y,
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
#' @return Updated mu_mat and theta_g
#' @keywords internal
refit_nb_gam_for_spatial_genes <- function(expr, coord, lib_factor, spatial_genes, mu_mat, theta_g) {
  for (g in spatial_genes) {
    dat <- data.frame(
      y       = as.numeric(expr[g, ]),
      row     = coord$row,
      col     = coord$col,
      log_lib = log(lib_factor)
    )
    gam_fit <- gam(y ~ offset(log_lib) +
                     s(row, col, bs = "tp", k = 50),
                   family = nb(link = "log"), data = dat, method = "REML")
    mu_mat[, g] <- fitted(gam_fit)                    # spot-wise means
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
#' @return List of mean_param, gene_relationship, and cov_gene
#' @keywords internal
prepare_input <- function(expr, Sigma, mu_mat, gene_names) {
  mean_param <- list(CT = mu_mat)
  
  gene_relationship <- as.data.table(
    matrix(1, nrow = nrow(expr), ncol = 1,
           dimnames = list(rownames(expr), "CT"))
  )
  
  cov_gene <- array(NA, dim = c(nrow(expr), nrow(expr), 1),
                    dimnames = list(rownames(expr), rownames(expr), "CT"))
  cov_gene[,,"CT"] <- array(Sigma,
                    dim       = c(nrow(expr), nrow(expr),1),
                    dimnames  = list(rownames(expr), rownames(expr), "CT"))
  
  return(list(mean_param = mean_param, gene_relationship = gene_relationship, cov_gene = cov_gene))
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
#' @return List of estimated parameters
#' @export
get_real_param <- function(expr, coord, gene_names, moran_thr = 0.15, cor_thr = 0.05, block_size = 1024, n_cores = 1) {
  # Step 1: Calculate negative binomial parameters
  nb_params <- calculate_nb_params(expr)
  mu_g <- nb_params$mu
  theta_g <- nb_params$theta
  
  # Step 2: Generate z-scores
  Z <- generate_z_matrix(expr, mu_g, theta_g, block = block_size, n_cores = n_cores)
  
  # Step 3: Calculate the correlation matrix
  Sigma <- calculate_correlation_matrix(Z, thr = cor_thr)
  
  # Step 4: Create spatial gene detection
  coords_mat <- as.matrix(coord[, c("row", "col")])
  Wnb <- knearneigh(coords_mat, k = 6)
  W <- nb2listw(knn2nb(Wnb))
  spatial_genes <- detect_spatial_genes(expr, W, gene_names, moran_thr)
  
  # Step 5: Prepare mu_mat
  lib_size <- colSums(expr)
  lib_factor <- lib_size / mean(lib_size)
  mu_mat <- prepare_mu_mat(expr, mu_g, lib_factor)
  
  # Step 6: Refit NB-GAM for spatial genes
  gam_results <- refit_nb_gam_for_spatial_genes(expr, coord, lib_factor, spatial_genes, mu_mat, theta_g)
  mu_mat <- gam_results$mu_mat
  theta_g <- gam_results$theta_g
  
  # Step 7: Prepare YZ input objects
  real_param <- prepare_input(expr, Sigma, mu_mat, gene_names)
  
  return(real_param)
}