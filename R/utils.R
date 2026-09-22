#' Calculate pairwise distances between two sets of points
#' 
#' @param X Matrix of coordinates for first set of points
#' @param Y Matrix of coordinates for second set of points
#' @return Distance matrix between X and Y points
#' @keywords internal
cdist <- function(X, Y) {
  as.matrix(dist(rbind(X, Y)))[1:nrow(X), (nrow(X)+1):nrow(rbind(X,Y))]
}

#' Make a correlation matrix positive definite
#'
#' @param rho Correlation matrix to adjust.
#' @param epsilon Minimum eigenvalue after repair.
#' @param return_diagnostics Return the repaired matrix and perturbation
#'   diagnostics in a list.
#' @return A positive-definite correlation matrix, or a list containing the
#'   matrix and diagnostics when `return_diagnostics = TRUE`.
#' @keywords internal
makespd <- function(rho, epsilon = 1e-6, return_diagnostics = FALSE) {
  if (!is.matrix(rho) || nrow(rho) != ncol(rho)) {
    stop("rho must be a square matrix.")
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L ||
      !is.finite(epsilon) || epsilon <= 0) {
    stop("epsilon must be one finite positive number.")
  }

  requested <- (rho + t(rho)) / 2
  diag(requested) <- 1
  eig <- eigen(requested, symmetric = TRUE)
  repaired_values <- pmax(eig$values, epsilon)
  repaired <- eig$vectors %*%
    diag(repaired_values, nrow = length(repaired_values)) %*%
    t(eig$vectors)
  repaired <- stats::cov2cor(repaired)
  repaired <- (repaired + t(repaired)) / 2
  diag(repaired) <- 1
  dimnames(repaired) <- dimnames(rho)

  diagnostics <- list(
    min_eigen_before = min(eig$values),
    min_eigen_after = min(eigen(
      repaired, symmetric = TRUE, only.values = TRUE
    )$values),
    minimum_eigenvalue_lift = max(epsilon - min(eig$values), 0),
    n_eigenvalues_modified = sum(eig$values < epsilon),
    relative_frobenius_perturbation =
      sqrt(sum((repaired - requested)^2)) /
      max(sqrt(sum(requested^2)), .Machine$double.eps),
    maximum_entrywise_change = max(abs(repaired - requested)),
    epsilon = epsilon
  )

  if (return_diagnostics) {
    return(list(rho = repaired, diagnostics = diagnostics))
  }
  attr(repaired, "repair_diagnostics") <- diagnostics
  repaired
}

#' Sample from a gaussian copula
#' 
#' @param cov_matrix Correlation matrix for gaussian copula
#' @param nSpots Number of samples to draw
#' @return Matrix of samples from gaussian copula
#' @keywords internal
randcop <- function(cov_matrix, nSpots, epsilon = 1e-6) {
  rho <- makespd(cov_matrix, epsilon = epsilon)
  upper_chol <- chol(rho)
  nGenes <- nrow(rho)
  latent <- t(upper_chol) %*%
    matrix(stats::rnorm(nGenes * nSpots), nrow = nGenes, ncol = nSpots)
  stats::pnorm(latent)
}
