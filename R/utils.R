#' Calculate distances from points to centers
#'
#' Computes Euclidean distances directly, without forming the full pairwise
#' distance matrix, so memory grows linearly with the number of points.
#' Coordinates are read from `x_coord`/`y_coord` columns when present and from
#' the first two columns otherwise.
#'
#' @param X Coordinates of the points.
#' @param Y Coordinates of the centers.
#' @return Distance vector when `Y` has one row; otherwise a points-by-centers
#'   distance matrix.
#' @keywords internal
cdist <- function(X, Y) {
  xy <- function(Z) {
    columns <- colnames(Z)
    if (!is.null(columns) && all(c("x_coord", "y_coord") %in% columns)) {
      if (is.matrix(Z)) return(cbind(as.numeric(Z[, "x_coord"]), as.numeric(Z[, "y_coord"])))
      return(cbind(as.numeric(Z[["x_coord"]]), as.numeric(Z[["y_coord"]])))
    }
    Z <- as.matrix(Z)
    cbind(as.numeric(Z[, 1L]), as.numeric(Z[, 2L]))
  }
  points <- xy(X)
  centers <- xy(Y)
  out <- vapply(seq_len(nrow(centers)), function(k) {
    sqrt((points[, 1L] - centers[k, 1L])^2 + (points[, 2L] - centers[k, 2L])^2)
  }, numeric(nrow(points)))
  if (nrow(centers) == 1L) drop(out) else out
}

#' Make a correlation matrix positive definite
#'
#' @param rho Correlation matrix to adjust.
#' @param epsilon Numerical eigenvalue floor used for repair. Eigenvalues below
#'   this value are individually clipped to the floor before unit-diagonal
#'   renormalization.
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
  eigenvalue_lifts <- pmax(epsilon - eig$values, 0)
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
    minimum_eigenvalue_lift = max(eigenvalue_lifts),
    total_eigenvalue_lift = sum(eigenvalue_lifts),
    n_eigenvalues_below_floor = sum(eig$values < epsilon),
    n_eigenvalues_modified = sum(eigenvalue_lifts > 0),
    relative_frobenius_perturbation =
      sqrt(sum((repaired - requested)^2)) /
      max(sqrt(sum(requested^2)), .Machine$double.eps),
    maximum_entrywise_change = max(abs(repaired - requested)),
    epsilon = epsilon,
    repair_method = "individual_eigenvalue_clipping"
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
