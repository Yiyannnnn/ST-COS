#' Calculate pairwise distances between two sets of points
#' 
#' @param X Matrix of coordinates for first set of points
#' @param Y Matrix of coordinates for second set of points
#' @return Distance matrix between X and Y points
#' @keywords internal
cdist <- function(X, Y) {
  as.matrix(dist(rbind(X, Y)))[1:nrow(X), (nrow(X)+1):nrow(rbind(X,Y))]
}

#' Make a matrix positive definite
#' 
#' @param rho Correlation matrix to adjust
#' @return Positive definite correlation matrix
#' @keywords internal
makespd <- function(rho) {
  er = eigen(rho)
  if(min(er$values)<0){
    oldsum = sum(er$values)
    er$values = er$values - min(er$values) + 1e-6
    newsum = sum(er$values)
    er$values = er$values/newsum*oldsum 
    rhocop = er$vectors %*% diag(er$values) %*% t(er$vectors)
    rhocop = rhocop%*%diag(1/diag(rhocop))
    return(rhocop)
  }
  else{
    return(rho)
  }
}

#' Sample from a gaussian copula
#' 
#' @param cov_matrix Correlation matrix for gaussian copula
#' @param nSpots Number of samples to draw
#' @return Matrix of samples from gaussian copula
#' @keywords internal
randcop <- function(cov_matrix, nSpots) {
  Col = chol(cov_matrix)
  nGenes = nrow(cov_matrix)
  copular = matrix(rnorm(nGenes*nSpots), ncol = nSpots)
  copular = t(Col) %*% copular
  copular = pnorm(copular)
  return(copular)
}