#' Generate random latent covariance blocks for genes
#'
#' For each cell type, the genes with non-missing `gene_relationship` entries
#' form one block. The block is the matrix exponential, [expm::expm()], of
#' \eqn{A + A^\top} with \eqn{A_{ij} \sim N(0, 1)}, which is positive
#' definite. Correlations can be positive or negative. The diagonal of the
#' remaining genes is \eqn{|N(0, 1)|}.
#'
#' @param gene_relationship Data frame with one row per gene and one column
#'   per cell type.
#' @return Named list of gene-by-gene covariance matrices, one per cell type,
#'   accepted directly by [generate_true_count()]. To impose specific latent
#'   correlation values, supply correlation matrices to
#'   [generate_true_count()] instead.
#' @export
generate_covariance <- function(gene_relationship) {
  cov_matrix = list()
  gene_name = rownames(gene_relationship)
  for(ct in colnames(gene_relationship)) {
    cov0 = matrix(rep(0,nrow(gene_relationship)^2),nrow=nrow(gene_relationship))
    diag(cov0) = abs(rnorm(nrow(gene_relationship)))
    correlated_gene = which(!is.na(gene_relationship[[ct]]))
    random_matrix = matrix(rnorm(length(correlated_gene)^2),nrow=length(correlated_gene))
    random_matrix = random_matrix+t(random_matrix)
    random_matrix = expm(random_matrix)
    cov0[correlated_gene,correlated_gene]=random_matrix
    cov_matrix[[ct]] = cov0
  }
  return(cov_matrix)
}

#' Generate group-based latent covariance blocks
#'
#' For each cell type and each group label in `gene_relationship`, the genes
#' carrying that label form one block: the matrix exponential, [expm::expm()],
#' of \eqn{A + A^\top} with \eqn{A_{ij} = |Z_{ij}|} and
#' \eqn{Z_{ij} \sim N(0, 1)}. Blocks are positive definite with positive
#' entries; genes in different groups are uncorrelated. Because the matrix
#' exponential is dominated by its leading eigenvector, the implied latent
#' correlations are close to one (median about 0.997 for five-gene blocks),
#' and blocks of ten or more genes are numerically singular before the
#' eigenvalue repair in [generate_true_count()]. To impose specific latent
#' correlation values, supply correlation matrices to [generate_true_count()]
#' instead.
#'
#' @param gene_relationship Data frame with one row per gene and one column
#'   per cell type; entries are group labels, with `NA` for genes outside any
#'   block.
#' @return Gene-by-gene-by-cell-type covariance array.
#' @export
generate_covariance_group <- function(gene_relationship) {
  cov_matrix = array(rep(0,ncol(gene_relationship)*nrow(gene_relationship)*nrow(gene_relationship)),
                     dim = c(nrow(gene_relationship),nrow(gene_relationship),ncol(gene_relationship)))
  dimnames(cov_matrix) <- list(rownames(gene_relationship),rownames(gene_relationship),colnames(gene_relationship))
  gene_name = rownames(gene_relationship)
  for(ct in colnames(gene_relationship)) {
    cov0 = matrix(rep(0,nrow(gene_relationship)^2),nrow=nrow(gene_relationship))
    diag(cov0) = abs(rnorm(nrow(gene_relationship)))
    for(gp in unique(unlist(gene_relationship))) {
      if (is.na(gp)) {
        next
      } else {
        correlated_gene = which(gene_relationship[[ct]] == gp)
        random_matrix = matrix(abs(rnorm(length(correlated_gene)^2)),nrow=length(correlated_gene))
        random_matrix = random_matrix+t(random_matrix)
        random_matrix = expm(random_matrix)
        if (length(correlated_gene)) {
          cov0[correlated_gene,correlated_gene]=random_matrix
        }
      }
    }
    cov_matrix[,,ct] = cov0
  }
  return(cov_matrix)
}
