#' Generate covariance matrix for genes
#' 
#' @param gene_relationship Matrix indicating gene relationships
#' @return List of covariance matrices per cell type
#' @export
generate_covariance <- function(gene_relationship) {
  cov_matrix = list()
  gene_name = rownames(gene_relationship)
  for(ct in colnames(gene_relationship)) {
    cov0 = matrix(rep(0,nrow(gene_relationship)^2),nrow=nrow(gene_relationship))
    diag(cov0) = abs(rnorm(nrow(gene_relationship)))
    correlated_gene = which(unlist(gene_relationship%>%select(ct)))
    random_matrix = matrix(rnorm(length(correlated_gene)^2),nrow=length(correlated_gene))
    random_matrix = random_matrix+t(random_matrix)
    random_matrix = expm(random_matrix)
    cov0[correlated_gene,correlated_gene]=random_matrix
    cov_matrix[[ct]] = cov0
  }
  return(cov_matrix)
}

#' Generate group-based covariance matrix
#' 
#' @param gene_relationship Matrix indicating gene group relationships
#' @return Array of covariance matrices
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
        non_correlated_gene = which(unlist(gene_relationship |> select(ct)) == gp)
        random_matrix = matrix(rep(0,length(non_correlated_gene)^2),nrow=length(non_correlated_gene))
      }
      else {
        correlated_gene = which(unlist(gene_relationship |> select(ct)) == gp)
        random_matrix = matrix(abs(rnorm(length(correlated_gene)^2)),nrow=length(correlated_gene))
        random_matrix = random_matrix+t(random_matrix)
        random_matrix = expm(random_matrix)
        cov0[correlated_gene,correlated_gene]=random_matrix
      }
    }
    cov_matrix[,,ct] = cov0
  }
  return(cov_matrix)
}