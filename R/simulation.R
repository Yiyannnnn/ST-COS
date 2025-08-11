#' Generate true count matrix
#' 
#' @param coord Coordinate matrix
#' @param mean_param List of mean parameters
#' @param gene_relationship Gene relationship matrix
#' @param cov_gene Gene covariance array
#' @param lib_size Library sizes
#' @param num_Core Number of cores for parallel processing
#' @return List of count matrices per cell type
#' @export
generate_true_count <- function(coord, mean_param, 
                               gene_relationship = NULL,
                               cov_gene = NULL,
                               lib_size = NULL,
                               num_Core = NULL) {
  
  if(is.null(lib_size)){
    lib_size <- lapply(mean_param,function(m) abs(round(rnorm(nrow(m),1e3,50))))
  }
  
  count = list()
  if(sum(gene_relationship,na.rm = TRUE)>0){
    mulpo<-function(m,c,lib,corrgene){
      
      count = rnbinom(length(m),size  = lib,mu = m)
      names(count) <- names(c)
      count[corrgene] = qnbinom(c[corrgene],size = lib, mu = m[corrgene])
      return(count)
    }
    for (i in names(mean_param)){
      mu = mean_param[[i]]
      nGenes = ncol(mu)
      nSpots = nrow(coord)
      colnames(mu) <- paste0("gene", 1:nGenes)
      lib_ct <- lib_size[[i]]
      rho <- makespd(cov2cor(cov_gene[,,i]))
      copular = randcop(rho,nSpots = nrow(coord))
      
      corrgene <- which(!is.na(gene_relationship[,..i]))
      
      count_ct <- matrix(NA,nrow = nSpots,ncol = nGenes)
      colnames(count_ct) <- paste0("gene_",1:nGenes)
      rownames(count_ct) <- paste0("spot_",1:nSpots)
      for (j in 1:length(lib_ct)){
        count_ct[j,] <- mulpo(m = mu[j,],c = copular[,j],lib = lib_ct[j],corrgene = corrgene)
        count[[i]] <- count_ct
      }
    }
  }
  else{
    for (i in names(mean_param)){
      nGenes = ncol(mean_param[[i]])
      nSpots = nrow(mean_param[[i]])
      count[[i]] <- apply(mean_param[[i]],1,function(m) rnbinom(length(m),size  = lib_size[[i]],mu = m))
      colnames(count[[i]]) <- paste0("spot_",1:nSpots)
      rownames(count[[i]]) <- paste0("gene_",1:nGenes)
    }
  }
  return(count)
}

#' Generate spot-level count matrix
#' 
#' @param true_count List of true count matrices
#' @param nGenes Number of genes
#' @param cell_number Number of cells per spot
#' @param celltype_proportion Cell type proportions
#' @return Matrix of spot-level counts
#' @export
generate_spot_count <- function(true_count, nGenes, cell_number, celltype_proportion) {
  spot_count <- matrix(0, nrow = nrow(celltype_proportion), ncol = nGenes)
  for (i in colnames(celltype_proportion)){
    count_attr_cells <- apply(as.matrix(true_count[[i]]), 2,
                             function(gene) gene*cell_number*celltype_proportion[,i])
    spot_count <- spot_count + as.matrix(count_attr_cells)
  }
  return(spot_count)
}

#' Generate batch effects in counts
#' 
#' @param true_count True count matrices
#' @param p0 Dropout probability
#' @param celltype_proportion Cell type proportions
#' @param seq_depth Sequencing depth
#' @return List of count matrices with batch effects
#' @export
generate_batch_count <- function(true_count, p0, celltype_proportion, seq_depth) {
  #capturing
  count <- lapply(true_count,function(ct) apply(ct,c(1,2),function(n) rbinom(1,n,1-p0)))
  
  #library size
  count_failed <- lapply(count,function(ct) 
    sapply(seq_len(nrow(ct)), function(n) ifelse(sum(ct[n,])>seq_depth*celltype_proportion[n], TRUE, FALSE)))
  
  count_obs = count
  for (ct in names(true_count)){
    count_obs[[ct]][which(count_failed[[ct]]),] = round(proportions(count[[ct]][which(count_failed[[ct]]),],1)*seq_depth*celltype_proportion[which(count_failed[[ct]]),ct])
  }
  
  return(count_obs)
}
