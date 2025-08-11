library(ggplot2)
library(data.table)
library(cowplot)
library(ggsci)
library(egg)
library(dplyr)
library(expm)
library(scuttle)
library(scater)
library(energy)
library(gamlss)
library(gamlss.dist)
library(copula)
library(parallel)
library(Matrix)      # for sparse handling
library(mgcv)        # GAM NB fits
library(spdep)       # Moran's I
library(data.table)



#scale_color_npg()


# utils
cdist <- function(X, Y) {
  as.matrix(dist(rbind(X, Y)))[1:nrow(X), (nrow(X)+1):nrow(rbind(X,Y))]
}

#generate_cov_matrix() <- function(n){
  
#}

# generate coordinates for each spot 
# pattern in c("grid", "hex", "random")
generate_coordinates = function(x_len,y_len,pattern,spot_distance=NULL,n_spot=NULL,seed=123){
  set.seed(seed)
  if(pattern == "grid"){
    x = seq(0,x_len,spot_distance)
    y = seq(0,y_len,spot_distance)
    return(data.table(x_coord = rep(x,times=length(y)),
                      y_coord = rep(y,each=length(x))))
  }else if(pattern == "hex"){
    x = seq(0,x_len,spot_distance)
    y = seq(0,y_len,spot_distance*sqrt(3)/2)
    if(length(y)%%2==0){
      x_coord = rep(c(x,x+spot_distance/2),times=length(y)/2)
    }else{
      x_coord = c(rep(c(x,x+spot_distance/2),times=length(y)%/%2),x)
    }
    return(data.table(x_coord = x_coord,
                      y_coord = rep(y,each=length(x))) %>% filter(x_coord<=x_len))
    
  }else if(pattern == "random"){
    x <- runif(n_spot, min = 0, max = x_len)
    y <- runif(n_spot, min = 0, max = y_len)
    return(data.table(x_coord = rep(x,times=length(y)),
                      y_coord = rep(y,each=length(x)))[sample(seq(1,n_spot^2),n_spot),])
  }else{
    return('Not implemented')
  }
}

#generate param for different gene pattern
#pattern in c("hot_spot", "streak", "layer")
generate_pattern = function(coord,pattern,
                            center_coord=NULL,center_param=NULL,background_param=NULL,hot_spot_size=NULL,
                            streak_x=NULL,streak_x_param=NULL,streak_size=NULL,
                            layer=NULL,layer_param=NULL,seed=123){
  set.seed(seed)
  # hot spot
  if(pattern == "hot_spot"){
    #center_coord = data.table(x_coord = center_coord[1],y_coord = center_coord[2])
    param = rep(background_param,nrow(coord))
    for( i in 1:nrow(center_coord)){
      distance_to_center = cdist(coord,center_coord[1,])
      change_param  = (center_param[1]-background_param)*(1-distance_to_center/hot_spot_size[1])
      change_param[distance_to_center>hot_spot_size[1]] = 0
      param = param + change_param
      center_coord = center_coord[-1,]
      center_param = center_param[-1]
      hot_spot_size = hot_spot_size[-1]
    }
    return(param)
  }else if(pattern == "streak"){
    streak_size=streak_size/2
    param = rep(background_param,nrow(coord))
    for( i in 1:length(streak_x_param)){
      distance_to_center = abs(unlist(coord[,1])- streak_x[1])
      change_param  = (streak_x_param[i]-background_param)*(1-distance_to_center/streak_size[1])
      change_param[distance_to_center>streak_size[1]] = 0
      param = param + change_param
      streak_x = streak_x[-1]
      streak_x_param = streak_x_param[-1]
      streak_size = streak_size[-1]
    }
    return(param)
  }else if(pattern == "layer"){
    layer_factor = as.numeric(as.factor(layer))
    param = rep(-1,nrow(coord))
    for(i in unique(layer_factor)){
      param[layer_factor==i]=layer_param[i]
    }
    return(param)
  }else{
    return('Not implemented')
  }
}




#generate gene covariance matrix
#diag element: positive rnorm
#off-diag element: relationship(T/F) |> rnorm
generate_covariance = function(gene_relationship){
  cov_matrix = list()
  gene_name = rownames(gene_relationship)
  for(ct in colnames(gene_relationship)){
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

#based on gene group (multiple groups)
generate_covariance_group = function(gene_relationship){
  cov_matrix = array(rep(0,ncol(gene_relationship)*nrow(gene_relationship)*nrow(gene_relationship)),
                     dim = c(nrow(gene_relationship),nrow(gene_relationship),ncol(gene_relationship)))
  dimnames(cov_matrix) <- list(rownames(gene_relationship),rownames(gene_relationship),colnames(gene_relationship))
  gene_name = rownames(gene_relationship)
  for(ct in colnames(gene_relationship)){
    cov0 = matrix(rep(0,nrow(gene_relationship)^2),nrow=nrow(gene_relationship))
    diag(cov0) = abs(rnorm(nrow(gene_relationship)))
    for(gp in unique(unlist(gene_relationship))){
      if (is.na(gp)){
        non_correlated_gene = which(unlist(gene_relationship |> select(ct)) == gp)
        random_matrix = matrix(rep(0,length(non_correlated_gene)^2),nrow=length(non_correlated_gene))
      }
      else{
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


generate_strap_layer <- function(coord, x_start, x_end, y_start, y_end, layer_id) {
  coord_layer <- coord |> filter(x_coord >x_start-0.01 & x_coord <= x_end & 
                             y_coord > y_start & y_coord <= y_end) |>
    mutate(layer = layer_id)
  
  return(coord_layer)
}


#generate mean(mu) matrix 
generate_gene_pattern = function(coord,gene_pattern,
                                 center_coord=NULL,center_param=NULL,background_param=NULL,hot_spot_size=NULL,
                                 streak_x=NULL,streak_x_param=NULL,streak_size=NULL,
                                 layer=NULL,layer_param=NULL,
                                 all_lower_mean=10,all_upper_mean=150,seed=1){
  set.seed(seed)
  gene_pattern = gene_pattern %>% as.data.frame(check.names=F)
  mean_matrix_all = list()
  gene_name = rownames(gene_pattern)
  for(ct in colnames(gene_pattern)){
    
    mean_matrix = c()
    for(i in 1:nGenes){
      if(gene_pattern[i,ct]=='hotspot'){
        pos <- which(gene_pattern[,ct]=="hotspot")
        if(is.null(center_coord)){
          center_coord = data.table(x_coord = mean(coord$x_coord),y_coord = mean(coord$y_coord))
        }
        if(is.null(center_param)){
          center_param = runif(1,50,80)
        }
        if(is.null(background_param)){
          background_param = center_param + sample(c(-1, 1), 1, replace = TRUE)*runif(1,20,50)
        }
        if(is.null(hot_spot_size)){
          hot_spot_size = 0.45*min(max(coord$x_coord)-min(coord$x_coord),max(coord$y_coord)-min(coord$y_coord))
        }
        param = generate_pattern(coord,'hot_spot',
                                 center_coord=center_coord[[ct]][which(pos == i),],center_param=center_param,
                                 background_param=background_param,hot_spot_size=hot_spot_size)
        
      }else if(gene_pattern[i,ct]=="streak"){
        pos <- which(gene_pattern[,ct]=="streak")
        
        if(!any(colnames(streak_x) == ct)){
          streak_x[,ct] = rep(mean(coord$x_coord),length(pos))
        }
        if(!any(colnames(streak_x_param) == ct)){
          streak_x_param = runif(1,50,80)
        }
        if(!any(colnames(background_param) == ct)){
          background_param = streak_x_param + sample(c(-1, 1), 1, replace = TRUE)*runif(1,20,50)
        }
        if(!any(colnames(streak_size) == ct)){
          streak_size = 0.5*(max(coord$x_coord)-min(coord$x_coord))
        }
        param = generate_pattern(coord,'streak',
                                 streak_x=streak_x[which(pos == i),ct],streak_x_param=streak_x_param,
                                 background_param=background_param,streak_size=streak_size)
        
      }else if(gene_pattern[i,ct]=='layer'){
        
        #if(is.null(layer)){
      
          layer_num = sample(1:(max(coord$x_coord)/10),1)
          layer_df <- data.frame()
          for (i in 0:layer_num-1){
            layer_i <- generate_strap_layer(coord,
                                            i*max(coord$x_coord)/layer_num, 
                                            (i+1)*max(coord$x_coord)/layer_num,
                                            0-0.01, 100+0.01,
                                            i+1)
            layer_df <- rbind(layer_df,layer_i)
          }
          
          layer_df <- layer_df[!duplicated(layer_df[,c(1,2)]), ] |> 
                      mutate(layer = ifelse(layer == 0,1, layer))
          layer = merge(coord,layer_df, by = c("x_coord","y_coord"),sort = FALSE)$layer
          #}
        #if(is.null(layer_param)){
          set.seed(NULL)
          layer_param =  runif(length(unique(layer)),0,150)
          #}
    
        param = generate_pattern(coord,'layer',
                                 layer = layer,layer_param = layer_param)
        
      }else if(gene_pattern[i,ct]=='no_pattern'){
         param = rep(runif(1,all_lower_mean,all_upper_mean),nrow(coord)) 
      }
      mean_matrix = cbind(mean_matrix,param)
    }
    colnames(mean_matrix) = rownames(gene_pattern)
    mean_matrix_all[[ct]] = mean_matrix
  }
  return(mean_matrix_all)
}



#' Sample from a gaussian copula 
#' @param cov_matrix the correlation matrix in the gaussian copula
#' @param nSpots the number of samples to draw
#' @return matrix of nCell rows, where each row is a sample from the gaussian copula
randcop <-function(cov_matrix, nSpots){
  Col = chol(cov_matrix)
  nGenes = nrow(cov_matrix)
  copular = matrix(rnorm(nGenes*nSpots), ncol = nSpots)
  copular = t(Col) %*% copular
  copular = pnorm(copular)
  return(copular)
}

makespd<-function(rho){
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

#generate true gene expression count for each cell type
#coord,lib_size = NULL,mean_param,gene_relationship,cov_gene,num_Core = NULL(later for parallel)
generate_true_count = function(coord,mean_param, 
                               gene_relationship = NULL,cov_gene = NULL,lib_size = NULL,num_Core = NULL){
  
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





generate_spot_count = function(true_count, nGenes, cell_number, celltype_proportion){
  
  spot_count <- matrix(0, nrow = nrow(celltype_proportion), ncol = nGenes)
  for (i in colnames(celltype_proportion)){
    
    count_attr_cells <- apply(as.matrix(true_count[[i]]), 1,
                              function(gene) gene*cell_number*celltype_proportion[,i])
    spot_count <- spot_count+as.matrix(count_attr_cells)
  }
  
  return(spot_count)
}




generate_strap_layer <- function(coord, x_start, x_end, y_start, y_end, layer_id) {
  coord <- coord |> filter(x_coord > x_start & x_coord <= x_end & 
                             y_coord > y_start & y_coord <= y_end) |>
    mutate(layer = layer_id)
  
  return(coord)
}


generate_layer <- function(coord, xmin, xmax, ymin, ymax){
  if (length(xmin) != length(xmax) | length(ymin) != length(ymax) | length(xmin) != length(ymin)){
    print("layer boundry length not match.")
  }
  x_min <- c(xmin[1]-10e-6, xmin[2]-10e-6, xmin[3]-10e-6, xmin[4]-10e-6)
  x_max <- c(xmax[1]+10e-6, xmax[2]+10e-6, xmax[3]+10e-6, xmax[4]+10e-6)
  y_min <- c(ymin[1]-10e-6,ymin[2]-10e-6,ymin[3]-10e-6,ymin[4]-10e-6)
  y_max <- c(ymax[1]+10e-6,ymax[2]+10e-6,ymax[3]+10e-6,ymax[4]+10e-6)
  layer_df <- data.frame()
  for (i in 1:length(xmin)){
    layer_i <- generate_strap_layer(coord,
                                    x_min[i], 
                                    x_max[i],
                                    y_min[i], 
                                    y_max[i],
                                    i)
    layer_df <- rbind(layer_df,layer_i)
  }
  
  layer_df <- layer_df[!duplicated(layer_df[,c(1,2)]), ] |> 
    mutate(layer = ifelse(layer == 0,1, layer))
  
  layer_coord = merge(coord,layer_df, by = c("x_coord","y_coord"),sort = FALSE)
  return(layer_coord)
}

generate_cell_prop <- function(layer_coord,layer_param){
  layer_factor = as.numeric(as.factor(layer_coord$layer))
  param = list()
  for (i in unique(layer_factor)){
    for (j in colnames(layer_param)){
      param[[j]][layer_factor==i] = rnorm(sum(layer_factor==i),layer_param[i,j],0.0001)
      param[[j]][param[[j]]<0] <- 0
    }
  }
  param <- do.call(cbind,param)
  celltype_proportion <- param/rowSums(param)
  return(celltype_proportion)
}







#Batch effect: amplifying & library size (Symsim: https://www.nature.com/articles/s41467-019-10500-w)


generate_batch_count <- function(true_count, p0, celltype_proportion, seq_depth){
  
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




## ---- Function to fit NB moments --------------------------------------------------
fit_nb_fast <- function(y) {
  mu <- mean(y); s2 <- var(y)
  th <- if (s2 > mu) mu^2 / (s2 - mu) else 1e6
  return(c(mu = mu, theta = th))
}

# Function to calculate the negative binomial parameters for the whole expression table
calculate_nb_params <- function(expr) {
  nb_par <- t(apply(expr, 1, fit_nb_fast))
  mu_g <- nb_par[, "mu"]
  theta_g <- nb_par[, "theta"]
  return(list(mu = mu_g, theta = theta_g))
}

## ---- Function to calculate z-scores -------------------------------------------------
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

# Function to generate the Z matrix (gene expression z-scores)
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

## ---- Function to calculate correlation matrix ---------------------------------
calculate_correlation_matrix <- function(Z, thr = 0.05) {
  Sigma <- cor(Z, use = "pairwise.complete.obs")  # handle constant-variance genes
  Sigma[is.na(Sigma)] <- 0                          # undefined → no corr
  diag(Sigma) <- 1                                  # proper 1’s on the diagonal
  Sigma[abs(Sigma) < thr] <- 0                      # optional sparsity
  return(makespd(Sigma))                         # repair PD if needed
}

#' Calculate Moran's I flag
#' 
#' @param y Vector of values
#' @param listw Spatial weights list
#' @param thr Threshold for significance
#' @return Logical indicating if spatially autocorrelated
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
#' @param W Spatial weights matrix
#' @param gene_names Vector of gene names
#' @param moran_thr Threshold for Moran's I
#' @return Vector of spatial gene names
#' @keywords internal
detect_spatial_genes <- function(expr, W, gene_names, moran_thr = 0.15) {
  spatial_genes <- gene_names[sapply(gene_names, function(g)
    moran_flag(expr[g, ], W, thr = moran_thr))]
  cat(length(spatial_genes), "spatial genes detected\n")
  return(spatial_genes)
}

#' Refit negative binomial GAM for spatial genes
#' 
#' @param expr Expression matrix
#' @param coord Coordinates data frame
#' @param lib_factor Library size factors
#' @param spatial_genes Vector of spatial gene names
#' @param mu_mat Matrix of means
#' @param theta_g Vector of dispersion parameters
#' @return List with updated mu_mat and theta_g
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

## ---- Function to prepare mu_mat -------------------------------------------------
prepare_mu_mat <- function(expr, mu_g, lib_factor) {
  mu_mat <- outer(lib_factor, mu_g)                   # spots × genes
  dimnames(mu_mat) <- list(colnames(expr), rownames(expr))
  return(mu_mat)
}

## ---- Function to refit NB-GAM for spatial genes ---------------------------------
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

## ---- Function to prepare YZ input objects -------------------------------------
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

## ---- Main function -------------------------------------------------------------
# Input:
# expr: The gene expression matrix (rows = genes, columns = spots)
# coord: Coordinates matrix with columns row and col for spatial coordinates
# gene_names: Vector of gene names
# moran_thr: Threshold for Moran's I (default is 0.15)
# cor_thr: Correlation threshold for the correlation matrix (default is 0.05)
# block_size: Block size for processing (default is 1024)
# n_cores: Number of cores for parallel processing (default is 1)

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