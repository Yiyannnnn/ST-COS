
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
  n_genes <- nrow(gene_pattern)

  select_value <- function(value, cell_type, position, default) {
    if (is.null(value)) return(default())
    if (is.list(value) && !is.data.frame(value)) value <- value[[cell_type]]
    if (is.data.frame(value) && cell_type %in% colnames(value)) value <- value[[cell_type]]
    if (length(value) > 1L) value <- value[min(position, length(value))]
    as.numeric(value)
  }

  for(ct in colnames(gene_pattern)){
    mean_matrix <- matrix(NA_real_, nrow = nrow(coord), ncol = n_genes)
    for(i in seq_len(n_genes)){
      if(gene_pattern[i,ct]=='hotspot'){
        pos <- which(gene_pattern[,ct]=="hotspot")
        position <- which(pos == i)
        center <- if (is.list(center_coord) && ct %in% names(center_coord)) {
          center_coord[[ct]][position, , drop = FALSE]
        } else {
          data.frame(x_coord = mean(coord$x_coord), y_coord = mean(coord$y_coord))
        }
        center_value <- select_value(center_param, ct, position, function() runif(1, 50, 80))
        background_value <- select_value(
          background_param, ct, position,
          function() center_value + sample(c(-1, 1), 1) * runif(1, 20, 50)
        )
        radius <- select_value(
          hot_spot_size, ct, position,
          function() 0.45 * min(diff(range(coord$x_coord)), diff(range(coord$y_coord)))
        )
        param = generate_pattern(coord,'hot_spot',
                                 center_coord=center,center_param=center_value,
                                 background_param=background_value,hot_spot_size=radius)
        
      }else if(gene_pattern[i,ct]=="streak"){
        pos <- which(gene_pattern[,ct]=="streak")
        position <- which(pos == i)
        streak_position <- select_value(streak_x, ct, position, function() mean(coord$x_coord))
        streak_value <- select_value(streak_x_param, ct, position, function() runif(1, 50, 80))
        background_value <- select_value(
          background_param, ct, position,
          function() streak_value + sample(c(-1, 1), 1) * runif(1, 20, 50)
        )
        width <- select_value(
          streak_size, ct, position,
          function() 0.5 * diff(range(coord$x_coord))
        )
        param = generate_pattern(coord,'streak',
                                 streak_x=streak_position,streak_x_param=streak_value,
                                 background_param=background_value,streak_size=width)
        
      }else if(gene_pattern[i,ct]=='layer'){
        layer_values <- layer
        if (is.list(layer_values)) layer_values <- layer_values[[ct]]
        if(is.null(layer_values)){
          layer_num = sample(1:(max(coord$x_coord)/10),1)
          layer_df <- data.frame()
          for (layer_i_index in 0:(layer_num-1)){
            layer_i <- generate_strap_layer(coord,
                                            layer_i_index*max(coord$x_coord)/layer_num, 
                                            (layer_i_index+1)*max(coord$x_coord)/layer_num,
                                            0-0.01, 100+0.01,
                                            layer_i_index+1)
            layer_df <- rbind(layer_df,layer_i)
          }
          layer_df <- layer_df[!duplicated(layer_df[,c(1,2)]), ] |> 
                      mutate(layer = ifelse(layer == 0,1, layer))
          layer_values = merge(coord,layer_df, by = c("x_coord","y_coord"),sort = FALSE)$layer
        }
        layer_levels <- as.integer(as.factor(layer_values))
        values <- layer_param
        if (is.list(values)) values <- values[[ct]]
        if (is.data.frame(values) && ct %in% colnames(values)) values <- values[[ct]]
        if (is.null(values)) values <- runif(length(unique(layer_levels)),0,150)
        param = generate_pattern(coord,'layer',
                                 layer = layer_levels,layer_param = values)
        
      }else if(gene_pattern[i,ct]=='no_pattern'){
         param = rep(runif(1,all_lower_mean,all_upper_mean),nrow(coord)) 
      }else{
        stop(sprintf("Unknown pattern '%s' for gene %s.", gene_pattern[i,ct], rownames(gene_pattern)[i]))
      }
      mean_matrix[,i] = param
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
randcop <-function(cov_matrix, nSpots, epsilon = 1e-6){
  rho <- makespd(cov_matrix, epsilon = epsilon)
  Col = chol(rho)
  nGenes = nrow(cov_matrix)
  copular = matrix(rnorm(nGenes*nSpots), ncol = nSpots)
  copular = t(Col) %*% copular
  copular = pnorm(copular)
  return(copular)
}

makespd<-function(rho, epsilon = 1e-6){
  requested <- (rho + t(rho)) / 2
  diag(requested) <- 1
  er <- eigen(requested, symmetric = TRUE)
  uniform_shift <- max(0, epsilon - min(er$values))
  repaired_values <- er$values + uniform_shift
  rhocop <- er$vectors %*% diag(repaired_values) %*% t(er$vectors)
  rhocop <- cov2cor(rhocop)
  rhocop <- (rhocop + t(rhocop)) / 2
  diag(rhocop) <- 1
  return(rhocop)
}

#generate true gene expression count for each cell type
#coord,lib_size = NULL,mean_param,gene_relationship,cov_gene,num_Core = NULL(later for parallel)
# theta_g: named list (one element per cell type) of gene-specific NB dispersion vectors.
#          If NULL (default), falls back to original behaviour (lib_size used as dispersion),
#          preserving backward compatibility with synthetic scenarios (S1-S8).
generate_true_count = function(coord,mean_param,
                               gene_relationship = NULL,cov_gene = NULL,
                               lib_size = NULL, theta_g = NULL, num_Core = NULL){

  if(is.null(lib_size)){
    lib_size <- lapply(mean_param,function(m) abs(round(rnorm(nrow(m),1e3,50))))
  }

  count = list()
  if(sum(gene_relationship,na.rm = TRUE)>0){
    mulpo <- function(m, c, lib, corrgene, theta_vec = NULL){
      nb_size <- if (!is.null(theta_vec)) theta_vec else rep(lib, length(m))
      count = rnbinom(length(m), size = nb_size, mu = m)
      names(count) <- names(c)
      count[corrgene] = qnbinom(c[corrgene], size = nb_size[corrgene], mu = m[corrgene])
      return(count)
    }
    for (i in names(mean_param)){
      mu = mean_param[[i]]
      nGenes = ncol(mu)
      nSpots = nrow(coord)
      colnames(mu) <- paste0("gene", 1:nGenes)
      lib_ct <- lib_size[[i]]
      theta_ct <- if (!is.null(theta_g)) {
        tv <- theta_g[[i]]; if (is.null(tv)) theta_g[[1]] else tv
      } else NULL
      rho <- makespd(cov2cor(cov_gene[,,i]))
      copular = randcop(rho,nSpots = nrow(coord))

      corrgene <- which(!is.na(gene_relationship[,..i]))

      count_ct <- matrix(NA,nrow = nSpots,ncol = nGenes)
      colnames(count_ct) <- paste0("gene_",1:nGenes)
      rownames(count_ct) <- paste0("spot_",1:nSpots)
        for (j in 1:length(lib_ct)){
          count_ct[j,] <- mulpo(m = mu[j,], c = copular[,j],
                                lib = lib_ct[j], corrgene = corrgene,
                                theta_vec = theta_ct)
          count[[i]] <- count_ct
        }
      }
  }
  else{
    for (i in names(mean_param)){
      nGenes = ncol(mean_param[[i]])
      nSpots = nrow(mean_param[[i]])
      theta_ct <- if (!is.null(theta_g)) {
        tv <- theta_g[[i]]; if (is.null(tv)) theta_g[[1]] else tv
      } else NULL
      if (!is.null(theta_ct)) {
        count[[i]] <- apply(mean_param[[i]], 1,
                            function(m) rnbinom(length(m), size = theta_ct, mu = m))
      } else {
        count[[i]] <- apply(mean_param[[i]], 1,
                            function(m) rnbinom(length(m), size = lib_size[[i]], mu = m))
      }
      colnames(count[[i]]) <- paste0("spot_",1:nSpots)
      rownames(count[[i]]) <- paste0("gene_",1:nGenes)
    }
  }
  return(count)
}





generate_spot_count_abundance_weighted = function(true_count, nGenes, cell_number, celltype_proportion){
  proportions <- as.matrix(celltype_proportion)
  n_spots <- nrow(proportions)
  if (length(cell_number) == 1L) cell_number <- rep(cell_number, n_spots)
  spot_count <- matrix(0, nrow = n_spots, ncol = nGenes)
  for (cell_type in colnames(proportions)) {
    count_ct <- as.matrix(true_count[[cell_type]])
    if (nrow(count_ct) != n_spots && ncol(count_ct) == n_spots) {
      count_ct <- t(count_ct)
    }
    if (any(dim(count_ct) != c(n_spots, nGenes))) {
      stop("Each true_count matrix must contain one row per spot and one column per gene.")
    }
    abundance <- cell_number * proportions[, cell_type]
    spot_count <- spot_count + sweep(count_ct, 1L, abundance, `*`)
  }
  colnames(spot_count) <- colnames(count_ct)
  return(spot_count)
}

# Legacy aggregation above scales one cell-type draw by its expected abundance.
# This function instead allocates an integer number of cells to each type and
# sums independently sampled cell-level expression vectors.
generate_spot_count_independent <- function(
    mean_param, cell_number, celltype_proportion,
    gene_relationship = NULL, cov_gene = NULL,
    lib_size = NULL, theta_g = NULL,
    p0 = 0, seq_depth = NULL, seed = NULL,
    return_details = FALSE) {

  if (!is.null(seed)) set.seed(seed)
  if (length(mean_param) == 0L || is.null(names(mean_param))) {
    stop("mean_param must be a named list with one matrix per cell type.")
  }

  prop <- as.matrix(celltype_proportion)
  ct_names <- intersect(names(mean_param), colnames(prop))
  if (length(ct_names) == 0L) {
    stop("No cell types are shared by mean_param and celltype_proportion.")
  }
  prop <- prop[, ct_names, drop = FALSE]
  if (any(!is.finite(prop)) || any(prop < 0)) {
    stop("celltype_proportion must contain finite, non-negative values.")
  }
  prop_sum <- rowSums(prop)
  if (any(prop_sum <= 0)) stop("Every spot must have a positive proportion sum.")
  prop <- prop / prop_sum

  n_spots <- nrow(prop)
  if (length(cell_number) == 1L) cell_number <- rep(cell_number, n_spots)
  if (length(cell_number) != n_spots || any(!is.finite(cell_number)) ||
      any(cell_number < 0) || any(abs(cell_number - round(cell_number)) > 1e-8)) {
    stop("cell_number must supply one non-negative integer total per spot.")
  }
  cell_number <- as.integer(round(cell_number))

  first_mu <- as.matrix(mean_param[[ct_names[1]]])
  if (nrow(first_mu) != n_spots) {
    stop("Each mean_param matrix must have one row per spot.")
  }
  n_genes <- ncol(first_mu)
  gene_names <- colnames(first_mu)
  if (is.null(gene_names)) gene_names <- paste0("gene_", seq_len(n_genes))

  get_ct_item <- function(x, ct, ct_index) {
    if (is.null(x)) return(NULL)
    if (is.list(x) && !is.array(x)) {
      ans <- x[[ct]]
      if (is.null(ans) && length(x) >= ct_index) ans <- x[[ct_index]]
      return(ans)
    }
    x
  }

  get_covariance <- function(ct, ct_index) {
    if (is.null(cov_gene)) return(NULL)
    if (length(dim(cov_gene)) == 2L) return(as.matrix(cov_gene))
    third_names <- dimnames(cov_gene)[[3]]
    slice <- if (!is.null(third_names) && ct %in% third_names) ct else ct_index
    as.matrix(cov_gene[, , slice])
  }

  get_correlated_genes <- function(ct, ct_index, has_covariance) {
    if (!has_covariance) return(integer())
    if (is.null(gene_relationship)) return(seq_len(n_genes))
    rel <- if (is.data.frame(gene_relationship) || data.table::is.data.table(gene_relationship)) {
      if (ct %in% colnames(gene_relationship)) gene_relationship[[ct]] else gene_relationship[[ct_index]]
    } else if (is.matrix(gene_relationship)) {
      if (!is.null(colnames(gene_relationship)) && ct %in% colnames(gene_relationship)) {
        gene_relationship[, ct]
      } else {
        gene_relationship[, ct_index]
      }
    } else {
      gene_relationship
    }
    which(!is.na(rel))
  }

  copula_info <- vector("list", length(ct_names))
  names(copula_info) <- ct_names
  for (ct_index in seq_along(ct_names)) {
    ct <- ct_names[ct_index]
    mu_ct <- as.matrix(mean_param[[ct]])
    if (!identical(dim(mu_ct), c(n_spots, n_genes))) {
      stop(sprintf("mean_param for %s has inconsistent dimensions.", ct))
    }
    sigma_ct <- get_covariance(ct, ct_index)
    corr_idx <- get_correlated_genes(ct, ct_index, !is.null(sigma_ct))
    lower_chol <- NULL
    if (length(corr_idx) > 1L) {
      rho_ct <- makespd(cov2cor(sigma_ct))
      if (min(eigen(rho_ct, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
        rho_ct <- as.matrix(Matrix::nearPD(rho_ct, corr = TRUE)$mat)
      }
      # R's chol() returns U with t(U) %*% U = rho; L = t(U).
      lower_chol <- t(chol(rho_ct))
    }
    copula_info[[ct]] <- list(corr_idx = corr_idx, lower_chol = lower_chol)
  }

  resolve_size <- function(ct, ct_index, spot_index) {
    theta_ct <- get_ct_item(theta_g, ct, ct_index)
    if (!is.null(theta_ct)) {
      if (length(theta_ct) == 1L) return(rep(theta_ct, n_genes))
      if (length(theta_ct) != n_genes) {
        stop(sprintf("theta_g for %s must have length 1 or %d.", ct, n_genes))
      }
      return(as.numeric(theta_ct))
    }

    lib_ct <- get_ct_item(lib_size, ct, ct_index)
    if (is.null(lib_ct)) return(rep(1e3, n_genes))
    if (length(lib_ct) == n_spots) return(rep(lib_ct[spot_index], n_genes))
    if (length(lib_ct) == n_genes) return(as.numeric(lib_ct))
    if (length(lib_ct) == 1L) return(rep(lib_ct, n_genes))
    stop(sprintf("lib_size for %s has an unsupported length.", ct))
  }

  sample_cells <- function(mu, size, n_cells, info) {
    if (n_cells == 0L) return(rep(0, n_genes))
    mu_rep <- rep(as.numeric(mu), times = n_cells)
    size_rep <- rep(pmax(as.numeric(size), 1e-8), times = n_cells)
    counts <- matrix(
      rnbinom(n_genes * n_cells, size = size_rep, mu = mu_rep),
      nrow = n_genes, ncol = n_cells
    )

    corr_idx <- info$corr_idx
    if (length(corr_idx) > 1L) {
      latent <- info$lower_chol %*% matrix(
        rnorm(n_genes * n_cells), nrow = n_genes, ncol = n_cells
      )
      uniforms <- pmin(pmax(pnorm(latent), 1e-12), 1 - 1e-12)
      for (g in corr_idx) {
        counts[g, ] <- qnbinom(uniforms[g, ], size = size[g], mu = mu[g])
      }
    }
    rowSums(counts)
  }

  cell_allocations <- matrix(
    0L, nrow = n_spots, ncol = length(ct_names),
    dimnames = list(rownames(prop), ct_names)
  )
  for (j in seq_len(n_spots)) {
    cell_allocations[j, ] <- as.vector(rmultinom(1L, cell_number[j], prop[j, ]))
  }

  celltype_spot_counts <- lapply(ct_names, function(ct) {
    matrix(0, nrow = n_spots, ncol = n_genes,
           dimnames = list(rownames(prop), gene_names))
  })
  names(celltype_spot_counts) <- ct_names

  for (ct_index in seq_along(ct_names)) {
    ct <- ct_names[ct_index]
    mu_ct <- as.matrix(mean_param[[ct]])
    for (j in seq_len(n_spots)) {
      ct_sum <- sample_cells(
        mu = mu_ct[j, ],
        size = resolve_size(ct, ct_index, j),
        n_cells = cell_allocations[j, ct],
        info = copula_info[[ct]]
      )
      if (p0 > 0) {
        ct_sum <- rbinom(n_genes, size = round(ct_sum), prob = 1 - p0)
      }
      if (!is.null(seq_depth)) {
        depth_limit <- max(0L, as.integer(round(seq_depth * prop[j, ct])))
        if (sum(ct_sum) > depth_limit && depth_limit > 0L) {
          ct_sum <- as.vector(rmultinom(1L, depth_limit, prob = ct_sum))
        } else if (depth_limit == 0L) {
          ct_sum[] <- 0
        }
      }
      celltype_spot_counts[[ct]][j, ] <- ct_sum
    }
  }

  spot_count <- Reduce(`+`, celltype_spot_counts)
  storage.mode(spot_count) <- "integer"
  if (return_details) {
    return(list(
      spot_count = spot_count,
      cell_allocations = cell_allocations,
      celltype_spot_counts = celltype_spot_counts
    ))
  }
  spot_count
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
  proportions <- as.matrix(celltype_proportion)
  if (!all(names(true_count) %in% colnames(proportions))) {
    stop("celltype_proportion must contain every cell type in true_count.")
  }
  if (length(seq_depth) == 1L) seq_depth <- rep(seq_depth, nrow(proportions))
  if (length(seq_depth) != nrow(proportions)) {
    stop("seq_depth must have length one or one value per spot.")
  }

  count_obs <- vector("list", length(true_count))
  names(count_obs) <- names(true_count)
  for (ct in names(true_count)) {
    count_ct <- as.matrix(true_count[[ct]])
    thinned <- matrix(
      rbinom(length(count_ct), size = round(count_ct), prob = 1 - p0),
      nrow = nrow(count_ct), ncol = ncol(count_ct),
      dimnames = dimnames(count_ct)
    )
    depth_limit <- pmax(0L, as.integer(round(seq_depth * proportions[, ct])))
    for (j in seq_len(nrow(thinned))) {
      if (sum(thinned[j, ]) > depth_limit[j]) {
        thinned[j, ] <- if (depth_limit[j] == 0L) {
          0
        } else {
          as.vector(rmultinom(1L, depth_limit[j], thinned[j, ] / sum(thinned[j, ])))
        }
      }
    }
    count_obs[[ct]] <- thinned
  }
  count_obs
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

## ---- Function for Moran's I flag -------------------------------------------------
moran_flag <- function(y, listw, thr = 0.15) {
  S0  <- Szero(listw)                          # sum of all weights
  res <- moran(x  = y,
               listw = listw,
               n    = length(y),
               S0   = S0,
               zero.policy = TRUE)$I
  return(res > thr)
}

# Function to detect spatial genes based on Moran's I
detect_spatial_genes <- function(expr, W, gene_names, moran_thr = 0.15) {
  spatial_genes <- gene_names[sapply(gene_names, function(g)
    moran_flag(expr[g, ], W, thr = moran_thr))]
  cat(length(spatial_genes), "spatial genes detected\n")
  return(spatial_genes)
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
