#' Generate expression pattern
#' 
#' @param coord Coordinates data.table
#' @param pattern Pattern type ("hot_spot", "streak", "layer")
#' @param center_coord Center coordinates for hotspot
#' @param center_param Center parameters
#' @param background_param Background parameters 
#' @param hot_spot_size Size of hotspots
#' @param streak_x X positions for streaks
#' @param streak_x_param Streak parameters
#' @param streak_size Size of streaks
#' @param layer Layer assignments
#' @param layer_param Layer parameters
#' @param seed Random seed
#' @return Vector of expression values
#' @export
generate_pattern <- function(coord, pattern, 
                           center_coord=NULL, center_param=NULL,
                           background_param=NULL, hot_spot_size=NULL,
                           streak_x=NULL, streak_x_param=NULL, streak_size=NULL,
                           layer=NULL, layer_param=NULL, seed=123) {
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

#' Generate gene expression patterns
#' 
#' @param coord Coordinates data.table
#' @param gene_pattern Gene pattern specifications
#' @param center_coord Center coordinates for hotspots
#' @param center_param Center parameters
#' @param background_param Background parameters
#' @param hot_spot_size Hotspot sizes
#' @param streak_x Streak x positions
#' @param streak_x_param Streak parameters
#' @param streak_size Streak sizes
#' @param layer Layer assignments
#' @param layer_param Layer parameters
#' @param all_lower_mean Lower bound for mean
#' @param all_upper_mean Upper bound for mean
#' @param seed Random seed
#' @return List of mean matrices per cell type
#' @export
generate_gene_pattern <- function(coord, gene_pattern,
                                center_coord=NULL, center_param=NULL,
                                background_param=NULL, hot_spot_size=NULL,
                                streak_x=NULL, streak_x_param=NULL, streak_size=NULL,
                                layer=NULL, layer_param=NULL,
                                all_lower_mean=10, all_upper_mean=150, seed=1) {
  set.seed(seed)
  # Initialize mean matrix
  mean_matrix <- matrix(0, nrow = nrow(coord), ncol = length(gene_pattern))
  colnames(mean_matrix) <- paste0("Gene", 1:length(gene_pattern))
  
  for (g in 1:length(gene_pattern)) {
    pattern <- gene_pattern[[g]]
    
    # Hot spot parameters
    if(!is.null(pattern$hot_spot)){
      hot_spot_size = pattern$hot_spot$size
      center_coord = pattern$hot_spot$center
      center_param = pattern$hot_spot$param
      background_param = pattern$hot_spot$background_param
      mean_matrix[,g] <- generate_pattern(coord, "hot_spot", 
                                         center_coord, center_param,
                                         background_param, hot_spot_size)
    }
    
    # Streak parameters
    if(!is.null(pattern$streak)){
      streak_x = pattern$streak$x
      streak_x_param = pattern$streak$param
      streak_size = pattern$streak$size
      mean_matrix[,g] <- generate_pattern(coord, "streak", 
                                         streak_x = streak_x, 
                                         streak_x_param = streak_x_param, 
                                         streak_size = streak_size)
    }
    
    # Layer parameters
    if(!is.null(pattern$layer)){
      layer = pattern$layer$assignment
      layer_param = pattern$layer$param
      mean_matrix[,g] <- generate_pattern(coord, "layer", 
                                         layer = layer, 
                                         layer_param = layer_param)
    }
    
    # Apply global mean adjustment
    global_mean = runif(1, all_lower_mean, all_upper_mean)
    mean_matrix[,g] <- mean_matrix[,g] + global_mean
  }
  
  return(mean_matrix)
}