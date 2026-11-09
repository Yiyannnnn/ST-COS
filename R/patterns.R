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
#' @param layer_assignments Optional vector of layer labels shared by cell
#'   types, or a named list with one vector per cell type.
#' @return Named list of spot-by-gene mean matrices, one per cell type.
#' @export
generate_gene_pattern <- function(coord, gene_pattern,
                                center_coord=NULL, center_param=NULL,
                                background_param=NULL, hot_spot_size=NULL,
                                streak_x=NULL, streak_x_param=NULL, streak_size=NULL,
                                layer=NULL, layer_param=NULL,
                                all_lower_mean=10, all_upper_mean=150, seed=1,
                                layer_assignments = NULL) {
  set.seed(seed)
  gene_pattern <- as.data.frame(gene_pattern, check.names = FALSE)
  n_genes <- nrow(gene_pattern)
  gene_names <- rownames(gene_pattern)
  if (is.null(gene_names)) gene_names <- paste0("gene_", seq_len(n_genes))
  mean_matrix_all <- vector("list", ncol(gene_pattern))
  names(mean_matrix_all) <- colnames(gene_pattern)

  for (cell_type in colnames(gene_pattern)) {
    mean_matrix <- matrix(
      NA_real_, nrow = nrow(coord), ncol = n_genes,
      dimnames = list(rownames(coord), gene_names)
    )
    for (g in seq_len(n_genes)) {
      pattern <- gene_pattern[g, cell_type]
      if (pattern %in% c("hotspot", "hot_spot")) {
        positions <- which(gene_pattern[, cell_type] %in% c("hotspot", "hot_spot"))
        center <- if (is.list(center_coord) && cell_type %in% names(center_coord)) {
          center_coord[[cell_type]][which(positions == g), , drop = FALSE]
        } else {
          data.frame(
            x_coord = mean(coord$x_coord),
            y_coord = mean(coord$y_coord)
          )
        }
        center_value <- if (is.null(center_param)) runif(1, 50, 80) else center_param
        background_value <- if (is.null(background_param)) {
          center_value + sample(c(-1, 1), 1) * runif(1, 20, 50)
        } else {
          background_param
        }
        radius <- if (is.null(hot_spot_size)) {
          0.45 * min(diff(range(coord$x_coord)), diff(range(coord$y_coord)))
        } else {
          hot_spot_size
        }
        param <- generate_pattern(
          coord, "hot_spot", center_coord = center,
          center_param = center_value, background_param = background_value,
          hot_spot_size = radius
        )
      } else if (pattern == "streak") {
        x_value <- if (is.null(streak_x)) mean(coord$x_coord) else streak_x
        streak_value <- if (is.null(streak_x_param)) runif(1, 50, 80) else streak_x_param
        background_value <- if (is.null(background_param)) {
          streak_value + sample(c(-1, 1), 1) * runif(1, 20, 50)
        } else {
          background_param
        }
        width <- if (is.null(streak_size)) diff(range(coord$x_coord)) / 2 else streak_size
        param <- generate_pattern(
          coord, "streak", streak_x = x_value,
          streak_x_param = streak_value,
          background_param = background_value, streak_size = width
        )
      } else if (pattern == "layer") {
        layer_values <- layer
        if (!is.null(layer_assignments)) {
          if (is.list(layer_assignments)) {
            layer_values <- layer_assignments[[cell_type]]
            if (is.null(layer_values)) layer_values <- layer_assignments[["default"]]
          } else {
            layer_values <- layer_assignments
          }
        }
        if (is.null(layer_values) || length(layer_values) != nrow(coord)) {
          stop("Layer-pattern genes require one layer assignment per spot.")
        }
        layer_levels <- as.integer(factor(layer_values))
        values <- layer_param
        if (is.null(values)) values <- runif(length(unique(layer_levels)), 0, 150)
        param <- generate_pattern(
          coord, "layer", layer = layer_levels, layer_param = values
        )
      } else if (pattern == "no_pattern") {
        param <- rep(runif(1, all_lower_mean, all_upper_mean), nrow(coord))
      } else {
        stop(sprintf("Unknown pattern '%s' for gene %s.", pattern, gene_names[g]))
      }
      mean_matrix[, g] <- param
    }
    mean_matrix_all[[cell_type]] <- mean_matrix
  }

  mean_matrix_all
}
