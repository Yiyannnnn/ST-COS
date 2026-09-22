#' Generate a single layer stripe
#' 
#' @param coord Coordinate matrix
#' @param x_start Start x position
#' @param x_end End x position
#' @param y_start Start y position 
#' @param y_end End y position
#' @param layer_id Layer identifier
#' @return Data frame with layer assignments
#' @keywords internal
generate_strap_layer <- function(coord, x_start, x_end, y_start, y_end, layer_id) {
  keep <- coord$x_coord > x_start & coord$x_coord <= x_end &
    coord$y_coord > y_start & coord$y_coord <= y_end
  coord_layer <- coord[keep, , drop = FALSE]
  coord_layer$layer <- layer_id
  coord_layer
}

#' Generate layer coordinates
#' 
#' @param coord Base coordinates
#' @param xmin Vector of minimum x values
#' @param xmax Vector of maximum x values
#' @param ymin Vector of minimum y values
#' @param ymax Vector of maximum y values
#' @return Data frame with layer assignments
#' @export
generate_layer <- function(coord, xmin, xmax, ymin, ymax) {
  if (length(xmin) != length(xmax) | length(ymin) != length(ymax) | length(xmin) != length(ymin)) {
    print("layer boundry length not match.")
  }
  x_min <- xmin - 10e-6
  x_max <- xmax + 10e-6
  y_min <- ymin - 10e-6
  y_max <- ymax + 10e-6
  
  layer_df <- data.frame()
  for (i in 1:length(xmin)) {
    layer_i <- generate_strap_layer(coord,
                                  x_min[i], 
                                  x_max[i],
                                  y_min[i], 
                                  y_max[i],
                                  i)
    layer_df <- rbind(layer_df, layer_i)
  }
  
  layer_df <- layer_df[!duplicated(layer_df[, c("x_coord", "y_coord")]), ]
  layer_df$layer[layer_df$layer == 0] <- 1
  
  layer_coord = merge(coord, layer_df, by = c("x_coord","y_coord"), sort = FALSE)
  return(layer_coord)
}

#' Generate cell proportion layers
#' 
#' @param layer_coord Coordinate matrix with layer assignments
#' @param layer_param Layer parameters matrix
#' @return Matrix of cell type proportions
#' @export
generate_cell_prop <- function(layer_coord, layer_param) {
  layer_factor = as.numeric(as.factor(layer_coord$layer))
  param = list()
  for (i in unique(layer_factor)) {
    for (j in colnames(layer_param)) {
      param[[j]][layer_factor==i] = rnorm(sum(layer_factor==i), layer_param[i,j], 0.0001)
      param[[j]][param[[j]]<0] <- 0
    }
  }
  param <- do.call(cbind, param)
  celltype_proportion <- param/rowSums(param)
  return(celltype_proportion)
}
