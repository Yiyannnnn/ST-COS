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
#' @param layer Layer assignments, one per spot.
#' @param layer_param Mean for each layer, in the order of the sorted layer
#'   labels.
#' @param seed Retained for backward compatibility. `generate_pattern()` is
#'   deterministic once its pattern parameters are supplied and does not alter
#'   the random-number generator state.
#' @return Vector of expression values
#' @export
generate_pattern <- function(coord, pattern, 
                           center_coord=NULL, center_param=NULL,
                           background_param=NULL, hot_spot_size=NULL,
                           streak_x=NULL, streak_x_param=NULL, streak_size=NULL,
                           layer=NULL, layer_param=NULL, seed=123) {
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
    parameter_lengths <- c(
      streak_x = length(streak_x),
      streak_x_param = length(streak_x_param),
      streak_size = length(streak_size)
    )
    if (any(parameter_lengths == 0L)) {
      stop("streak_x, streak_x_param, and streak_size must be supplied.")
    }
    n_streaks <- max(parameter_lengths)
    if (any(!parameter_lengths %in% c(1L, n_streaks))) {
      stop("Streak parameters must be scalars or have one value per streak.")
    }

    streak_x <- rep(streak_x, length.out = n_streaks)
    streak_x_param <- rep(streak_x_param, length.out = n_streaks)
    streak_size <- rep(streak_size, length.out = n_streaks)
    if (any(!is.finite(streak_size)) || any(streak_size <= 0)) {
      stop("streak_size must contain positive finite values.")
    }

    half_width <- streak_size / 2
    coord_x <- as.numeric(if ("x_coord" %in% names(coord)) coord$x_coord else coord[[1L]])
    param <- rep(background_param, nrow(coord))
    for (i in seq_len(n_streaks)) {
      distance_to_center <- abs(coord_x - streak_x[i])
      change_param <- (streak_x_param[i] - background_param) *
        (1 - distance_to_center / half_width[i])
      change_param[distance_to_center > half_width[i]] <- 0
      param <- param + change_param
    }
    return(param)
  }else if(pattern == "layer"){
    layer_factor = as.numeric(as.factor(layer))
    if (length(layer_param) < max(layer_factor)) {
      stop("layer_param must contain one mean per layer (", max(layer_factor),
           " layers, ", length(layer_param), " values).")
    }
    param = layer_param[layer_factor]
    return(param)
  }else{
    return('Not implemented')
  }
}

#' Generate gene expression patterns
#'
#' Builds the spatial mean surface \eqn{\mu_{g,c}(x, y)} for every gene and
#' cell type. Each gene-cell-type entry of `gene_pattern` is one of
#' `"hotspot"` (or `"hot_spot"`), `"streak"`, `"layer"`, or `"no_pattern"`.
#'
#' Pattern parameters can be given as one scalar for all genes and cell types,
#' a vector with one value per gene (or per gene carrying that pattern), a data
#' frame or matrix with one column per cell type, or a named list with one
#' entry per cell type. Omitted levels are sampled separately for each
#' gene-cell-type combination.
#'
#' @param coord Coordinates data.table with `x_coord` and `y_coord`.
#' @param gene_pattern Data frame with one row per gene and one column per
#'   cell type.
#' @param center_coord Named list with one data frame of hotspot centers per
#'   cell type (one row per hotspot gene); defaults to the domain center.
#' @param center_param Hotspot center level.
#' @param background_param Background level for hotspots and streaks.
#' @param hot_spot_size Hotspot radius \eqn{r}.
#' @param streak_x Streak position \eqn{s} on the x axis.
#' @param streak_x_param Streak level.
#' @param streak_size Streak width \eqn{w} (half-width \eqn{w/2}).
#' @param layer Layer labels shared by all cell types (see also
#'   `layer_assignments`).
#' @param layer_param Layer means \eqn{\ell_{i,g,c}}: a vector with one mean
#'   per layer shared by the layer-pattern genes of a cell type, or, inside a
#'   named list by cell type, a matrix with one row per gene (or per
#'   layer-pattern gene) and one column per layer. Sampled from
#'   \eqn{U(0, 150)} per gene when omitted.
#' @param all_lower_mean,all_upper_mean Bounds \eqn{[\alpha_{g,c},
#'   \beta_{g,c}]} of the uniform distribution used for `"no_pattern"` genes;
#'   scalars or per-gene and per-cell-type values as described above.
#' @param seed Optional seed passed to [set.seed()] before any draw. `NULL`
#'   (the default) uses the current random-number state. Supplying a seed
#'   resets the global random-number generator.
#' @param layer_assignments Optional vector of layer labels shared by cell
#'   types, or a named list with one vector per cell type.
#' @return Named list of spot-by-gene mean matrices, one per cell type.
#' @export
generate_gene_pattern <- function(coord, gene_pattern,
                                center_coord=NULL, center_param=NULL,
                                background_param=NULL, hot_spot_size=NULL,
                                streak_x=NULL, streak_x_param=NULL, streak_size=NULL,
                                layer=NULL, layer_param=NULL,
                                all_lower_mean=10, all_upper_mean=150, seed=NULL,
                                layer_assignments = NULL) {
  if (!is.null(seed)) set.seed(seed)
  gene_pattern <- as.data.frame(gene_pattern, check.names = FALSE)
  n_genes <- nrow(gene_pattern)
  gene_names <- rownames(gene_pattern)
  if (is.null(gene_names)) gene_names <- paste0("gene_", seq_len(n_genes))

  cell_type_value <- function(value, cell_type) {
    if (is.null(value)) return(NULL)
    if (is.list(value) && !is.data.frame(value)) {
      selected <- value[[cell_type]]
      if (!is.null(selected)) return(selected)
    }
    if ((is.data.frame(value) || is.matrix(value)) &&
        cell_type %in% colnames(value)) {
      return(value[, cell_type])
    }
    value
  }

  gene_pattern_value <- function(value, cell_type, gene_index,
                                 pattern_positions) {
    selected <- cell_type_value(value, cell_type)
    if (is.null(selected)) return(NULL)
    selected <- unlist(selected, use.names = FALSE)
    if (length(selected) == 1L) return(selected)
    if (length(selected) == n_genes) return(selected[gene_index])
    local_index <- match(gene_index, pattern_positions)
    if (!is.na(local_index) && local_index <= length(selected)) {
      return(selected[local_index])
    }
    stop("A pattern parameter could not be matched to its gene and cell type.")
  }

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
        center_value <- gene_pattern_value(
          center_param, cell_type, g, positions
        )
        if (is.null(center_value)) center_value <- runif(1, 50, 80)
        background_value <- if (is.null(background_param)) {
          center_value + sample(c(-1, 1), 1) * runif(1, 20, 50)
        } else {
          gene_pattern_value(background_param, cell_type, g, positions)
        }
        radius <- gene_pattern_value(hot_spot_size, cell_type, g, positions)
        if (is.null(radius)) {
          radius <- 0.45 * min(
            diff(range(coord$x_coord)), diff(range(coord$y_coord))
          )
        }
        param <- generate_pattern(
          coord, "hot_spot", center_coord = center,
          center_param = center_value, background_param = background_value,
          hot_spot_size = radius
        )
      } else if (pattern == "streak") {
        positions <- which(gene_pattern[, cell_type] == "streak")
        x_value <- gene_pattern_value(streak_x, cell_type, g, positions)
        if (is.null(x_value)) x_value <- mean(coord$x_coord)
        streak_value <- gene_pattern_value(
          streak_x_param, cell_type, g, positions
        )
        if (is.null(streak_value)) streak_value <- runif(1, 50, 80)
        background_value <- if (is.null(background_param)) {
          streak_value + sample(c(-1, 1), 1) * runif(1, 20, 50)
        } else {
          gene_pattern_value(background_param, cell_type, g, positions)
        }
        width <- gene_pattern_value(streak_size, cell_type, g, positions)
        if (is.null(width)) width <- diff(range(coord$x_coord)) / 2
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
        values <- cell_type_value(layer_param, cell_type)
        if (is.matrix(values) || is.data.frame(values)) {
          values <- as.matrix(values)
          layer_positions <- which(gene_pattern[, cell_type] == "layer")
          row <- if (nrow(values) == n_genes) g else match(g, layer_positions)
          if (is.na(row) || row > nrow(values)) {
            stop("Layer means for gene ", gene_names[g], " and cell type ",
                 cell_type, " could not be matched.")
          }
          values <- values[row, ]
        }
        if (is.null(values)) values <- runif(length(unique(layer_levels)), 0, 150)
        param <- generate_pattern(
          coord, "layer", layer = layer_levels, layer_param = values
        )
      } else if (pattern == "no_pattern") {
        positions <- which(gene_pattern[, cell_type] == "no_pattern")
        lower <- gene_pattern_value(all_lower_mean, cell_type, g, positions)
        upper <- gene_pattern_value(all_upper_mean, cell_type, g, positions)
        if (!is.finite(lower) || !is.finite(upper) || lower < 0 || lower > upper) {
          stop("Uniform mean bounds for gene ", gene_names[g], " and cell type ",
               cell_type, " must satisfy 0 <= lower <= upper.")
        }
        param <- rep(runif(1, lower, upper), nrow(coord))
      } else {
        stop(sprintf("Unknown pattern '%s' for gene %s.", pattern, gene_names[g]))
      }
      mean_matrix[, g] <- param
    }
    mean_matrix_all[[cell_type]] <- mean_matrix
  }

  mean_matrix_all
}
