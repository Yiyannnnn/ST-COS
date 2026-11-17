#' Assign spots to rectangular layers
#'
#' Spot \eqn{j} belongs to layer \eqn{i} when
#' \eqn{x_{\min,i} < x_j \le x_{\max,i}} and
#' \eqn{y_{\min,i} < y_j \le y_{\max,i}} (with a tolerance of
#' \eqn{10^{-5}}). When rectangles overlap, the first listed layer wins.
#' Every spot must fall inside at least one rectangle.
#'
#' @param coord Coordinates with `x_coord` and `y_coord`.
#' @param xmin,xmax,ymin,ymax Rectangle bounds, one value per layer.
#' @return `coord` with an integer `layer` column, in the original row order.
#' @export
generate_layer <- function(coord, xmin, xmax, ymin, ymax) {
  n_layers <- length(xmin)
  if (n_layers == 0L || any(c(length(xmax), length(ymin), length(ymax)) != n_layers)) {
    stop("xmin, xmax, ymin, and ymax must be non-empty and have the same length.")
  }
  tolerance <- 1e-5
  x <- coord$x_coord
  y <- coord$y_coord
  layer <- rep(NA_integer_, length(x))
  for (i in seq_len(n_layers)) {
    inside <- is.na(layer) &
      x > xmin[i] - tolerance & x <= xmax[i] + tolerance &
      y > ymin[i] - tolerance & y <= ymax[i] + tolerance
    layer[inside] <- i
  }
  if (anyNA(layer)) {
    stop(sum(is.na(layer)), " of ", length(layer),
         " spots fall outside every layer rectangle.")
  }
  if (data.table::is.data.table(coord)) {
    out <- data.table::copy(coord)
    data.table::set(out, j = "layer", value = layer)
  } else {
    out <- coord
    out$layer <- layer
  }
  out
}

#' Generate cell proportion layers
#'
#' Raw abundance parameters \eqn{\alpha_{i,c} \sim N(\gamma_{i,c}, \tau^2)}
#' are drawn once per layer \eqn{i} and cell type \eqn{c}, exponentiated, and
#' normalized across cell types, so every spot \eqn{j} in layer \eqn{i(j)}
#' receives \eqn{p_{jc} = \exp(\alpha_{i(j),c}) / \sum_{c'}
#' \exp(\alpha_{i(j),c'})}.
#'
#' @param layer_coord Coordinate matrix with layer assignments
#' @param layer_param Finite layer-by-cell-type mean logits
#'   \eqn{\gamma_{i,c}}.
#' @param noise_sd Standard deviation \eqn{\tau} of the per-layer Gaussian
#'   perturbation of each logit.
#' @return Matrix of cell type proportions
#' @export
generate_cell_prop <- function(layer_coord, layer_param, noise_sd = 1e-4) {
  if (!is.numeric(noise_sd) || length(noise_sd) != 1L ||
      !is.finite(noise_sd) || noise_sd < 0) {
    stop("noise_sd must be one nonnegative finite number.")
  }
  layer_param <- as.matrix(layer_param)
  if (!is.numeric(layer_param) || any(!is.finite(layer_param))) {
    stop("layer_param must contain finite mean logits.")
  }
  if (is.null(colnames(layer_param))) {
    stop("layer_param must have cell-type column names.")
  }

  layer_factor <- as.numeric(as.factor(layer_coord$layer))
  if (nrow(layer_param) < max(layer_factor)) {
    stop("layer_param must contain one row for every layer.")
  }
  logits <- matrix(
    NA_real_, nrow = length(layer_factor), ncol = ncol(layer_param),
    dimnames = list(NULL, colnames(layer_param))
  )
  for (i in unique(layer_factor)) {
    for (j in seq_len(ncol(layer_param))) {
      logits[layer_factor == i, j] <- stats::rnorm(1L, layer_param[i, j], noise_sd)
    }
  }
  logits <- logits - apply(logits, 1L, max)
  weights <- exp(logits)
  weights / rowSums(weights)
}
