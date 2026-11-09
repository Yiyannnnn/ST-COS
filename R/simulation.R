#' Generate true count matrix
#' 
#' @param coord Coordinate matrix
#' @param mean_param List of mean parameters
#' @param gene_relationship Gene relationship matrix
#' @param cov_gene Gene covariance array
#' @param lib_size Library sizes
#' @param theta_g Optional named list of gene-specific negative-binomial size
#'   vectors, one per cell type. When omitted, `lib_size` is retained as the
#'   legacy size parameter.
#' @param eigen_floor Minimum eigenvalue used to repair latent correlation
#'   matrices before Cholesky factorization.
#' @param num_Core Number of cores for parallel processing
#' @return List of count matrices per cell type
#' @export
generate_true_count <- function(coord, mean_param, 
                               gene_relationship = NULL,
                               cov_gene = NULL,
                               lib_size = NULL,
                               theta_g = NULL,
                               eigen_floor = 1e-6,
                               num_Core = NULL) {
  if (length(mean_param) == 0L || is.null(names(mean_param)) ||
      any(names(mean_param) == "")) {
    stop("mean_param must be a named list with one matrix per cell type.")
  }
  
  if(is.null(lib_size)){
    lib_size <- lapply(
      mean_param,
      function(m) pmax(1, abs(round(stats::rnorm(nrow(m), 1e3, 50))))
    )
  }
  
  get_item <- function(x, cell_type, index) {
    if (is.null(x)) return(NULL)
    if (is.list(x) && !is.array(x)) {
      value <- x[[cell_type]]
      if (is.null(value) && length(x) >= index) value <- x[[index]]
      return(value)
    }
    x
  }

  relationship_for <- function(cell_type, index, n_genes) {
    if (is.null(gene_relationship)) return(integer())
    if (is.data.frame(gene_relationship) ||
        data.table::is.data.table(gene_relationship)) {
      rel <- if (cell_type %in% colnames(gene_relationship)) {
        gene_relationship[[cell_type]]
      } else {
        gene_relationship[[index]]
      }
    } else if (is.matrix(gene_relationship)) {
      rel <- if (!is.null(colnames(gene_relationship)) &&
                 cell_type %in% colnames(gene_relationship)) {
        gene_relationship[, cell_type]
      } else {
        gene_relationship[, index]
      }
    } else {
      rel <- gene_relationship
    }
    if (length(rel) != n_genes) {
      stop("gene_relationship must have one entry per gene.")
    }
    which(!is.na(rel))
  }

  covariance_for <- function(cell_type, index) {
    if (is.null(cov_gene)) return(NULL)
    if (length(dim(cov_gene)) == 2L) return(as.matrix(cov_gene))
    third_names <- dimnames(cov_gene)[[3]]
    slice <- if (!is.null(third_names) && cell_type %in% third_names) {
      cell_type
    } else {
      index
    }
    as.matrix(cov_gene[, , slice])
  }

  count <- vector("list", length(mean_param))
  names(count) <- names(mean_param)
  for (index in seq_along(mean_param)) {
    cell_type <- names(mean_param)[index]
    mu <- as.matrix(mean_param[[index]])
    n_spots <- nrow(mu)
    n_genes <- ncol(mu)
    if (n_spots != nrow(coord)) {
      stop("Each mean_param matrix must have one row per coordinate.")
    }
    gene_names <- colnames(mu)
    if (is.null(gene_names)) gene_names <- paste0("gene_", seq_len(n_genes))
    spot_names <- rownames(mu)
    if (is.null(spot_names)) spot_names <- paste0("spot_", seq_len(n_spots))

    theta_ct <- get_item(theta_g, cell_type, index)
    lib_ct <- get_item(lib_size, cell_type, index)
    if (is.null(lib_ct)) lib_ct <- rep(1e3, n_spots)
    if (length(lib_ct) == 1L) lib_ct <- rep(lib_ct, n_spots)
    if (length(lib_ct) != n_spots) {
      stop("lib_size must have length one or one value per spot.")
    }
    if (!is.null(theta_ct) && !length(theta_ct) %in% c(1L, n_genes)) {
      stop("Each theta_g entry must have length one or one value per gene.")
    }

    correlated_genes <- relationship_for(cell_type, index, n_genes)
    uniforms <- NULL
    if (length(correlated_genes) > 1L) {
      covariance <- covariance_for(cell_type, index)
      if (is.null(covariance)) {
        stop("cov_gene is required when correlated genes are specified.")
      }
      rho <- stats::cov2cor(covariance)
      uniforms <- randcop(rho, nSpots = n_spots, epsilon = eigen_floor)
    }

    count_ct <- matrix(
      0, nrow = n_spots, ncol = n_genes,
      dimnames = list(spot_names, gene_names)
    )
    for (j in seq_len(n_spots)) {
      size <- if (is.null(theta_ct)) lib_ct[j] else theta_ct
      count_ct[j, ] <- stats::rnbinom(n_genes, size = size, mu = mu[j, ])
      if (length(correlated_genes) > 1L) {
        correlated_size <- if (length(size) == 1L) {
          size
        } else {
          size[correlated_genes]
        }
        count_ct[j, correlated_genes] <- stats::qnbinom(
          uniforms[correlated_genes, j],
          size = correlated_size,
          mu = mu[j, correlated_genes]
        )
      }
    }
    count[[cell_type]] <- count_ct
  }
  count
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
  proportions <- as.matrix(celltype_proportion)
  n_spots <- nrow(proportions)
  if (length(cell_number) == 1L) cell_number <- rep(cell_number, n_spots)
  if (length(cell_number) != n_spots) {
    stop("cell_number must have length one or one value per spot.")
  }
  if (any(!is.finite(proportions)) || any(proportions < 0)) {
    stop("celltype_proportion must contain finite non-negative values.")
  }
  if (any(abs(rowSums(proportions) - 1) > 1e-8)) {
    stop("Each row of celltype_proportion must sum to one.")
  }

  gene_names <- colnames(as.matrix(true_count[[colnames(proportions)[1L]]]))
  spot_count <- matrix(
    0, nrow = n_spots, ncol = nGenes,
    dimnames = list(rownames(proportions), gene_names)
  )
  for (cell_type in colnames(proportions)) {
    count_ct <- as.matrix(true_count[[cell_type]])
    if (length(dim(count_ct)) != 2L ||
        any(dim(count_ct) != c(n_spots, nGenes))) {
      stop("Each true_count matrix must have spots in rows and genes in columns.")
    }
    abundance <- cell_number * proportions[, cell_type]
    spot_count <- spot_count + sweep(count_ct, 1L, abundance, `*`)
  }
  spot_count
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
  for (cell_type in names(true_count)) {
    count_ct <- as.matrix(true_count[[cell_type]])
    thinned <- matrix(
      stats::rbinom(length(count_ct), size = round(count_ct), prob = 1 - p0),
      nrow = nrow(count_ct), ncol = ncol(count_ct),
      dimnames = dimnames(count_ct)
    )
    depth_limit <- pmax(
      0L, as.integer(round(seq_depth * proportions[, cell_type]))
    )
    for (j in seq_len(nrow(thinned))) {
      if (sum(thinned[j, ]) > depth_limit[j]) {
        thinned[j, ] <- if (depth_limit[j] == 0L) {
          0
        } else {
          as.vector(stats::rmultinom(
            1L, size = depth_limit[j], prob = thinned[j, ]
          ))
        }
      }
    }
    count_obs[[cell_type]] <- thinned
  }
  count_obs
}
