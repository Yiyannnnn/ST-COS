#' Generate cell-type-specific spot-level expression profiles
#'
#' Draws one negative-binomial expression profile per cell type and spot.
#' Genes flagged in `gene_relationship` share a Gaussian copula: their latent
#' correlation matrix is repaired with [makespd()] (individual eigenvalue
#' clipping), factorized with [chol()], and the resulting uniforms are mapped
#' through the negative-binomial inverse CDF. Other genes are drawn
#' independently.
#'
#' @param coord Coordinate matrix or data frame with one row per spot.
#' @param mean_param Named list of non-negative spot-by-gene mean matrices,
#'   one per cell type.
#' @param gene_relationship Data frame or matrix with one row per gene and one
#'   column per cell type; non-missing entries mark the genes that enter the
#'   copula for that cell type.
#' @param cov_gene Latent covariance or correlation input for correlated genes:
#'   one gene-by-gene matrix shared by all cell types, a gene-by-gene-by-cell
#'   type array, or a named list of per-cell-type matrices such as the output
#'   of [generate_covariance()]. It is converted to a correlation matrix with
#'   [stats::cov2cor()].
#' @param lib_size Legacy spot-specific negative-binomial size \eqn{s_j}, used
#'   for every gene when `theta_g` is omitted. It is a dispersion parameter,
#'   not a library size. Supply one vector shared by all cell types or a named
#'   list with one vector per cell type. When both `lib_size` and `theta_g` are
#'   omitted, \eqn{s_j = \max\{1, |\mathrm{round}(s_j^*)|\}} with
#'   \eqn{s_j^* \sim N(1000, 50^2)}.
#' @param theta_g Optional named list of gene-specific negative-binomial size
#'   vectors, one per cell type. When supplied, it is used instead of
#'   `lib_size`.
#' @param library_size Optional observed library sizes \eqn{L_j}, one value
#'   per spot (or a named list with one vector per cell type). Each spot's
#'   means are multiplied by \eqn{f_j = L_j / \bar L}, as for the
#'   library-normalized means returned by [get_real_param()].
#' @param eigen_floor Minimum eigenvalue used to repair latent correlation
#'   matrices before Cholesky factorization.
#' @param num_Core Unused; retained for backward compatibility.
#' @return Named list of spot-by-gene expression matrices, one per cell type.
#'   Each matrix represents a latent cell-type-specific average expression
#'   profile at each spot, not counts from independently sampled individual
#'   cells.
#' @export
generate_true_count <- function(coord, mean_param, 
                               gene_relationship = NULL,
                               cov_gene = NULL,
                               lib_size = NULL,
                               theta_g = NULL,
                               library_size = NULL,
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
    if (is.list(cov_gene) && !is.data.frame(cov_gene) && !is.array(cov_gene)) {
      value <- cov_gene[[cell_type]]
      if (is.null(value) && is.null(names(cov_gene)) && length(cov_gene) >= index) {
        value <- cov_gene[[index]]
      }
      return(if (is.null(value)) NULL else as.matrix(value))
    }
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
    if (any(!is.finite(mu)) || any(mu < 0)) {
      stop("mean_param for cell type '", cell_type,
           "' contains negative or non-finite means.")
    }
    library_ct <- get_item(library_size, cell_type, index)
    if (!is.null(library_ct)) {
      if (length(library_ct) != n_spots || any(!is.finite(library_ct)) ||
          any(library_ct < 0)) {
        stop("library_size must contain one non-negative value per spot.")
      }
      if (mean(library_ct) > 0) mu <- mu * (library_ct / mean(library_ct))
    } else if (isTRUE(attr(mean_param[[index]], "library_normalized"))) {
      warning("mean_param for cell type '", cell_type, "' is library-normalized; ",
              "supply library_size to restore library-size variation.")
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
      if (any(dim(covariance) != n_genes)) {
        stop("cov_gene for cell type '", cell_type,
             "' must be a ", n_genes, " x ", n_genes, " matrix.")
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

#' Construct an abundance-weighted spot-level expression matrix
#'
#' For spot \eqn{j}, gene \eqn{g}, and cell type \eqn{c}, constructs
#' \deqn{X_{jg} = \sum_c N_j p_{jc} Y_{jg}^{(c)},}
#' which propagates variation among cell-type-specific spot profiles and does
#' not simulate or sum independent individual cells within a spot. When
#' `seq_depth` is supplied, any spot whose total exceeds the depth threshold
#' \eqn{\Delta} is resampled to exactly \eqn{\Delta} counts with a
#' multinomial draw whose probabilities are proportional to that spot's
#' profile; spots at or below \eqn{\Delta} are unchanged.
#'
#' @param true_count Named list of spot-by-gene cell-type-specific expression
#'   matrices. Each row is a latent average profile for that cell type at one
#'   spot, not an individual-cell observation.
#' @param nGenes Number of genes
#' @param cell_number Total cell-equivalent abundance per spot, supplied as one
#'   value or one value per spot.
#' @param celltype_proportion Spot-by-cell-type proportion matrix.
#' @param seq_depth Optional depth threshold \eqn{\Delta}, one value or one
#'   value per spot.
#' @return Abundance-weighted spot-by-gene expression matrix.
#' @export
generate_spot_count <- function(true_count, nGenes, cell_number, celltype_proportion,
                                seq_depth = NULL) {
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

  if (!is.null(seq_depth)) {
    if (length(seq_depth) == 1L) seq_depth <- rep(seq_depth, n_spots)
    if (length(seq_depth) != n_spots || any(!is.finite(seq_depth)) ||
        any(seq_depth < 0)) {
      stop("seq_depth must be one non-negative value or one value per spot.")
    }
    depth <- round(seq_depth)
    for (j in which(rowSums(spot_count) > depth)) {
      spot_count[j, ] <- if (depth[j] == 0) {
        0
      } else {
        as.vector(stats::rmultinom(1L, size = depth[j], prob = spot_count[j, ]))
      }
    }
  }
  spot_count
}

#' Apply binomial thinning to cell-type-specific profiles
#'
#' Replaces each count \eqn{Y} of every cell-type-specific profile returned by
#' [generate_true_count()] with \eqn{Y' \sim \mathrm{Binomial}(Y, 1 - p_0)},
#' modelling capture efficiency. The sequencing-depth threshold is applied to
#' the mixed spot profile by [generate_spot_count()].
#'
#' @param true_count Named list of spot-by-gene cell-type profiles.
#' @param p0 Per-molecule binomial thinning probability. This operation can
#'   increase the observed zero proportion but is not a zero-inflated mixture
#'   model.
#' @param celltype_proportion Unused; retained for backward compatibility.
#' @param seq_depth Not supported here; supply the depth threshold to
#'   [generate_spot_count()].
#' @return Named list of thinned cell-type profiles.
#' @export
generate_batch_count <- function(true_count, p0, celltype_proportion = NULL,
                                 seq_depth = NULL) {
  if (!is.null(seq_depth)) {
    stop("The depth threshold applies to the mixed spot profile; ",
         "pass seq_depth to generate_spot_count().")
  }
  if (!is.numeric(p0) || length(p0) != 1L || !is.finite(p0) ||
      p0 < 0 || p0 > 1) {
    stop("p0 must be one finite probability between zero and one.")
  }
  lapply(true_count, function(count_ct) {
    count_ct <- as.matrix(count_ct)
    matrix(
      stats::rbinom(length(count_ct), size = round(count_ct), prob = 1 - p0),
      nrow = nrow(count_ct), ncol = ncol(count_ct),
      dimnames = dimnames(count_ct)
    )
  })
}
