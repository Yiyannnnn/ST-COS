#!/usr/bin/env Rscript

# Image-driven spatial simulator (minimal deps)
# - Segments a histology image into K domains via k-means on RGB(+XY)
# - Maps a hex grid of spots to image domains
# - Assigns domain-specific cell-type proportions and domain-specific gene means
# - Simulates counts per spot using simple NB sampling

suppressPackageStartupMessages({
  # No hard deps beyond base R and (jpeg or png) for image I/O
  # Optional: data.table (for fwrite) and scales (palette). We fall back to base otherwise.
})

# Robust column-wise scaling: center and divide by sd; if sd==0 or NA, return zeros.
robust_scale <- function(M) {
  if (!is.matrix(M)) M <- as.matrix(M)
  p <- ncol(M)
  N <- nrow(M)
  if (N == 0L || p == 0L) return(M)
  for (j in seq_len(p)) {
    col <- M[, j]
    mu <- mean(col, na.rm = TRUE)
    sdv <- stats::sd(col, na.rm = TRUE)
    if (is.na(sdv) || sdv == 0) {
      M[, j] <- 0
    } else {
      M[, j] <- (col - mu) / sdv
      M[is.na(M[, j]), j] <- 0
    }
  }
  M
}

# Try to read JPEG first, then PNG as fallback
read_image_any <- function(path) {
  img <- tryCatch({
    suppressPackageStartupMessages(library(jpeg))
    jpeg::readJPEG(path)
  }, error = function(e) NULL)
  if (is.null(img)) {
    img <- tryCatch({
      suppressPackageStartupMessages(library(png))
      png::readPNG(path)
    }, error = function(e) NULL)
  }
  if (is.null(img)) stop("Failed to read image. Please install 'jpeg' or 'png' package and check file path.")
  # Coerce to 3-channel RGB if grayscale
  if (length(dim(img)) == 2L) {
    img <- array(rep(img, 3L), dim = c(dim(img), 3L))
  }
  return(img)
}

downscale_image <- function(img, max_side = 512) {
  h <- dim(img)[1]; w <- dim(img)[2]
  s <- max(h, w)
  if (s <= max_side) return(img)
  # simple nearest-neighbor downscale to avoid heavy deps
  ratio <- max_side / s
  new_h <- max(1L, round(h * ratio))
  new_w <- max(1L, round(w * ratio))
  ys <- unique(round(seq(1, h, length.out = new_h)))
  xs <- unique(round(seq(1, w, length.out = new_w)))
  img[ys, xs, , drop = FALSE]
}

toLab_image <- function(img) {
  # Convert 3-channel sRGB [0,1] image to CIE Lab using grDevices::convertColor
  h <- dim(img)[1]; w <- dim(img)[2]
  RGB <- cbind(as.vector(img[,,1]), as.vector(img[,,2]), as.vector(img[,,3]))
  Lab <- tryCatch({
    grDevices::convertColor(RGB, from = "sRGB", to = "Lab")
  }, error = function(e) NULL)
  if (is.null(Lab)) return(NULL)
  L <- matrix(Lab[,1], nrow = h, ncol = w)
  a <- matrix(Lab[,2], nrow = h, ncol = w)
  b <- matrix(Lab[,3], nrow = h, ncol = w)
  list(L = L, a = a, b = b)
}

compute_gradients <- function(img) {
  # Compute gradient magnitude on luminance using central differences
  h <- dim(img)[1]; w <- dim(img)[2]
  lum <- (img[,,1] + img[,,2] + img[,,3]) / 3
  dx <- matrix(0, nrow = h, ncol = w)
  dy <- matrix(0, nrow = h, ncol = w)
  if (w > 2) dx[, 2:(w-1)] <- (lum[, 3:w] - lum[, 1:(w-2)]) / 2
  if (h > 2) dy[2:(h-1), ] <- (lum[3:h, ] - lum[1:(h-2), ]) / 2
  sqrt(dx^2 + dy^2)
}

compute_tissue_mask <- function(img, white_thr = 0.92) {
  # img values assumed in [0,1]
  R <- img[,,1]; G <- img[,,2]; B <- img[,,3]
  lum <- (R + G + B) / 3
  tissue <- lum < white_thr
  tissue
}

prepare_segmentation_features <- function(img,
                                          xy_weight = 0.15,
                                          mask_background = TRUE,
                                          white_thr = 0.92,
                                          use_lab = TRUE) {
  h <- dim(img)[1]; w <- dim(img)[2]
  tissue_mask <- if (mask_background) compute_tissue_mask(img, white_thr) else matrix(TRUE, nrow = h, ncol = w)
  mask_vec <- as.vector(tissue_mask)
  if (isTRUE(use_lab)) {
    lab <- toLab_image(img)
    if (!is.null(lab)) {
      R <- as.vector(lab$L)[mask_vec]
      G <- as.vector(lab$a)[mask_vec]
      B <- as.vector(lab$b)[mask_vec]
    } else {
      R <- as.vector(img[,,1])[mask_vec]
      G <- as.vector(img[,,2])[mask_vec]
      B <- as.vector(img[,,3])[mask_vec]
    }
  } else {
    R <- as.vector(img[,,1])[mask_vec]
    G <- as.vector(img[,,2])[mask_vec]
    B <- as.vector(img[,,3])[mask_vec]
  }
  yy_full <- rep(seq_len(h), times = w)
  xx_full <- rep(seq_len(w), each = h)
  yy <- yy_full[mask_vec]
  xx <- xx_full[mask_vec]
  X <- cbind(R, G, B,
             x = (xx - 1) / max(1, (w - 1)) * xy_weight,
             y = (yy - 1) / max(1, (h - 1)) * xy_weight)
  X <- robust_scale(X)
  index_map <- matrix(NA_integer_, nrow = h, ncol = w)
  index_map[cbind(yy, xx)] <- seq_len(nrow(X))
  list(
    X = X,
    tissue_mask = tissue_mask,
    coords = cbind(xx = xx, yy = yy),
    dims = c(h = h, w = w),
    index_map = index_map
  )
}

mode_filter_labels <- function(seg, iterations = 1L) {
  if (!is.matrix(seg) || iterations <= 0L) return(seg)
  mode_once <- function(A) {
    h <- nrow(A); w <- ncol(A)
    if (h < 3 || w < 3) return(A)
    B <- A
    for (i in 2:(h-1)) {
      for (j in 2:(w-1)) {
        if (!is.na(A[i,j])) {
          nb <- as.integer(c(A[(i-1):(i+1), (j-1):(j+1)]))
          nb <- nb[!is.na(nb)]
          if (length(nb) > 0) {
            tab <- table(nb)
            B[i,j] <- as.integer(names(tab)[which.max(tab)])
          }
        }
      }
    }
    B
  }
  out <- seg
  for (t in seq_len(iterations)) out <- mode_once(out)
  out
}

merge_segments_to_K <- function(seg, centers, targetK, sizes = NULL) {
  if (is.null(targetK) || is.na(targetK)) {
    return(list(seg = seg, centers = centers, sizes = sizes))
  }
  targetK <- max(1L, as.integer(targetK))
  merge_step <- function(seg, centers, sizes) {
    labs <- sort(unique(seg[!is.na(seg)]))
    L <- length(labs)
    if (L <= targetK) return(list(seg = seg, centers = centers, sizes = sizes, done = TRUE))
    adj <- matrix(FALSE, nrow = L, ncol = L)
    lab_to_idx <- setNames(seq_len(L), labs)
    # Right neighbors
    if (ncol(seg) > 1) {
      rr <- seg[, -ncol(seg), drop = FALSE]
      rl <- seg[, -1, drop = FALSE]
      mask_rr <- !is.na(rr) & !is.na(rl) & rr != rl
      if (any(mask_rr)) {
        a <- lab_to_idx[as.character(rr[mask_rr])]
        b <- lab_to_idx[as.character(rl[mask_rr])]
        adj[cbind(pmin(a,b), pmax(a,b))] <- TRUE
      }
    }
    # Down neighbors
    if (nrow(seg) > 1) {
      dr <- seg[-nrow(seg), , drop = FALSE]
      dl <- seg[-1, , drop = FALSE]
      mask_dd <- !is.na(dr) & !is.na(dl) & dr != dl
      if (any(mask_dd)) {
        a <- lab_to_idx[as.character(dr[mask_dd])]
        b <- lab_to_idx[as.character(dl[mask_dd])]
        adj[cbind(pmin(a,b), pmax(a,b))] <- TRUE
      }
    }
    adj <- adj | t(adj)
    if (!any(adj)) return(list(seg = seg, centers = centers, sizes = sizes, done = TRUE))
    labs_char <- as.character(labs)
    C <- centers[labs_char, , drop = FALSE]
    best_i <- NA_integer_; best_j <- NA_integer_; best_d <- Inf
    for (i in seq_len(L-1)) {
      js <- which(adj[i, (i+1):L]) + i
      if (length(js)) {
        d <- apply(C[js, , drop = FALSE], 1, function(v) sum((v - C[i, ])^2))
        j <- js[which.min(d)]
        if (min(d) < best_d) {
          best_d <- min(d); best_i <- i; best_j <- j
        }
      }
    }
    if (is.infinite(best_d)) return(list(seg = seg, centers = centers, sizes = sizes, done = TRUE))
    li <- labs[best_i]; lj <- labs[best_j]
    seg[seg == lj] <- li
    si <- sizes[[as.character(li)]]; sj <- sizes[[as.character(lj)]]
    if (is.null(si) || is.na(si)) si <- sum(seg == li, na.rm = TRUE)
    if (is.null(sj) || is.na(sj)) sj <- sum(seg == lj, na.rm = TRUE)
    ci <- centers[as.character(li), , drop = FALSE]
    cj <- centers[as.character(lj), , drop = FALSE]
    cn <- as.numeric((ci * si + cj * sj) / pmax(1, (si + sj)))
    centers[as.character(li), ] <- cn
    centers <- centers[setdiff(rownames(centers), as.character(lj)), , drop = FALSE]
    sizes[as.character(li)] <- si + sj
    sizes <- sizes[setdiff(names(sizes), as.character(lj))]
    list(seg = seg, centers = centers, sizes = sizes, done = (length(unique(seg[!is.na(seg)])) <= targetK))
  }
  centers <- as.matrix(centers)
  if (is.null(rownames(centers))) {
    rownames(centers) <- as.character(seq_len(nrow(centers)))
  }
  if (is.null(sizes)) {
    labs <- sort(unique(seg[!is.na(seg)]))
    sizes <- setNames(as.numeric(table(seg[!is.na(seg)])), as.character(labs))
  }
  repeat {
    if (length(unique(seg[!is.na(seg)])) <= targetK) break
    ms <- merge_step(seg, centers, sizes)
    seg <- ms$seg; centers <- ms$centers; sizes <- ms$sizes
    if (isTRUE(ms$done)) break
  }
  list(seg = seg, centers = centers, sizes = sizes)
}

segment_image_kmeans <- function(img, K = 4, xy_weight = 0.15, seed = 123,
                                 mask_background = TRUE, white_thr = 0.92,
                                 use_lab = TRUE, apply_mode_filter = TRUE, mode_iters = 1,
                                 merge_to_K = NULL, K_init = NULL) {
  set.seed(seed)
  seg_feats <- prepare_segmentation_features(
    img = img,
    xy_weight = xy_weight,
    mask_background = mask_background,
    white_thr = white_thr,
    use_lab = use_lab
  )
  X <- seg_feats$X
  tissue_mask <- seg_feats$tissue_mask
  dims <- seg_feats$dims
  # Determine initial K for over-segmentation if requested
  nd <- nrow(unique(X))
  K0_req <- if (is.null(K_init)) K else K_init
  K0 <- max(1L, min(K0_req, nd))
  if (K0 < K0_req) warning(sprintf("Adjusted initial K from %d to %d due to limited distinct pixels", K0_req, K0))
  if (K0 == 1L) {
    km_res <- list(cluster = rep(1L, nrow(X)), centers = matrix(colMeans(X), nrow = 1), size = nrow(X))
  } else {
    km_res <- kmeans(X, centers = K0, iter.max = 50, nstart = 5)
  }
  km_cluster <- km_res$cluster
  seg <- matrix(NA_integer_, nrow = dims[['h']], ncol = dims[['w']])
  seg[as.vector(tissue_mask)] <- km_cluster
  centers <- km_res$centers
  # assign rownames to centers to match cluster labels
  rownames(centers) <- as.character(sort(unique(km_cluster)))
  sizes <- if (!is.null(km_res$size)) km_res$size else as.numeric(table(km_cluster))
  names(sizes) <- rownames(centers)

  # Optional hierarchical merging to target K
  if (!is.null(merge_to_K)) {
    merged <- merge_segments_to_K(seg, centers, merge_to_K, sizes)
    seg <- merged$seg
    centers <- merged$centers
    sizes <- merged$sizes
  }
  # optional 3x3 mode filter to smooth labels
  if (isTRUE(apply_mode_filter) && mode_iters > 0) {
    seg <- mode_filter_labels(seg, iterations = mode_iters)
  }
  list(map = seg, centers = centers)
}

segment_image_gmm <- function(img, K = 4, xy_weight = 0.15, seed = 123,
                              mask_background = TRUE, white_thr = 0.92,
                              use_lab = TRUE, apply_mode_filter = TRUE, mode_iters = 1,
                              merge_to_K = NULL) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    warning("Package 'mclust' not available; falling back to k-means segmentation.")
    return(segment_image_kmeans(img, K, xy_weight, seed, mask_background, white_thr,
                                use_lab, apply_mode_filter, mode_iters, merge_to_K))
  }
  set.seed(seed)
  seg_feats <- prepare_segmentation_features(
    img = img,
    xy_weight = xy_weight,
    mask_background = mask_background,
    white_thr = white_thr,
    use_lab = use_lab
  )
  X <- seg_feats$X
  tissue_mask <- seg_feats$tissue_mask
  dims <- seg_feats$dims
  if (K <= 1) {
    seg <- matrix(1L, nrow = dims[['h']], ncol = dims[['w']])
    seg[!tissue_mask] <- NA_integer_
    return(list(map = seg, centers = matrix(colMeans(X), nrow = 1)))
  }
  gmm <- tryCatch({
    mclust::Mclust(X, G = K, modelNames = "VVV", verbose = FALSE)
  }, error = function(e) {
    warning("GMM segmentation failed (", conditionMessage(e), "); falling back to k-means.")
    NULL
  })
  if (is.null(gmm)) {
    return(segment_image_kmeans(img, K, xy_weight, seed, mask_background, white_thr,
                                use_lab, apply_mode_filter, mode_iters, merge_to_K))
  }
  clusters <- as.integer(gmm$classification)
  seg <- matrix(NA_integer_, nrow = dims[['h']], ncol = dims[['w']])
  seg[as.vector(tissue_mask)] <- clusters
  centers <- t(gmm$parameters$mean)
  rownames(centers) <- as.character(seq_len(nrow(centers)))
  sizes <- as.numeric(table(clusters))
  names(sizes) <- as.character(seq_along(sizes))
  if (!is.null(merge_to_K)) {
    merged <- merge_segments_to_K(seg, centers, merge_to_K, sizes)
    seg <- merged$seg
    centers <- merged$centers
    sizes <- merged$sizes
  }
  if (isTRUE(apply_mode_filter) && mode_iters > 0) {
    seg <- mode_filter_labels(seg, iterations = mode_iters)
  }
  list(map = seg, centers = centers)
}

segment_image_slic <- function(img, K = 4, xy_weight = 0.15, seed = 123,
                               mask_background = TRUE, white_thr = 0.92,
                               use_lab = TRUE, apply_mode_filter = TRUE, mode_iters = 1,
                               merge_to_K = NULL, superpixel_count = 1200,
                               compactness = 20) {
  if (!requireNamespace("OpenImageR", quietly = TRUE)) {
    warning("Package 'OpenImageR' not available; falling back to k-means segmentation.")
    return(segment_image_kmeans(img, K, xy_weight, seed, mask_background, white_thr,
                                use_lab, apply_mode_filter, mode_iters, merge_to_K))
  }
  set.seed(seed)
  seg_feats <- prepare_segmentation_features(
    img = img,
    xy_weight = xy_weight,
    mask_background = mask_background,
    white_thr = white_thr,
    use_lab = use_lab
  )
  tissue_mask <- seg_feats$tissue_mask
  dims <- seg_feats$dims
  h <- dims[['h']]; w <- dims[['w']]
  img255 <- img
  img255[img255 < 0] <- 0
  img255[img255 > 1] <- 1
  img255 <- img255 * 255
  sp_res <- tryCatch({
    OpenImageR::superpixels(
      input_image = img255,
      method = "slic",
      superpixel = min(superpixel_count, h * w),
      compactness = compactness,
      return_slic_data = TRUE
    )
  }, error = function(e) {
    warning("SLIC superpixel segmentation failed (", conditionMessage(e), "); falling back to k-means.")
    NULL
  })
  if (is.null(sp_res) || is.null(sp_res$labels)) {
    return(segment_image_kmeans(img, K, xy_weight, seed, mask_background, white_thr,
                                use_lab, apply_mode_filter, mode_iters, merge_to_K))
  }
  labels <- matrix(sp_res$labels, nrow = h, ncol = w, byrow = TRUE)
  if (!all(dim(labels) == c(h, w))) {
    labels <- t(matrix(sp_res$labels, ncol = h, nrow = w))
  }
  sp_ids <- sort(unique(as.integer(labels)))
  sp_ids <- sp_ids[sp_ids > 0]
  n_sp <- length(sp_ids)
  feats <- matrix(NA_real_, nrow = n_sp, ncol = ncol(seg_feats$X))
  colnames(feats) <- colnames(seg_feats$X)
  for (i in seq_along(sp_ids)) {
    id <- sp_ids[i]
    mask <- labels == id & tissue_mask
    if (!any(mask)) {
      mask <- labels == id
    }
    idx <- seg_feats$index_map[mask]
    idx <- idx[!is.na(idx)]
    if (!length(idx)) {
      idx <- seg_feats$index_map[labels == id]
      idx <- idx[!is.na(idx)]
    }
    if (!length(idx)) {
      feats[i, ] <- colMeans(seg_feats$X)
    } else if (length(idx) == 1) {
      feats[i, ] <- seg_feats$X[idx, ]
    } else {
      feats[i, ] <- colMeans(seg_feats$X[idx, , drop = FALSE])
    }
  }
  if (n_sp <= K) {
    cluster_ids <- seq_len(n_sp)
  } else {
    cluster_ids <- kmeans(feats, centers = K, iter.max = 50, nstart = 5)$cluster
  }
  cluster_levels <- sort(unique(cluster_ids))
  cluster_map <- setNames(seq_along(cluster_levels), cluster_levels)
  cluster_ids <- cluster_map[as.character(cluster_ids)]
  sp_to_cluster <- setNames(as.integer(cluster_ids), sp_ids)
  seg <- matrix(NA_integer_, nrow = h, ncol = w)
  seg[tissue_mask] <- sp_to_cluster[as.character(labels[tissue_mask])]
  seg[is.na(seg) & labels > 0] <- sp_to_cluster[as.character(labels[is.na(seg) & labels > 0])]
  K_used <- length(cluster_levels)
  centers <- t(sapply(seq_len(K_used), function(k) colMeans(feats[cluster_ids == k, , drop = FALSE])))
  rownames(centers) <- as.character(seq_len(nrow(centers)))
  if (!is.null(merge_to_K)) {
    merged <- merge_segments_to_K(seg, centers, merge_to_K)
    seg <- merged$seg
    centers <- merged$centers
  }
  if (isTRUE(apply_mode_filter) && mode_iters > 0) {
    seg <- mode_filter_labels(seg, iterations = mode_iters)
  }
  list(map = seg, centers = centers)
}

# Build hex-grid coordinates aligned to image pixel space (base R)
generate_hex_coordinates <- function(width_px, height_px, spot_distance = 6) {
  x_len <- max(0, width_px - 1)
  y_len <- max(0, height_px - 1)
  # Hex pattern generation similar to simulator_YZ but in base R
  x <- seq(0, x_len, spot_distance)
  y <- seq(0, y_len, spot_distance * sqrt(3) / 2)
  if ((length(y) %% 2) == 0) {
    x_coord <- rep(c(x, x + spot_distance / 2), times = length(y) / 2)
  } else {
    x_coord <- c(rep(c(x, x + spot_distance / 2), times = length(y) %/% 2), x)
  }
  coord <- data.frame(x_coord = x_coord, y_coord = rep(y, each = length(x)))
  coord <- coord[coord$x_coord >= 0 & coord$x_coord <= x_len & coord$y_coord >= 0 & coord$y_coord <= y_len, , drop = FALSE]
  rownames(coord) <- NULL
  coord
}

# Map spot coordinates to nearest pixel and read domain label
assign_domains_to_spots <- function(coord, domain_map, search_radius = 3L) {
  h <- nrow(domain_map); w <- ncol(domain_map)
  xi <- pmin(w, pmax(1, round(coord$x_coord) + 1L))
  yi <- pmin(h, pmax(1, round(coord$y_coord) + 1L))
  domain <- domain_map[cbind(yi, xi)]
  # if NA (background), look around within a small window for nearest labeled pixel
  na_idx <- which(is.na(domain))
  if (length(na_idx) > 0) {
    for (k in na_idx) {
      y0 <- yi[k]; x0 <- xi[k]
      found <- NA_integer_
      rr <- search_radius
      y_min <- max(1L, y0 - rr); y_max <- min(h, y0 + rr)
      x_min <- max(1L, x0 - rr); x_max <- min(w, x0 + rr)
      window <- domain_map[y_min:y_max, x_min:x_max]
      vals <- unique(as.integer(window[!is.na(window)]))
      if (length(vals) > 0) {
        found <- vals[1]
      } else {
        found <- 0L
      }
      domain[k] <- found
    }
  }
  cbind(coord, layer = as.integer(domain))
}

# Sample per-domain cell-type proportions via Dirichlet using Gamma draws
dirichlet_gamma <- function(alpha) {
  r <- rgamma(length(alpha), shape = alpha, rate = 1)
  r / sum(r)
}

generate_layer_param <- function(n_layers, n_ct, alpha = 1) {
  mat <- t(replicate(n_layers, dirichlet_gamma(rep(alpha, n_ct))))
  colnames(mat) <- paste0('ct', seq_len(n_ct))
  as.data.frame(mat, check.names = FALSE)
}

# Minimal layer-based pattern helper (no dplyr)
generate_pattern_layer <- function(layer_vec, layer_param) {
  f <- as.integer(factor(layer_vec))
  out <- numeric(length(f))
  u <- sort(unique(f))
  for (i in u) out[f == i] <- layer_param[i]
  out
}

# Build mean_param list using domain-driven means for a subset of genes per CT
build_mean_param_from_domains <- function(coord, layer_vec, nGenes = 150, n_ct = 3,
                                          svg_per_ct = 5,
                                          mean_lo = 10, mean_hi = 150,
                                          seed = 123) {
  set.seed(seed)
  n_spots <- nrow(coord)
  genes <- paste0('gene_', seq_len(nGenes))
  layers <- sort(unique(as.integer(layer_vec)))
  n_layers <- length(layers)

  mean_param <- vector("list", n_ct)
  names(mean_param) <- paste0('ct', seq_len(n_ct))
  for (ct in seq_len(n_ct)) {
    mean_mat <- matrix(NA_real_, nrow = n_spots, ncol = nGenes)
    for (g in seq_len(nGenes)) {
      if (g <= svg_per_ct) {
        layer_means <- runif(n_layers, min = mean_lo, max = mean_hi)
        param <- generate_pattern_layer(layer_vec, layer_means)
      } else {
        param <- rep(runif(1, mean_lo, mean_hi), n_spots)
      }
      mean_mat[, g] <- param
    }
    colnames(mean_mat) <- genes
    mean_param[[ct]] <- mean_mat
  }
  return(mean_param)
}

# Simple negative binomial generator per CT returning gene x spot matrix
generate_true_count_simple <- function(mean_param, size = 100) {
  counts <- list()
  for (nm in names(mean_param)) {
    M <- mean_param[[nm]]             # spots x genes
    n_spots <- nrow(M); n_genes <- ncol(M)
    C <- t(apply(M, 1, function(mu_row) rnbinom(n_genes, size = size, mu = mu_row)))
    # Now C is spots x genes; transpose to genes x spots
    C <- t(C)
    rownames(C) <- paste0('gene_', seq_len(n_genes))
    colnames(C) <- paste0('spot_', seq_len(n_spots))
    counts[[nm]] <- C
  }
  counts
}

# No-op batch (keeps structure list[ct] = gene x spot)
generate_batch_count_simple <- function(true_count) true_count

save_segmentation_figure <- function(
  out_png,
  img,
  domain_map = NULL,
  coord = NULL,
  labels = NULL,
  spot_radius = 6,
  class_labels = NULL,
  class_palette = NULL,
  draw_boundaries = TRUE,
  white_thr = 0.92,
  title = NULL,
  overlay_alpha = 0.35
) {
  h <- dim(img)[1]; w <- dim(img)[2]
  # Determine K
  if (!is.null(domain_map)) {
    domain_vec <- as.integer(domain_map)
  } else if (!is.null(labels)) {
    domain_vec <- as.integer(labels)
  } else {
    stop("Either domain_map or labels must be provided")
  }
  labs <- sort(unique(domain_vec[!is.na(domain_vec) & domain_vec != 0]))
  K <- length(labs)
  if (is.null(class_palette)) {
    class_palette <- tryCatch({
      suppressPackageStartupMessages(library(scales))
      scales::hue_pal()(K)
    }, error = function(e) grDevices::rainbow(K))
  }
  if (is.null(class_labels)) {
    class_labels <- paste("Class", seq_len(K))
  } else if (length(class_labels) != K) {
    class_labels <- class_labels[seq_len(K)]
  }
  pal_map <- setNames(class_palette, labs)

  png(out_png, width = w, height = h)
  par(mar = c(2,2,2,2), xaxs = 'i', yaxs = 'i')
  plot.new(); plot.window(xlim = c(0, w), ylim = c(0, h))
  # background image
  try({ rasterImage(grDevices::as.raster(img), 0, 0, w, h, interpolate = TRUE) }, silent = TRUE)

  if (!is.null(domain_map)) {
    # Draw boundaries per class
    if (isTRUE(draw_boundaries)) {
      hz <- nrow(domain_map); wz <- ncol(domain_map)
      sx <- ifelse(wz > 0, w / wz, 1)
      sy <- ifelse(hz > 0, h / hz, 1)
      xs <- (0:(wz - 1)) * sx
      ys <- (0:(hz - 1)) * sy
      ok <- try({
        for (i in seq_along(labs)) {
          k <- labs[i]
          z <- ifelse(domain_map == k, 1, 0)
          # Flip rows to align with rasterImage orientation
          zf <- z[hz:1, , drop = FALSE]
          contour(x = xs, y = ys, z = zf,
                  levels = 0.5, add = TRUE, drawlabels = FALSE,
                  col = class_palette[i], lwd = 3)
        }
        # Tissue outline (non-NA region)
        mask <- ifelse(is.na(domain_map), 0, 1)
        mf <- mask[hz:1, , drop = FALSE]
        contour(x = xs, y = ys, z = mf,
                levels = 0.5, add = TRUE, drawlabels = FALSE,
                col = "red", lwd = 4, lty = 1)
        TRUE
      }, silent = TRUE)
      if (inherits(ok, "try-error") || isFALSE(ok)) {
        # Fallback: draw label raster overlay
        lab_mat <- matrix(NA_character_, nrow = hz, ncol = wz)
        for (i in seq_along(labs)) {
          lab_mat[domain_map == labs[i]] <- class_palette[i]
        }
        lab_mat[is.na(domain_map)] <- NA
        # Apply alpha to palette colors
        add_alpha <- function(col, a) {
          cc <- col2rgb(col, TRUE)/255
          rgb(cc[1,], cc[2,], cc[3,], alpha = a)
        }
        lab_col <- ifelse(is.na(lab_mat), NA, add_alpha(lab_mat, overlay_alpha))
        lab_col <- lab_col[hz:1, , drop = FALSE]
        rasterImage(grDevices::as.raster(lab_col), 0, 0, w, h, interpolate = FALSE)
      }
    }
    # Optional: faint spot centers
    if (!is.null(coord)) {
      points(coord$x_coord + 1, h - (coord$y_coord + 1), pch = 16, col = "#00000040", cex = 0.3)
    }
  } else if (!is.null(coord) && !is.null(labels)) {
    # Spots-only visualization: draw colored disks
    cols <- class_palette[as.integer(labels)]
    x_plot <- coord$x_coord + 1
    y_plot <- h - (coord$y_coord + 1)
    symbols(x_plot, y_plot, circles = rep(spot_radius, length(x_plot)), inches = FALSE,
            bg = cols, fg = NA, add = TRUE)
  }

  if (!is.null(title)) title(main = title, cex.main = 1.2)
  legend("topright", legend = class_labels, fill = class_palette, cex = 0.8, bg = "white")
  dev.off()
}

simulator_image <- function(
  image_path,
  out_prefix = "./Image_feature/simulator_image",
  K_domains = 7,
  spot_distance = 6,
  nGenes = 150,
  n_celltypes = 3,
  svg_per_ct = 5,
  seq_depth = 1e6,
  xy_weight = 0.15,
  seed = 123,
  # background control
  mask_background = TRUE,
  white_thr = 0.92,
  # mapping control
  use_mapping = FALSE,
  sample = NULL,
  mapping_path = NULL,
  mapping_dir = "./data_external/HER2/polaris_processed",
  mapping_cols = c('X','x','y','new_x','new_y','pixel_x','pixel_y','selected','new_spotID','Row.names','layer','pixel_spotID'),
  mapping_selected_col = 'selected',
  mapping_x_col = 'pixel_x',
  mapping_y_col = 'pixel_y',
  mapping_truth_layer_col = 'layer',
  # spots-only segmentation
  spots_only = FALSE,
  spot_radius = 6,
  # segmentation-only mode (no simulation or proportions)
  segment_only = FALSE,
  segmentation_method = c("kmeans", "gmm", "slic"),
  slic_superpixel_count = 1200,
  slic_compactness = 20,
  # quality controls
  use_lab = TRUE,
  add_texture = TRUE,
  smooth_labels = TRUE,
  smooth_k = 8,
  smooth_iter = 1
) {
  set.seed(seed)

  # 1) Read and (optionally) downscale image
  img0 <- read_image_any(image_path)
  img  <- downscale_image(img0, max_side = 512)
  h <- dim(img)[1]; w <- dim(img)[2]
  h0 <- dim(img0)[1]; w0 <- dim(img0)[2]
  segmentation_method <- match.arg(segmentation_method)

  # 2) Build coordinates in pixel space (mapping or hex grid)
  if (use_mapping) {
    if (is.null(mapping_path)) {
      if (is.null(sample)) stop("Provide either mapping_path or sample")
      mapping_path <- file.path(mapping_dir, paste0(sample, "_mapping.tsv"))
    }
    m1 <- tryCatch(read.table(mapping_path, sep = '\t', header = TRUE, check.names = FALSE), error = function(e) NULL)
    if (is.null(m1)) stop("Failed to read mapping file: ", mapping_path)
    if (length(colnames(m1)) == length(mapping_cols)) colnames(m1) <- mapping_cols
    # Filter to selected tissue spots
    sel <- m1[[mapping_selected_col]]
    if (is.logical(sel)) {
      m1 <- m1[sel, , drop = FALSE]
    } else {
      # treat as 0/1 or string 'TRUE'/'FALSE'
      sel2 <- suppressWarnings(as.numeric(sel))
      if (all(is.na(sel2))) sel2 <- tolower(as.character(sel)) %in% c('true','t','yes','y','selected','1')
      m1 <- m1[as.logical(sel2), , drop = FALSE]
    }
    # Build coord from pixel columns; ensure within image bounds
    cx0 <- as.numeric(m1[[mapping_x_col]])
    cy0 <- as.numeric(m1[[mapping_y_col]])
    # scale original pixel coords to downscaled image size
    cx <- cx0 * (w / w0)
    cy <- cy0 * (h / h0)
    coord <- data.frame(x_coord = pmin(w - 1, pmax(0, cx)),
                        y_coord = pmin(h - 1, pmax(0, cy)))
    rownames(coord) <- if ('Row.names' %in% colnames(m1)) as.character(m1[['Row.names']]) else NULL
  } else {
    coord <- generate_hex_coordinates(width_px = w, height_px = h, spot_distance = spot_distance)
  }

  # 3) Segment into K domains: either per-pixel or spots-only
  domain_map <- NULL
  K_vec <- unique(as.integer(K_domains))
  if (spots_only) {
    # Extract features for each spot region and cluster spots only
    compute_spot_features <- function(img, coord, radius = 6, xy_weight = 0.15, white_thr = 0.92,
                                      use_lab = TRUE, add_texture = TRUE) {
      h <- dim(img)[1]; w <- dim(img)[2]
      n <- nrow(coord)
      # Precompute helpers
      lab <- if (isTRUE(use_lab)) toLab_image(img) else NULL
      grad <- if (isTRUE(add_texture)) compute_gradients(img) else NULL
      feats <- NULL
      # Build features per spot with preallocation after first fv
      for (i in seq_len(n)) {
        x0 <- min(w, max(1L, round(coord$x_coord[i]) + 1L))
        y0 <- min(h, max(1L, round(coord$y_coord[i]) + 1L))
        xr <- max(1L, x0 - radius):min(w, x0 + radius)
        yr <- max(1L, y0 - radius):min(h, y0 + radius)
        patch <- img[yr, xr, , drop = FALSE]
        if (!is.null(lab)) {
          Lp <- as.vector(lab$L[yr, xr]); ap <- as.vector(lab$a[yr, xr]); bp <- as.vector(lab$b[yr, xr])
          lum <- Lp / max(1, max(lab$L, na.rm = TRUE))
        } else {
          Rp <- as.vector(patch[,,1]); Gp <- as.vector(patch[,,2]); Bp <- as.vector(patch[,,3])
          lum <- (Rp + Gp + Bp) / 3
        }
        keep <- lum < white_thr
        if (!is.null(lab)) {
          Lm <- if (any(keep)) mean(Lp[keep]) else mean(Lp)
          am <- if (any(keep)) mean(ap[keep]) else mean(ap)
          bm <- if (any(keep)) mean(bp[keep]) else mean(bp)
          fv <- c(Lm, am, bm)
        } else {
          Rp <- as.vector(patch[,,1]); Gp <- as.vector(patch[,,2]); Bp <- as.vector(patch[,,3])
          Rm <- if (any(keep)) mean(Rp[keep]) else mean(Rp)
          Gm <- if (any(keep)) mean(Gp[keep]) else mean(Gp)
          Bm <- if (any(keep)) mean(Bp[keep]) else mean(Bp)
          fv <- c(Rm, Gm, Bm)
        }
        # add texture features
        if (!is.null(grad)) {
          gp <- as.vector(grad[yr, xr])
          gm <- if (any(keep)) mean(gp[keep]) else mean(gp)
          lv <- if (any(keep)) var(lum[keep]) else var(lum)
          if (is.na(lv)) lv <- 0
          fv <- c(fv, gm, lv)
        }
        # add XY
        fv <- c(fv,
                (x0 - 1) / max(1, (w - 1)) * xy_weight,
                (y0 - 1) / max(1, (h - 1)) * xy_weight)
        if (is.null(feats)) {
          feats <- matrix(NA_real_, nrow = n, ncol = length(fv))
          colnames(feats) <- c(if (!is.null(lab)) c('L','a','b') else c('R','G','B'),
                               if (!is.null(grad)) c('gmag','lvar') else NULL,
                               c('x','y'))
        }
        feats[i, ] <- fv
      }
      # standardize safely
      feats <- robust_scale(feats)
      feats
    }
    Xs <- compute_spot_features(img, coord, radius = spot_radius, xy_weight = xy_weight, white_thr = white_thr,
                                use_lab = use_lab, add_texture = add_texture)
    lc_list <- list()
    cluster_spot_features <- function(X, K, method, seed, label) {
      nd <- nrow(unique(X))
      k_eff <- max(1L, min(K, nd))
      if (k_eff < K) warning(sprintf("[%s] Adjusted K from %d to %d due to limited distinct spots", label, K, k_eff))
      if (k_eff <= 1L) return(rep(1L, nrow(X)))
      set.seed(seed)
      if (method == "gmm") {
        if (!requireNamespace("mclust", quietly = TRUE)) {
          warning("Package 'mclust' not available; falling back to k-means for spot clustering.")
        } else {
          gmm <- tryCatch({
            mclust::Mclust(X, G = k_eff, modelNames = "VVV", verbose = FALSE)
          }, error = function(e) NULL)
          if (!is.null(gmm)) return(as.integer(gmm$classification))
          warning("Spot-level GMM failed; falling back to k-means.")
        }
        method <- "kmeans"
      }
      if (method == "slic") {
        warning("Spot-level SLIC is not supported; using k-means on spot features instead.")
        method <- "kmeans"
      }
      kmeans(X, centers = k_eff, iter.max = 50, nstart = 5)$cluster
    }
    for (k in K_vec) {
      cl <- cluster_spot_features(Xs, k, segmentation_method, seed, paste0("K", k))
      lc <- cbind(coord, layer_seg = as.integer(cl))
      if (use_mapping && exists("m1") && mapping_truth_layer_col %in% colnames(m1)) {
        lc$layer_mapping <- m1[[mapping_truth_layer_col]]
      }
      # simple neighbor smoothing by majority vote
      if (isTRUE(smooth_labels) && smooth_iter > 0) {
        pts <- as.matrix(coord[, c("x_coord", "y_coord")])
        D <- as.matrix(dist(pts))
        for (t in seq_len(smooth_iter)) {
          labs <- lc$layer_seg
          new_labs <- labs
          for (i in seq_len(nrow(pts))) {
            ord <- order(D[i, ])
            nn <- ord[ord != i][1:min(smooth_k, nrow(pts) - 1)]
            nb <- labs[nn]
            tab <- table(nb)
            new_labs[i] <- as.integer(names(tab)[which.max(tab)])
          }
          lc$layer_seg <- new_labs
        }
      }
      lc_list[[as.character(k)]] <- lc
    }
    # If multiple K and segmentation-only requested, write per-K results and return
    if (isTRUE(segment_only) && length(K_vec) > 1) {
      out_dir <- dirname(out_prefix)
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      # write coordinates once
      write.csv(coord, paste0(out_prefix, "_coord.csv"), row.names = FALSE)
      # write per-K layers and PNGs
      for (k in K_vec) {
        lc <- lc_list[[as.character(k)]]
        write.csv(lc, paste0(out_prefix, "_layer_K", k, ".csv"), row.names = TRUE)
        # PNG overlay
        out_png <- paste0(out_prefix, "_segmentation_K", k, ".png")
        h <- dim(img)[1]; w <- dim(img)[2]
        domain_vec <- lc$layer_seg
        KK <- length(unique(domain_vec))
        pal <- tryCatch({
          suppressPackageStartupMessages(library(scales))
          scales::hue_pal()(KK)
        }, error = function(e) grDevices::rainbow(KK))
        cols <- pal[as.integer(domain_vec)]
        png(out_png, width = w, height = h)
        par(mar = c(0,0,0,0), xaxs = 'i', yaxs = 'i')
        plot.new(); plot.window(xlim = c(0, w), ylim = c(0, h))
        try({ rasterImage(grDevices::as.raster(img), 0, 0, w, h, interpolate = TRUE) }, silent = TRUE)
        x_plot <- lc$x_coord + 1
        y_plot <- h - (lc$y_coord + 1)
        symbols(x_plot, y_plot, circles = rep(spot_radius, length(x_plot)), inches = FALSE,
                bg = cols, fg = NA, add = TRUE)
        dev.off()
      }
      return(invisible(lc_list))
    }
    # single K case: proceed with one set
    layer_coord <- lc_list[[as.character(K_vec[1])]]
  } else {
    # Per-pixel segmentation then map spots to nearest pixel domain
    lc_list <- list(); dm_list <- list()
    for (k in K_vec) {
      over_k <- max(k + 4L, k * 2L)
      seg <- switch(segmentation_method,
                    kmeans = segment_image_kmeans(
                      img, K = k, xy_weight = xy_weight, seed = seed,
                      mask_background = mask_background, white_thr = white_thr,
                      use_lab = use_lab, apply_mode_filter = TRUE, mode_iters = 1,
                      merge_to_K = k, K_init = over_k
                    ),
                    gmm = segment_image_gmm(
                      img, K = k, xy_weight = xy_weight, seed = seed,
                      mask_background = mask_background, white_thr = white_thr,
                      use_lab = use_lab, apply_mode_filter = TRUE, mode_iters = 1,
                      merge_to_K = k
                    ),
                    slic = segment_image_slic(
                      img, K = k, xy_weight = xy_weight, seed = seed,
                      mask_background = mask_background, white_thr = white_thr,
                      use_lab = use_lab, apply_mode_filter = TRUE, mode_iters = 1,
                      merge_to_K = k, superpixel_count = slic_superpixel_count,
                      compactness = slic_compactness
                    ))
      dm <- seg$map
      lc <- assign_domains_to_spots(coord, dm)
      if (use_mapping && exists("m1") && mapping_truth_layer_col %in% colnames(m1)) {
        lc$layer_mapping <- m1[[mapping_truth_layer_col]]
        names(lc)[names(lc) == 'layer'] <- 'layer_seg'
      }
      lc_list[[as.character(k)]] <- lc
      dm_list[[as.character(k)]] <- dm
    }
    if (isTRUE(segment_only) && length(K_vec) > 1) {
      out_dir <- dirname(out_prefix)
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      write.csv(coord, paste0(out_prefix, "_coord.csv"), row.names = FALSE)
      for (k in K_vec) {
        lc <- lc_list[[as.character(k)]]
        write.csv(lc, paste0(out_prefix, "_layer_K", k, ".csv"), row.names = TRUE)
        # prefer spots overlay for clarity
        out_png <- paste0(out_prefix, "_segmentation_K", k, ".png")
        h <- dim(img)[1]; w <- dim(img)[2]
        domain_vec <- if ("layer_seg" %in% names(lc)) lc$layer_seg else lc$layer
        KK <- length(unique(domain_vec))
        pal <- tryCatch({
          suppressPackageStartupMessages(library(scales))
          scales::hue_pal()(KK)
        }, error = function(e) grDevices::rainbow(KK))
        cols <- pal[as.integer(domain_vec)]
        png(out_png, width = w, height = h)
        par(mar = c(0,0,0,0), xaxs = 'i', yaxs = 'i')
        plot.new(); plot.window(xlim = c(0, w), ylim = c(0, h))
        try({ rasterImage(grDevices::as.raster(img), 0, 0, w, h, interpolate = TRUE) }, silent = TRUE)
        x_plot <- lc$x_coord + 1
        y_plot <- h - (lc$y_coord + 1)
        symbols(x_plot, y_plot, circles = rep(spot_radius, length(x_plot)), inches = FALSE,
                bg = cols, fg = NA, add = TRUE)
        dev.off()
      }
      return(invisible(lc_list))
    }
    # single K case
    layer_coord <- lc_list[[as.character(K_vec[1])]]
    domain_map  <- dm_list[[as.character(K_vec[1])]]
  }

  # If only segmentation is requested, write minimal outputs and return
  if (isTRUE(segment_only)) {
    out_dir <- dirname(out_prefix)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    # Minimal outputs: coordinates, segmentation labels, and preview PNG
    write.csv(coord,            paste0(out_prefix, "_coord.csv"), row.names = FALSE)
    write.csv(layer_coord,      paste0(out_prefix, "_layer.csv"), row.names = TRUE)
    # Preview image with boundaries and legend
    Kcur <- length(unique((if ("layer_seg" %in% names(layer_coord)) layer_coord$layer_seg else layer_coord$layer)))
    class_labels <- paste0("cluster ", seq_len(Kcur))
    tryCatch({
      if (!is.null(domain_map)) {
        save_segmentation_figure(
          out_png = paste0(out_prefix, "_segmentation.png"),
          img = img,
          domain_map = domain_map,
          coord = coord,
          labels = NULL,
          spot_radius = spot_radius,
          class_labels = class_labels,
          class_palette = NULL,
          draw_boundaries = TRUE,
          white_thr = white_thr,
          title = NULL
        )
      } else {
        layer_col <- if ("layer_seg" %in% names(layer_coord)) "layer_seg" else "layer"
        save_segmentation_figure(
          out_png = paste0(out_prefix, "_segmentation.png"),
          img = img,
          domain_map = NULL,
          coord = coord,
          labels = layer_coord[[layer_col]],
          spot_radius = spot_radius,
          class_labels = class_labels,
          class_palette = NULL,
          draw_boundaries = FALSE,
          white_thr = white_thr,
          title = NULL
        )
      }
    }, error = function(e) message("Warning: could not write segmentation preview PNG: ", e$message))
    return(invisible(list(coord = coord, layer = layer_coord)))
  }

  # 4) Domain-specific cell-type proportions
  layer_col <- if ("layer_seg" %in% names(layer_coord)) "layer_seg" else "layer"
  domain_vec <- layer_coord[[layer_col]]
  n_layers <- length(unique(domain_vec))
  layer_param <- generate_layer_param(n_layers = n_layers, n_ct = n_celltypes, alpha = 1)
  # Generate celltype proportions per spot based on layer
  generate_cell_prop <- function(domain_vec, layer_param) {
    f <- as.integer(factor(domain_vec))
    n_spots <- length(domain_vec)
    n_ct <- ncol(layer_param)
    param <- matrix(0, nrow = n_spots, ncol = n_ct)
    for (i in sort(unique(f))) {
      idx <- which(f == i)
      for (j in seq_len(n_ct)) {
        mu <- layer_param[i, j]
        # small noise around domain mean, truncated at 0
        param[idx, j] <- pmax(0, rnorm(length(idx), mean = mu, sd = 1e-4))
      }
    }
    rs <- rowSums(param)
    rs[rs == 0] <- 1
    param / rs
  }
  celltype_proportion <- generate_cell_prop(domain_vec, layer_param)

  # 5) Build mean parameters using domain info for a subset of genes
  mean_param <- build_mean_param_from_domains(
    coord = coord,
    layer_vec = domain_vec,
    nGenes = nGenes,
    n_ct = n_celltypes,
    svg_per_ct = svg_per_ct,
    seed = seed
  )

  # 6) Simulate one latent spot-level profile per cell type and combine the
  # profiles using their continuous cell-equivalent abundances.
  generate_spot_count_abundance_simple <- function(
      mean_param, size, cell_number, celltype_proportion) {
    prop <- as.matrix(celltype_proportion)
    prop <- prop / rowSums(prop)
    n_spots <- nrow(prop)
    if (length(cell_number) == 1L) cell_number <- rep(cell_number, n_spots)
    spot_count <- matrix(0, nrow = n_spots, ncol = nGenes)
    for (i in seq_along(mean_param)) {
      profile <- matrix(
        rnbinom(n_spots * nGenes, size = size, mu = as.numeric(mean_param[[i]])),
        nrow = n_spots, ncol = nGenes
      )
      abundance <- cell_number * prop[, i]
      spot_count <- spot_count + sweep(profile, 1L, abundance, `*`)
    }
    spot_count
  }
  spot_count <- generate_spot_count_abundance_simple(
    mean_param, size = 100, cell_number = rep(10L, nrow(coord)),
    celltype_proportion = celltype_proportion
  )

  # 7) Save outputs
  out_dir <- dirname(out_prefix)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  coord_dt <- coord
  rownames(spot_count) <- paste0('spot_', seq_len(nrow(spot_count)))
  colnames(spot_count) <- paste0('gene_', seq_len(ncol(spot_count)))
  rownames(layer_coord) <- paste0('spot_', seq_len(nrow(layer_coord)))
  
  # Write outputs using base R to avoid extra deps
  write.csv(coord_dt,      paste0(out_prefix, "_coord.csv"), row.names = FALSE)
  write.csv(layer_coord,   paste0(out_prefix, "_layer.csv"), row.names = TRUE)
  write.csv(spot_count,    paste0(out_prefix, "_expression.csv"))
  write.csv(as.data.frame(celltype_proportion), paste0(out_prefix, "_celltype_proportion.csv"), row.names = FALSE)
  # Previews with boundaries and legend
  tryCatch({
    Kcur <- length(unique(domain_vec))
    class_labels <- paste0("cluster ", seq_len(Kcur))
    if (!is.null(domain_map)) {
      save_segmentation_figure(
        out_png = paste0(out_prefix, "_segmentation.png"),
        img = img,
        domain_map = domain_map,
        coord = coord,
        labels = NULL,
        spot_radius = spot_radius,
        class_labels = class_labels,
        class_palette = NULL,
        draw_boundaries = TRUE,
        white_thr = white_thr,
        title = NULL
      )
    } else {
      save_segmentation_figure(
        out_png = paste0(out_prefix, "_segmentation.png"),
        img = img,
        domain_map = NULL,
        coord = coord,
        labels = domain_vec,
        spot_radius = spot_radius,
        class_labels = class_labels,
        class_palette = NULL,
        draw_boundaries = FALSE,
        white_thr = white_thr,
        title = NULL
      )
    }
  }, error = function(e) { message("Warning: could not write segmentation preview PNG: ", e$message) })

  invisible(list(
    coord = coord_dt,
    layer = layer_coord,
    celltype_proportion = celltype_proportion,
    mean_param = mean_param,
    counts = list(true = true_count, seq = seq_count, spot = spot_count),
    segmentation = domain_map
  ))
}

if (sys.nframe() == 0) {
  # Default run when executed as a script
  image_path <- "./data_external/HER2/9769_C1_HE_small.jpg"
  out_prefix <- "./Image_feature/simulator_image"
  message("Running simulator_image with default parameters...\n  Image: ", image_path, "\n  Output prefix: ", out_prefix)
  res <- simulator_image(
    image_path = image_path,
    out_prefix = out_prefix,
    K_domains = 7,
    spot_distance = 6,
    nGenes = 150,
    n_celltypes = 3,
    svg_per_ct = 5,
    seq_depth = 1e6,
    xy_weight = 0.15,
    seed = 123
  )
  message("Done. Outputs written with prefix: ", out_prefix)
}
