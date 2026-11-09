#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(STCOS)
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
rep_id <- if (length(args) >= 1L) suppressWarnings(as.integer(args[1])) else 1L
if (is.na(rep_id) || !rep_id %in% 1:5) {
  stop("replicate must be an integer from 1 to 5")
}

# Load STCOS helpers --------------------------------------------------------

normalize_vec <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (any(!is.finite(rng)) || diff(rng) < 1e-8) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

normalize_mix <- function(x, ct_names) {
  out <- setNames(rep(0, length(ct_names)), ct_names)
  out[names(x)] <- x
  out[out < 0] <- 0
  s <- sum(out)
  if (s <= 0) {
    out["CAF"] <- 1
    return(out)
  }
  out / s
}

ellipse_masks <- function(coord, center, a, b, theta, core_frac = 0.65) {
  dx <- coord$x_coord - center[1]
  dy <- coord$y_coord - center[2]
  xr <- dx * cos(theta) + dy * sin(theta)
  yr <- -dx * sin(theta) + dy * cos(theta)
  scaled <- sqrt((xr / a)^2 + (yr / b)^2)
  core <- scaled <= core_frac
  boundary <- scaled > core_frac & scaled <= 1
  list(core = core, boundary = boundary)
}

stripe_masks <- function(coord, p1, p2, core_width, boundary_width) {
  x0 <- p1[1]; y0 <- p1[2]; x1 <- p2[1]; y1 <- p2[2]
  A <- y0 - y1
  B <- x1 - x0
  C <- x0 * y1 - x1 * y0
  denom <- sqrt(A^2 + B^2)
  if (denom == 0) {
    return(list(core = rep(FALSE, nrow(coord)), boundary = rep(FALSE, nrow(coord))))
  }
  dist <- abs(A * coord$x_coord + B * coord$y_coord + C) / denom
  dx <- x1 - x0
  dy <- y1 - y0
  seg_len_sq <- dx^2 + dy^2
  t_param <- ((coord$x_coord - x0) * dx + (coord$y_coord - y0) * dy) / seg_len_sq
  on_segment <- t_param >= -0.2 & t_param <= 1.2
  core <- dist <= core_width & on_segment
  boundary <- dist > core_width & dist <= boundary_width & on_segment
  list(core = core, boundary = boundary)
}

apply_assignment <- function(prop_mat, layer_vec, region_vec, priority_vec,
                             mask, mix, layer_value, region_label, priority) {
  idx <- which(mask & priority >= priority_vec)
  if (length(idx) == 0) {
    return(list(prop_mat = prop_mat,
                layer_vec = layer_vec,
                region_vec = region_vec,
                priority_vec = priority_vec))
  }
  prop_mat[idx, ] <- matrix(rep(mix, each = length(idx)),
                            nrow = length(idx), byrow = FALSE)
  layer_vec[idx] <- layer_value
  region_vec[idx] <- region_label
  priority_vec[idx] <- priority
  list(prop_mat = prop_mat,
       layer_vec = layer_vec,
       region_vec = region_vec,
       priority_vec = priority_vec)
}

aggregate_spot_counts_abundance_weighted <- function(count_list, cell_numbers, celltype_prop) {
  n_spots <- nrow(celltype_prop)
  if (length(count_list) == 0) {
    stop("Empty count list provided.")
  }
  first_mat <- as.matrix(count_list[[1]])
  if (nrow(first_mat) != n_spots && ncol(first_mat) == n_spots) {
    first_mat <- t(first_mat)
  }
  if (nrow(first_mat) != n_spots) {
    stop("Spot dimension mismatch in count list.")
  }
  n_genes <- ncol(first_mat)
  spot_count <- matrix(0, nrow = n_spots, ncol = n_genes)
  for (ct in names(count_list)) {
    if (!ct %in% colnames(celltype_prop)) {
      next
    }
    mat <- as.matrix(count_list[[ct]])
    if (nrow(mat) != n_spots && ncol(mat) == n_spots) {
      mat <- t(mat)
    }
    if (nrow(mat) != n_spots) {
      stop(sprintf("Spot dimension mismatch for cell type %s.", ct))
    }
    weights <- cell_numbers * celltype_prop[, ct]
    spot_count <- spot_count + sweep(mat, 1, weights, `*`)
  }
  colnames(spot_count) <- paste0("gene_", seq_len(n_genes))
  spot_count
}

make_mean_param <- function(coord, prop_mat, ct_names, n_genes, svg_gene_sets,
                            marker_gene_sets, rep_id, seed) {
  set.seed(seed + 5000L)
  n_spots <- nrow(coord)
  center_x <- mean(coord$x_coord)
  center_y <- mean(coord$y_coord)
  grad_x <- normalize_vec(coord$x_coord)
  grad_y <- normalize_vec(coord$y_coord)
  radial <- normalize_vec((coord$x_coord - center_x)^2 + (coord$y_coord - center_y)^2)
  wave_a <- normalize_vec(sin(coord$x_coord * 0.21 + rep_id * 0.7) +
                            cos(coord$y_coord * 0.17 + rep_id * 0.3))
  wave_b <- normalize_vec(cos(coord$x_coord * 0.19 - rep_id * 0.4) -
                            sin(coord$y_coord * 0.23 + rep_id * 0.5))
  modifiers <- list(grad_x, grad_y, 1 - radial, wave_a, wave_b)

  mean_param <- setNames(vector("list", length(ct_names)), ct_names)
  for (ct in ct_names) {
    base_gene <- runif(n_genes, 3, 7)
    mat <- matrix(rep(base_gene, each = n_spots), nrow = n_spots)
    mat <- mat + matrix(rnorm(n_spots * n_genes, mean = 0, sd = 0.35),
                        nrow = n_spots, ncol = n_genes)

    ct_prop <- prop_mat[, ct]
    core_signal <- normalize_vec(ct_prop^1.25)
    boundary_signal <- normalize_vec(sqrt(pmax(ct_prop * (1 - ct_prop), 0)))

    svg_idx <- svg_gene_sets[[ct]]
    core_amp <- runif(5, 24, 34)
    edge_amp <- runif(5, 20, 30)
    for (j in 1:5) {
      signal <- normalize_vec(0.78 * core_signal + 0.22 * modifiers[[j]])
      mat[, svg_idx[j]] <- mat[, svg_idx[j]] + core_amp[j] * signal
    }
    for (j in 1:5) {
      mod_idx <- ((j + 2L - 1L) %% 5L) + 1L
      signal <- normalize_vec(0.70 * boundary_signal + 0.30 * modifiers[[mod_idx]])
      mat[, svg_idx[5 + j]] <- mat[, svg_idx[5 + j]] + edge_amp[j] * signal
    }

    marker_idx <- marker_gene_sets[[ct]]
    marker_boost <- runif(length(marker_idx), 12, 20)
    mat[, marker_idx] <- sweep(mat[, marker_idx, drop = FALSE], 2, marker_boost, `+`)

    mat[mat < 0.1] <- 0.1
    colnames(mat) <- paste0("gene_", seq_len(n_genes))
    rownames(mat) <- paste0("spot_", seq_len(n_spots))
    mean_param[[ct]] <- mat
  }
  mean_param
}

# Global simulation settings -------------------------------------------------
ct_names <- c("Cancer", "Tcell", "Endothelial", "CAF")
nGenes <- 200
spot_distance <- 1.5
gene_ids <- paste0("gene_", seq_len(nGenes))

svg_gene_sets <- list(
  Cancer = 1:10,
  Tcell = 11:20,
  Endothelial = 21:30,
  CAF = 31:40
)

marker_gene_sets <- list(
  Cancer = 41:50,
  Tcell = 51:60,
  Endothelial = 61:70,
  CAF = 71:80
)

background_gene_sets <- list(
  shared_a = 81:90,
  shared_b = 91:100
)

base_gene_relationship <- data.table(
  matrix(NA_integer_, nrow = nGenes, ncol = length(ct_names))
)
setnames(base_gene_relationship, ct_names)
setattr(base_gene_relationship, "row.names", paste0("gene", seq_len(nGenes)))

for (ct in ct_names) {
  svg_idx <- svg_gene_sets[[ct]]
  marker_idx <- marker_gene_sets[[ct]]
  base_gene_relationship[svg_idx[1:5], (ct) := 1L]
  base_gene_relationship[svg_idx[6:10], (ct) := 2L]
  base_gene_relationship[marker_idx[1:5], (ct) := 3L]
  base_gene_relationship[marker_idx[6:10], (ct) := 4L]
  base_gene_relationship[background_gene_sets$shared_a, (ct) := 5L]
  base_gene_relationship[background_gene_sets$shared_b, (ct) := 6L]
}

replicate_designs <- list(
  list(
    id = 1,
    seed = 7101,
    label = "central_tumor_ring",
    background = c(Cancer = 0.05, Tcell = 0.10, Endothelial = 0.15, CAF = 0.70),
    tumors = list(
      list(center = c(25, 25), a = 11.5, b = 8.5, theta = pi / 10, core_frac = 0.63),
      list(center = c(33, 31), a = 5.8, b = 4.0, theta = -pi / 8, core_frac = 0.60)
    ),
    immune = list(
      list(center = c(16, 27), a = 6.2, b = 4.5, theta = pi / 9, core_frac = 0.62),
      list(center = c(28, 37), a = 5.0, b = 4.0, theta = pi / 4, core_frac = 0.60)
    ),
    vessels = list(
      list(p1 = c(4, 12), p2 = c(48, 31), core_width = 1.0, boundary_width = 2.0),
      list(p1 = c(7, 43), p2 = c(45, 19), core_width = 0.9, boundary_width = 1.8)
    ),
    caf_islands = list(
      list(center = c(9, 8), a = 4.8, b = 3.5, theta = pi / 6, core_frac = 0.65),
      list(center = c(42, 43), a = 4.5, b = 3.2, theta = -pi / 6, core_frac = 0.65)
    ),
    cell_per_spot = 9,
    p0 = 0.05,
    seq_depth = 5.5e5
  ),
  list(
    id = 2,
    seed = 7102,
    label = "bi_lobed_tumor",
    background = c(Cancer = 0.08, Tcell = 0.12, Endothelial = 0.20, CAF = 0.60),
    tumors = list(
      list(center = c(17, 23), a = 8.8, b = 6.3, theta = pi / 7, core_frac = 0.60),
      list(center = c(33, 22), a = 8.6, b = 6.1, theta = -pi / 9, core_frac = 0.60),
      list(center = c(25, 24), a = 5.5, b = 3.0, theta = 0, core_frac = 0.58)
    ),
    immune = list(
      list(center = c(13, 34), a = 6.0, b = 4.8, theta = pi / 12, core_frac = 0.62),
      list(center = c(38, 34), a = 5.7, b = 4.0, theta = -pi / 5, core_frac = 0.62)
    ),
    vessels = list(
      list(p1 = c(5, 7), p2 = c(46, 42), core_width = 1.0, boundary_width = 2.1),
      list(p1 = c(6, 29), p2 = c(46, 16), core_width = 0.9, boundary_width = 1.9)
    ),
    caf_islands = list(
      list(center = c(8, 42), a = 4.2, b = 3.0, theta = pi / 8, core_frac = 0.65),
      list(center = c(43, 8), a = 5.0, b = 3.2, theta = -pi / 7, core_frac = 0.65)
    ),
    cell_per_spot = 9,
    p0 = 0.05,
    seq_depth = 5.8e5
  ),
  list(
    id = 3,
    seed = 7103,
    label = "upper_dome_tumor",
    background = c(Cancer = 0.07, Tcell = 0.13, Endothelial = 0.18, CAF = 0.62),
    tumors = list(
      list(center = c(24, 34), a = 11.0, b = 7.8, theta = -pi / 14, core_frac = 0.62),
      list(center = c(30, 26), a = 6.0, b = 4.8, theta = pi / 5, core_frac = 0.60)
    ),
    immune = list(
      list(center = c(18, 18), a = 7.5, b = 4.8, theta = pi / 10, core_frac = 0.62),
      list(center = c(33, 41), a = 4.8, b = 3.8, theta = -pi / 4, core_frac = 0.60)
    ),
    vessels = list(
      list(p1 = c(3, 37), p2 = c(47, 25), core_width = 0.95, boundary_width = 1.9),
      list(p1 = c(8, 10), p2 = c(40, 44), core_width = 0.95, boundary_width = 1.8)
    ),
    caf_islands = list(
      list(center = c(11, 7), a = 5.2, b = 3.4, theta = pi / 7, core_frac = 0.65),
      list(center = c(39, 9), a = 4.9, b = 3.0, theta = -pi / 5, core_frac = 0.65)
    ),
    cell_per_spot = 10,
    p0 = 0.05,
    seq_depth = 6.0e5
  ),
  list(
    id = 4,
    seed = 7104,
    label = "split_tumor_vessel_axis",
    background = c(Cancer = 0.06, Tcell = 0.14, Endothelial = 0.22, CAF = 0.58),
    tumors = list(
      list(center = c(20, 24), a = 8.2, b = 6.0, theta = pi / 9, core_frac = 0.60),
      list(center = c(31, 26), a = 7.8, b = 5.8, theta = -pi / 8, core_frac = 0.60),
      list(center = c(25, 16), a = 5.0, b = 3.3, theta = 0, core_frac = 0.58)
    ),
    immune = list(
      list(center = c(14, 34), a = 6.8, b = 4.5, theta = pi / 8, core_frac = 0.62),
      list(center = c(36, 34), a = 6.2, b = 4.2, theta = -pi / 6, core_frac = 0.62)
    ),
    vessels = list(
      list(p1 = c(24, 3), p2 = c(25, 47), core_width = 1.0, boundary_width = 2.2),
      list(p1 = c(6, 20), p2 = c(44, 29), core_width = 0.9, boundary_width = 1.9)
    ),
    caf_islands = list(
      list(center = c(7, 12), a = 4.5, b = 3.3, theta = pi / 10, core_frac = 0.65),
      list(center = c(43, 40), a = 5.4, b = 3.2, theta = -pi / 9, core_frac = 0.65)
    ),
    cell_per_spot = 10,
    p0 = 0.05,
    seq_depth = 6.2e5
  ),
  list(
    id = 5,
    seed = 7105,
    label = "multifocal_tumor",
    background = c(Cancer = 0.06, Tcell = 0.16, Endothelial = 0.18, CAF = 0.60),
    tumors = list(
      list(center = c(16, 20), a = 6.2, b = 4.8, theta = pi / 7, core_frac = 0.60),
      list(center = c(25, 30), a = 7.0, b = 5.6, theta = -pi / 8, core_frac = 0.60),
      list(center = c(35, 20), a = 6.5, b = 5.0, theta = pi / 9, core_frac = 0.60),
      list(center = c(28, 12), a = 4.8, b = 3.2, theta = pi / 4, core_frac = 0.58)
    ),
    immune = list(
      list(center = c(11, 33), a = 5.8, b = 4.4, theta = pi / 10, core_frac = 0.62),
      list(center = c(39, 33), a = 5.5, b = 4.2, theta = -pi / 5, core_frac = 0.62)
    ),
    vessels = list(
      list(p1 = c(4, 8), p2 = c(46, 38), core_width = 0.95, boundary_width = 1.9),
      list(p1 = c(7, 42), p2 = c(44, 15), core_width = 0.95, boundary_width = 2.0)
    ),
    caf_islands = list(
      list(center = c(8, 42), a = 5.1, b = 3.4, theta = pi / 6, core_frac = 0.65),
      list(center = c(43, 8), a = 5.1, b = 3.2, theta = -pi / 8, core_frac = 0.65)
    ),
    cell_per_spot = 10,
    p0 = 0.05,
    seq_depth = 6.3e5
  )
)

mix_tumor_core <- c(Cancer = 1, Tcell = 0, Endothelial = 0, CAF = 0)
mix_tumor_boundary <- c(Cancer = 0.62, Tcell = 0.20, Endothelial = 0.10, CAF = 0.08)
mix_immune_core <- c(Cancer = 0, Tcell = 1, Endothelial = 0, CAF = 0)
mix_immune_boundary <- c(Cancer = 0.22, Tcell = 0.56, Endothelial = 0.14, CAF = 0.08)
mix_vessel_core <- c(Cancer = 0.10, Tcell = 0.10, Endothelial = 0.80, CAF = 0.00)
mix_vessel_boundary <- c(Cancer = 0.24, Tcell = 0.22, Endothelial = 0.40, CAF = 0.14)
mix_caf_core <- c(Cancer = 0.05, Tcell = 0.05, Endothelial = 0.10, CAF = 0.80)
mix_caf_boundary <- c(Cancer = 0.20, Tcell = 0.12, Endothelial = 0.14, CAF = 0.54)

marker_registry <- list()
coexp_registry <- list()

# Simulation driver ---------------------------------------------------------
output_dir <- "./Data/S7_gene/simulated"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
replicate_designs <- replicate_designs[vapply(
  replicate_designs, function(x) x$id == rep_id, logical(1)
)]
if (length(replicate_designs) != 1L) {
  stop("Exactly one S7 replicate design must be selected.")
}

for (cfg in replicate_designs) {
  set.seed(cfg$seed)
  coord <- generate_coordinates(
    x_len = 50,
    y_len = 50,
    pattern = "hex",
    spot_distance = spot_distance,
    seed = cfg$seed
  )
  n_spots <- nrow(coord)

  base_mix <- normalize_mix(cfg$background, ct_names)
  prop_mat <- matrix(rep(base_mix, each = n_spots), nrow = n_spots)
  colnames(prop_mat) <- ct_names
  layer_vec <- rep(3L, n_spots)
  region_vec <- rep("background", n_spots)
  priority_vec <- rep(1L, n_spots)

  for (idx in seq_along(cfg$tumors)) {
    tm <- cfg$tumors[[idx]]
    masks <- ellipse_masks(coord, tm$center, tm$a, tm$b, tm$theta, tm$core_frac)
    ans <- apply_assignment(
      prop_mat, layer_vec, region_vec, priority_vec,
      masks$core, normalize_mix(mix_tumor_core, ct_names),
      layer_value = 1L, region_label = paste0("tumor_core_", idx), priority = 100L
    )
    prop_mat <- ans$prop_mat; layer_vec <- ans$layer_vec
    region_vec <- ans$region_vec; priority_vec <- ans$priority_vec
    ans <- apply_assignment(
      prop_mat, layer_vec, region_vec, priority_vec,
      masks$boundary, normalize_mix(mix_tumor_boundary, ct_names),
      layer_value = 2L, region_label = paste0("tumor_boundary_", idx), priority = 80L
    )
    prop_mat <- ans$prop_mat; layer_vec <- ans$layer_vec
    region_vec <- ans$region_vec; priority_vec <- ans$priority_vec
  }

  for (idx in seq_along(cfg$immune)) {
    im <- cfg$immune[[idx]]
    masks <- ellipse_masks(coord, im$center, im$a, im$b, im$theta, im$core_frac)
    ans <- apply_assignment(
      prop_mat, layer_vec, region_vec, priority_vec,
      masks$core, normalize_mix(mix_immune_core, ct_names),
      layer_value = 2L, region_label = paste0("immune_core_", idx), priority = 70L
    )
    prop_mat <- ans$prop_mat; layer_vec <- ans$layer_vec
    region_vec <- ans$region_vec; priority_vec <- ans$priority_vec
    ans <- apply_assignment(
      prop_mat, layer_vec, region_vec, priority_vec,
      masks$boundary, normalize_mix(mix_immune_boundary, ct_names),
      layer_value = 2L, region_label = paste0("immune_boundary_", idx), priority = 60L
    )
    prop_mat <- ans$prop_mat; layer_vec <- ans$layer_vec
    region_vec <- ans$region_vec; priority_vec <- ans$priority_vec
  }

  for (idx in seq_along(cfg$vessels)) {
    vs <- cfg$vessels[[idx]]
    masks <- stripe_masks(coord, vs$p1, vs$p2, vs$core_width, vs$boundary_width)
    ans <- apply_assignment(
      prop_mat, layer_vec, region_vec, priority_vec,
      masks$core, normalize_mix(mix_vessel_core, ct_names),
      layer_value = 2L, region_label = paste0("vessel_core_", idx), priority = 65L
    )
    prop_mat <- ans$prop_mat; layer_vec <- ans$layer_vec
    region_vec <- ans$region_vec; priority_vec <- ans$priority_vec
    ans <- apply_assignment(
      prop_mat, layer_vec, region_vec, priority_vec,
      masks$boundary, normalize_mix(mix_vessel_boundary, ct_names),
      layer_value = 2L, region_label = paste0("vessel_boundary_", idx), priority = 55L
    )
    prop_mat <- ans$prop_mat; layer_vec <- ans$layer_vec
    region_vec <- ans$region_vec; priority_vec <- ans$priority_vec
  }

  for (idx in seq_along(cfg$caf_islands)) {
    cf <- cfg$caf_islands[[idx]]
    masks <- ellipse_masks(coord, cf$center, cf$a, cf$b, cf$theta, cf$core_frac)
    ans <- apply_assignment(
      prop_mat, layer_vec, region_vec, priority_vec,
      masks$core, normalize_mix(mix_caf_core, ct_names),
      layer_value = 3L, region_label = paste0("stroma_core_", idx), priority = 45L
    )
    prop_mat <- ans$prop_mat; layer_vec <- ans$layer_vec
    region_vec <- ans$region_vec; priority_vec <- ans$priority_vec
    ans <- apply_assignment(
      prop_mat, layer_vec, region_vec, priority_vec,
      masks$boundary, normalize_mix(mix_caf_boundary, ct_names),
      layer_value = 3L, region_label = paste0("stroma_boundary_", idx), priority = 35L
    )
    prop_mat <- ans$prop_mat; layer_vec <- ans$layer_vec
    region_vec <- ans$region_vec; priority_vec <- ans$priority_vec
  }

  prop_mat[prop_mat < 0] <- 0
  row_tot <- rowSums(prop_mat)
  zero_rows <- which(row_tot <= 0)
  if (length(zero_rows) > 0) {
    prop_mat[zero_rows, ] <- matrix(rep(normalize_mix(c(CAF = 1), ct_names), each = length(zero_rows)),
                                    nrow = length(zero_rows), byrow = FALSE)
    row_tot <- rowSums(prop_mat)
  }
  prop_mat <- prop_mat / row_tot

  dominant <- ct_names[max.col(prop_mat, ties.method = "first")]
  layer_vec[layer_vec == 3L & dominant == "Cancer" & prop_mat[, "Cancer"] > 0.75] <- 1L
  layer_vec[layer_vec == 3L &
              (prop_mat[, "Tcell"] + prop_mat[, "Endothelial"]) > 0.40] <- 2L

  mean_param <- make_mean_param(
    coord = coord,
    prop_mat = prop_mat,
    ct_names = ct_names,
    n_genes = nGenes,
    svg_gene_sets = svg_gene_sets,
    marker_gene_sets = marker_gene_sets,
    rep_id = cfg$id,
    seed = cfg$seed
  )

  gene_relationship_s <- copy(base_gene_relationship)
  cov_gene <- generate_covariance_group(gene_relationship_s)
  true_count <- generate_true_count(
    coord,
    mean_param,
    gene_relationship = gene_relationship_s,
    cov_gene = cov_gene
  )
  true_count <- generate_batch_count(
    true_count, p0 = cfg$p0, celltype_proportion = prop_mat,
    seq_depth = cfg$seq_depth
  )
  spot_count <- aggregate_spot_counts_abundance_weighted(
    count_list = true_count,
    cell_numbers = rep(as.integer(cfg$cell_per_spot), n_spots),
    celltype_prop = prop_mat
  )

  coord_dt <- as.data.table(coord)
  expr_dt <- as.data.table(spot_count)
  region_dt <- cbind(
    coord_dt,
    data.table(region = region_vec),
    as.data.table(prop_mat)
  )
  layer_dt <- data.table(
    x_coord = coord$x_coord,
    y_coord = coord$y_coord,
    layer = layer_vec
  )

  sc_ref <- rbindlist(lapply(names(true_count), function(ct) {
    mat <- as.matrix(true_count[[ct]])
    if (nrow(mat) != n_spots && ncol(mat) == n_spots) {
      mat <- t(mat)
    }
    dt <- as.data.table(mat)
    dt[, Cell_type := ct]
    dt
  }), fill = TRUE)
  setcolorder(sc_ref, c("Cell_type", setdiff(colnames(sc_ref), "Cell_type")))
  if ((ncol(sc_ref) - 1L) == length(gene_ids)) {
    colnames(sc_ref)[-1] <- gene_ids
  } else {
    colnames(sc_ref)[-1] <- paste0("gene_", seq_len(ncol(sc_ref) - 1L))
  }

  fwrite(coord_dt, file = file.path(output_dir, sprintf("S7_coord_%d.csv", cfg$id)))
  fwrite(expr_dt, file = file.path(output_dir, sprintf("S7_gene_expression_%d.csv", cfg$id)))
  fwrite(region_dt, file = file.path(output_dir, sprintf("S7_region_%d.csv", cfg$id)))
  fwrite(layer_dt, file = file.path(output_dir, sprintf("S7_layer_%d.csv", cfg$id)))
  fwrite(sc_ref, file = file.path(output_dir, sprintf("S7_scref_%d.csv", cfg$id)))

  marker_registry[[as.character(cfg$id)]] <- rbindlist(lapply(names(svg_gene_sets), function(ct) {
    idx <- svg_gene_sets[[ct]]
    data.table(
      scenario = cfg$id,
      cell_type = ct,
      gene_index = idx,
      gene_id = gene_ids[idx]
    )
  }))

  coexp_ct <- rbindlist(lapply(ct_names, function(ct) {
    svg_idx <- svg_gene_sets[[ct]]
    marker_idx <- marker_gene_sets[[ct]]
    module_base <- (match(ct, ct_names) - 1L) * 10L
    rbindlist(list(
      data.table(
        scenario = cfg$id,
        cell_type = ct,
        module_id = module_base + 1L,
        module_label = paste0(ct, "_svg_core"),
        gene_index = svg_idx[1:5],
        gene_id = gene_ids[svg_idx[1:5]]
      ),
      data.table(
        scenario = cfg$id,
        cell_type = ct,
        module_id = module_base + 2L,
        module_label = paste0(ct, "_svg_boundary"),
        gene_index = svg_idx[6:10],
        gene_id = gene_ids[svg_idx[6:10]]
      ),
      data.table(
        scenario = cfg$id,
        cell_type = ct,
        module_id = module_base + 3L,
        module_label = paste0(ct, "_marker_a"),
        gene_index = marker_idx[1:5],
        gene_id = gene_ids[marker_idx[1:5]]
      ),
      data.table(
        scenario = cfg$id,
        cell_type = ct,
        module_id = module_base + 4L,
        module_label = paste0(ct, "_marker_b"),
        gene_index = marker_idx[6:10],
        gene_id = gene_ids[marker_idx[6:10]]
      )
    ))
  }))

  coexp_bg <- rbindlist(list(
    data.table(
      scenario = cfg$id,
      cell_type = "Background",
      module_id = 901L,
      module_label = "background_shared_a",
      gene_index = background_gene_sets$shared_a,
      gene_id = gene_ids[background_gene_sets$shared_a]
    ),
    data.table(
      scenario = cfg$id,
      cell_type = "Background",
      module_id = 902L,
      module_label = "background_shared_b",
      gene_index = background_gene_sets$shared_b,
      gene_id = gene_ids[background_gene_sets$shared_b]
    )
  ))
  coexp_registry[[as.character(cfg$id)]] <- rbindlist(list(coexp_ct, coexp_bg), fill = TRUE)

  message(sprintf(
    "Replicate %d (%s) completed with %d spots.",
    cfg$id, cfg$label, n_spots
  ))
}

if (length(marker_registry) > 0) {
  marker_summary <- rbindlist(marker_registry, fill = TRUE)
  fwrite(marker_summary, file = file.path(output_dir, "S7_marker_gene_truth.csv"))
}

if (length(coexp_registry) > 0) {
  coexp_summary <- rbindlist(coexp_registry, fill = TRUE)
  fwrite(coexp_summary, file = file.path(output_dir, "S7_coexpression_truth.csv"))
}
