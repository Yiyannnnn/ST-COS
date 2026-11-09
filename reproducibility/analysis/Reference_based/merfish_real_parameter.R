#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Estimate MERFISH parameters and simulate with ST-COS
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(mgcv)
  library(spdep)
  library(parallel)
  library(ggplot2)
  library(tidyr)
})
## Optional simulators
suppressPackageStartupMessages({
  library(SRTsim)
  library(scDesign3)
  library(SingleCellExperiment)
})

source("./Data/simulator_YZ.R")

count_path <- "./data_external/MERFISH/count/X0_04_count.csv"
mapping_path <- "./data_external/MERFISH/mapping/X0_04_mapping.csv"

output_prefix <- "./Reference_based/MERFISH_X0_04"
block_size <- 1024
n_cores <- 1
spatial_check_n <- 300
spatial_moran_thr <- 0.05
knn_k <- 10
gene_subset_n <- NULL

message("[INFO] Loading MERFISH counts ...")
counts_dt <- fread(count_path)
spot_ids <- counts_dt[[1]]
gene_names <- colnames(counts_dt)[-1]
expr_mat_spot_gene <- as.matrix(counts_dt[, -1, with = FALSE])
expr_mat <- t(expr_mat_spot_gene)
rownames(expr_mat) <- gene_names
colnames(expr_mat) <- spot_ids

message("[INFO] Loading coordinates ...")
mapping_dt <- fread(mapping_path)
coords_mat <- as.matrix(mapping_dt[match(spot_ids, mapping_dt$row_name), c("x", "y"), with = FALSE])
if (any(is.na(coords_mat))) {
  stop("Missing coordinates for some spots after matching mapping table.")
}
rownames(coords_mat) <- spot_ids

fit_nb_fast <- function(y) {
  mu <- mean(y)
  s2 <- var(y)
  th <- if (s2 > mu) mu^2 / (s2 - mu) else 1e6
  c(mu = mu, theta = th)
}

message("[INFO] Estimating NB moments ...")
nb_par <- t(apply(expr_mat, 1, fit_nb_fast))
mu_g <- nb_par[, "mu"]
theta_g <- nb_par[, "theta"]

counts2z_block <- function(idx, Y, mu, th, rng_block, eps = 1e-6) {
  out <- matrix(0, nrow = length(idx), ncol = ncol(Y))
  for (ii in seq_along(idx)) {
    g <- idx[ii]
    y <- Y[g, ]
    u0 <- pnbinom(pmax(y - 1, 0), size = th[g], mu = mu[g])
    u1 <- pnbinom(y, size = th[g], mu = mu[g])
    u <- u0 + rng_block[ii, ] * (u1 - u0)
    u <- pmin(pmax(u, eps), 1 - eps)
    out[ii, ] <- qnorm(u)
  }
  out
}

counts2z_mat <- function(Y, mu, th, block = 1024, n_cores = 1, eps = 1e-6) {
  g_n <- nrow(Y)
  idx_ls <- split(seq_len(g_n), ceiling(seq_len(g_n) / block))
  rng <- matrix(runif(g_n * ncol(Y)), nrow = g_n)
  worker <- function(idx) {
    counts2z_block(idx, Y, mu, th, rng[idx, , drop = FALSE], eps)
  }
  if (n_cores > 1) {
    z_parts <- mclapply(idx_ls, worker, mc.cores = n_cores)
  } else {
    z_parts <- lapply(idx_ls, worker)
  }
  do.call(rbind, z_parts)
}

message("[INFO] Building correlation matrix ...")
Z <- t(counts2z_mat(expr_mat, mu_g, theta_g, block = block_size, n_cores = n_cores))
Sigma <- cor(Z, use = "pairwise.complete.obs")
Sigma[is.na(Sigma)] <- 0
diag(Sigma) <- 1
#Sigma[abs(Sigma) < 0.05] <- 0
Sigma <- makespd(Sigma)
if (min(eigen(Sigma, only.values = TRUE)$values) <= 0) {
  Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = TRUE)$mat)
}
rownames(Sigma) <- gene_names
colnames(Sigma) <- gene_names

message("[INFO] Computing library sizes ...")
lib_size <- colSums(expr_mat)
lib_factor <- lib_size / mean(lib_size)

message("[INFO] Detecting spatial genes (subset) ...")
Wnb <- knearneigh(coords_mat, k = knn_k)
W <- nb2listw(knn2nb(Wnb))

moran_flag <- function(y, listw, thr = 0.15) {
  S0 <- Szero(listw)
  res <- moran(x = y, listw = listw, n = length(y), S0 = S0, zero.policy = TRUE)$I
  res > thr
}

spatial_candidates <- gene_names[seq_len(min(spatial_check_n, length(gene_names)))]
spatial_genes <- spatial_candidates[sapply(spatial_candidates, function(g) {
  moran_flag(expr_mat[g, ], W, thr = spatial_moran_thr)
})]
message(sprintf("[INFO] %d spatial genes detected (subset)", length(spatial_genes)))

message("[INFO] Building mean matrix ...")
mu_mat <- outer(lib_factor, mu_g)
dimnames(mu_mat) <- list(colnames(expr_mat), rownames(expr_mat))

if (length(spatial_genes) > 0) {
  message("[INFO] Refitting NB-GAM for spatial genes ...")
  for (g in spatial_genes) {
    dat <- data.frame(
      y = as.numeric(expr_mat[g, ]),
      row = coords_mat[, 1],
      col = coords_mat[, 2],
      log_lib = log(lib_factor)
    )
    gam_fit <- gam(
      y ~ offset(log_lib) + s(row, col, bs = "tp", k = 50),
      family = nb(link = "log"),
      data = dat,
      method = "REML"
    )
    mu_mat[, g] <- fitted(gam_fit)
    theta_g[g] <- gam_fit$family$getTheta(TRUE)
  }
}

message("[INFO] Preparing ST-COS inputs ...")
mean_param <- list(CT = mu_mat)
gene_relationship <- as.data.table(
  matrix(1, nrow = nrow(expr_mat), ncol = 1,
         dimnames = list(rownames(expr_mat), "CT"))
)
cov_gene <- array(NA, dim = c(nrow(expr_mat), nrow(expr_mat), 1),
                  dimnames = list(rownames(expr_mat), rownames(expr_mat), "CT"))
cov_gene[, , "CT"] <- array(Sigma,
                            dim = c(nrow(expr_mat), nrow(expr_mat), 1),
                            dimnames = list(rownames(expr_mat), rownames(expr_mat), "CT"))

message("[INFO] Running ST-COS simulation ...")
sim_list <- generate_true_count(
  coord = coords_mat,
  mean_param = mean_param,
  gene_relationship = gene_relationship,
  cov_gene = cov_gene,
  lib_size = list(CT = lib_size)
)

sim_spot_counts <- sim_list[["CT"]]
rownames(sim_spot_counts) <- spot_ids
colnames(sim_spot_counts) <- gene_names

message("[INFO] Writing outputs ...")
counts_out <- paste0(output_prefix, "_simulated_counts.csv")
coords_out <- paste0(output_prefix, "_simulated_coordinates.csv")
write.csv(sim_spot_counts, file = counts_out, row.names = TRUE)
write.csv(
  data.frame(s1 = coords_mat[, 1], s2 = coords_mat[, 2], row.names = spot_ids),
  file = coords_out,
  row.names = TRUE
)

## ---- SRTsim -------------------------------------------------------------
message("[INFO] Running SRTsim simulation ...")
gene_use <- gene_names
if (!is.null(gene_subset_n)) {
  gene_use <- gene_names[seq_len(min(gene_subset_n, length(gene_names)))]
}
srt_counts <- expr_mat[gene_use, ]
srt_counts <- as.matrix(srt_counts)
srt_loc <- data.frame(
  y = coords_mat[, 1],
  x = coords_mat[, 2],
  label = rep("CT", nrow(coords_mat)),
  row.names = colnames(srt_counts)
)
simSRT <- createSRT(count_in = srt_counts, loc_in = srt_loc)
simSRT <- srtsim_fit(simSRT, sim_schem = "tissue")
simSRT <- srtsim_count(simSRT)
srt_sim_counts <- simCounts(simSRT)
srt_out <- paste0(output_prefix, "_SRTsim_counts.csv")
write.csv(as.matrix(srt_sim_counts), file = srt_out, row.names = TRUE)

## ---- scDesign3 ----------------------------------------------------------
message("[INFO] Running scDesign3 simulation ...")
sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(srt_counts)),
  colData = DataFrame(
    spatial1 = coords_mat[, 1],
    spatial2 = coords_mat[, 2],
    cell_type = rep("CT", nrow(coords_mat))
  )
)
scd3_fit <- scdesign3(
  sce = sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 200)",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = n_cores,
  usebam = FALSE,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = FALSE,
  parallelization = "pbmcapply"
)
scd3_counts <- t(scd3_fit$new_count)
scd3_out <- paste0(output_prefix, "_scDesign3_counts.csv")
write.csv(scd3_counts, file = scd3_out, row.names = TRUE)

message("[INFO] Done.")

## ---- Optional scCube results -------------------------------------------
scCube_expr <- NULL
scCube_path <- "./Reference_based/scCube_MERFISH_04.csv"
scCube_meta_path <- "./Reference_based/scCube_MERFISH_04_meta.csv"
scCube_train_expr <- NULL
scCube_train_path <- "./Reference_based/scCube_train_MERFISH_X0_04/MERFISH_X0_04_sim_data.csv.gz"
scCube_train_meta_path <- "./Reference_based/scCube_train_MERFISH_X0_04/MERFISH_X0_04_sim_meta.csv.gz"

map_scube_rows <- function(sc_expr, coords, spot_ids, meta_dt) {
  meta_dt$key_xy <- paste0(round(meta_dt$x), "x", round(meta_dt$y))
  row_ids <- rownames(sc_expr)
  spot_key_xy <- paste0(round(coords[, 1]), "x", round(coords[, 2]))
  spot_key_yx <- paste0(round(coords[, 2]), "x", round(coords[, 1]))
  spot_map_xy <- setNames(spot_ids, spot_key_xy)
  spot_map_yx <- setNames(spot_ids, spot_key_yx)
  mapped_xy <- spot_map_xy[row_ids]
  mapped_yx <- spot_map_yx[row_ids]
  n_xy <- sum(!is.na(mapped_xy))
  n_yx <- sum(!is.na(mapped_yx))
  mapped <- if (n_yx > n_xy) mapped_yx else mapped_xy
  rownames(sc_expr) <- ifelse(is.na(mapped), rownames(sc_expr), mapped)
  sc_expr
}

map_scube_cols <- function(sc_expr, coords, spot_ids, meta_dt) {
  meta_dt$key_xy <- paste0(round(meta_dt$x), "x", round(meta_dt$y))
  col_ids <- colnames(sc_expr)
  spot_key_xy <- paste0(round(coords[, 1]), "x", round(coords[, 2]))
  spot_key_yx <- paste0(round(coords[, 2]), "x", round(coords[, 1]))
  spot_map_xy <- setNames(spot_ids, spot_key_xy)
  spot_map_yx <- setNames(spot_ids, spot_key_yx)
  cell_keys <- meta_dt$key_xy[match(col_ids, meta_dt$Cell)]
  mapped_xy <- spot_map_xy[cell_keys]
  mapped_yx <- spot_map_yx[cell_keys]
  n_xy <- sum(!is.na(mapped_xy))
  n_yx <- sum(!is.na(mapped_yx))
  mapped <- if (n_yx > n_xy) mapped_yx else mapped_xy
  colnames(sc_expr) <- ifelse(is.na(mapped), colnames(sc_expr), mapped)
  sc_expr
}

orient_by_names <- function(mat, gene_ids, spot_ids) {
  gene_in_rows <- length(intersect(rownames(mat), gene_ids))
  gene_in_cols <- length(intersect(colnames(mat), gene_ids))
  spot_in_rows <- length(intersect(rownames(mat), spot_ids))
  spot_in_cols <- length(intersect(colnames(mat), spot_ids))
  if (gene_in_cols > gene_in_rows && spot_in_rows > spot_in_cols) {
    mat <- t(mat)
  }
  mat
}

if (file.exists(scCube_path)) {
  scCube_expr <- read.csv(scCube_path, row.names = 1, check.names = FALSE)
  scCube_expr <- as.matrix(scCube_expr)
  if (file.exists(scCube_meta_path)) {
    scCube_meta <- read.csv(scCube_meta_path, check.names = FALSE)
    if (!("key_xy" %in% colnames(scCube_meta))) {
      scCube_meta$key_xy <- paste0(round(scCube_meta$x), "x", round(scCube_meta$y))
    }
    if (length(intersect(rownames(scCube_expr), scCube_meta$key_xy)) > 0) {
      scCube_expr <- map_scube_rows(scCube_expr, coords_mat, spot_ids, scCube_meta)
    } else if (length(intersect(rownames(scCube_expr), scCube_meta$Cell)) > 0) {
      spot_key <- paste0(round(coords_mat[, 1]), "x", round(coords_mat[, 2]))
      spot_key_map <- setNames(spot_ids, spot_key)
      cell_keys <- scCube_meta$key_xy[match(rownames(scCube_expr), scCube_meta$Cell)]
      mapped <- spot_key_map[cell_keys]
      rownames(scCube_expr) <- ifelse(is.na(mapped), rownames(scCube_expr), mapped)
    }
  } else {
    message(sprintf("[WARN] scCube meta not found at %s; skipping coordinate mapping.", scCube_meta_path))
  }
  scCube_expr <- orient_by_names(scCube_expr, gene_names, spot_ids)
} else {
  message(sprintf("[WARN] scCube results not found at %s; skipping scCube plots.", scCube_path))
}

if (file.exists(scCube_train_path)) {
  scCube_train_dt <- fread(cmd = paste("zcat", shQuote(scCube_train_path)))
  rownames(scCube_train_dt) <- scCube_train_dt[[1]]
  scCube_train_expr <- as.matrix(scCube_train_dt[, -1, with = FALSE])
  if (file.exists(scCube_train_meta_path)) {
    scCube_train_meta <- fread(cmd = paste("zcat", shQuote(scCube_train_meta_path)))
    if (length(intersect(colnames(scCube_train_expr), scCube_train_meta$Cell)) > 0) {
      scCube_train_expr <- map_scube_cols(scCube_train_expr, coords_mat, spot_ids, scCube_train_meta)
    } else if (length(intersect(rownames(scCube_train_expr), scCube_train_meta$Cell)) > 0) {
      scCube_train_expr <- map_scube_rows(scCube_train_expr, coords_mat, spot_ids, scCube_train_meta)
    }
  } else {
    message(sprintf("[WARN] scCube self-trained meta not found at %s; skipping coordinate mapping.", scCube_train_meta_path))
  }
  scCube_train_expr <- orient_by_names(scCube_train_expr, gene_names, spot_ids)
} else {
  message(sprintf("[WARN] scCube self-trained results not found at %s; skipping scCube_train plots.", scCube_train_path))
}

## ---- Optional other simulated results -----------------------------------
stcos_path <- paste0(output_prefix, "_simulated_counts.csv")
srt_path <- paste0(output_prefix, "_SRTsim_counts.csv")
scd3_path <- paste0(output_prefix, "_scDesign3_counts.csv")

if (file.exists(stcos_path)) {
  sim_spot_counts <- read.csv(stcos_path, row.names = 1, check.names = FALSE)
  sim_spot_counts <- as.matrix(sim_spot_counts)
}
if (file.exists(srt_path)) {
  srt_sim_counts <- read.csv(srt_path, row.names = 1, check.names = FALSE)
  srt_sim_counts <- as.matrix(srt_sim_counts)
}
if (file.exists(scd3_path)) {
  scd3_counts <- read.csv(scd3_path, row.names = 1, check.names = FALSE)
  scd3_counts <- as.matrix(scd3_counts)
}

## ---- Harmonize genes across methods -------------------------------------
subset_genes_matrix <- function(mat, genes) {
  if (length(intersect(rownames(mat), genes)) > 0) {
    return(mat[intersect(rownames(mat), genes), , drop = FALSE])
  }
  if (length(intersect(colnames(mat), genes)) > 0) {
    return(mat[, intersect(colnames(mat), genes), drop = FALSE])
  }
  mat
}

common_genes <- Reduce(intersect, list(
  rownames(expr_mat),
  colnames(sim_spot_counts),
  rownames(srt_sim_counts),
  colnames(scd3_counts),
  if (!is.null(scCube_expr)) rownames(scCube_expr)[-which(apply(scCube_expr, 1, function(x) all(x == 0)))] else rownames(expr_mat),
  if (!is.null(scCube_train_expr)) rownames(scCube_train_expr)[-which(apply(scCube_train_expr, 1, function(x) all(x == 0)))] else rownames(expr_mat)
))
gene_names <- common_genes
expr_mat <- expr_mat[gene_names, , drop = FALSE]
sim_spot_counts <- subset_genes_matrix(sim_spot_counts, gene_names)
srt_sim_counts <- subset_genes_matrix(srt_sim_counts, gene_names)
scd3_counts <- subset_genes_matrix(scd3_counts, gene_names)
Sigma_sub <- Sigma[gene_names, gene_names, drop = FALSE]
if (!is.null(scCube_expr)) {
  scCube_expr <- subset_genes_matrix(scCube_expr, gene_names)
}
if (!is.null(scCube_train_expr)) {
  scCube_train_expr <- subset_genes_matrix(scCube_train_expr, gene_names)
}

library(corrplot)
cor_pal <- colorRampPalette(c("#3B4CC0", "#BFD3E6", "#FFFFFF", "#F4B6C2", "#B40426"))(200)

plot_corr_dual <- function(base_mat, method_mat, title, ord = NULL) {
  if (!is.null(ord)) {
    base_mat <- base_mat[ord, ord, drop = FALSE]
    method_mat <- method_mat[ord, ord, drop = FALSE]
  }
  corrplot::corrplot(
    base_mat,
    method = "color",
    col = cor_pal,
    type = "lower",
    order = "original",
    addrect = 4,
    tl.pos = "n",
    cl.pos = "n",
    diag = FALSE,
    mar = c(0, 0, 2, 0),
    title = title
  )
  corrplot::corrplot(
    method_mat,
    method = "color",
    col = cor_pal,
    type = "upper",
    order = "original",
    add = TRUE,
    tl.pos = "n",
    cl.cex = 0.7,
    diag = FALSE
  )
}

base_corr <- cor(t(expr_mat))
base_ord <- corrplot::corrMatOrder(base_corr, order = "hclust")

pdf(paste0(output_prefix, "simulated_cor_comparison.pdf"), width = 8, height = 8)
par(mfrow = c(2, 2))
plot_corr_dual(base_corr, Sigma_sub, "ST-COS vs True", ord = base_ord)
#plot_corr_dual(base_corr, Sigma, "Estimated vs True (Sigma)", ord = base_ord)
#plot_corr_dual(base_corr, cor(as.matrix(sim_spot_counts)), "ST-COS vs True", ord = base_ord)
plot_corr_dual(base_corr, cor(as.matrix(t(srt_sim_counts))), "SRTsim vs True", ord = base_ord)
plot_corr_dual(base_corr, cor(as.matrix(scd3_counts)), "scDesign3 vs True", ord = base_ord)
if (!is.null(scCube_expr)) {
  plot_corr_dual(base_corr, cor(as.matrix(scCube_expr)), "scCube vs True", ord = base_ord)
}
if (!is.null(scCube_train_expr)) {
  plot_corr_dual(base_corr, cor(as.matrix(scCube_train_expr)), "scCube_train vs True", ord = base_ord)
}
dev.off()

## ---- Spatial expression plots -------------------------------------------
message("[INFO] Plotting spatial gene expression ...")
library(viridis)

message("[INFO] Computing Moran's I for all genes ...")
all_moran <- sapply(gene_names, function(g) {
  moran(x = expr_mat[g, ], listw = W, n = length(spot_ids), S0 = Szero(W), zero.policy = TRUE)$I
})
plot_genes <- names(sort(all_moran, decreasing = TRUE))[seq_len(min(50, length(all_moran)))]
plot_genes <- plot_genes[seq_len(min(160, length(plot_genes)))]

as_gene_by_spot <- function(mat, gene_ids, spot_ids) {
  mat <- as.matrix(mat)
  if (length(intersect(gene_ids, colnames(mat))) > 0 &&
      length(intersect(spot_ids, rownames(mat))) > 0) {
    mat <- t(mat)
  }
  gene_keep <- intersect(gene_ids, rownames(mat))
  spot_keep <- intersect(spot_ids, colnames(mat))
  if (length(gene_keep) == 0 || length(spot_keep) == 0) {
    stop("No overlapping genes or spots found for expression plotting.")
  }
  mat <- mat[gene_keep, spot_keep, drop = FALSE]
  mat
}

coords_df <- data.frame(
  spot = spot_ids,
  x = coords_mat[, 1],
  y = coords_mat[, 2]
)

expr_true <- as_gene_by_spot(expr_mat, plot_genes, spot_ids)
expr_stcos <- as_gene_by_spot(sim_spot_counts, plot_genes, spot_ids)
expr_srt <- as_gene_by_spot(srt_sim_counts, plot_genes, spot_ids)
expr_scd3 <- as_gene_by_spot(scd3_counts, plot_genes, spot_ids)
expr_scube <- NULL
if (!is.null(scCube_expr)) {
  expr_scube <- as_gene_by_spot(scCube_expr, plot_genes, spot_ids)
}
expr_scube_train <- NULL
if (!is.null(scCube_train_expr)) {
  expr_scube_train <- as_gene_by_spot(scCube_train_expr, plot_genes, spot_ids)
}

make_long <- function(expr, label) {
  dt <- as.data.table(t(expr))
  dt[, spot := spot_ids]
  dt <- as.data.table(merge(dt, coords_df, by = "spot", all.x = TRUE))
  dt <- melt(
    dt,
    id.vars = c("spot", "x", "y"),
    measure.vars = plot_genes,
    variable.name = "gene",
    value.name = "count"
  )
  dt <- as.data.table(dt)
  dt[, source := label]
  dt
}

plot_dt <- rbind(
  make_long(expr_true, "True"),
  make_long(expr_stcos, "ST-COS"),
  make_long(expr_srt, "SRTsim"),
  make_long(expr_scd3, "scDesign3")
)
if (!is.null(expr_scube)) {
  plot_dt <- rbind(plot_dt, make_long(expr_scube, "scCube"))
}
if (!is.null(expr_scube_train)) {
  plot_dt <- rbind(plot_dt, make_long(expr_scube_train, "scCube_train"))
}

plot_dt$source <- factor(plot_dt$source, levels = c("True", "ST-COS", "SRTsim", "scDesign3", "scCube", "scCube_train"))
# plot_dt <- plot_dt[, {
#   gmin <- min(count, na.rm = TRUE)
#   gmax <- max(count, na.rm = TRUE)
#   scaled <- if (gmax > gmin) (count - gmin) / (gmax - gmin) else 0
#   .(spot, x, y, gene, source, count, count_scaled = scaled,
#     gene_label = sprintf("%s (min=%.2f, max=%.2f)", gene[1], gmin, gmax))
# }, by = gene]

plot_dt_df <- as.data.frame(plot_dt)
names(plot_dt_df) <- make.unique(names(plot_dt_df))

expr_plot <- #ggplot(plot_dt_df, aes(x = x, y = y, color = count_scaled)) +
  ggplot(plot_dt_df, aes(x = x, y = y, color = count)) +
  geom_point(size = 0.2) +
  scale_color_gradientn(
    colors = viridis(7, direction = -1),
    values = c(0, seq(0.01, 1, length.out = 6))
  ) +
  coord_fixed() +
  #facet_grid(gene_label ~ source,) +
  facet_grid(gene ~ source) +
  theme_void() +
  theme(strip.text = element_text(size = 8))

genes_per_page <- 5
gene_chunks <- split(plot_genes, ceiling(seq_along(plot_genes) / genes_per_page))
for (i in seq_along(gene_chunks)) {
  expr_plot_page <- expr_plot %+% plot_dt_df[plot_dt_df$gene %in% gene_chunks[[i]], , drop = FALSE]
  expr_out <- paste0(output_prefix, "_expression_spatial_p", i, ".pdf")
  plot_sources <- levels(plot_dt_df$source)
  width_in <- max(10, length(plot_sources) * 3)
  height_in <- max(6, length(gene_chunks[[i]]) * 2)
  ggsave(expr_out, expr_plot_page, width = width_in, height = height_in, units = "in")
}

## ---- Metrics comparison plots ------------------------------------------
gene_metrics <- plot_dt[, .(
  GeneMean = mean(count, na.rm = TRUE),
  GeneVar = var(count, na.rm = TRUE),
  GeneCV = sd(count, na.rm = TRUE) / mean(count, na.rm = TRUE),
  GeneZeroProp = mean(count == 0, na.rm = TRUE)
), by = .(gene, source)]

location_metrics <- plot_dt[, .(
  LocZeroProp = mean(count == 0, na.rm = TRUE),
  LocLibSize = sum(count, na.rm = TRUE)
), by = .(spot, source)]

gene_metrics_plot <- gene_metrics %>%
  pivot_longer(cols = c(GeneMean, GeneVar, GeneCV, GeneZeroProp),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = source, y = value, fill = source)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene-level Metrics Comparison")

location_metrics_plot <- location_metrics %>%
  pivot_longer(cols = c(LocZeroProp, LocLibSize),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = source, y = value, fill = source)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Location-level Metrics Comparison")

metrics_out <- paste0(output_prefix, "_metrics_comparison.pdf")
pdf(metrics_out, width = 12, height = 8)
print(gene_metrics_plot)
print(location_metrics_plot)
dev.off()

## ---- Pearson Correlation Coefficient (PCC) Analysis --------------------
message("[INFO] Computing PCC for each method vs Truth ...")

# Helper function to compute PCC between two matrices
compute_pcc <- function(true_mat, sim_mat, gene_subset = NULL, spot_subset = NULL) {
  if (!is.null(gene_subset)) {
    common_genes <- intersect(rownames(true_mat), gene_subset)
    true_mat <- true_mat[common_genes, , drop = FALSE]
    sim_mat <- sim_mat[common_genes, , drop = FALSE]
  }
  
  if (!is.null(spot_subset)) {
    common_spots <- intersect(colnames(true_mat), spot_subset)
    true_mat <- true_mat[, common_spots, drop = FALSE]
    sim_mat <- sim_mat[, common_spots, drop = FALSE]
  }
  
  # Ensure matrices have same dimensions and ordering
  common_genes <- intersect(rownames(true_mat), rownames(sim_mat))
  common_spots <- intersect(colnames(true_mat), colnames(sim_mat))
  
  if (length(common_genes) == 0 || length(common_spots) == 0) {
    return(NA_real_)
  }

  true_mat <- true_mat[common_genes, common_spots, drop = FALSE]
  sim_mat <- sim_mat[common_genes, common_spots, drop = FALSE]
  
  # Flatten matrices and compute correlation
  true_vec <- as.vector(true_mat)
  sim_vec <- as.vector(sim_mat)
  
  cor(true_vec, sim_vec, use = "complete.obs")
}

collapse_dup_cols <- function(mat) {
  cn <- colnames(mat)
  if (is.null(cn) || length(cn) == length(unique(cn))) {
    return(mat)
  }
  keep <- !is.na(cn)
  mat <- mat[, keep, drop = FALSE]
  cn <- cn[keep]
  split_idx <- split(seq_along(cn), cn)
  mat_out <- do.call(cbind, lapply(split_idx, function(idx) {
    if (length(idx) == 1) {
      mat[, idx, drop = FALSE]
    } else {
      rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
    }
  }))
  colnames(mat_out) <- names(split_idx)
  mat_out
}

sanitize_matrix <- function(mat) {
  mat <- as.matrix(mat)
  if (!is.null(colnames(mat))) {
    mat <- mat[, !is.na(colnames(mat)), drop = FALSE]
  }
  if (!is.null(rownames(mat))) {
    mat <- mat[!is.na(rownames(mat)), , drop = FALSE]
  }
  mat <- collapse_dup_cols(mat)
  mat
}

# Helper function to orient gene x spot matrices
orient_gene_spot <- function(mat, gene_names, spot_ids) {
  mat <- as.matrix(mat)
  gene_in_rows <- length(intersect(rownames(mat), gene_names))
  gene_in_cols <- length(intersect(colnames(mat), gene_names))
  spot_in_rows <- length(intersect(rownames(mat), spot_ids))
  spot_in_cols <- length(intersect(colnames(mat), spot_ids))
  if (gene_in_cols > gene_in_rows && spot_in_rows > spot_in_cols) {
    mat <- t(mat)
  }
  mat
}

# Prepare matrices for comparison
true_expr <- sanitize_matrix(expr_mat)
stcos_expr <- sanitize_matrix(t(sim_spot_counts))
srt_expr <- sanitize_matrix(srt_sim_counts)
scd3_expr <- sanitize_matrix(t(scd3_counts))
if (is.null(scCube_expr)) {
  scCube_path <- "./Reference_based/scCube_MERFISH_04.csv"
  if (file.exists(scCube_path)) {
    scCube_expr <- read.csv(scCube_path, row.names = 1, check.names = FALSE)
    scCube_meta_path <- "./Reference_based/scCube_MERFISH_04_meta.csv"
    if (file.exists(scCube_meta_path)) {
      scCube_meta <- read.csv(scCube_meta_path, check.names = FALSE)
      scCube_meta$key_xy <- paste0(round(scCube_meta$x), "x", round(scCube_meta$y))
      row_ids <- rownames(scCube_expr)
      spot_key_xy <- paste0(round(coords_mat[, 1]), "x", round(coords_mat[, 2]))
      spot_key_yx <- paste0(round(coords_mat[, 2]), "x", round(coords_mat[, 1]))
      spot_map_xy <- setNames(spot_ids, spot_key_xy)
      spot_map_yx <- setNames(spot_ids, spot_key_yx)
      if (length(intersect(row_ids, names(spot_map_xy))) > 0) {
        mapped_xy <- spot_map_xy[row_ids]
        mapped_yx <- spot_map_yx[row_ids]
        n_xy <- sum(!is.na(mapped_xy))
        n_yx <- sum(!is.na(mapped_yx))
        mapped_spots <- if (n_yx > n_xy) mapped_yx else mapped_xy
      } else if (length(intersect(row_ids, scCube_meta$key_xy)) > 0) {
        mapped_xy <- spot_map_xy[row_ids]
        mapped_yx <- spot_map_yx[row_ids]
        n_xy <- sum(!is.na(mapped_xy))
        n_yx <- sum(!is.na(mapped_yx))
        mapped_spots <- if (n_yx > n_xy) mapped_yx else mapped_xy
      } else {
        cell_keys <- scCube_meta$key_xy[match(row_ids, scCube_meta$Cell)]
        mapped_spots <- spot_map_xy[cell_keys]
      }
      rownames(scCube_expr) <- ifelse(is.na(mapped_spots), rownames(scCube_expr), mapped_spots)
    }
    scCube_expr <- orient_gene_spot(scCube_expr, rownames(true_expr), colnames(true_expr))
  } else {
    message(sprintf("[WARN] scCube results not found at %s; skipping scCube PCC.", scCube_path))
  }
}
if (!is.null(scCube_expr)) {
  scCube_expr <- sanitize_matrix(scCube_expr)
  gene_in_rows <- length(intersect(rownames(scCube_expr), rownames(true_expr)))
  gene_in_cols <- length(intersect(colnames(scCube_expr), rownames(true_expr)))
  spot_in_rows <- length(intersect(rownames(scCube_expr), colnames(true_expr)))
  spot_in_cols <- length(intersect(colnames(scCube_expr), colnames(true_expr)))
  if (gene_in_cols > gene_in_rows && spot_in_rows > spot_in_cols) {
    scCube_expr <- t(scCube_expr)
  }
}
if (is.null(scCube_train_expr) && file.exists(scCube_train_path)) {
  scCube_train_dt <- fread(cmd = paste("zcat", shQuote(scCube_train_path)))
  rownames(scCube_train_dt) <- scCube_train_dt[[1]]
  scCube_train_expr <- as.matrix(scCube_train_dt[, -1, with = FALSE])
  if (file.exists(scCube_train_meta_path)) {
    scCube_train_meta <- fread(cmd = paste("zcat", shQuote(scCube_train_meta_path)))
    if (length(intersect(colnames(scCube_train_expr), scCube_train_meta$Cell)) > 0) {
      scCube_train_expr <- map_scube_cols(scCube_train_expr, coords_mat, spot_ids, scCube_train_meta)
    } else if (length(intersect(rownames(scCube_train_expr), scCube_train_meta$Cell)) > 0) {
      scCube_train_expr <- map_scube_rows(scCube_train_expr, coords_mat, spot_ids, scCube_train_meta)
    }
  }
  scCube_train_expr <- orient_gene_spot(scCube_train_expr, rownames(true_expr), colnames(true_expr))
}
if (!is.null(scCube_train_expr)) {
  scCube_train_expr <- sanitize_matrix(scCube_train_expr)
  gene_in_rows <- length(intersect(rownames(scCube_train_expr), rownames(true_expr)))
  gene_in_cols <- length(intersect(colnames(scCube_train_expr), rownames(true_expr)))
  spot_in_rows <- length(intersect(rownames(scCube_train_expr), colnames(true_expr)))
  spot_in_cols <- length(intersect(colnames(scCube_train_expr), colnames(true_expr)))
  if (gene_in_cols > gene_in_rows && spot_in_rows > spot_in_cols) {
    scCube_train_expr <- t(scCube_train_expr)
  }
}

# Compute overall PCC for each method
pcc_stcos <- compute_pcc(true_expr, stcos_expr)
pcc_srt <- compute_pcc(true_expr, srt_expr)
pcc_scd3 <- compute_pcc(true_expr, scd3_expr)
pcc_scube <- if (!is.null(scCube_expr)) compute_pcc(true_expr, scCube_expr) else NA_real_
pcc_scube_train <- if (!is.null(scCube_train_expr)) compute_pcc(true_expr, scCube_train_expr) else NA_real_

# Compute gene-wise PCC
gene_pcc_results <- data.frame(
  gene = character(0),
  ST_COS = numeric(0),
  SRTsim = numeric(0),
  scDesign3 = numeric(0),
  scCube = numeric(0),
  scCube_train = numeric(0),
  ST_COS_MSE = numeric(0),
  SRTsim_MSE = numeric(0),
  scDesign3_MSE = numeric(0),
  scCube_MSE = numeric(0),
  scCube_train_MSE = numeric(0),
  stringsAsFactors = FALSE
)

common_genes_all <- Reduce(intersect, list(
  rownames(true_expr),
  rownames(stcos_expr),
  rownames(srt_expr),
  rownames(scd3_expr),
  if (!is.null(scCube_expr)) rownames(scCube_expr) else rownames(true_expr),
  if (!is.null(scCube_train_expr)) rownames(scCube_train_expr) else rownames(true_expr)
))

for (g in common_genes_all) {
  pcc_g_stcos <- cor(true_expr[g, ], stcos_expr[g, ], use = "complete.obs")
  pcc_g_srt <- cor(true_expr[g, ], srt_expr[g, ], use = "complete.obs")
  pcc_g_scd3 <- cor(true_expr[g, ], scd3_expr[g, ], use = "complete.obs")
  pcc_g_scube <- if (!is.null(scCube_expr)) cor(true_expr[g, ], scCube_expr[g, ], use = "complete.obs") else NA_real_
  pcc_g_scube_train <- if (!is.null(scCube_train_expr)) cor(true_expr[g, ], scCube_train_expr[g, ], use = "complete.obs") else NA_real_
  mse_stcos <- mean((true_expr[g, ] - stcos_expr[g, ])^2, na.rm = TRUE)
  mse_srt <- mean((true_expr[g, ] - srt_expr[g, ])^2, na.rm = TRUE)
  mse_scd3 <- mean((true_expr[g, ] - scd3_expr[g, ])^2, na.rm = TRUE)
  mse_scube <- if (!is.null(scCube_expr)) mean((true_expr[g, ] - scCube_expr[g, ])^2, na.rm = TRUE) else NA_real_
  mse_scube_train <- if (!is.null(scCube_train_expr)) mean((true_expr[g, ] - scCube_train_expr[g, ])^2, na.rm = TRUE) else NA_real_
  
  gene_pcc_results <- rbind(gene_pcc_results, data.frame(
    gene = g,
    ST_COS = pcc_g_stcos,
    SRTsim = pcc_g_srt,
    scDesign3 = pcc_g_scd3,
    scCube = pcc_g_scube,
    scCube_train = pcc_g_scube_train,
    ST_COS_MSE = mse_stcos,
    SRTsim_MSE = mse_srt,
    scDesign3_MSE = mse_scd3,
    scCube_MSE = mse_scube,
    scCube_train_MSE = mse_scube_train,
    stringsAsFactors = FALSE
  ))
}

# Print overall results
message(sprintf("[RESULTS] Overall PCC - ST-COS: %.4f", pcc_stcos))
message(sprintf("[RESULTS] Overall PCC - SRTsim: %.4f", pcc_srt))
message(sprintf("[RESULTS] Overall PCC - scDesign3: %.4f", pcc_scd3))
if (!is.null(scCube_expr)) {
  message(sprintf("[RESULTS] Overall PCC - scCube: %.4f", pcc_scube))
}
if (!is.null(scCube_train_expr)) {
  message(sprintf("[RESULTS] Overall PCC - scCube_train: %.4f", pcc_scube_train))
}

# Create summary table
pcc_summary <- data.frame(
  Method = c("ST-COS", "SRTsim", "scDesign3", "scCube", "scCube_train"),
  Overall_PCC = c(pcc_stcos, pcc_srt, pcc_scd3, pcc_scube, pcc_scube_train),
  Mean_Gene_PCC = c(
    mean(gene_pcc_results$ST_COS, na.rm = TRUE),
    mean(gene_pcc_results$SRTsim, na.rm = TRUE),
    mean(gene_pcc_results$scDesign3, na.rm = TRUE),
    mean(gene_pcc_results$scCube, na.rm = TRUE),
    mean(gene_pcc_results$scCube_train, na.rm = TRUE)
  ),
  Median_Gene_PCC = c(
    median(gene_pcc_results$ST_COS, na.rm = TRUE),
    median(gene_pcc_results$SRTsim, na.rm = TRUE),
    median(gene_pcc_results$scDesign3, na.rm = TRUE),
    median(gene_pcc_results$scCube, na.rm = TRUE),
    median(gene_pcc_results$scCube_train, na.rm = TRUE)
  )
)

# Save results
pcc_out <- paste0(output_prefix, "_PCC_results.csv")
write.csv(pcc_summary, file = pcc_out, row.names = FALSE)

gene_pcc_out <- paste0(output_prefix, "_gene_PCC_results.csv")
write.csv(gene_pcc_results, file = gene_pcc_out, row.names = FALSE)

gene_metrics_out <- paste0(output_prefix, "_gene_corr_mse_results.csv")
write.csv(gene_pcc_results, file = gene_metrics_out, row.names = FALSE)

message(sprintf("[INFO] PCC results saved to %s and %s", pcc_out, gene_pcc_out))
message(sprintf("[INFO] Gene-level correlation/MSE results saved to %s", gene_metrics_out))

# Create PCC comparison plot
library(reshape2)
gene_pcc_long <- melt(gene_pcc_results, 
                      id.vars = "gene", 
                      variable.name = "Method", 
                      value.name = "PCC")

pcc_plot <- ggplot(gene_pcc_long, aes(x = Method, y = PCC, fill = Method)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Gene-wise PCC Distribution by Method",
       y = "Pearson Correlation Coefficient",
       x = "Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pcc_plot_out <- paste0(output_prefix, "_PCC_comparison.pdf")
ggsave(pcc_plot_out, pcc_plot, width = 8, height = 6)

message(sprintf("[INFO] PCC comparison plot saved to %s", pcc_plot_out))
