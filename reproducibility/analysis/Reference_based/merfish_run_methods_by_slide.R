#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Run ST-COS, SRTsim, scDesign3 for a given MERFISH slide
## Usage: Rscript merfish_run_methods_by_slide.R 04
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(mgcv)
  library(spdep)
  library(parallel)
})

suppressPackageStartupMessages({
  library(SRTsim)
  library(scDesign3)
  library(SingleCellExperiment)
})

source("./Data/simulator_YZ.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide slide number, e.g. 04")
}
slide_id <- sprintf("%02d", as.integer(args[1]))

count_path <- file.path("./data_external/MERFISH/count",
                        paste0("X0_", slide_id, "_count.csv"))
mapping_path <- file.path("./data_external/MERFISH/mapping",
                          paste0("X0_", slide_id, "_mapping.csv"))

output_prefix <- file.path("./Reference_based",
                           paste0("MERFISH_X0_", slide_id))
block_size <- 1024
n_cores <- 1
spatial_check_n <- 300
spatial_moran_thr <- 0
knn_k <- 10
gene_subset_n <- NULL

message(sprintf("[INFO] Slide %s: loading counts", slide_id))
counts_dt <- fread(count_path)
spot_ids <- counts_dt[[1]]
gene_names <- colnames(counts_dt)[-1]
expr_mat_spot_gene <- as.matrix(counts_dt[, -1, with = FALSE])
expr_mat <- t(expr_mat_spot_gene)
rownames(expr_mat) <- gene_names
colnames(expr_mat) <- spot_ids

message(sprintf("[INFO] Slide %s: loading coordinates", slide_id))
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

message("[INFO] Writing ST-COS outputs ...")
counts_out <- paste0(output_prefix, "_simulated_counts.csv")
coords_out <- paste0(output_prefix, "_simulated_coordinates.csv")
write.csv(sim_spot_counts, file = counts_out, row.names = TRUE)
write.csv(
  data.frame(s1 = coords_mat[, 1], s2 = coords_mat[, 2], row.names = spot_ids),
  file = coords_out,
  row.names = TRUE
)

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
