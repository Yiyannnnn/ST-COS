
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

count_path <- "./data_external/STEREOSEQ/RNA_counts.tsv.gz"
mapping_path <- "./data_external/STEREOSEQ/position.tsv"

output_prefix <- "./Reference_based/MouseOlfactoryBulb_99"
block_size <- 1024
n_cores <- 1
spatial_check_n <- 8500
spatial_moran_thr <- 0
knn_k <- 10
gene_subset_n <- NULL

message("[INFO] Loading Mouse Olfactory Bulb counts ...")
counts_dt <- fread(count_path, header = TRUE, check.names = FALSE)
gene_col <- if ("geneID" %in% names(counts_dt)) "geneID" else names(counts_dt)[1]
gene_names <- counts_dt[[gene_col]]
spot_ids <- setdiff(names(counts_dt), gene_col)
expr_mat <- as.matrix(counts_dt[, ..spot_ids])
rownames(expr_mat) <- gene_names
colnames(expr_mat) <- spot_ids

message("[INFO] Filtering genes with zero proportion > 99% ...")
zero_props <- rowMeans(expr_mat == 0)
keep_genes <- zero_props <= 0.99
expr_mat <- expr_mat[keep_genes, ]
message(sprintf("[INFO] Retained %d genes after zero proportion filtering", nrow(expr_mat)))
gene_names <- rownames(expr_mat)

message("[INFO] Computing library sizes ...")
lib_size <- colSums(expr_mat)
lib_factor <- lib_size / mean(lib_size)
message("[INFO] Filtering spots with zero library size ...")
zero_lib_spots <- lib_size == 0
if (any(zero_lib_spots)) {
  message(sprintf("[INFO] Removing %d spots with zero library size", sum(zero_lib_spots)))
  keep_spots <- !zero_lib_spots
  expr_mat <- expr_mat[, keep_spots, drop = FALSE]
  lib_size <- lib_size[keep_spots]
  lib_factor <- lib_factor[keep_spots]
  spot_ids <- spot_ids[keep_spots]
  message(sprintf("[INFO] Retained %d spots after zero library size filtering", length(spot_ids)))
}

message("[INFO] Loading coordinates ...")
mapping_dt <- fread(mapping_path)
spot_ids_chr <- as.character(spot_ids)
spot_ids_num <- suppressWarnings(as.integer(spot_ids_chr))
mapping_dt[, label := suppressWarnings(as.integer(label))]

coords_dt <- mapping_dt[.(spot_ids_num), on = .(label)]
coords_mat <- as.matrix(coords_dt[, .(x, y)])
if (any(is.na(coords_mat))) {
  missing_ids <- spot_ids_chr[is.na(coords_mat[, 1]) | is.na(coords_mat[, 2])]
  stop(sprintf("Missing coordinates for %d spots after matching mapping table. Example IDs: %s",
               length(missing_ids), paste(head(missing_ids, 5), collapse = ", ")))
}
rownames(coords_mat) <- spot_ids_chr

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

message("[INFO] Detecting spatial genes (subset) ...")
Wnb <- knearneigh(coords_mat, k = knn_k)
W <- nb2listw(knn2nb(Wnb))
#W <- readRDS("./Reference_based/MouseOlfactoryBulb_99_spatial_weights.rds")
moran_flag <- function(y, listw, thr = 0.15) {
  S0 <- Szero(listw)
  res <- moran(x = y, listw = listw, n = length(y), S0 = S0, zero.policy = TRUE)$I
  res >= thr
}
all_moran <- sapply(gene_names, function(g) {
  moran(x = expr_mat[g, ], listw = W, n = length(spot_ids), S0 = Szero(W), zero.policy = TRUE)$I
})
spatial_candidates <- gene_names[seq_len(min(spatial_check_n, length(gene_names)))]
spatial_genes <- spatial_candidates[sapply(spatial_candidates, function(g) {
  moran_flag(expr_mat[g, ], W, thr = summary(all_moran)[["Min."]])
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
  lib_size = list(CT = lib_size),
  theta_g = list(CT = theta_g)
)

sim_spot_counts <- sim_list[["CT"]]
rownames(sim_spot_counts) <- spot_ids
colnames(sim_spot_counts) <- gene_names

message("[INFO] Writing outputs ...")
counts_out <- paste0(output_prefix, "_simulated_counts_all_gene.csv")
coords_out <- paste0(output_prefix, "_simulated_coordinates_all_gene.csv")
write.csv(sim_spot_counts, file = counts_out, row.names = TRUE)
write.csv(
  data.frame(s1 = coords_mat[, 1], s2 = coords_mat[, 2], row.names = spot_ids),
  file = coords_out,
  row.names = TRUE
)

# message("[INFO] Saving spatial weights matrix ...")
# weights_out <- paste0(output_prefix, "_spatial_weights.rds")
# saveRDS(W, file = weights_out)