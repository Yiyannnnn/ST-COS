#!/usr/bin/env Rscript
## ---------------------------------------------------------------
## Run ST-COS, SRTsim, and scDesign3 for a single DLPFC (HumanPilot10X) slide
## Usage: Rscript DLPFC_run_methods_by_slide.R <slide_id>
## Example: Rscript DLPFC_run_methods_by_slide.R 151673
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(mgcv)
  library(spdep)
  library(parallel)
  library(SRTsim)
  library(scDesign3)
  library(SingleCellExperiment)
})

source("./Data/simulator_YZ.R")

## ---- Command-line args -----------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript DLPFC_run_methods_by_slide.R <slide_id> [method]")
slide_id   <- as.character(args[1])
run_method <- if (length(args) >= 2) tolower(args[2]) else "all"
# run_method: "all" | "stcos" | "srtsim" | "scdesign3"
message(sprintf("[%s] ===== Starting simulation pipeline =====", slide_id))

## ---- Paths -----------------------------------------------------------
data_dir   <- "./data_external/DLPFC/"
out_dir    <- "./Reference_based/DLPFC"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

output_prefix <- file.path(out_dir, slide_id)
h5_file       <- paste0(slide_id, "_filtered_feature_bc_matrix.h5")
layer_file    <- file.path(data_dir, "layers", paste0(slide_id, "_layers.txt"))

## ---- Tuning parameters -----------------------------------------------
block_size        <- 1024
n_cores           <- 1
spatial_moran_thr <- 0.15
knn_k             <- 6
corr_thr          <- 0.05    # sparsify Sigma: correlations below this -> 0
min_spot_pct      <- 0.05    # gene must be expressed in >= 5% of annotated spots

## ---- 1. Load 10X Spatial Data ----------------------------------------
message(sprintf("[%s] Loading 10X Spatial h5 data ...", slide_id))
srt <- Seurat::Load10X_Spatial(
  data.dir      = file.path(data_dir, slide_id, ""),
  filename      = h5_file,
  assay         = "Spatial",
  slice         = "slice1",
  filter.matrix = TRUE,
  to.upper      = FALSE
)

## ---- 2. Coordinates --------------------------------------------------
coord_df <- as.data.frame(GetTissueCoordinates(srt, image = "slice1"))
colnames(coord_df) <- c("row", "col", "cellnames")

## ---- 3. Layer annotations (used for spot subsetting and SRTsim labels)
if (!file.exists(layer_file)) {
  warning(sprintf("[%s] Layer file not found: %s — using all spots", slide_id, layer_file))
  layer_coord        <- coord_df
  layer_coord$label  <- "CT"
} else {
  layer_dt    <- fread(layer_file, header = FALSE, col.names = c("cellnames", "label"))
  layer_coord <- merge(coord_df, layer_dt, by = "cellnames", all = FALSE)
  message(sprintf("[%s] %d / %d spots with layer annotation",
                  slide_id, nrow(layer_coord), nrow(coord_df)))
}

spot_ids   <- layer_coord$cellnames
coords_mat <- as.matrix(layer_coord[, c("row", "col")])
rownames(coords_mat) <- spot_ids
n_spots    <- length(spot_ids)

## ---- 4. Gene Filtering -----------------------------------------------
message(sprintf("[%s] Filtering genes (non-MT, >= %.0f%% of %d spots) ...",
                slide_id, min_spot_pct * 100, n_spots))
all_counts <- srt@assays$Spatial$counts[, spot_ids]
min_spots_thr <- max(1, round(min_spot_pct * n_spots))

not_mt     <- !startsWith(rownames(all_counts), "MT-")
expr_n_spots <- Matrix::rowSums(all_counts > 0)
keep_genes <- not_mt & (expr_n_spots >= min_spots_thr)

expr_mat   <- as.matrix(all_counts[keep_genes, ])
gene_names <- rownames(expr_mat)
n_genes    <- length(gene_names)
message(sprintf("[%s] %d genes retained", slide_id, n_genes))

## ---- 5. NB Moments ---------------------------------------------------
fit_nb_fast <- function(y) {
  mu <- mean(y); s2 <- var(y)
  th <- if (s2 > mu) mu^2 / (s2 - mu) else 1e6
  c(mu = mu, theta = th)
}

message(sprintf("[%s] Estimating NB moments ...", slide_id))
nb_par  <- t(apply(expr_mat, 1, fit_nb_fast))
mu_g    <- nb_par[, "mu"]
theta_g <- nb_par[, "theta"]

## ---- 6. PIT / Z Matrix -----------------------------------------------
counts2z_block <- function(idx, Y, mu, th, rng_block, eps = 1e-6) {
  out <- matrix(0, nrow = length(idx), ncol = ncol(Y))
  for (ii in seq_along(idx)) {
    g  <- idx[ii]; y <- Y[g, ]
    u0 <- pnbinom(pmax(y - 1, 0), size = th[g], mu = mu[g])
    u1 <- pnbinom(y,              size = th[g], mu = mu[g])
    u  <- u0 + rng_block[ii, ] * (u1 - u0)
    u  <- pmin(pmax(u, eps), 1 - eps)
    out[ii, ] <- qnorm(u)
  }
  out
}

counts2z_mat <- function(Y, mu, th, block = 1024, n_cores = 1, eps = 1e-6) {
  G      <- nrow(Y)
  idx_ls <- split(seq_len(G), ceiling(seq_len(G) / block))
  rng    <- matrix(runif(G * ncol(Y)), nrow = G)
  worker <- function(idx) counts2z_block(idx, Y, mu, th, rng[idx, , drop = FALSE], eps)
  if (n_cores > 1) {
    do.call(rbind, mclapply(idx_ls, worker, mc.cores = n_cores))
  } else {
    do.call(rbind, lapply(idx_ls, worker))
  }
}

message(sprintf("[%s] Computing PIT Z matrix (%d genes x %d spots) ...",
                slide_id, n_genes, n_spots))
Z <- t(counts2z_mat(expr_mat, mu_g, theta_g, block = block_size, n_cores = n_cores))
# Z is spots x genes

## ---- 7. Correlation Matrix Sigma (ST-COS) ----------------------------
message(sprintf("[%s] Building gene-gene correlation matrix ...", slide_id))
Sigma <- cor(Z, use = "pairwise.complete.obs")
Sigma[is.na(Sigma)] <- 0
diag(Sigma) <- 1
Sigma[abs(Sigma) < corr_thr] <- 0   # sparsify
Sigma <- makespd(Sigma)
if (min(eigen(Sigma, only.values = TRUE)$values) <= 0) {
  message(sprintf("[%s] Repairing positive-definiteness via nearPD ...", slide_id))
  Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = TRUE)$mat)
}
rownames(Sigma) <- gene_names
colnames(Sigma) <- gene_names

## ---- 8. Library Sizes ------------------------------------------------
lib_size   <- colSums(expr_mat)
lib_factor <- lib_size / mean(lib_size)

## ---- 9. Spatial Gene Detection (Moran's I) ---------------------------
message(sprintf("[%s] Detecting spatially variable genes ...", slide_id))
Wnb <- knearneigh(coords_mat, k = knn_k)
W   <- nb2listw(knn2nb(Wnb))

moran_flag <- function(y, listw, thr = 0.15) {
  S0 <- Szero(listw)
  moran(x = y, listw = listw, n = length(y), S0 = S0, zero.policy = TRUE)$I > thr
}

spatial_genes <- gene_names[sapply(gene_names, function(g)
  moran_flag(expr_mat[g, ], W, thr = spatial_moran_thr))]
message(sprintf("[%s] %d / %d spatial genes detected", slide_id, length(spatial_genes), n_genes))

## ---- 10. Mean Matrix with NB-GAM for Spatial Genes ------------------
mu_mat <- outer(lib_factor, mu_g)
dimnames(mu_mat) <- list(spot_ids, gene_names)

if (length(spatial_genes) > 0) {
  message(sprintf("[%s] Refitting NB-GAM for %d spatial genes ...", slide_id, length(spatial_genes)))
  for (g in spatial_genes) {
    dat <- data.frame(
      y       = as.numeric(expr_mat[g, ]),
      row     = coords_mat[, 1],
      col     = coords_mat[, 2],
      log_lib = log(pmax(lib_factor, 1e-6))
    )
    gam_fit <- tryCatch(
      gam(y ~ offset(log_lib) + s(row, col, bs = "tp", k = 50),
          family = nb(link = "log"), data = dat, method = "REML"),
      error = function(e) {
        message(sprintf("  [WARN] GAM failed for gene %s: %s", g, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(gam_fit)) {
      mu_mat[, g] <- fitted(gam_fit)
      theta_g[g]  <- gam_fit$family$getTheta(TRUE)
    }
  }
}

## ===========================================================================
##  METHOD 1: ST-COS
## ===========================================================================
if (run_method %in% c("all", "stcos")) {
message(sprintf("[%s] --- Running ST-COS ---", slide_id))
t0_stcos <- proc.time()

mean_param   <- list(CT = mu_mat)
gene_rel     <- as.data.table(
  matrix(1, nrow = n_genes, ncol = 1, dimnames = list(gene_names, "CT"))
)
cov_gene     <- array(NA, dim = c(n_genes, n_genes, 1),
                      dimnames = list(gene_names, gene_names, "CT"))
cov_gene[, , "CT"] <- array(Sigma, dim = c(n_genes, n_genes, 1),
                             dimnames = list(gene_names, gene_names, "CT"))

sim_list     <- generate_true_count(
  coord            = coords_mat,
  mean_param       = mean_param,
  gene_relationship = gene_rel,
  cov_gene         = cov_gene,
  lib_size         = list(CT = lib_size),
  theta_g          = list(CT = theta_g)
)
stcos_counts <- sim_list[["CT"]]          # spots x genes
rownames(stcos_counts) <- spot_ids
colnames(stcos_counts) <- gene_names

stcos_time <- (proc.time() - t0_stcos)[["elapsed"]]
message(sprintf("[%s] ST-COS completed in %.1f sec", slide_id, stcos_time))

stcos_out  <- paste0(output_prefix, "_STCOS_simulated_counts.csv")
coords_out <- paste0(output_prefix, "_simulated_coordinates.csv")
write.csv(stcos_counts,
          file = stcos_out, row.names = TRUE)
write.csv(data.frame(row = coords_mat[, 1], col = coords_mat[, 2],
                     label = layer_coord$label, row.names = spot_ids),
          file = coords_out, row.names = TRUE)
message(sprintf("[%s] ST-COS counts saved: %s", slide_id, stcos_out))
} # end stcos

## ===========================================================================
##  METHOD 2: SRTsim
## ===========================================================================
if (run_method %in% c("all", "srtsim")) {
message(sprintf("[%s] --- Running SRTsim ---", slide_id))
t0_srt <- proc.time()

srt_loc <- data.frame(
  y     = coords_mat[, 1],
  x     = coords_mat[, 2],
  label = layer_coord$label,
  row.names = spot_ids
)
simSRT <- createSRT(count_in = expr_mat, loc_in = srt_loc)
simSRT <- srtsim_fit(simSRT, sim_schem = "tissue")
simSRT <- srtsim_count(simSRT)
srt_sim_counts <- simCounts(simSRT)

srt_time <- (proc.time() - t0_srt)[["elapsed"]]
message(sprintf("[%s] SRTsim completed in %.1f sec", slide_id, srt_time))

srt_out <- paste0(output_prefix, "_SRTsim_counts.csv")
write.csv(as.matrix(srt_sim_counts), file = srt_out, row.names = TRUE)
message(sprintf("[%s] SRTsim counts saved: %s", slide_id, srt_out))
} # end srtsim

## ===========================================================================
##  METHOD 3: scDesign3
## ===========================================================================
if (run_method %in% c("all", "scdesign3")) {
message(sprintf("[%s] --- Running scDesign3 ---", slide_id))
t0_scd3 <- proc.time()

sce <- SingleCellExperiment(
  assays  = list(counts = expr_mat),
  colData = DataFrame(
    spatial1  = coords_mat[, 1],
    spatial2  = coords_mat[, 2],
    cell_type = layer_coord$label
  )
)
scd3_fit <- scdesign3(
  sce              = sce,
  assay_use        = "counts",
  celltype         = "cell_type",
  pseudotime       = NULL,
  spatial          = c("spatial1", "spatial2"),
  other_covariates = NULL,
  mu_formula       = "s(spatial1, spatial2, bs = 'gp', k = 200)",
  sigma_formula    = "1",
  family_use       = "nb",
  n_cores          = n_cores,
  usebam           = FALSE,
  corr_formula     = "1",
  copula           = "gaussian",
  DT               = TRUE,
  pseudo_obs       = FALSE,
  return_model     = FALSE,
  nonzerovar       = FALSE,
  parallelization  = "pbmcapply"
)
scd3_counts <- t(scd3_fit$new_count)   # spots x genes

scd3_time <- (proc.time() - t0_scd3)[["elapsed"]]
message(sprintf("[%s] scDesign3 completed in %.1f sec", slide_id, scd3_time))

scd3_out <- paste0(output_prefix, "_scDesign3_counts.csv")
write.csv(scd3_counts, file = scd3_out, row.names = TRUE)
message(sprintf("[%s] scDesign3 counts saved: %s", slide_id, scd3_out))
} # end scdesign3

## ===========================================================================
##  Summary
## ===========================================================================
timing_df <- data.frame(
  slide    = slide_id,
  n_genes  = n_genes,
  n_spots  = n_spots,
  STCOS_sec     = if (exists("stcos_time"))   stcos_time   else NA,
  SRTsim_sec    = if (exists("srt_time"))     srt_time     else NA,
  scDesign3_sec = if (exists("scd3_time"))    scd3_time    else NA
)
timing_out <- paste0(output_prefix, "_timing.csv")
write.csv(timing_df, file = timing_out, row.names = FALSE)

message(sprintf("[%s] ===== All methods complete =====", slide_id))
message(sprintf("[%s]   n_genes=%d  n_spots=%d", slide_id, n_genes, n_spots))
if (exists("stcos_time"))   message(sprintf("[%s]   ST-COS:    %.1f sec", slide_id, stcos_time))
if (exists("srt_time"))     message(sprintf("[%s]   SRTsim:    %.1f sec", slide_id, srt_time))
if (exists("scd3_time"))    message(sprintf("[%s]   scDesign3: %.1f sec", slide_id, scd3_time))
message(sprintf("[%s] Timing saved: %s", slide_id, timing_out))
