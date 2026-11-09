#!/usr/bin/env Rscript
## ---------------------------------------------------------------
## Gene-gene correlation evaluation for one DLPFC slide.
## Computes Pearson correlation of log1p(counts) on a random
## subsample of genes, then saves:
##   - <slide>_gene_corr_summary.csv   (RMSE, MAE, Pearson r per method)
##   - <slide>_gene_corr_pairs.csv.gz  (sampled gene pairs for scatter)
## Usage: Rscript DLPFC_gene_corr_by_slide.R <slide_id>
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript DLPFC_gene_corr_by_slide.R <slide_id>")
slide_id <- as.character(args[1])
message(sprintf("[%s] Gene-gene correlation evaluation ...", slide_id))

data_dir  <- "./data_external/DLPFC/"
out_dir   <- "./Reference_based/DLPFC"
pfx       <- file.path(out_dir, slide_id)

n_genes_sub <- 1000   # subsample for correlation (→ ~500K pairs)
set.seed(42)

## ---- helpers ---------------------------------------------------------
read_mat <- function(path, orient = c("gene_x_spot", "spot_x_gene")) {
  orient <- match.arg(orient)
  dt  <- fread(path)
  rn  <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- rn
  if (orient == "spot_x_gene") mat <- t(mat)
  mat
}

## ---- Load simulated matrices -----------------------------------------
message(sprintf("[%s] Loading simulated matrices ...", slide_id))
stcos_mat <- read_mat(paste0(pfx, "_STCOS_simulated_counts.csv"),  orient = "spot_x_gene")
srt_mat   <- read_mat(paste0(pfx, "_SRTsim_counts.csv"),           orient = "gene_x_spot")
scd3_mat  <- read_mat(paste0(pfx, "_scDesign3_counts.csv"),        orient = "spot_x_gene")

## ---- Load real data --------------------------------------------------
message(sprintf("[%s] Loading truth ...", slide_id))
srt_obj   <- Seurat::Load10X_Spatial(
  data.dir      = file.path(data_dir, slide_id, ""),
  filename      = paste0(slide_id, "_filtered_feature_bc_matrix.h5"),
  assay         = "Spatial", slice = "slice1",
  filter.matrix = TRUE, to.upper = FALSE
)
truth_mat <- as.matrix(srt_obj@assays$Spatial$counts)

## ---- Load scCube (after truth so real_lib is available) --------------
sccube_lognorm_path <- paste0(pfx, "_scCube_lognorm.csv.gz")
sccube_meta_path    <- paste0(pfx, "_scCube_meta.csv")
sccube_mat <- NULL
sccube_lognorm_mat <- NULL
if (file.exists(sccube_lognorm_path) && file.exists(sccube_meta_path)) {
  sccube_dt   <- fread(sccube_lognorm_path)
  sccube_gene <- sccube_dt[[1]]
  sccube_raw  <- as.matrix(sccube_dt[, -1, with = FALSE])
  rownames(sccube_raw) <- sccube_gene
  sccube_meta_df <- fread(sccube_meta_path)

  spatial_dir <- file.path("./data_external/DLPFC",
                           slide_id, "spatial")
  for (f in c("tissue_positions.txt","tissue_positions_list.txt")) {
    tp_path <- file.path(spatial_dir, f)
    if (file.exists(tp_path)) {
      tp <- tryCatch(fread(tp_path), error = function(e) NULL)
      if (!is.null(tp) && !("barcode" %in% names(tp)))
        setnames(tp, c("barcode","in_tissue","array_row","array_col",
                       "pxl_row_in_fullres","pxl_col_in_fullres"))
      if (!is.null(tp)) {
        if ("in_tissue" %in% names(tp)) tp <- tp[in_tissue == 1]
        tp[, xy_key := paste0(array_col, "_", array_row)]
        sccube_meta_df[, xy_key := paste0(x, "_", y)]
        xy_to_bc <- tp$barcode; names(xy_to_bc) <- tp$xy_key
        matched_bc <- xy_to_bc[sccube_meta_df$xy_key]
        valid <- !is.na(matched_bc)
        real_lib <- colSums(truth_mat)
        matched_lib <- real_lib[matched_bc[valid]]
        matched_lib[is.na(matched_lib)] <- median(real_lib, na.rm = TRUE)
        sc_sub <- sccube_raw[, valid, drop = FALSE]
        sccube_mat <- round(sweep(expm1(sc_sub), 2, matched_lib / 1e4, "*"))
        colnames(sccube_mat) <- matched_bc[valid]
        # keep log-norm matrix for correlation (avoids double log1p from back-transform)
        sccube_lognorm_mat <- sc_sub
        colnames(sccube_lognorm_mat) <- matched_bc[valid]
        message(sprintf("[%s]   scCube: %d genes x %d spots",
                        slide_id, nrow(sccube_mat), ncol(sccube_mat)))
      }
      break
    }
  }
}

## ---- Intersect genes & spots -----------------------------------------
gene_sets  <- list(rownames(stcos_mat), rownames(srt_mat), rownames(scd3_mat), rownames(truth_mat))
spot_sets  <- list(colnames(stcos_mat), colnames(truth_mat))
if (!is.null(sccube_mat)) {
  gene_sets <- c(gene_sets, list(rownames(sccube_mat)))
  spot_sets <- c(spot_sets, list(colnames(sccube_mat)))
}
common_genes <- Reduce(intersect, gene_sets)
common_spots <- Reduce(intersect, spot_sets)

truth_mat <- truth_mat[common_genes, common_spots]
stcos_mat <- stcos_mat[common_genes, common_spots]
srt_mat   <- srt_mat[common_genes, ]    # SRTsim spots may differ — use all
scd3_mat  <- scd3_mat[common_genes, common_spots]
if (!is.null(sccube_mat))        sccube_mat        <- sccube_mat[common_genes, common_spots]
if (!is.null(sccube_lognorm_mat)) sccube_lognorm_mat <- sccube_lognorm_mat[common_genes, common_spots]

n_common <- length(common_genes)
message(sprintf("[%s] %d common genes, %d spots", slide_id, n_common, length(common_spots)))

## ---- Subsample genes -------------------------------------------------
sub_idx   <- sort(sample(seq_len(n_common), min(n_genes_sub, n_common)))
sub_genes <- common_genes[sub_idx]

truth_sub <- truth_mat[sub_genes, ]
stcos_sub <- stcos_mat[sub_genes, ]
srt_sub   <- srt_mat[sub_genes, ]
scd3_sub  <- scd3_mat[sub_genes, ]

## ---- log1p transform -------------------------------------------------
log1p_cor <- function(mat) {
  cor(t(log1p(mat)), method = "pearson", use = "pairwise.complete.obs")
}

## ---- Estimate Sigma via PIT on truth (ST-COS input correlation) ------
fit_nb_moment <- function(x) {
  mu <- mean(x); v <- var(x)
  theta <- if (v > mu) mu^2 / (v - mu) else 1e4
  c(mu = mu, theta = min(theta, 1e6))
}

pit_to_z <- function(y, mu, theta, eps = 1e-6) {
  u0 <- pnbinom(pmax(y - 1L, 0L), size = theta, mu = mu)
  u1 <- pnbinom(y,                 size = theta, mu = mu)
  u  <- u0 + runif(length(y)) * (u1 - u0)
  qnorm(pmin(pmax(u, eps), 1 - eps))
}

message(sprintf("[%s] Computing Sigma (PIT of truth, %d genes) ...", slide_id, length(sub_genes)))
set.seed(42)
Z_sub <- t(apply(truth_sub, 1, function(y) {
  par <- fit_nb_moment(y)
  pit_to_z(y, par["mu"], par["theta"])
}))
# Z_sub: genes x spots  →  cor() expects spots x genes
C_sigma <- cor(t(Z_sub), method = "pearson", use = "pairwise.complete.obs")
C_sigma[is.na(C_sigma)] <- 0; diag(C_sigma) <- 1

## ---- Subsample scCube to same genes (if available) -------------------
# Use log-norm matrix directly — already in log1p(norm×1e4) space, no back-transform needed
sccube_lognorm_sub <- if (!is.null(sccube_lognorm_mat)) sccube_lognorm_mat[sub_genes, ] else NULL

message(sprintf("[%s] Computing correlation matrices (%d genes) ...", slide_id, length(sub_genes)))
C_truth  <- log1p_cor(truth_sub)
C_stcos  <- log1p_cor(stcos_sub)
C_srt    <- log1p_cor(srt_sub)
C_scd3   <- log1p_cor(scd3_sub)
C_sccube <- if (!is.null(sccube_lognorm_sub))
              cor(t(sccube_lognorm_sub), method = "pearson", use = "pairwise.complete.obs") else NULL
message(sprintf("[%s] Correlation matrices done.", slide_id))

## ---- Upper triangle indices ------------------------------------------
ut <- which(upper.tri(C_truth))

truth_vec  <- C_truth[ut]
stcos_vec  <- C_stcos[ut]
srt_vec    <- C_srt[ut]
scd3_vec   <- C_scd3[ut]
sccube_vec <- if (!is.null(C_sccube)) C_sccube[ut] else rep(NA_real_, length(ut))
sigma_vec  <- C_sigma[ut]

## ---- Summary metrics per method --------------------------------------
eval_metrics <- function(sim_vec, truth_vec, label) {
  valid <- !is.na(sim_vec)
  diff  <- sim_vec[valid] - truth_vec[valid]
  data.table(
    method = label,
    RMSE   = sqrt(mean(diff^2)),
    MAE    = mean(abs(diff)),
    PearsonR = cor(sim_vec[valid], truth_vec[valid])
  )
}

rows <- list(
  eval_metrics(sigma_vec,  truth_vec, "STCOS_Sigma"),
  eval_metrics(stcos_vec,  truth_vec, "STCOS"),
  eval_metrics(srt_vec,    truth_vec, "SRTsim"),
  eval_metrics(scd3_vec,   truth_vec, "scDesign3")
)
if (!is.null(C_sccube)) rows <- c(rows, list(eval_metrics(sccube_vec, truth_vec, "scCube")))

summary_dt <- rbindlist(rows)
summary_dt[, slide := slide_id]
setcolorder(summary_dt, c("slide", "method", "RMSE", "MAE", "PearsonR"))

summary_out <- paste0(pfx, "_gene_corr_summary.csv")
fwrite(summary_dt, summary_out)
message(sprintf("[%s] Summary saved: %s", slide_id, summary_out))

## ---- Sample pairs for scatter (max 50K per method) -------------------
n_pairs  <- length(truth_vec)
n_sample <- min(50000, n_pairs)
samp_idx <- sample(n_pairs, n_sample)

pairs_dt <- data.table(
  slide       = slide_id,
  truth       = truth_vec[samp_idx],
  STCOS_Sigma = sigma_vec[samp_idx],
  STCOS       = stcos_vec[samp_idx],
  SRTsim      = srt_vec[samp_idx],
  scDesign3   = scd3_vec[samp_idx],
  scCube      = sccube_vec[samp_idx]
)

pairs_out <- paste0(pfx, "_gene_corr_pairs.csv.gz")
fwrite(pairs_dt, pairs_out, compress = "gzip")
message(sprintf("[%s] Pairs saved: %s  (n=%d)", slide_id, pairs_out, n_sample))
message(sprintf("[%s] Done.", slide_id))
