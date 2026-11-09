
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
spatial_check_n <- 300
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


## ---- SRTsim -------------------------------------------------------------
message("[INFO] Running SRTsim simulation ...")
gene_use <- gene_names
if (!is.null(gene_subset_n)) {
  gene_use <- gene_names[seq_len(min(gene_subset_n, length(gene_names)))]
}
srt_counts <- expr_mat#[gene_use, ]
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