#!/usr/bin/env Rscript

.libPaths(unique(c(Sys.getenv("R_LIBS_USER", unset = .libPaths()[1]), .Library, .libPaths())))

suppressPackageStartupMessages({
  library(data.table)
  library(spMOCA)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please provide the two arguments: S and suffix (_N).")
}

scenario <- args[1]
suffix <- args[2]
rep_id <- sub("^_", "", suffix)

base_dir <- file.path("./Data", paste0(scenario, "_gene"))
sim_dir <- file.path(base_dir, "simulated")
output_dir <- "./coexpression/spMOCA_Agent"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(sim_dir, paste0(scenario, "_gene_expression", suffix, ".csv"))
coord_file <- file.path(sim_dir, paste0(scenario, "_coord", suffix, ".csv"))

for (path in c(expr_file, coord_file)) {
  if (!file.exists(path)) {
    stop("Missing input file: ", path)
  }
}

cat("Starting spMOCA for", paste0(scenario, suffix), "\n")
cat("Expression file:", expr_file, "\n")
cat("Coordinate file:", coord_file, "\n")

expr_dt <- fread(expr_file)
coord_dt <- fread(coord_file)

count_mat <- t(as.matrix(expr_dt))
storage.mode(count_mat) <- "numeric"
gene_labels <- colnames(expr_dt)
spot_labels <- paste0("spot_", seq_len(nrow(expr_dt)))
rownames(count_mat) <- gene_labels
colnames(count_mat) <- spot_labels

coords <- as.data.frame(coord_dt[, 1:2, drop = FALSE])
colnames(coords) <- c("x_coord", "y_coord")
rownames(coords) <- spot_labels

feature_list <- gene_labels

obj <- createSPMOCAObject(
  spatial_countMat = count_mat,
  spatial_location = coords,
  feature_list = feature_list
)
obj <- normalizeSPMOCA(obj)
obj <- featureSelectionSPMOCA(
  object = obj,
  feature.selection = NULL,
  feature.list = feature_list
)
obj <- createSpatialKernel(obj, kernelType = "Gaussian", kernelBw = "Rule-Of-Thumb")
obj <- estimatespMOCA(obj)

corr_est <- as.matrix(obj@corr_est)
storage.mode(corr_est) <- "numeric"

write.csv(
  corr_est,
  file = file.path(output_dir, paste0(scenario, "_", rep_id, "_correlation_matrix.csv")),
  row.names = FALSE
)
fwrite(
  data.table(gene_name = gene_labels),
  file.path(output_dir, paste0(scenario, "_", rep_id, "_node_labels.csv"))
)
saveRDS(
  obj,
  file.path(output_dir, paste0(scenario, "_", rep_id, "_spMOCA_object.rds"))
)

summary_dt <- data.table(
  scenario = scenario,
  replicate = as.integer(rep_id),
  n_spots = ncol(count_mat),
  n_genes = nrow(count_mat),
  kernel_type = "Gaussian",
  kernel_bandwidth = "Rule-Of-Thumb",
  min_corr = min(corr_est, na.rm = TRUE),
  max_corr = max(corr_est, na.rm = TRUE)
)
fwrite(summary_dt, file.path(output_dir, paste0(scenario, "_", rep_id, "_summary_stats.csv")))

cat("spMOCA completed for", paste0(scenario, suffix), "\n")
print(summary_dt)
