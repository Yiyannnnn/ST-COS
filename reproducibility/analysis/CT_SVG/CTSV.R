suppressPackageStartupMessages({
  # Keep the module library first, but include the user library so batch
  # reruns can still load locally installed method packages.
  .libPaths(unique(c(.Library, Sys.getenv("R_LIBS_USER", unset = .libPaths()[1]), .libPaths())))
  library(CTSV)
  library(SpatialExperiment)
  library(spacexr)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please provide the two arguments: S and suffix (_N).")
}

dataset_id <- args[1]
suffix <- args[2]
base_dir <- file.path("./Data", paste0(dataset_id, "_gene"))
sim_dir <- file.path(base_dir, "simulated")

sc_ref_file <- file.path(sim_dir, paste0(dataset_id, "_scref", suffix, ".csv"))
gene_expression_file <- file.path(sim_dir, paste0(dataset_id, "_gene_expression", suffix, ".csv"))
coord_file <- file.path(sim_dir, paste0(dataset_id, "_coord", suffix, ".csv"))
region_file <- file.path(sim_dir, paste0(dataset_id, "_region", suffix, ".csv"))
layer_file <- file.path(sim_dir, paste0(dataset_id, "_layer", suffix, ".csv"))

required_inputs <- c(sc_ref_file, gene_expression_file, coord_file)
if (!all(file.exists(required_inputs))) {
  missing <- required_inputs[!file.exists(required_inputs)]
  stop(sprintf("Missing required input(s): %s", paste(missing, collapse = ", ")))
}

coord_dt <- fread(coord_file)
gene_expr_dt <- fread(gene_expression_file)
sc_ref_dt <- fread(sc_ref_file)

if (!all(c("x_coord", "y_coord") %in% names(coord_dt))) {
  stop("Coordinate file must contain columns named x_coord and y_coord.")
}

spot_ids <- paste0("spot_", seq_len(nrow(coord_dt)))
spot_gene_names <- colnames(gene_expr_dt)
sc_gene_names <- setdiff(colnames(sc_ref_dt), "Cell_type")
common_genes <- intersect(spot_gene_names, sc_gene_names)
if (length(common_genes) == 0) {
  stop("No overlapping genes between sc reference and spatial expression.")
}
gene_names <- common_genes

spot_counts <- t(round(as.matrix(gene_expr_dt[, ..common_genes])))
rownames(spot_counts) <- gene_names
colnames(spot_counts) <- spot_ids

coords_mat <- as.matrix(coord_dt[, .(x_coord, y_coord)])
rownames(coords_mat) <- spot_ids

nUMI_spots <- colSums(spot_counts)
puck <- SpatialRNA(
  coords = data.frame(coords_mat, row.names = spot_ids),
  counts = spot_counts,
  nUMI = nUMI_spots
)

cell_ids <- paste0("cell_", seq_len(nrow(sc_ref_dt)))
cell_types <- sc_ref_dt$Cell_type
names(cell_types) <- cell_ids

sc_counts <- t(round(as.matrix(sc_ref_dt[, ..common_genes])))
rownames(sc_counts) <- gene_names
colnames(sc_counts) <- cell_ids
nUMI_cells <- colSums(sc_counts)

reference <- Reference(sc_counts, as.factor(cell_types), nUMI_cells)

build_truth <- function() {
  if (file.exists(region_file)) {
    region_dt <- fread(region_file)
    region_ct_cols <- setdiff(colnames(region_dt), c("x_coord", "y_coord", "region"))
    if (length(region_ct_cols) == 0) {
      stop("Region file does not contain cell type proportion columns.")
    }
    truth_mat <- as.matrix(region_dt[, ..region_ct_cols])
    rownames(truth_mat) <- spot_ids
    colnames(truth_mat) <- region_ct_cols
    return(truth_mat)
  }

  if (!file.exists(layer_file)) {
    stop("Neither region file nor layer file exists; cannot build truth proportions.")
  }
  layer_dt <- fread(layer_file)
  if (!("layer" %in% colnames(layer_dt))) {
    stop("Layer file must contain a 'layer' column.")
  }
  if (nrow(layer_dt) != length(spot_ids)) {
    stop("Layer file row count does not match number of spots.")
  }

  if (dataset_id %in% c("S1", "S4")) {
    truth_df <- data.frame(
      ct1 = rep(1 / 3, nrow(layer_dt)),
      ct2 = rep(1 / 3, nrow(layer_dt)),
      ct3 = rep(1 / 3, nrow(layer_dt))
    )
  } else {
    truth_df <- data.frame(
      ct1 = ifelse(layer_dt$layer == 1, 1 / 2, ifelse(layer_dt$layer == 2, 1 / 4, 0)),
      ct2 = ifelse(layer_dt$layer == 1, 1 / 4, ifelse(layer_dt$layer == 2, 3 / 4, 0)),
      ct3 = ifelse(layer_dt$layer == 1, 1 / 4, ifelse(layer_dt$layer == 2, 0, 1))
    )
  }
  truth_mat <- as.matrix(truth_df)
  rownames(truth_mat) <- spot_ids
  colnames(truth_mat) <- colnames(truth_df)
  truth_mat
}

truth_mat <- build_truth()
region_ct_cols <- colnames(truth_mat)

ct_sums <- colSums(truth_mat)
keep_ct <- names(ct_sums[ct_sums > 0])
if (length(keep_ct) == 0) {
  stop("All cell type proportions are zero.")
}
if (length(keep_ct) < length(region_ct_cols)) {
  message(sprintf(
    "Dropping %d cell type(s) with zero total proportion: %s",
    length(setdiff(region_ct_cols, keep_ct)),
    paste(setdiff(region_ct_cols, keep_ct), collapse = ", ")
  ))
}
truth_mat <- truth_mat[, keep_ct, drop = FALSE]
region_ct_cols <- keep_ct

spot_spe <- SpatialExperiment(
  assays = list(counts = spot_counts),
  colData = DataFrame(spot_id = spot_ids, row.names = spot_ids),
  spatialCoords = coords_mat
)

extract_significant <- function(ctsv_result, cell_types) {
  if (is.null(ctsv_result$qval)) {
    return(setNames(vector("list", length(cell_types)), cell_types))
  }
  sig <- svGene(ctsv_result$qval, 0.05)
  out <- setNames(vector("list", length(cell_types)), cell_types)
  if (!is.null(sig$SVGene)) {
    for (nm in intersect(names(sig$SVGene), cell_types)) {
      out[[nm]] <- sig$SVGene[[nm]]
    }
  }
  out
}

empty_sig <- setNames(vector("list", length(region_ct_cols)), region_ct_cols)
rctd_error <- NULL
myRCTD <- tryCatch(
  {
    # `spacexr` fails on the PSOCK worker path here for this environment.
    # Running serially avoids the missing-function worker error.
    rctd_obj <- create.RCTD(puck, reference, max_cores = 1)
    run.RCTD(rctd_obj, doublet_mode = "full")
  },
  error = function(e) {
    rctd_error <<- conditionMessage(e)
    message(sprintf("RCTD failed for %s%s: %s", dataset_id, suffix, rctd_error))
    NULL
  }
)

common_ct <- region_ct_cols
if (!is.null(myRCTD)) {
  weights_est <- as.matrix(myRCTD@results$weights)
  common_ct <- intersect(colnames(weights_est), region_ct_cols)
  if (length(common_ct) == 0) {
    stop("No overlapping cell types between estimated weights and truth.")
  }
  weights_est <- weights_est[spot_ids, common_ct, drop = FALSE]
  result_est <- CTSV(spot_spe, weights_est, num_core = 1)
  sig_est <- extract_significant(result_est, common_ct)
} else {
  sig_est <- empty_sig
}
saveRDS(sig_est, file.path(base_dir, paste0("CTSV", suffix, ".RDS")))

truth_mat <- truth_mat[, common_ct, drop = FALSE]
result_truth <- CTSV(spot_spe, truth_mat, num_core = 1)
sig_truth <- extract_significant(result_truth, common_ct)
saveRDS(sig_truth, file.path(base_dir, paste0("CTSV_truth", suffix, ".RDS")))

if (!is.null(rctd_error)) {
  writeLines(
    c(
      sprintf("dataset=%s", dataset_id),
      sprintf("suffix=%s", suffix),
      sprintf("rctd_error=%s", rctd_error)
    ),
    file.path(base_dir, paste0("CTSV_debug", suffix, ".txt"))
  )
}
