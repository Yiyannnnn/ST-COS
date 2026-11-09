suppressPackageStartupMessages({
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

default_ct <- sort(unique(as.character(cell_types)))
empty_sig <- setNames(vector("list", length(default_ct)), default_ct)

myRCTD <- tryCatch(
  create.RCTD(puck, reference, max_cores = 2),
  error = function(e) {
    message("create.RCTD failed: ", e$message)
    NULL
  }
)
if (is.null(myRCTD)) {
  saveRDS(empty_sig, file.path(base_dir, paste0("CSIDE", suffix, ".RDS")))
  saveRDS(empty_sig, file.path(base_dir, paste0("CSIDE_trueprop", suffix, ".RDS")))
  quit(save = "no", status = 0)
}

myRCTD <- tryCatch(
  run.RCTD(myRCTD, doublet_mode = "full"),
  error = function(e) {
    message("run.RCTD failed: ", e$message)
    NULL
  }
)
if (is.null(myRCTD)) {
  saveRDS(empty_sig, file.path(base_dir, paste0("CSIDE", suffix, ".RDS")))
  saveRDS(empty_sig, file.path(base_dir, paste0("CSIDE_trueprop", suffix, ".RDS")))
  quit(save = "no", status = 0)
}

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

denom <- max(coords_mat[, "x_coord"])
if (isTRUE(all.equal(denom, 0))) {
  denom <- 1
}
explanatory.variable <- coords_mat[, "x_coord"] / denom
names(explanatory.variable) <- spot_ids

extract_sig_genes <- function(rctd_obj, cell_types) {
  sig_list <- rctd_obj@de_results$sig_gene_list
  out <- setNames(vector("list", length(cell_types)), cell_types)
  if (length(sig_list) == 0) {
    return(out)
  }
  for (nm in intersect(names(sig_list), cell_types)) {
    df <- sig_list[[nm]]
    if (is.data.frame(df) && nrow(df) > 0) {
      out[[nm]] <- rownames(df)
    }
  }
  out
}

safe_run_cside <- function(rctd_obj) {
  tryCatch({
    run.CSIDE.single(
      rctd_obj,
      gene_threshold = 0,
      cell_type_threshold = 0,
      explanatory.variable = explanatory.variable,
      doublet_mode = FALSE
    )
  }, error = function(e) {
    message("Error in CSIDE analysis: ", e$message)
    invisible(NULL)
  })
}

myRCTD_est <- safe_run_cside(myRCTD)
if (!is.null(myRCTD_est)) {
  myRCTD <- myRCTD_est
}
sig_est <- extract_sig_genes(myRCTD, region_ct_cols)
saveRDS(sig_est, file.path(base_dir, paste0("CSIDE", suffix, ".RDS")))

truth_sparse <- as(as.matrix(truth_mat), "dgCMatrix")
rownames(truth_sparse) <- spot_ids
myRCTD@results$weights <- truth_sparse

myRCTD_truth <- safe_run_cside(myRCTD)
if (!is.null(myRCTD_truth)) {
  myRCTD <- myRCTD_truth
}
sig_truth <- extract_sig_genes(myRCTD, region_ct_cols)
saveRDS(sig_truth, file.path(base_dir, paste0("CSIDE_trueprop", suffix, ".RDS")))
