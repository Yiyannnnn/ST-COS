suppressPackageStartupMessages({
  library(CELINA)
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

build_truth <- function(spot_ids) {
  if (file.exists(region_file)) {
    region_dt <- fread(region_file)
    ct_cols <- setdiff(colnames(region_dt), c("x_coord", "y_coord", "region"))
    if (length(ct_cols) > 0) {
      truth <- as.data.frame(region_dt[, ..ct_cols])
      if (nrow(truth) != length(spot_ids)) {
        stop("Region-derived truth row count does not match spot count.")
      }
      rownames(truth) <- spot_ids
      return(truth)
    }
  }

  if (!file.exists(layer_file)) {
    stop("Neither region file nor layer file exists; cannot derive truth proportions.")
  }
  layer_dt <- fread(layer_file)
  if (!("layer" %in% colnames(layer_dt))) {
    stop("Layer file must contain a 'layer' column.")
  }
  if (nrow(layer_dt) != length(spot_ids)) {
    stop("Layer-derived truth row count does not match spot count.")
  }

  if (dataset_id %in% c("S1", "S4")) {
    truth <- data.frame(
      ct1 = rep(1 / 3, nrow(layer_dt)),
      ct2 = rep(1 / 3, nrow(layer_dt)),
      ct3 = rep(1 / 3, nrow(layer_dt))
    )
  } else {
    truth <- data.frame(
      ct1 = ifelse(layer_dt$layer == 1, 1 / 2, ifelse(layer_dt$layer == 2, 1 / 4, 0)),
      ct2 = ifelse(layer_dt$layer == 1, 1 / 4, ifelse(layer_dt$layer == 2, 3 / 4, 0)),
      ct3 = ifelse(layer_dt$layer == 1, 1 / 4, ifelse(layer_dt$layer == 2, 0, 1))
    )
  }
  rownames(truth) <- spot_ids
  truth
}

sc_ref_dt <- fread(sc_ref_file)
if (!("Cell_type" %in% colnames(sc_ref_dt))) {
  stop("sc_ref file must contain 'Cell_type' column.")
}

spatial_dt <- fread(gene_expression_file)
coord_dt <- fread(coord_file)
if (!all(c("x_coord", "y_coord") %in% colnames(coord_dt))) {
  stop("Coordinate file must contain columns x_coord and y_coord.")
}

spot_ids <- paste0("spot_", seq_len(nrow(coord_dt)))
location <- as.data.frame(coord_dt[, .(x_coord, y_coord)])
colnames(location) <- c("x", "y")
rownames(location) <- spot_ids

sc_gene_names <- setdiff(colnames(sc_ref_dt), "Cell_type")
spatial_gene_names <- colnames(spatial_dt)
common_genes <- intersect(sc_gene_names, spatial_gene_names)
if (length(common_genes) == 0) {
  stop("No overlapping genes between scRNA reference and spatial expression.")
}

sc_meta <- data.frame(
  cellType = sc_ref_dt$Cell_type,
  cellID = paste0("cell_", seq_len(nrow(sc_ref_dt))),
  sampleInfo = 1,
  stringsAsFactors = FALSE
)
rownames(sc_meta) <- sc_meta$cellID

sc_counts <- t(as.matrix(sc_ref_dt[, ..common_genes]))
colnames(sc_counts) <- sc_meta$cellID
rownames(sc_counts) <- common_genes

spatial_counts <- t(round(as.matrix(spatial_dt[, ..common_genes])))
spatial_counts <- abs(spatial_counts) / 100
rownames(spatial_counts) <- common_genes
colnames(spatial_counts) <- spot_ids

truth_prop <- build_truth(spot_ids)
cell_types_to_test <- colnames(truth_prop)
cell_types_to_test <- cell_types_to_test[!is.na(cell_types_to_test) & nzchar(cell_types_to_test)]
if (length(cell_types_to_test) == 0) {
  stop("No valid cell-type columns in truth proportions.")
}

obj <- Create_Celina_Object(
  celltype_mat = t(as.matrix(truth_prop)),
  gene_expression_mat = as.matrix(spatial_counts),
  location = as.matrix(location),
  covariates = NULL,
  project = paste0(dataset_id, "_trueprop")
)

obj <- preprocess_input(
  obj,
  cell_types_to_test = cell_types_to_test,
  scRNA_count = as.matrix(sc_counts),
  sc_cell_type_labels = sc_meta$cellType,
  threshold = 0
)

obj <- Calculate_Kernel(obj, approximation = FALSE)
obj <- tryCatch(
  Testing_interaction_all(obj, num_cores = 1),
  error = function(e) {
    message("Testing_interaction_all failed: ", e$message)
    NULL
  }
)

significant_genes <- setNames(vector("list", length(cell_types_to_test)), cell_types_to_test)
if (!is.null(obj) && length(obj@result) > 0) {
  for (ct in names(obj@result)) {
    df <- obj@result[[ct]]
    if (is.data.frame(df) && "CombinedPvals" %in% colnames(df)) {
      sig_rows <- df[df$CombinedPvals <= 0.05, , drop = FALSE]
      significant_genes[[ct]] <- rownames(sig_rows)
    }
  }
}

out_file <- file.path(base_dir, paste0("CELINA_trueprop", suffix, ".RDS"))
saveRDS(significant_genes, out_file)

for (ct in names(significant_genes)) {
  cat(ct, ":", length(significant_genes[[ct]]), "significant genes\n")
}
cat("Saved:", out_file, "\n")
