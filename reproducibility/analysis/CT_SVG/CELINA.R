#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CELINA)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: CELINA.R SCENARIO SUFFIX [both|estimated|trueprop]")
}

dataset_id <- args[1L]
suffix <- args[2L]
mode <- if (length(args) == 3L) tolower(args[3L]) else "both"
if (!mode %in% c("both", "estimated", "trueprop")) {
  stop("Mode must be one of: both, estimated, trueprop.")
}

base_dir <- file.path(".", "Data", paste0(dataset_id, "_gene"))
sim_dir <- file.path(base_dir, "simulated")
sc_ref_file <- file.path(sim_dir, paste0(dataset_id, "_scref", suffix, ".csv"))
expression_file <- file.path(sim_dir, paste0(dataset_id, "_gene_expression", suffix, ".csv"))
coord_file <- file.path(sim_dir, paste0(dataset_id, "_coord", suffix, ".csv"))
region_file <- file.path(sim_dir, paste0(dataset_id, "_region", suffix, ".csv"))
layer_file <- file.path(sim_dir, paste0(dataset_id, "_layer", suffix, ".csv"))

required_inputs <- c(sc_ref_file, expression_file, coord_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "))
}

build_truth <- function(spot_ids) {
  if (file.exists(region_file)) {
    region_dt <- fread(region_file)
    ct_cols <- setdiff(names(region_dt), c("x_coord", "y_coord", "region"))
    if (length(ct_cols) > 0L) {
      truth <- as.data.frame(region_dt[, ..ct_cols])
      if (nrow(truth) != length(spot_ids)) {
        stop("Region-derived truth row count does not match spot count.")
      }
      rownames(truth) <- spot_ids
      return(truth)
    }
  }

  if (!file.exists(layer_file)) {
    stop("Neither region nor layer input is available for truth proportions.")
  }
  layer_dt <- fread(layer_file)
  if (!"layer" %in% names(layer_dt)) {
    stop("Layer input must contain a 'layer' column.")
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
    # Row-normalized S2/S3/S5/S6 generator weights:
    # ct1 = (1/2, 0, 0), ct2 = (1/4, 3/4, 0), ct3 = (1/3, 0, 1).
    truth <- data.frame(
      ct1 = ifelse(layer_dt$layer == 1, 6 / 13, 0),
      ct2 = ifelse(layer_dt$layer == 1, 3 / 13, ifelse(layer_dt$layer == 2, 1, 0)),
      ct3 = ifelse(layer_dt$layer == 1, 4 / 13, ifelse(layer_dt$layer == 3, 1, 0))
    )
  }
  rownames(truth) <- spot_ids
  truth
}

# CELINA reports unadjusted ACAT-combined P values; genes are called at
# BH-adjusted FDR 0.05 within each cell type (Supplementary Methods S2).
extract_significant_genes <- function(obj, cell_types) {
  result <- setNames(vector("list", length(cell_types)), cell_types)
  for (cell_type in intersect(names(obj@result), cell_types)) {
    table <- obj@result[[cell_type]]
    if (is.data.frame(table) && "CombinedPvals" %in% names(table)) {
      adjusted <- stats::p.adjust(table$CombinedPvals, method = "BH")
      keep <- !is.na(adjusted) & adjusted <= 0.05
      result[[cell_type]] <- rownames(table)[keep]
    }
  }
  result
}

run_celina <- function(proportions, project, output_file) {
  proportions <- as.data.frame(proportions)
  if (nrow(proportions) != ncol(spatial_counts)) {
    stop(project, " proportion row count does not match the spot count.")
  }
  rownames(proportions) <- spot_ids
  cell_types <- names(proportions)
  cell_types <- cell_types[!is.na(cell_types) & nzchar(cell_types)]
  if (length(cell_types) == 0L) {
    stop("No valid cell types were available for ", project, ".")
  }

  obj <- Create_Celina_Object(
    celltype_mat = t(as.matrix(proportions[, cell_types, drop = FALSE])),
    gene_expression_mat = spatial_counts,
    location = as.matrix(location),
    covariates = NULL,
    project = project
  )
  obj <- preprocess_input(
    obj,
    cell_types_to_test = cell_types,
    scRNA_count = sc_counts,
    sc_cell_type_labels = sc_meta$cellType,
    threshold = 0
  )
  obj <- Calculate_Kernel(obj, approximation = FALSE)
  obj <- tryCatch(
    Testing_interaction_all(obj, num_cores = 1),
    error = function(error) {
      message(project, " failed: ", conditionMessage(error))
      NULL
    }
  )
  if (is.null(obj) || length(obj@result) == 0L) {
    message(project, " produced no test results; no output was written.")
    return(FALSE)
  }

  saveRDS(obj@result, sub("\\.RDS$", "_pvalues.RDS", output_file))
  significant_genes <- extract_significant_genes(obj, cell_types)
  saveRDS(significant_genes, output_file)
  for (cell_type in names(significant_genes)) {
    cat(project, cell_type, ":", length(significant_genes[[cell_type]]),
        "significant genes\n")
  }
  cat("Saved:", output_file, "\n")
  TRUE
}

sc_ref_dt <- fread(sc_ref_file)
if (!"Cell_type" %in% names(sc_ref_dt)) {
  stop("Single-cell reference input must contain a 'Cell_type' column.")
}
spatial_dt <- fread(expression_file)
coord_dt <- fread(coord_file)
if (!all(c("x_coord", "y_coord") %in% names(coord_dt))) {
  stop("Coordinate input must contain x_coord and y_coord columns.")
}

spot_ids <- paste0("spot_", seq_len(nrow(coord_dt)))
location <- as.data.frame(coord_dt[, .(x_coord, y_coord)])
names(location) <- c("x", "y")
rownames(location) <- spot_ids

common_genes <- intersect(
  setdiff(names(sc_ref_dt), "Cell_type"),
  names(spatial_dt)
)
if (length(common_genes) == 0L) {
  stop("No overlapping genes between the reference and spatial inputs.")
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

spatial_counts <- t(abs(round(as.matrix(spatial_dt[, ..common_genes])))) / 100
rownames(spatial_counts) <- common_genes
colnames(spatial_counts) <- spot_ids
truth_prop <- build_truth(spot_ids)

failed <- character(0)

if (mode %in% c("both", "trueprop")) {
  ok <- run_celina(
    proportions = truth_prop,
    project = paste0(dataset_id, "_trueprop"),
    output_file = file.path(base_dir, paste0("CELINA_trueprop", suffix, ".RDS"))
  )
  if (!ok) failed <- c(failed, "CELINA_trueprop")
}

if (mode %in% c("both", "estimated")) {
  card_obj <- CARD::createCARDObject(
    sc_count = sc_counts,
    sc_meta = sc_meta,
    spatial_count = spatial_counts,
    spatial_location = location,
    ct.varname = "cellType",
    ct.select = unique(sc_meta$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 0,
    minCountSpot = 0
  )
  card_obj <- tryCatch(
    CARD::CARD_deconvolution(CARD_object = card_obj),
    error = function(error) {
      message("CARD deconvolution failed: ", conditionMessage(error))
      NULL
    }
  )
  # Default CELINA requires estimated proportions. If CARD fails, the run is
  # recorded as failed rather than substituting (perturbed) true proportions.
  if (is.null(card_obj)) {
    failed <- c(failed, "CELINA (CARD deconvolution failed)")
  } else {
    ok <- run_celina(
      proportions = as.data.frame(card_obj@Proportion_CARD),
      project = paste0(dataset_id, "_estimated"),
      output_file = file.path(base_dir, paste0("CELINA", suffix, ".RDS"))
    )
    if (!ok) failed <- c(failed, "CELINA")
  }
}

if (length(failed) > 0L) {
  stop("CELINA runs failed for ", dataset_id, suffix, ": ", paste(failed, collapse = "; "))
}
