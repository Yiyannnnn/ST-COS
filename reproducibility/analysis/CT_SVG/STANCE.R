suppressPackageStartupMessages({
  library(STANCE)
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

gene_expression_file <- file.path(sim_dir, paste0(dataset_id, "_gene_expression", suffix, ".csv"))
coord_file <- file.path(sim_dir, paste0(dataset_id, "_coord", suffix, ".csv"))
region_file <- file.path(sim_dir, paste0(dataset_id, "_region", suffix, ".csv"))
layer_file <- file.path(sim_dir, paste0(dataset_id, "_layer", suffix, ".csv"))

required_inputs <- c(gene_expression_file, coord_file)
if (!all(file.exists(required_inputs))) {
  missing <- required_inputs[!file.exists(required_inputs)]
  stop(sprintf("Missing required input(s): %s", paste(missing, collapse = ", ")))
}

coord_dt <- fread(coord_file)
gene_expr_dt <- fread(gene_expression_file)

spot_ids <- paste0("spot_", seq_len(nrow(coord_dt)))
gene_names <- colnames(gene_expr_dt)

counts <- t(round(as.matrix(gene_expr_dt)))
rownames(counts) <- gene_names
colnames(counts) <- spot_ids

coord_mat <- as.data.frame(coord_dt[, .(x_coord, y_coord)])
rownames(coord_mat) <- spot_ids

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

truth <- build_truth()
region_ct_cols <- colnames(truth)

mySTANCE <- creatSTANCEobject(
  counts = counts,
  pos = coord_mat,
  prop = truth,
  covariates = NULL
)

mySTANCE <- data_preprocess(
  object = mySTANCE,
  gene.threshold = 0,
  spot.threshold = 0,
  normalized = FALSE
)
mySTANCE <- build_kernelMatrix(mySTANCE)

mySTANCE <- runTest1(
  object = mySTANCE,
  correction = FALSE,
  pv.adjust = "BY"
)

genes_list <- rownames(mySTANCE@gene_expression)
sig_genes_overall <- genes_list[mySTANCE@Test_1$p_value_adj < 0.05]

mySTANCE <- runTest2(
  object = mySTANCE,
  Genes_to_test = sig_genes_overall,
  Cell_types_to_test = NULL,
  correction = FALSE,
  ncores = 1
)

mySTANCE <- varcomp_est(
  object = mySTANCE,
  Genes_to_est = sig_genes_overall,
  ncores = 1
)

saveRDS(mySTANCE, file.path(base_dir, paste0("STANCE", suffix, ".RDS")))

sig_by_ct <- setNames(vector("list", length(region_ct_cols)), region_ct_cols)
if (length(mySTANCE@Test_2)) {
  for (ct in intersect(names(mySTANCE@Test_2), region_ct_cols)) {
    df <- mySTANCE@Test_2[[ct]]
    if (is.data.frame(df)) {
      if ("p_value_adj" %in% colnames(df)) {
        keep <- !is.na(df$p_value_adj) & df$p_value_adj < 0.05
        sig_by_ct[[ct]] <- rownames(df)[keep]
      } else if ("p_value" %in% colnames(df)) {
        p_adj <- p.adjust(df$p_value, method = "BH")
        keep <- !is.na(p_adj) & p_adj < 0.05
        sig_by_ct[[ct]] <- rownames(df)[keep]
      }
    }
  }
}

saveRDS(sig_by_ct, file.path(base_dir, paste0("STANCE_siggenes", suffix, ".RDS")))
