suppressPackageStartupMessages({
  .libPaths(unique(c(.Library, Sys.getenv("R_LIBS_USER", unset = .libPaths()[1]), .libPaths())))
  library(spVC)
  library(sp)
  library(BPST)
  library(Triangulation)
  library(MGLM)
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
layer_file <- file.path(sim_dir, paste0(dataset_id, "_layer", suffix, ".csv"))
region_file <- file.path(sim_dir, paste0(dataset_id, "_region", suffix, ".csv"))

required_inputs <- c(gene_expression_file, coord_file)
if (!all(file.exists(required_inputs))) {
  missing <- required_inputs[!file.exists(required_inputs)]
  stop(sprintf("Missing required input(s): %s", paste(missing, collapse = ", ")))
}

sample_dt <- fread(gene_expression_file)
coord_dt <- fread(coord_file)
if (!all(c("x_coord", "y_coord") %in% colnames(coord_dt))) {
  stop("Coordinate file must contain x_coord and y_coord columns.")
}

Y <- abs(t(round(as.matrix(sample_dt))))
rownames(Y) <- colnames(sample_dt)
colnames(Y) <- paste0("spot_", seq_len(nrow(sample_dt)))

S <- as.matrix(coord_dt[, .(x_coord, y_coord)])
rownames(S) <- colnames(Y)

build_truth <- function() {
  if (file.exists(region_file)) {
    region_dt <- fread(region_file)
    ct_cols <- setdiff(colnames(region_dt), c("x_coord", "y_coord", "region"))
    if (length(ct_cols) > 0) {
      X <- as.matrix(region_dt[, ..ct_cols])
      colnames(X) <- ct_cols
      rownames(X) <- colnames(Y)
      return(X)
    }
  }

  if (!file.exists(layer_file)) {
    stop("Neither region file nor layer file exists; cannot build truth proportions.")
  }
  layer_dt <- fread(layer_file)
  if (!("layer" %in% colnames(layer_dt))) {
    stop("Layer file must contain a 'layer' column.")
  }
  if (nrow(layer_dt) != ncol(Y)) {
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
  X <- as.matrix(truth_df)
  rownames(X) <- colnames(Y)
  colnames(X) <- colnames(truth_df)
  X
}

X <- build_truth()

x_min <- min(S[, 1], na.rm = TRUE)
x_max <- max(S[, 1], na.rm = TRUE)
y_min <- min(S[, 2], na.rm = TRUE)
y_max <- max(S[, 2], na.rm = TRUE)

boundary <- data.frame(
  x = c(x_min, x_max, x_max, x_min),
  y = c(y_min, y_min, y_max, y_max)
)

Tr.cell <- TriMesh(boundary, n = 2)
V <- as.matrix(Tr.cell$V)
Tr <- as.matrix(Tr.cell$Tr)

simu.results <- test.spVC(
  Y = Y,
  X = X,
  S = S,
  V = V,
  Tr = Tr,
  para.cores = 2,
  filter.min.nonzero = 0,
  filter.spot.counts = 0
)

ct_names <- colnames(X)
significant_genes <- setNames(vector("list", length(ct_names)), ct_names)

extract_constant_pvals <- function(result_entry, ct_names) {
  out <- setNames(rep(NA_real_, length(ct_names) + 1L), c(ct_names, "gamma_0"))
  if (is.null(result_entry) || is.null(result_entry$p.value)) {
    return(out)
  }

  pvals <- suppressWarnings(as.numeric(result_entry$p.value))
  pval_names <- names(result_entry$p.value)

  if (!is.null(pval_names)) {
    gamma_idx <- match("gamma_0", pval_names)
    if (!is.na(gamma_idx)) {
      out[["gamma_0"]] <- pvals[gamma_idx]
    }
    for (idx in seq_along(ct_names)) {
      candidate_names <- c(paste0("beta_", ct_names[[idx]]), paste0("beta_X", idx))
      match_idx <- match(candidate_names, pval_names)
      match_idx <- match_idx[!is.na(match_idx)]
      if (length(match_idx) > 0) {
        out[[ct_names[[idx]]]] <- pvals[match_idx[[1L]]]
      }
    }
  }

  expected_len <- length(ct_names) + 2L
  if (all(is.na(out)) && length(pvals) >= expected_len) {
    out[ct_names] <- pvals[2:(length(ct_names) + 1L)]
    out[["gamma_0"]] <- pvals[length(ct_names) + 2L]
  }

  out
}

extract_spvc_sig_genes <- function(spvc_result, ct_names, p_adjust_method = "BH", p_adjust_thresh = 0.05) {
  out <- setNames(vector("list", length(ct_names)), ct_names)
  constant_results <- spvc_result$results.constant
  varying_results <- spvc_result$results.varying

  if (is.null(constant_results) || length(constant_results) == 0 || is.null(varying_results) || length(varying_results) == 0) {
    return(out)
  }

  constant_pvals <- lapply(constant_results, extract_constant_pvals, ct_names = ct_names)
  pval_mat <- do.call(rbind, constant_pvals)
  rownames(pval_mat) <- names(constant_results)

  beta_adj <- vapply(
    ct_names,
    FUN.VALUE = numeric(nrow(pval_mat)),
    FUN = function(ct_name) p.adjust(pval_mat[, ct_name], method = p_adjust_method)
  )
  beta_adj <- as.matrix(beta_adj)
  colnames(beta_adj) <- ct_names
  rownames(beta_adj) <- rownames(pval_mat)

  gamma_adj <- p.adjust(pval_mat[, "gamma_0"], method = p_adjust_method)
  names(gamma_adj) <- rownames(pval_mat)

  varying_gene_ids <- intersect(names(varying_results), names(gamma_adj))
  if (length(varying_gene_ids) == 0) {
    return(out)
  }

  for (gene_id in varying_gene_ids) {
    if (is.na(gamma_adj[[gene_id]]) || gamma_adj[[gene_id]] >= p_adjust_thresh) {
      next
    }
    ct_hits <- ct_names[!is.na(beta_adj[gene_id, ]) & beta_adj[gene_id, ] < p_adjust_thresh]
    for (ct_name in ct_hits) {
      out[[ct_name]] <- c(out[[ct_name]], gene_id)
    }
  }

  lapply(out, unique)
}

significant_genes <- extract_spvc_sig_genes(simu.results, ct_names)

saveRDS(significant_genes, file.path(base_dir, paste0("spVC_trueprop", suffix, ".RDS")))
