#!/usr/bin/env Rscript
## ---------------------------------------------------------------
## Compute gene- and location-level performance metrics for
## ST-COS, SRTsim, scDesign3 outputs in the DLPFC run,
## mirroring plot_precomputed_simulation_metrics.R
##
## Usage: Rscript DLPFC_metrics_by_slide.R <slide_id>
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript DLPFC_metrics_by_slide.R <slide_id>")
slide_id <- as.character(args[1])
message(sprintf("[%s] Computing performance metrics ...", slide_id))

data_dir  <- "./data_external/DLPFC/"
out_dir   <- "./Reference_based/DLPFC"
pfx       <- file.path(out_dir, slide_id)

## ---- helpers ---------------------------------------------------------
# Read a CSV and return a gene x spot matrix.
# orient: "gene_x_spot" (rows=genes) or "spot_x_gene" (rows=spots, needs transpose)
read_mat <- function(path, orient = c("gene_x_spot", "spot_x_gene")) {
  orient <- match.arg(orient)
  dt  <- fread(path)
  rn  <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- rn
  if (orient == "spot_x_gene") mat <- t(mat)   # -> gene x spot
  mat
}

intersect_genes_spots <- function(mat_list, genes, spots) {
  genes <- intersect(genes, Reduce(intersect, lapply(mat_list, rownames)))
  spots <- intersect(spots, Reduce(intersect, lapply(mat_list, colnames)))
  lapply(mat_list, function(m) m[genes, spots, drop = FALSE])
}

## ---- 1. Load simulated matrices ---------------------------------------
stcos_path   <- paste0(pfx, "_STCOS_simulated_counts.csv")
srt_path     <- paste0(pfx, "_SRTsim_counts.csv")
scd3_path    <- paste0(pfx, "_scDesign3_counts.csv")
coord_path   <- paste0(pfx, "_simulated_coordinates.csv")
sccube_path  <- paste0(pfx, "_scCube_lognorm.csv.gz")
sccube_meta_path <- paste0(pfx, "_scCube_meta.csv")

for (p in c(stcos_path, srt_path, scd3_path)) {
  if (!file.exists(p)) stop(sprintf("Missing file: %s", p))
}

message(sprintf("[%s] Loading simulated matrices ...", slide_id))
# STCOS saved as spots x genes -> transpose to gene x spot
stcos_mat <- read_mat(stcos_path, orient = "spot_x_gene")
# SRTsim saved as genes x spots (simCounts output)
srt_mat   <- read_mat(srt_path,  orient = "gene_x_spot")
# scDesign3 saved as spots x genes -> transpose
scd3_mat  <- read_mat(scd3_path, orient = "spot_x_gene")

message(sprintf("[%s]   STCOS:    %d genes x %d spots", slide_id, nrow(stcos_mat), ncol(stcos_mat)))
message(sprintf("[%s]   SRTsim:   %d genes x %d spots", slide_id, nrow(srt_mat),   ncol(srt_mat)))
message(sprintf("[%s]   scDesign3:%d genes x %d spots", slide_id, nrow(scd3_mat),  ncol(scd3_mat)))

## ---- 2. Load real data -----------------------------------------------
message(sprintf("[%s] Loading real 10X Spatial data ...", slide_id))
srt_obj <- Seurat::Load10X_Spatial(
  data.dir      = file.path(data_dir, slide_id, ""),
  filename      = paste0(slide_id, "_filtered_feature_bc_matrix.h5"),
  assay         = "Spatial",
  slice         = "slice1",
  filter.matrix = TRUE,
  to.upper      = FALSE
)
truth_mat <- as.matrix(srt_obj@assays$Spatial$counts)   # gene x spot
message(sprintf("[%s]   Truth:    %d genes x %d spots", slide_id, nrow(truth_mat), ncol(truth_mat)))

## ---- 2b. Load scCube (optional) ----------------------------------------
sccube_mat <- NULL
has_sccube <- file.exists(sccube_path) && file.exists(sccube_meta_path)

if (has_sccube) {
  message(sprintf("[%s] Loading scCube ...", slide_id))
  sccube_dt   <- fread(sccube_path)
  sccube_gene <- sccube_dt[[1]]
  sccube_raw  <- as.matrix(sccube_dt[, -1, with = FALSE])
  rownames(sccube_raw) <- sccube_gene

  sccube_meta_df <- fread(sccube_meta_path)

  ## Back-transform log1p(norm×1e4) to pseudo-raw counts using real lib sizes.
  ## scCube meta has x=array_col, y=array_row matching reference spots.
  ## Build (x,y) -> barcode lookup from tissue_positions file.
  spatial_dir <- file.path(data_dir, slide_id, "spatial")
  tp_candidates <- c("tissue_positions.txt", "tissue_positions_list.txt")
  tp_path <- NULL
  for (f in tp_candidates) {
    p <- file.path(spatial_dir, f)
    if (file.exists(p)) { tp_path <- p; break }
  }

  if (!is.null(tp_path)) {
    tp <- tryCatch(fread(tp_path), error = function(e) NULL)
    if (!is.null(tp) && ncol(tp) >= 6 && !("barcode" %in% names(tp))) {
      setnames(tp, c("barcode","in_tissue","array_row","array_col",
                     "pxl_row_in_fullres","pxl_col_in_fullres"))
    }
    if (!is.null(tp) && "in_tissue" %in% names(tp)) tp <- tp[in_tissue == 1]

    if (!is.null(tp) && all(c("barcode","array_col","array_row") %in% names(tp))) {
      tp[, xy_key := paste0(array_col, "_", array_row)]
      sccube_meta_df[, xy_key := paste0(x, "_", y)]

      xy_to_bc       <- tp$barcode; names(xy_to_bc) <- tp$xy_key
      matched_barcodes <- xy_to_bc[sccube_meta_df$xy_key]
      valid_idx      <- !is.na(matched_barcodes)

      real_lib   <- colSums(truth_mat)
      matched_lib <- real_lib[matched_barcodes[valid_idx]]
      matched_lib[is.na(matched_lib)] <- median(real_lib, na.rm = TRUE)

      # keep only matched spots
      sccube_sub  <- sccube_raw[, valid_idx, drop = FALSE]
      # expm1(log1p(expr/1e4)) = expr/1e4 (normalized proportion × 1e4)
      # scale: proportion × real_lib_size = pseudo-raw
      sccube_norm <- expm1(sccube_sub)   # genes × spots, scale 0–1e4
      sccube_mat  <- round(sweep(sccube_norm, 2, matched_lib / 1e4, "*"))
      colnames(sccube_mat) <- matched_barcodes[valid_idx]
      message(sprintf("[%s]   scCube:   %d genes x %d spots (pseudo-raw)",
                      slide_id, nrow(sccube_mat), ncol(sccube_mat)))
    }
  }
  if (is.null(sccube_mat)) {
    message(sprintf("[%s]   scCube: tissue_positions not found, skipping", slide_id))
    has_sccube <- FALSE
  }
}

## ---- 3. Intersect genes & spots across matrices ----------------------
all_genes <- rownames(stcos_mat)
all_spots <- colnames(stcos_mat)

base_list <- list(Truth = truth_mat, STCOS = stcos_mat, SRTsim = srt_mat, scDesign3 = scd3_mat)
if (has_sccube) base_list[["scCube"]] <- sccube_mat

mat_list <- intersect_genes_spots(base_list, genes = all_genes, spots = all_spots)

n_genes <- nrow(mat_list[[1]]); n_spots <- ncol(mat_list[[1]])
message(sprintf("[%s] Shared: %d genes x %d spots", slide_id, n_genes, n_spots))

## ---- 4. Long-format assembly -----------------------------------------
combine_long <- function(mat, method) {
  dt      <- as.data.table(t(mat))          # spots x genes
  dt[, cellnames := colnames(mat)]
  melt(dt, id.vars = "cellnames",
       variable.name = "gene", value.name = "count")[, source := method]
}

method_levels <- c("Truth", "STCOS", "SRTsim", "scDesign3", "scCube")

message(sprintf("[%s] Building combined long table ...", slide_id))
combined <- rbindlist(Map(combine_long, mat_list, names(mat_list)), use.names = TRUE)
combined[, method := factor(source, levels = method_levels)]

## ---- 5. Gene-level metrics -------------------------------------------
gene_metrics <- combined[, .(
  GeneMean     = mean(count, na.rm = TRUE),
  GeneVar      = var(count,  na.rm = TRUE),
  GeneCV       = sd(count,   na.rm = TRUE) / mean(count, na.rm = TRUE),
  GeneZeroProp = mean(count == 0, na.rm = TRUE)
), by = .(gene, method)]

gene_summary <- gene_metrics[, .(
  Mean_GeneMean     = mean(GeneMean,     na.rm = TRUE),
  Mean_GeneVar      = mean(GeneVar,      na.rm = TRUE),
  Mean_GeneCV       = mean(GeneCV,       na.rm = TRUE),
  Mean_GeneZeroProp = mean(GeneZeroProp, na.rm = TRUE)
), by = method]

## ---- 6. Location-level metrics ---------------------------------------
location_metrics <- combined[, .(
  LocZeroProp = mean(count == 0, na.rm = TRUE),
  LocLibSize  = sum(count,       na.rm = TRUE)
), by = .(cellnames, method)]

location_summary <- location_metrics[, .(
  Mean_LocZeroProp = mean(LocZeroProp, na.rm = TRUE),
  Mean_LocLibSize  = mean(LocLibSize,  na.rm = TRUE)
), by = method]

## ---- 7. Print summaries ----------------------------------------------
cat("\n=== Gene-level metrics summary ===\n")
print(gene_summary)
cat("\n=== Location-level metrics summary ===\n")
print(location_summary)

## ---- 8. Save CSVs ----------------------------------------------------
fwrite(gene_metrics,     file = paste0(pfx, "_gene_metrics.csv"))
fwrite(location_metrics, file = paste0(pfx, "_location_metrics.csv"))
fwrite(gene_summary,     file = paste0(pfx, "_gene_summary.csv"))
fwrite(location_summary, file = paste0(pfx, "_location_summary.csv"))
message(sprintf("[%s] Metrics CSVs saved.", slide_id))

## ---- 9. Violin plots -------------------------------------------------
plot_violin <- function(dt, cols, title) {
  dt %>%
    pivot_longer(cols = all_of(cols), names_to = "metric", values_to = "value") %>%
    ggplot(aes(x = method, y = value, fill = method)) +
    geom_violin(trim = FALSE, alpha = 0.9) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.6, outlier.shape = NA) +
    facet_wrap(~metric, scales = "free_y") +
    theme_minimal(base_size = 18) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          strip.text   = element_text(size = 14)) +
    labs(title = title, x = "Method", y = "Value")
}

pdf_path <- paste0(pfx, "_metrics_comparison.pdf")
pdf(pdf_path, width = 14, height = 10)
print(plot_violin(as.data.frame(gene_metrics),
                  c("GeneMean", "GeneVar", "GeneCV", "GeneZeroProp"),
                  sprintf("Gene-level Metrics — DLPFC %s (%d genes)", slide_id, n_genes)))
print(plot_violin(as.data.frame(location_metrics),
                  c("LocZeroProp", "LocLibSize"),
                  sprintf("Location-level Metrics — DLPFC %s (%d spots)", slide_id, n_spots)))
dev.off()
message(sprintf("[%s] PDF saved: %s", slide_id, pdf_path))
message(sprintf("[%s] Done.", slide_id))
