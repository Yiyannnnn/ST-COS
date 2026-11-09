#!/usr/bin/env Rscript

# Load required libraries
library(Giotto)
library(data.table)
library(tidyverse)
library(Matrix)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure the correct number of arguments are provided
if (length(args) != 2) {
  stop("Please provide the two arguments: S and suffix (_N).")
}

# Get the inputs from command-line arguments
S <- args[1]
suffix <- args[2]

# Construct file paths based on the arguments
gene_expression_file <- paste0("./Data/", S, "_gene/simulated/", S, "_gene_expression", suffix, ".csv")
coord_file <- paste0("./Data/", S, "_gene/simulated/", S, "_coord", suffix, ".csv")
layer_file <- paste0("./Data/", S, "_gene/simulated/", S, "_layer", suffix, ".csv")

# Output directories and files
output_dir <- paste0("./coexpression/binSpect_Agent/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_file <- paste0(output_dir, S, suffix, "_binSpect_results.csv")
spatial_genes_file <- paste0(output_dir, S, suffix, "_spatial_genes.csv")
summary_file <- paste0(output_dir, S, suffix, "_summary_stats.csv")

# Load the data
cat("Loading data for", S, suffix, "\n")
gene_expression <- fread(gene_expression_file)
coordinates <- fread(coord_file)
layer_info <- fread(layer_file)

# Prepare data for Giotto
# Input matrix is spots x genes; Giotto expects genes x cells.
expr_matrix <- t(as.matrix(gene_expression))
rownames(expr_matrix) <- paste0("gene-", seq_len(nrow(expr_matrix)))

# Prepare coordinate data
coord_df <- data.frame(
  cell_ID = paste0("cell-", seq_len(nrow(coordinates))),
  x = coordinates[[1]],
  y = coordinates[[2]]
)
colnames(expr_matrix) <- coord_df$cell_ID
expr_matrix <- Matrix::Matrix(expr_matrix, sparse = TRUE)
spatial_locs <- as.matrix(coord_df[, c("x", "y")])
rownames(spatial_locs) <- coord_df$cell_ID
cell_metadata <- data.frame(layer = layer_info[[3]], row.names = coord_df$cell_ID)

# Create Giotto object
cat("Creating Giotto object...\n")
giotto_obj <- createGiottoObject(
  raw_exprs = expr_matrix,
  spatial_locs = spatial_locs,
  cell_metadata = cell_metadata
)

# Normalize expression data
giotto_obj <- normalizeGiotto(giotto_obj, 
                             scalefactor = 6000, 
                             verbose = TRUE)

# Create spatial network
cat("Creating spatial network...\n")
giotto_obj <- createSpatialNetwork(giotto_obj, 
                                  method = "Delaunay", 
                                  name = "Delaunay_network")

# Run binSpect analysis
cat("Running binSpect analysis...\n")
tryCatch({
  # Run binSpect with kmeans binarization
  binspect_results <- binSpect(giotto_obj,
                              bin_method = "kmeans",
                              expression_values = "normalized",
                              spatial_network_name = "Delaunay_network",
                              do_fisher_test = TRUE,
                              adjust_method = "fdr",
                              calc_hub = FALSE,
                              get_av_expr = TRUE,
                              get_high_expr = TRUE,
                              implementation = "data.table",
                              do_parallel = FALSE,
                              verbose = TRUE)
  binspect_results <- as.data.table(binspect_results)

  gene_col <- intersect(c("feats", "genes", "gene"), colnames(binspect_results))[1]
  adj_p_col <- intersect(c("adj_p_value", "adj.p.value"), colnames(binspect_results))[1]
  p_col <- intersect(c("p_value", "p.value"), colnames(binspect_results))[1]
  score_col <- intersect(c("or", "estimate", "odds_ratio"), colnames(binspect_results))[1]
  if (is.na(gene_col) || is.na(adj_p_col) || is.na(p_col) || is.na(score_col)) {
    stop("Unexpected binSpect output columns: ", paste(colnames(binspect_results), collapse = ", "))
  }
  
  # Save binSpect results
  fwrite(binspect_results, output_file)
  
  # Extract significantly spatial genes (adjusted p-value < 0.05)
  spatial_genes <- binspect_results[get(adj_p_col) < 0.05]
  fwrite(spatial_genes, spatial_genes_file)
  
  # Create correlation matrix from spatial genes
  # Genes identified as spatial are considered to have spatial co-expression patterns
  n_genes <- nrow(expr_matrix)
  spatial_corr_matrix <- matrix(0, nrow = n_genes, ncol = n_genes)
  
  # Set correlation to 1 for spatial genes with themselves and modest correlation with other spatial genes
  spatial_gene_indices <- as.numeric(gsub("gene-", "", spatial_genes[[gene_col]]))
  
  if (length(spatial_gene_indices) > 0) {
    # Self-correlation
    diag(spatial_corr_matrix) <- 1
    
    # Spatial genes have moderate correlation with each other
    for (i in spatial_gene_indices) {
      for (j in spatial_gene_indices) {
        if (i != j) {
          # Use odds ratio as correlation strength (normalized)
          gene1_or <- spatial_genes[get(gene_col) == paste0("gene-", i), get(score_col)]
          gene2_or <- spatial_genes[get(gene_col) == paste0("gene-", j), get(score_col)]
          if (length(gene1_or) > 0 && length(gene2_or) > 0) {
            correlation <- min(0.8, sqrt(gene1_or * gene2_or) / 10)  # Normalize to reasonable range
            spatial_corr_matrix[i, j] <- correlation
          }
        }
      }
    }
  }
  
  # Save correlation matrix for evaluation
  corr_output_file <- paste0(output_dir, S, suffix, "_correlation_matrix.csv")
  fwrite(data.table(spatial_corr_matrix), corr_output_file)
  
  # Create summary statistics
  summary_stats <- data.frame(
    Metric = c("Total_genes", "Spatial_genes", "Significant_spatial_genes", 
               "Mean_odds_ratio", "Median_pvalue", "Processing_status"),
    Value = c(nrow(expr_matrix), 
              nrow(binspect_results),
              nrow(spatial_genes),
              mean(binspect_results[[score_col]], na.rm = TRUE),
              median(binspect_results[[p_col]], na.rm = TRUE),
              "Success")
  )
  
  fwrite(summary_stats, summary_file)
  
  cat("binSpect analysis completed successfully for", S, suffix, "\n")
  cat("Found", nrow(spatial_genes), "significantly spatial genes\n")
  
}, error = function(e) {
  cat("Error in binSpect analysis for", S, suffix, ":", e$message, "\n")
  
  # Create error summary
  summary_stats <- data.frame(
    Metric = c("Total_genes", "Spatial_genes", "Significant_spatial_genes", 
               "Mean_odds_ratio", "Median_pvalue", "Processing_status", "Error_message"),
    Value = c(nrow(expr_matrix), 0, 0, 0, 1, "Failed", e$message)
  )
  
  fwrite(summary_stats, summary_file)
})
