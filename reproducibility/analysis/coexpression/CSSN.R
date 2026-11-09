#!/usr/bin/env Rscript

# CSSN: Cell-Specific Spatial Network Analysis Script
# Based on the CSSN GitHub tutorial: https://github.com/jingeyu/CSSN
#
# This script runs CSSN analysis on simulated spatial transcriptomics data
# to recover spatially-varying cell-specific gene co-expression networks
#
# Usage: Rscript CSSN.R S1 _1

# Set up R library path
.libPaths(c(Sys.getenv("R_LIBS_USER", unset = .libPaths()[1]), .libPaths()))

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(pheatmap)
})

# Install CSSN if not already installed
if (!require("CSSN", quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cran.r-project.org/")
  }
  devtools::install_github("jingeyu/CSSN")
  library(CSSN)
} else {
  library(CSSN)
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure the correct number of arguments are provided
if (length(args) != 2) {
  stop("Please provide the two arguments: S and suffix (_N).")
}

# Get the inputs from command-line arguments
S <- args[1]
suffix <- args[2]

# Create output directory
output_dir <- "./coexpression/CSSN"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Construct file paths based on the arguments
gene_expression_file <- paste0("./Data/", S, "_gene/simulated/", S, "_gene_expression", suffix, ".csv")
coord_file <- paste0("./Data/", S, "_gene/simulated/", S, "_coord", suffix, ".csv")
layer_file <- paste0("./Data/", S, "_gene/simulated/", S, "_layer", suffix, ".csv")

# Check if input files exist
for (file_path in c(gene_expression_file, coord_file, layer_file)) {
  if (!file.exists(file_path)) {
    stop(paste("Input file not found:", file_path))
  }
}

cat("Starting CSSN analysis for", paste0(S, suffix), "\n")
cat("Output directory:", output_dir, "\n")

# Load the data
cat("Loading data...\n")
spot_count_layers <- fread(gene_expression_file)
coord <- fread(coord_file)
layer_coord <- fread(layer_file)

# Print data dimensions for debugging
cat("Data dimensions:\n")
cat("Gene expression:", dim(spot_count_layers), "\n")
cat("Coordinates:", dim(coord), "\n")  
cat("Layers:", dim(layer_coord), "\n")

# Prepare data for CSSN
# CSSN expects: X (gene x cell matrix), cell.info (cell metadata with coordinates)
expr_matrix <- as.matrix(spot_count_layers)
n_coord_spots <- nrow(coord)

# Determine orientation using the coordinate table rather than guessing from shape.
# Simulated inputs are spots x genes, but keep a fallback for legacy files.
if (nrow(expr_matrix) == n_coord_spots) {
  # Rows are spots, columns are genes.
  X <- t(expr_matrix)
  spot_names <- rownames(expr_matrix)
  if (is.null(spot_names) || anyNA(spot_names) || any(spot_names == "")) {
    spot_names <- paste0("spot_", seq_len(ncol(X)))
  }
  gene_names <- colnames(expr_matrix)
} else if (ncol(expr_matrix) == n_coord_spots) {
  # Columns are spots, rows are genes.
  X <- expr_matrix
  spot_names <- colnames(expr_matrix)
  if (is.null(spot_names) || anyNA(spot_names) || any(spot_names == "")) {
    spot_names <- paste0("spot_", seq_len(ncol(X)))
  }
  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names) || anyNA(gene_names) || any(gene_names == "")) {
    gene_names <- paste0("gene_", seq_len(nrow(X)))
  }
} else {
  stop(
    sprintf(
      "Cannot align expression matrix to coordinates: expression=%s, coordinates=%d",
      paste(dim(expr_matrix), collapse = "x"),
      n_coord_spots
    )
  )
}

rownames(X) <- gene_names
colnames(X) <- spot_names

cat("Expression matrix dimensions (genes x cells):", dim(X), "\n")

# Prepare cell information (cell.info format: cell_type, X_coord, Y_coord)
n_spots <- ncol(X)
n_genes <- nrow(X)

# Use coordinates (assuming first two columns are x, y coordinates)
coords_matrix <- as.matrix(coord[1:n_spots, 1:2, drop = FALSE])
x_coords <- coords_matrix[, 1]
y_coords <- coords_matrix[, 2]

# Use layer information as cell types
if (nrow(layer_coord) >= n_spots) {
  if ("layer" %in% colnames(layer_coord)) {
    cell_types <- layer_coord$layer[1:n_spots]
  } else {
    cell_types <- layer_coord[[1]][1:n_spots]  # Use first column
  }
} else {
  cell_types <- rep(1, n_spots)  # Default to single cell type
}

# Create cell.info matrix (required format for CSSN)
cell.info <- data.frame(
  cell_type = cell_types,
  X = x_coords,
  Y = y_coords
)

cat("Cell info dimensions:", dim(cell.info), "\n")
cat("Number of cell types:", length(unique(cell.info$cell_type)), "\n")
cat("Cell type distribution:", table(cell.info$cell_type), "\n")

# Quality control - filter low expression genes and cells
cat("Performing quality control...\n")

# Filter genes with very low expression
gene_totals <- rowSums(X)
valid_genes <- gene_totals >= 50  # Keep genes with at least 50 total counts
X_filtered <- X[valid_genes, ]
gene_names_filtered <- gene_names[valid_genes]

# Filter cells with very low expression
cell_totals <- colSums(X_filtered)
valid_cells <- cell_totals >= 100  # Keep cells with at least 100 total counts
X_filtered <- X_filtered[, valid_cells]
cell.info_filtered <- cell.info[valid_cells, ]

cat("Filtered data dimensions (genes x cells):", dim(X_filtered), "\n")
cat("Genes filtered:", sum(valid_genes), "out of", length(valid_genes), "\n")
cat("Cells filtered:", sum(valid_cells), "out of", length(valid_cells), "\n")

# Update dimensions
G <- nrow(X_filtered)  # gene number
n <- ncol(X_filtered)  # cell number

cat("Final data for CSSN: G =", G, "genes, n =", n, "cells\n")

# Visualize spatial pattern
cat("Creating spatial visualization...\n")

# Set up colors for cell types
unique_types <- unique(cell.info_filtered$cell_type)
n_types <- length(unique_types)

if (n_types <= 2) {
  # Use predefined colors for 1-2 cell types
  pal <- c(rgb(221, 160, 221, maxColorValue = 255),
           rgb(0, 206, 209, maxColorValue = 255))
  pal <- pal[1:n_types]
} else {
  # Generate colors for more cell types
  pal <- rainbow(n_types)
}
names(pal) <- as.character(unique_types)

# Create spatial plot
tryCatch({
  gg <- ggplot(cell.info_filtered, aes(x = X, y = Y, col = as.factor(cell_type), shape = as.factor(cell_type)))
  pl <- gg + geom_point(size = 2.5) +
    scale_color_manual(values = pal) +
    theme_bw() +
    theme(legend.text = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold")) +
    labs(x = "X Coordinate", y = "Y Coordinate", title = paste0("Spatial Pattern - ", S, suffix)) +
    guides(color = guide_legend(title = "Cell Type",
                               title.theme = element_text(size = 16),
                               override.aes = list(size = 4)),
           shape = guide_legend(title = "Cell Type",
                               title.theme = element_text(size = 16),
                               override.aes = list(size = 4)))
  
  ggsave(file.path(output_dir, paste0(S, suffix, "_cell_spatial_pattern.pdf")), pl, width = 10, height = 8)
  cat("Spatial pattern plot saved\n")
}, error = function(e) {
  cat("Warning: Could not create spatial pattern plot:", e$message, "\n")
})

# Run CSSN Estimation
cat("Running CSSN estimation...\n")

tryCatch({
  # Set parameters for CSSN
  # nu: penalty parameter (recommended: 2*G for each cell)
  nu <- rep(2 * G, n)
  
  # d: bandwidth parameter for spatial kernel (adjust based on spatial resolution)
  # Start with a reasonable default and adjust if needed
  d_param <- 0.1
  
  # m.info: number of informative features (genes)
  m_info <- min(70, G)  # Use at most 70 genes or all available genes
  
  cat("CSSN parameters:\n")
  cat("- nu (penalty):", length(nu), "values (2*G per cell)\n")
  cat("- d (bandwidth):", d_param, "\n")
  cat("- m.info (informative genes):", m_info, "\n")
  
  # Run CSSN estimation for all cells
  cat("Running CSSNEst for all cells...\n")
  Result_all <- CSSNEst(X_filtered, cell.info_filtered, 
                        nu = nu, d = d_param, m.info = m_info, 
                        is.scale = TRUE, is.all = TRUE)
  
  cat("CSSN estimation completed for", length(Result_all), "cells\n")
  
  # Save the complete results
  saveRDS(Result_all, file.path(output_dir, paste0(S, suffix, "_CSSN_all_results.rds")))
  
  # Run CSSN estimation for a subset of cells with correlation output
  if (n >= 5) {
    # Select a subset of cells for detailed analysis
    indx.cell <- c(1, min(3, n), min(7, n), min(10, n), min(15, n))
    indx.cell <- indx.cell[indx.cell <= n]  # Remove indices beyond available cells
    indx.cell <- unique(indx.cell)  # Remove duplicates
    
    cat("Running CSSNEst for subset of cells:", paste(indx.cell, collapse = ", "), "\n")
    
    result_subset <- CSSNEst(X_filtered, cell.info_filtered, 
                            nu = nu, d = d_param, m.info = m_info, 
                            is.scale = TRUE, is.all = FALSE, 
                            indx.cell = indx.cell, output.corr = TRUE)
    
    # Save subset results
    saveRDS(result_subset, file.path(output_dir, paste0(S, suffix, "_CSSN_subset_results.rds")))
  }
  
}, error = function(e) {
  cat("Error in CSSN estimation:", e$message, "\n")
  cat("This might be due to data characteristics or parameter settings.\n")
  
  # Try with adjusted parameters
  cat("Trying with adjusted parameters...\n")
  tryCatch({
    # More lenient parameters
    nu_adj <- rep(G, n)  # Reduce penalty
    d_adj <- 0.2  # Increase bandwidth
    m_info_adj <- min(50, G)  # Reduce number of informative genes
    
    cat("Adjusted CSSN parameters:\n")
    cat("- nu (penalty):", length(nu_adj), "values (G per cell)\n")
    cat("- d (bandwidth):", d_adj, "\n")
    cat("- m.info (informative genes):", m_info_adj, "\n")
    
    Result_all <- CSSNEst(X_filtered, cell.info_filtered, 
                          nu = nu_adj, d = d_adj, m.info = m_info_adj, 
                          is.scale = TRUE, is.all = TRUE)
    
    cat("CSSN estimation completed with adjusted parameters for", length(Result_all), "cells\n")
    saveRDS(Result_all, file.path(output_dir, paste0(S, suffix, "_CSSN_all_results.rds")))
    
  }, error = function(e2) {
    cat("Error with adjusted parameters:", e2$message, "\n")
    Result_all <- NULL
  })
})

# Generate heatmaps for the first few networks (if estimation succeeded)
if (exists("Result_all") && !is.null(Result_all) && length(Result_all) > 0) {
  cat("Generating co-expression network heatmaps and saving networks...\n")
  
  # Save individual co-expression networks
  cat("Saving individual co-expression networks...\n")
  networks_dir <- file.path(output_dir, paste0(S, suffix, "_networks"))
  dir.create(networks_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save networks in multiple formats
  for (i in 1:length(Result_all)) {
    if (!is.null(Result_all[[i]]) && is.matrix(Result_all[[i]])) {
      network_matrix <- Result_all[[i]]
      
      # Save as CSV
      write.csv(network_matrix, 
                file.path(networks_dir, paste0("cell_", i, "_network.csv")), 
                row.names = TRUE)
      
      # Save as RDS (R format)
      saveRDS(network_matrix, 
              file.path(networks_dir, paste0("cell_", i, "_network.rds")))
      
      # Save edge list format (for network analysis tools)
      if (nrow(network_matrix) == ncol(network_matrix)) {
        # Create edge list from adjacency matrix
        edge_list <- data.frame(
          Gene1 = character(),
          Gene2 = character(),
          Correlation = numeric(),
          Cell_ID = integer(),
          stringsAsFactors = FALSE
        )
        
        # Extract upper triangle (avoid duplicates)
        for (row in 1:(nrow(network_matrix)-1)) {
          for (col in (row+1):ncol(network_matrix)) {
            correlation <- network_matrix[row, col]
            if (!is.na(correlation) && abs(correlation) > 0.01) {  # Only save significant correlations
              edge_list <- rbind(edge_list, data.frame(
                Gene1 = paste0("Gene_", row),
                Gene2 = paste0("Gene_", col),
                Correlation = correlation,
                Cell_ID = i
              ))
            }
          }
        }
        
        # Save edge list
        if (nrow(edge_list) > 0) {
          write.csv(edge_list, 
                    file.path(networks_dir, paste0("cell_", i, "_edges.csv")), 
                    row.names = FALSE)
        }
      }
    }
  }
  
  cat("Individual networks saved for", length(Result_all), "cells\n")
  
  # Create combined network summary
  cat("Creating combined network analysis...\n")
  
  # Aggregate networks - calculate mean correlation across all cells
  if (length(Result_all) > 1) {
    # Get dimensions from first valid network
    first_valid <- NULL
    for (i in 1:length(Result_all)) {
      if (!is.null(Result_all[[i]]) && is.matrix(Result_all[[i]])) {
        first_valid <- Result_all[[i]]
        break
      }
    }
    
    if (!is.null(first_valid)) {
      n_genes <- nrow(first_valid)
      mean_network <- array(0, dim = c(n_genes, n_genes))
      valid_count <- array(0, dim = c(n_genes, n_genes))
      
      # Sum all valid networks
      for (i in 1:length(Result_all)) {
        if (!is.null(Result_all[[i]]) && is.matrix(Result_all[[i]]) && 
            nrow(Result_all[[i]]) == n_genes && ncol(Result_all[[i]]) == n_genes) {
          mean_network <- mean_network + Result_all[[i]]
          valid_count <- valid_count + 1
        }
      }
      
      # Calculate mean (avoid division by zero)
      mean_network[valid_count > 0] <- mean_network[valid_count > 0] / valid_count[valid_count > 0]
      
      # Save mean network
      write.csv(mean_network, 
                file.path(output_dir, paste0(S, suffix, "_mean_coexpression_network.csv")), 
                row.names = TRUE)
      saveRDS(mean_network, 
              file.path(output_dir, paste0(S, suffix, "_mean_coexpression_network.rds")))
      
      cat("Mean co-expression network saved\n")
      
      # Create mean network heatmap
      tryCatch({
        colors_func <- colorRampPalette(c('blue', 'white', 'red'))
        colors <- colors_func(100)
        
        pdf(file.path(output_dir, paste0(S, suffix, "_mean_network_heatmap.pdf")), width = 8, height = 7)
        pheatmap(mean_network,
                 color = colors,
                 cluster_cols = TRUE, 
                 cluster_rows = TRUE,
                 show_rownames = FALSE, 
                 show_colnames = FALSE,
                 main = paste0("Mean Co-expression Network - ", S, suffix),
                 breaks = seq(-1, 1, length.out = 101))
        dev.off()
        
        cat("Mean network heatmap saved\n")
      }, error = function(e) {
        cat("Warning: Could not create mean network heatmap:", e$message, "\n")
      })
    }
  }
  
  # Create network connectivity analysis
  cat("Analyzing network connectivity...\n")
  
  connectivity_stats <- data.frame(
    Cell_ID = integer(),
    Total_Edges = integer(),
    Strong_Positive_Edges = integer(),
    Strong_Negative_Edges = integer(),
    Average_Connectivity = numeric(),
    Max_Connectivity = numeric(),
    Network_Density = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:min(length(Result_all), 50)) {  # Analyze first 50 cells
    if (!is.null(Result_all[[i]]) && is.matrix(Result_all[[i]])) {
      net <- Result_all[[i]]
      n_genes <- nrow(net)
      
      # Calculate connectivity metrics
      total_edges <- sum(abs(net[upper.tri(net)]) > 0.01)
      strong_pos <- sum(net[upper.tri(net)] > 0.5)
      strong_neg <- sum(net[upper.tri(net)] < -0.5)
      
      # Node connectivity (sum of absolute correlations for each gene)
      node_connectivity <- rowSums(abs(net)) - diag(abs(net))  # Exclude self-connections
      avg_connectivity <- mean(node_connectivity)
      max_connectivity <- max(node_connectivity)
      
      # Network density
      max_possible_edges <- n_genes * (n_genes - 1) / 2
      density <- total_edges / max_possible_edges
      
      connectivity_stats <- rbind(connectivity_stats, data.frame(
        Cell_ID = i,
        Total_Edges = total_edges,
        Strong_Positive_Edges = strong_pos,
        Strong_Negative_Edges = strong_neg,
        Average_Connectivity = avg_connectivity,
        Max_Connectivity = max_connectivity,
        Network_Density = density
      ))
    }
  }
  
  # Save connectivity analysis
  write.csv(connectivity_stats, 
            file.path(output_dir, paste0(S, suffix, "_network_connectivity.csv")), 
            row.names = FALSE)
  
  cat("Network connectivity analysis saved\n")
  
  tryCatch({
    # ...existing heatmap code...
    colors_func <- colorRampPalette(c('white', "red"))
    colors <- colors_func(100)
    
    for (i in 1:min(5, length(Result_all))) {
      if (!is.null(Result_all[[i]]) && is.matrix(Result_all[[i]])) {
        filename <- file.path(output_dir, paste0(S, suffix, "_Est_network_cell_", i, ".pdf"))
        
        # Create heatmap
        pdf(filename, width = 6, height = 5)
        pheatmap(abs(Result_all[[i]]),  # Use absolute values for better visualization
                 color = colors,
                 cluster_cols = FALSE, 
                 cluster_rows = FALSE,
                 show_rownames = FALSE, 
                 show_colnames = FALSE,
                 main = paste0("Cell ", i, " Co-expression Network"))
        dev.off()
        
        cat("Heatmap saved for cell", i, "\n")
      }
    }
  }, error = function(e) {
    cat("Warning: Could not create heatmaps:", e$message, "\n")
  })
  
  # ...existing network analysis summary code...
}

# Run CSSN Prediction (if estimation succeeded)
if (exists("Result_all") && !is.null(Result_all) && length(Result_all) > 0) {
  cat("Running CSSN prediction for missing locations...\n")
  
  tryCatch({
    # Generate random missing locations for prediction
    set.seed(123)
    miss.num <- 5
    miss.x <- runif(miss.num, min(cell.info_filtered$X), max(cell.info_filtered$X))
    miss.y <- runif(miss.num, min(cell.info_filtered$Y), max(cell.info_filtered$Y))
    miss.indx <- cbind(miss.x, miss.y)
    
    cat("Predicting networks at", miss.num, "missing locations\n")
    
    # Run prediction
    pre <- CSSNPredict(Result_all, cell.info_filtered, miss.indx)
    
    # Save prediction results
    saveRDS(pre, file.path(output_dir, paste0(S, suffix, "_CSSN_predictions.rds")))
    
    # Save missing locations
    missing_locations <- data.frame(
      Location_ID = 1:miss.num,
      X = miss.x,
      Y = miss.y
    )
    write.csv(missing_locations, file.path(output_dir, paste0(S, suffix, "_prediction_locations.csv")), row.names = FALSE)
    
    cat("Prediction completed and saved\n")
    
  }, error = function(e) {
    cat("Warning: Could not run prediction:", e$message, "\n")
  })
}

# Create summary report
cat("Creating summary report...\n")

summary_stats <- data.frame(
  Metric = c("Total_spots", "Total_genes", "Filtered_spots", "Filtered_genes", 
             "Cell_types", "Networks_estimated", "Predictions_made"),
  Value = c(
    ncol(X),
    nrow(X),
    if (exists("cell.info_filtered")) nrow(cell.info_filtered) else 0,
    if (exists("X_filtered")) nrow(X_filtered) else 0,
    if (exists("cell.info_filtered")) length(unique(cell.info_filtered$cell_type)) else 0,
    if (exists("Result_all") && !is.null(Result_all)) length(Result_all) else 0,
    if (exists("pre")) length(pre) else 0
  )
)

write.csv(summary_stats, file.path(output_dir, paste0(S, suffix, "_summary_stats.csv")), row.names = FALSE)

cat("CSSN analysis completed! All results saved to:", output_dir, "\n")
cat("Summary:\n")
print(summary_stats)
