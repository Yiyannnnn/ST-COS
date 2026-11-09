library(data.table)
library(tidyverse)
library(Seurat)
library(spruce)
library(mclust)
library(ggplot2)
library(maple)


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
output_file <- paste0("./Data/", S, "_gene/MAPLE", suffix, ".csv")

# Load the data using the provided file paths
spot_count_layers <- fread(gene_expression_file)
coord <- fread(coord_file)
layer_coord <- fread(layer_file)

# Create Seurat object
data <- CreateSeuratObject(t(spot_count_layers), assay = "Spatial")
data <- RenameCells(data, new.names = 1:nrow(coord))
colnames(coord) <- c("row", "col")
data@images$slice1 <- new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "slice1_",
  coordinates = coord
)

# Perform SCTransform on the data
data.t <- SCTransform(data, assay = "Spatial", verbose = FALSE)

# Run PCA on the transformed data
d1 <- RunPCA(data.t)

# Methods requiring a fixed cluster count receive the simulated truth count.
n_clusters <- data.table::uniqueN(layer_coord$layer, na.rm = TRUE)
if (n_clusters < 2L) {
  stop("Expected at least two spatial domains, found ", n_clusters, ".")
}
message("Using ", n_clusters, " spatial domains from the simulation truth.")
d1_fit_PCs <- fit_spruce(d1, K = n_clusters, emb = "PCs", MCAR = FALSE, CAR = FALSE)

# Combine clustering results with the layer coordinates
results <- cbind(cluster = d1_fit_PCs$z, layer_coord)

# Save the results to the specified file
fwrite(results, output_file)


