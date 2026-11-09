# module load r/4.1.0
library(SpaceX)
library(data.table)
library(tidyverse)

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
output_file <- paste0("./Data/", S, "_gene/SpaceX", suffix, ".csv")

# Load the data using the provided file paths
spot_count_layers <- abs(fread(gene_expression_file))
pos <- fread(coord_file)
layer_factor <- fread(layer_file)

# Run SpaceX
BC_fit <- SpaceX(spot_count_layers, pos, layer_factor[, 3], sPMM = FALSE, Post_process = TRUE, numCore = 1)

# Save the results
fwrite(BC_fit$Shared_network, output_file)

# pdf("test.pdf")
# corrplot(BC_fit$Cluster_network[,,1])
# corrplot(BC_fit$Cluster_network[,,2])
# corrplot(BC_fit$Cluster_network[,,3])
# dev.off()



