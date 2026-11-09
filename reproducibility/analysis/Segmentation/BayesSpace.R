library(BayesSpace)
library(SingleCellExperiment)
library(data.table)
library(tidyverse)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure the correct number of arguments are provided
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript BayesSpace.R <dataset> <suffix> [q]\n",
       "Example: Rscript BayesSpace.R S7 _1 4")
}

# Get the inputs from command-line arguments
S <- args[1]
suffix <- args[2]
if (length(args) == 3) {
  q_input <- suppressWarnings(as.integer(args[3]))
  if (is.na(q_input) || q_input <= 0) {
    stop("Third argument (q) must be a positive integer.")
  }
} else {
  q_input <- NULL
}

# Construct file paths based on the arguments
gene_expression_file <- paste0("./Data/", S, "_gene/simulated/", S, "_gene_expression", suffix, ".csv")
layer_file <- paste0("./Data/", S, "_gene/simulated/", S, "_layer", suffix, ".csv")
coord_file <- paste0("./Data/", S, "_gene/simulated/", S, "_coord", suffix, ".csv")
output_file <- paste0("./Data/", S, "_gene/BayesSpace", suffix, ".csv")

# Load the data using the provided file paths
spot_count_layers <- fread(gene_expression_file)
layer_coord <- fread(layer_file)
coord <- fread(coord_file)

coord_names <- tolower(colnames(coord))
x_idx <- which(coord_names %in% c("x", "x_coord", "row"))
y_idx <- which(coord_names %in% c("y", "y_coord", "col"))
if (length(x_idx) != 1L || length(y_idx) != 1L) {
  stop("Unable to identify x/y coordinate columns in: ", coord_file)
}

valid_idx <- !is.na(coord[[x_idx]]) & !is.na(coord[[y_idx]])
if (!all(valid_idx)) {
  spot_count_layers <- spot_count_layers[valid_idx, , drop = FALSE]
  layer_coord <- layer_coord[valid_idx, , drop = FALSE]
  coord <- coord[valid_idx, , drop = FALSE]
}

counts_mat <- t(as.matrix(spot_count_layers))
rownames(counts_mat) <- if (!is.null(colnames(spot_count_layers))) {
  colnames(spot_count_layers)
} else {
  paste0("gene_", seq_len(nrow(counts_mat)))
}
rownames(counts_mat) <- make.unique(rownames(counts_mat))
colnames(counts_mat) <- paste0("spot_", seq_len(ncol(counts_mat)))

# Create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = counts_mat))
if (ncol(sce) != nrow(layer_coord) || ncol(sce) != nrow(coord)) {
  stop("Spot counts, layer annotations, and coordinates have mismatched dimensions after filtering.")
}
colData(sce)$layer <- as.factor(layer_coord$layer)

colData(sce)$row <- coord[[x_idx]]
colData(sce)$col <- coord[[y_idx]]
rownames(colData(sce)) <- colnames(sce)

colData(sce)$array_row <- as.numeric(coord[[y_idx]])
colData(sce)$array_col <- as.numeric(coord[[x_idx]])
colData(sce)$spot.idx <- seq_len(ncol(sce))

# Methods requiring a fixed cluster count receive the simulated truth count
# unless an explicit value is supplied for a separate sensitivity analysis.
q <- if (is.null(q_input)) data.table::uniqueN(layer_coord$layer, na.rm = TRUE) else q_input
if (q < 2L) {
  stop("Expected at least two spatial domains, found ", q, ".")
}
message("Using ", q, " spatial domains from the simulation truth.")
d <- 15  # Number of PCs

# Preprocess spatial data
sce_BayesSpace <- spatialPreprocess(sce, skip.PCA = FALSE, platform = "Visium", log.normalize = TRUE)

# Perform spatial clustering
sce_BayesSpace <- spatialCluster(sce_BayesSpace, q = q, d = d, platform = "Visium", nrep = 50000, gamma = 3, save.chain = TRUE, model = "t")

# Extract BayesSpace labels
BayesSpace_labels <- as.factor(sce_BayesSpace$spatial.cluster)

# Optionally compute ARI (uncomment if needed)
# ari_value <- ARI(BayesSpace_labels, layer_coord$layer)
# print(paste("Adjusted Rand Index (ARI):", ari_value))

# Save the output to the specified file
fwrite(cbind(layer_coord, BayesSpace_labels), output_file)
