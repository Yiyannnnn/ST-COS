library(BASS)
library(data.table)
library(tidyverse)

#!/usr/bin/env Rscript

# Load necessary libraries
library(BASS)
library(data.table)
library(tidyverse)
library(aricode)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure the correct number of arguments are provided
if (length(args) < 2 || length(args) > 3) {
  stop("Please provide two required arguments (S and suffix) and an optional third argument for the number of clusters.")
}

# Get the inputs from command-line arguments
S1 <- args[1]
suffix <- args[2]
target_R <- NA_integer_
if (length(args) == 3) {
  arg3 <- tolower(args[3])
  if (arg3 %in% c("auto", "true", "default")) {
    target_R <- NA_integer_
  } else {
    target_R <- suppressWarnings(as.integer(args[3]))
    if (is.na(target_R) || target_R <= 0) {
      stop("Third argument must be a positive integer, 'auto', or omitted.")
    }
  }
}

# Construct file paths based on the arguments
gene_expression_file <- paste0("./Data/", S1, "_gene/simulated/", S1, "_gene_expression", suffix, ".csv")
layer_file <- paste0("./Data/", S1, "_gene/simulated/", S1, "_layer", suffix, ".csv")
coord_file <- paste0("./Data/", S1, "_gene/simulated/", S1, "_coord", suffix, ".csv")
cluster_suffix <- if (is.na(target_R)) "" else paste0("_K", target_R)
output_file <- paste0("./Data/", S1, "_gene/BASS", cluster_suffix, suffix, ".csv")


 # Load the data using the provided file paths
spot_count_layers <- fread(gene_expression_file)
layer_coord <- fread(layer_file)
coord <- fread(coord_file)

#   - Number of cell types C: 3
#   - Number of spatial domains R: 3/4
C = 3
if (is.na(target_R)) {
  R = length(unique(layer_coord$layer))
} else {
  R = target_R
}
beta_method = "SW"
burnin = 10000
doLogNormalize = TRUE
doPCA = TRUE
scaleFeature = TRUE
nPC = 20
doBatchCorrect = FALSE
  
  
  # Prepare the data for BASS analysis
cnts <- list(t(spot_count_layers))
colnames(cnts[[1]]) <- paste0("spot_", 1:ncol(cnts[[1]]))
xy <- list(as.matrix(coord))
rownames(xy[[1]]) <- paste0("spot_", 1:ncol(cnts[[1]]))
  
  # Create BASS object
BASS <- createBASSObject(cnts, xy, C = C, R = R, beta_method = beta_method, burnin = burnin)
  
  # List all hyperparameters (for debugging or checking)
listAllHyper(BASS)
  
  # Preprocess, run, and postprocess the BASS object
BASS <- BASS.preprocess(BASS, doLogNormalize = doLogNormalize, doPCA = doPCA, scaleFeature = scaleFeature, nPC = nPC, doBatchCorrect = doBatchCorrect)
BASS <- BASS.run(BASS)
BASS <- BASS.postprocess(BASS)
  
  # Extract the layer information and compute the ARI
BASS_layer <- as.factor(unlist(BASS@results$z))

fwrite(cbind(layer_coord, BASS_layer), output_file)




# spot_count_layers <- fread("./Data/S1_gene/S1_gene_expression_1.csv")
# layer_coord <- fread("./Data/S1_gene/S1_layer_1.csv")
# coord <- fread("./Data/S1_gene/S1_coord_1.csv")


# cnts <- list(t(spot_count_layers))
# colnames(cnts[[1]]) <- paste0("spot_",1:ncol(cnts[[1]]))
# xy <- list(as.matrix(coord))
# rownames(xy[[1]]) <- paste0("spot_",1:ncol(cnts[[1]]))
# BASS <- createBASSObject(cnts, xy, C = 3, R = 3, beta_method = "SW",burnin = 10000)

# listAllHyper(BASS)
# BASS <- BASS.preprocess(BASS,doLogNormalize = TRUE, doPCA = TRUE, scaleFeature = TRUE, nPC = 20,doBatchCorrect = FALSE)
# BASS <- BASS.run(BASS)
# BASS <- BASS.postprocess(BASS)

# BASS_layer <- as.factor(unlist(BASS@results$z))

# library(aricode)
# ARI(BASS_layer, layer_coord$layer)

# fwrite(cbind(layer_coord, BASS_layer), "./Data/S1_gene/BASS_1.csv")
