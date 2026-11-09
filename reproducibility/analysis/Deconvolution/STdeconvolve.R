library(STdeconvolve)
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
coord_file <- paste0("./Data/", S, "_gene/simulated/", S, "_coord", suffix, ".csv")
gene_expression_file <- paste0("./Data/", S, "_gene/simulated/", S, "_gene_expression", suffix, ".csv")
output_file <- paste0("./Data/", S, "_gene/STdeconvolve", suffix, ".csv")

# Load spatial coordinates and gene expression data
pos <- fread(coord_file)
colnames(pos) <- c("x", "y")
cd <- round(t(fread(gene_expression_file)))

# Remove pixels with too few genes
counts <- cleanCounts(cd, min.lib.size = 100)

# Feature select for genes
corpus <- restrictCorpus(counts, removeAbove = 2.0, removeBelow = 0.01)

# Choose optimal number of cell types (K)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(3, 3, by = 1))

# Get best model results
optLDA <- optimalModel(models = ldas, opt = "min")

# Extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta

# Save the deconvolved cell-type proportions
fwrite(deconProp, output_file)

