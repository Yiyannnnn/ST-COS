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
sc_ref_file <- paste0("./Data/", S, "_gene/simulated/", S, "_scref", suffix, ".csv")

# Load spatial coordinates and gene expression data
pos <- fread(coord_file)
colnames(pos) <- c("x", "y")
cd <- round(t(fread(gene_expression_file)))
colnames(cd) <- paste0("spot_", seq_len(ncol(cd)))

# STdeconvolve is reference-free; the reference is read only to set K to the
# scenario's number of simulated cell types (3 for S1-S6, 4 for S7).
n_cell_types <- length(unique(fread(sc_ref_file, select = "Cell_type")$Cell_type))

# Remove pixels with too few genes
counts <- cleanCounts(cd, min.lib.size = 100)

# Feature select for genes
corpus <- restrictCorpus(counts, removeAbove = 2.0, removeBelow = 0.01)

# Choose optimal number of cell types (K)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = n_cell_types)

# Get best model results
optLDA <- optimalModel(models = ldas, opt = "min")

# Extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta

# Save the deconvolved topic proportions with spot identifiers; topics are
# matched to cell types during evaluation.
fwrite(data.table(spot = rownames(deconProp), as.data.frame(deconProp)), output_file)

