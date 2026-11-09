library(Redeconve)
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
sc_ref_file <- paste0("./Data/", S, "_gene/simulated/", S, "_scref", suffix, ".csv")
gene_expression_file <- paste0("./Data/", S, "_gene/simulated/", S, "_gene_expression", suffix, ".csv")
output_file <- paste0("./Data/", S, "_gene/Redeconve", suffix, ".csv")

# Read in single-cell and spatial transcriptomics data
sc <- fread(sc_ref_file)
ref <- t(sc[,-1])
colnames(ref) <- as.vector(sc[,"Cell_type"])$Cell_type

st <- fread(gene_expression_file)
st <- t(st)
colnames(st) <- paste0("spot_", 1:ncol(st))

# Redeconve parallelizes the independent spot-level optimization when the
# scheduler assigns more than one CPU to a task.
ncores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (!is.finite(ncores) || ncores < 1L) ncores <- 1L
res <- deconvoluting(
  ref, st,
  genemode = "def", hpmode = "def",
  dopar = ncores > 1L, ncores = ncores
)

# Aggregate the cell-level coefficients by their actual cell-type labels.
# Capture labels before data.frame conversion makes duplicate row names unique.
res.prop.matrix <- to.proportion(res)
cell.type.labels <- rownames(res.prop.matrix)
res.prop <- as.data.frame(res.prop.matrix)
result_by_type <- aggregate(
  res.prop,
  by = list(cell_type = cell.type.labels),
  FUN = sum
)
result <- t(result_by_type[, -1, drop = FALSE])
colnames(result) <- result_by_type$cell_type

# Save the results
fwrite(result, output_file)
