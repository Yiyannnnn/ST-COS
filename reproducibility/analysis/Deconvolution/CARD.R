#r/4.3.1
#install.packages('RcppML',INSTALL_opts = '--no-lock')      
#devtools::install_github('YMa-lab/CARD')
#BiocManager::install("TOAST")
#devtools::install_github('xuranw/MuSiC')
library(CARD)
library(data.table)
library(tidyverse)
library(Matrix)
library(SummarizedExperiment)

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
coord_file <- paste0("./Data/", S, "_gene/simulated/", S, "_coord", suffix, ".csv")
output_file <- paste0("./Data/", S, "_gene/CARD", suffix, ".csv")

# Load the data using the provided file paths
sc_count <- fread(sc_ref_file)
sc_meta <- cbind(sc_count[,1])
rownames(sc_meta) <- paste0("cell_", 1:nrow(sc_meta))
sc_count <- as.matrix(t(sc_count[,-1]))
colnames(sc_count) <- paste0("cell_", 1:nrow(sc_meta))
rownames(sc_count) <- paste0("gene_", 1:nrow(sc_count))

spatial_count <- fread(gene_expression_file)
spatial_count <- as.matrix(t(spatial_count))
rownames(spatial_count) <- paste0("gene_", 1:nrow(spatial_count))
colnames(spatial_count) <- paste0("spot_", 1:ncol(spatial_count))

spatial_location <- fread(coord_file)
spatial_location <- as.data.frame(spatial_location)
rownames(spatial_location) <- paste0("spot_", 1:nrow(spatial_location))
colnames(spatial_location) <- c("x", "y")

# Create CARD object
CARD_object = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count*0.01, #colSums(spatial_countMat) <= 1e6)
  spatial_location = spatial_location,
  ct.varname = "Cell_type",
  ct.select = unique(sc_meta$Cell_type),
  sample.varname = NULL,
  minCountGene = 0,
  minCountSpot = 0
)

# Relaxed informative-gene selection to avoid empty sets on simulated data.
relaxed_selectInfo <- function(Basis, sc_eset, commonGene, ct.select, ct.varname) {
  gene1 <- lapply(ct.select, function(ict) {
    rest <- rowMeans(Basis[, colnames(Basis) != ict, drop = FALSE])
    FC <- log((Basis[, ict] + 1e-06)) - log((rest + 1e-06))
    rownames(Basis)[FC > 0 & Basis[, ict] > 0]
  })
  gene1 <- unique(unlist(gene1))
  gene1 <- intersect(gene1, commonGene)
  if (length(gene1) == 0) {
    return(commonGene)
  }
  counts <- assays(sc_eset)$counts
  counts <- counts[rownames(counts) %in% gene1, , drop = FALSE]
  if (nrow(counts) == 0) {
    return(commonGene)
  }
  ct.select <- names(table(colData(sc_eset)[, ct.varname]))[table(colData(sc_eset)[, ct.varname]) > 1]
  sd_within <- sapply(ct.select, function(ict) {
    temp <- counts[, colData(sc_eset)[, ct.varname] == ict, drop = FALSE]
    apply(temp, 1, var) / apply(temp, 1, mean)
  })
  if (is.null(dim(sd_within)) || nrow(sd_within) == 0) {
    return(rownames(counts))
  }
  gene2 <- rownames(sd_within)[apply(sd_within, 1, mean, na.rm = TRUE) <
                                 quantile(apply(sd_within, 1, mean, na.rm = TRUE), prob = 0.99, na.rm = TRUE)]
  if (length(gene2) == 0) commonGene else gene2
}
assignInNamespace("selectInfo", relaxed_selectInfo, ns = "CARD")

# Run CARD deconvolution
CARD_object = CARD_deconvolution(CARD_object)

# Save the results
fwrite(cbind(spot = rownames(CARD_object@Proportion_CARD), CARD_object@Proportion_CARD), file = paste0("./Data/", S, "_gene/CARD", suffix, ".csv"))

# pdf("test.pdf")
# p1 <- CARD.visualize.pie(
# 	proportion = CARD_object@Proportion_CARD,
# 	spatial_location = CARD_object@spatial_location, 
#  	colors = c("red", "blue", "green"), 
#   	radius = 1) ### You can choose radius = NULL or your own radius number
# print(p1)
# dev.off()
