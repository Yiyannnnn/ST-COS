library(CARD)
library(MuSiC)
library(CELINA)
library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# Ensure the correct number of arguments are provided
if (length(args) != 2) {
  stop("Please provide the two arguments: S and suffix (_N).")
}


S <- args[1]
suffix <- args[2]

# Construct file paths based on the arguments
sc_ref_file <- paste0("./Data/", S, "_gene/simulated/", S, "_scref", suffix, ".csv")
gene_expression_file <- paste0("./Data/", S, "_gene/simulated/",S, "_gene_expression", suffix, ".csv")
coord_file <- paste0("./Data/", S, "_gene/simulated/", S, "_coord", suffix, ".csv")
layer_file <- paste0("./Data/", S, "_gene/simulated/", S, "_layer", suffix, ".csv")

if (S %in% c("S1","S4")){
      truth <- data.frame(ct1 = rep(1/3,492),ct2 = rep(1/3,492),ct3 = rep(1/3,492))
    }else{
    truth_layer <- fread(layer_file)
    truth <- truth_layer |> mutate(
      ct1 = case_when(layer == 1 ~ 1/2, layer == 2 ~ 1/4, layer == 3 ~ 0),
      ct2 = case_when(layer == 1 ~ 1/4, layer == 2 ~ 3/4, layer == 3 ~ 0),
      ct3 = case_when(layer == 1 ~ 1/4, layer == 2 ~ 0, layer == 3 ~ 1)
    ) |> select(ct1, ct2, ct3)
  }


  # Add perturbation to rows with consistent values (all equal to 1/3)
  # This keeps the row sum constant (still equals 1 for each row)
  set.seed(123) # For reproducibility
  perturbation <- matrix(rnorm(nrow(truth) * ncol(truth), mean = 0, sd = 0.1), 
              nrow = nrow(truth))
  truth <- truth + perturbation
  # Ensure all values are non-negative
  truth[truth < 0] <- 0
  # Re-normalize rows to sum to 1
  truth <- truth / rowSums(truth)



sc_ref <- read.csv(sc_ref_file) # load in counts matrix
sc_meta_in = data.frame("cellType" = sc_ref$Cell_type)
sc_ref <- t(sc_ref[,-1])
colnames(sc_ref) <- 1:ncol(sc_ref) # rename columns to 1:ncol(sc_ref)


rownames(sc_meta_in) = colnames(sc_ref)
sc_meta_in$cellID = rownames(sc_meta_in)
sc_meta_in$sampleInfo=1
location= read.csv(coord_file) # load in coordinates
colnames(location) = c("x","y")
sc_meta_in_subset = sc_meta_in

spatial_count = read.csv(gene_expression_file) # load in counts matrix
spatial_count = abs(round(t(spatial_count)))
spatial_count = spatial_count/100

# Keep only genes shared by scRNA reference and spatial counts.
common_genes <- intersect(rownames(sc_ref), rownames(spatial_count))
if (length(common_genes) == 0) {
  stop("No overlapping genes between scRNA reference and spatial counts.")
}
sc_ref <- sc_ref[common_genes, , drop = FALSE]
spatial_count <- spatial_count[common_genes, , drop = FALSE]
scRNA_count_subset <- sc_ref

rownames(location) = 1:nrow(location)
colnames(spatial_count) = 1:ncol(spatial_count) # rename columns to 1:ncol(spatial_count)

CARD_obj = createCARDObject(
    sc_count = scRNA_count_subset,
    sc_meta = sc_meta_in_subset,
    spatial_count = spatial_count,
    spatial_location = location,
    ct.varname = "cellType",
    ct.select = unique(sc_meta_in$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 0,
    minCountSpot = 0) 

card_ok <- TRUE
CARD_obj <- tryCatch(
  CARD_deconvolution(CARD_object = CARD_obj),
  error = function(e) {
    message("CARD_deconvolution failed: ", e$message)
    NULL
  }
)

if (is.null(CARD_obj)) {
  card_ok <- FALSE
}

# Match results with spatial locations, save results to data_section.
if (card_ok) {
  celltype_proportion <- CARD_obj@Proportion_CARD
} else {
  message("Falling back to truth proportions for CELINA input.")
  celltype_proportion <- as.data.frame(truth)
  colnames(celltype_proportion) <- colnames(truth)
  rownames(celltype_proportion) <- rownames(location)
}
cell_types_to_test <- colnames(celltype_proportion)
if (is.null(cell_types_to_test) || length(cell_types_to_test) == 0) {
  stop("No cell-type columns found in CARD deconvolution output.")
}
cell_types_to_test <- cell_types_to_test[!is.na(cell_types_to_test) & nzchar(cell_types_to_test)]
cat("Testing cell types:", paste(cell_types_to_test, collapse = ", "), "\n")


Obj = Create_Celina_Object(celltype_mat = t(celltype_proportion), 
                              gene_expression_mat = as.matrix(spatial_count), 
                              location = as.matrix(location),
                              covariates = NULL,
                              project = "Simulator")

Obj = preprocess_input(Obj, 
                       # Celina object
                       cell_types_to_test = cell_types_to_test, 
                       # a vector of cell types to be used for testing 
                       scRNA_count = as.matrix(scRNA_count_subset), 
                       # a gene x cell expression matrix of reference scRNA-seq data
                       sc_cell_type_labels = sc_meta_in$cellType, 
                       # a vector of cell type labels for each cell in scRNA_count
                       threshold = 0)

Obj = Calculate_Kernel(Obj, approximation = F)
Obj <- tryCatch(
  Testing_interaction_all(Obj, num_cores = 1),
  error = function(e) {
    message("Testing_interaction_all failed: ", e$message)
    NULL
  }
)

# Filter results where CombinedPvals <= 0.05 and create a list of significant genes
significant_genes <- setNames(vector("list", length(cell_types_to_test)), cell_types_to_test)

if (!is.null(Obj) && length(Obj@result) > 0) {
  for (ct in names(Obj@result)) {
    sig_rows <- Obj@result[[ct]][Obj@result[[ct]]$CombinedPvals <= 0.05, ]
    significant_genes[[ct]] <- rownames(sig_rows)
  }
}

# Save the list of significant genes as RDS file
saveRDS(significant_genes, paste0("./Data/", S, "_gene/CELINA", suffix, ".RDS"))

# Print summary of significant genes found
for (ct in names(significant_genes)) {
  cat(ct, ":", length(significant_genes[[ct]]), "significant genes\n")
}

# # Truth
# Obj = Create_Celina_Object(celltype_mat = as.matrix(celltype_proportion), 
#                               gene_expression_mat = as.matrix(spatial_count), 
#                               location = as.matrix(location),
#                               covariates = NULL,
#                               project = "Simulator_truth")

# Obj = preprocess_input(Obj, 
#                        # Celina object
#                        cell_types_to_test = c('ct1', 'ct2','ct3'), 
#                        # a vector of cell types to be used for testing 
#                        scRNA_count = as.matrix(scRNA_count_subset), 
#                        # a gene x cell expression matrix of reference scRNA-seq data
#                        sc_cell_type_labels = sc_meta_in$cellType, 
#                        # a vector of cell type labels for each cell in scRNA_count
#                        threshold = 0)

# Obj = Calculate_Kernel(Obj, approximation = F)
# Obj = Testing_interaction_all(Obj,num_cores=20)

# # Filter results where CombinedPvals <= 0.05 and create a list of significant genes
# significant_genes <- list()

# # Loop through each cell type in the results
# for (ct in names(Obj@result)) {
#   # Filter rows where CombinedPvals <= 0.05
#   sig_rows <- Obj@result[[ct]][Obj@result[[ct]]$CombinedPvals <= 0.05, ]
  
#   # Extract rownames of significant genes
#   significant_genes[[ct]] <- rownames(sig_rows)
# }

# # Save the list of significant genes as RDS file
# saveRDS(significant_genes, paste0("./Data/",S,"_gene/CELINA_trueprop", suffix, ".RDS"))

# # Print summary of significant genes found
# for (ct in names(significant_genes)) {
#   cat(ct, ":", length(significant_genes[[ct]]), "significant genes\n")
# }
