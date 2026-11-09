#r/4.3.1
#remotes::install_github('RfastOfficial/Rfast')
#remotes::install_github("dmcable/RCTD")
library(spacexr)
library(Matrix)
library(tidyverse)
library(data.table)

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
coord_file <- paste0("./Data/", S, "_gene/simulated/", S, "_coord", suffix, ".csv")
output_file <- paste0("./Data/", S, "_gene/RCTD", suffix, ".csv")

##### Read in single cell and ST data ####
ref <- fread(sc_ref_file)
cluster <- cbind(cell_type = ref[, 1], cell = paste0("cell_", 1:nrow(ref)))
ref <- ref[, -1]
ref <- as.data.frame(round(t(ref)), check.names = FALSE)
colnames(ref) <- paste0("cell_", 1:ncol(ref))

spots <- fread(gene_expression_file)
spots <- as.data.frame(abs(round(t(spots))), check.names = FALSE)
colnames(spots) <- paste0("spot_",1:ncol(spots))
# Keep only overlapping genes with a minimum count > 0
spots <- spots[rowSums(spots) > 0, ]
overlap_gene <- intersect(rownames(spots), rownames(ref))
ref_overlap <- ref[overlap_gene, ]
spots <- spots[overlap_gene, ]

# Clean the cluster labels
cell_types <- cluster$cell_type
names(cell_types) <- cluster$cell
cell_types <- as.factor(cell_types)

nUMI <- colSums(ref_overlap)
names(nUMI) <- colnames(ref_overlap)

## Create reference data object
reference <- Reference(ref_overlap, cell_types, nUMI)

#### Spatial RNA files ####
coords <- tibble(barcodes = colnames(spots))

#### Read and process ST data coordinates
coord <- fread(coord_file)
coords <- coords %>%
  separate(barcodes, into = c("xcoord", "ycoord"), sep = "x", remove = FALSE)

coords <- as.data.frame(coords)
coords$barcodes <- NULL
coords$xcoord <- as.numeric(coord$x_coord)
coords$ycoord <- as.numeric(coord$y_coord)
rownames(coords) <- paste0("spot_",1:nrow(coords))
nUMI_spot <- colSums(spots)

### Create SpatialRNA object
puck <- SpatialRNA(coords, spots, nUMI_spot)

## Create RCTD object (relaxed DE-gene thresholds for simulated data)
myRCTD <- create.RCTD(
  puck, reference,
  max_cores = 1,
  CELL_MIN_INSTANCE = 0,
  gene_cutoff = 0,
  fc_cutoff = 0,
  gene_cutoff_reg = 0,
  fc_cutoff_reg = 0,
  counts_MIN = 0,
  UMI_min = 0
)

#### Run RCTD ####
myRCTD_full <- run.RCTD(myRCTD, doublet_mode = 'full')

#### Summarize results ####
RCTD_results_full <- myRCTD_full@results
RCTD_norm_weights <- as.matrix(sweep(RCTD_results_full$weights, 1, rowSums(RCTD_results_full$weights), '/'))
RCTD_norm_weights_tibb <- as_tibble(RCTD_norm_weights, rownames = "barcodes")
RCTD_norm_weights_tibb$barcodes <- paste0("spot_", gsub("[^0-9]", "", RCTD_norm_weights_tibb$barcodes))

# Save the results
fwrite(RCTD_norm_weights_tibb, output_file)


# # Load necessary C-SIDE packages and functions (if applicable)
# library(Rfast)
# library(spacexr)
# layer_file <- paste0("./Data/", S, "_gene/simulated/", S, "_layer", suffix, ".csv")
# region <- fread(layer_file)
# region_list <- list(layer_1 = paste0("spot_",which(region$layer == 1)),
# layer_2 = paste0("spot_",which(region$layer == 2)),
# layer_3 = paste0("spot_",which(region$layer == 3)))
# myRCTD_SVG <- run.CSIDE.regions(myRCTD_full, region_list,cell_types = cell_types, gene_threshold = .01, cell_type_threshold = 10, fdr = 0.25, doublet_mode = F, weight_threshold = 0.1) 

# saveRDS(myRCTD_SVG,paste0("./Data/", S, "_gene/CSIDE", suffix, ".RDS"))
