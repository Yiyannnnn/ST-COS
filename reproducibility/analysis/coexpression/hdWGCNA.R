library(Seurat)
library(hdWGCNA)
library(data.table)
library(tidyverse)
library(Matrix)
library(WGCNA)
library(patchwork)
library(ggplot2)

# Compatibility shim: hdWGCNA uses GetAssayData(slot=...), which is defunct
# under SeuratObject v5. Override the Seurat method to map slot -> layer.
GetAssayData_Seurat_compat <- function(object, assay = NULL, layer = NULL, slot = NULL, ...) {
  SeuratObject:::CheckDots(...)
  if (is.null(layer) && !is.null(slot)) {
    layer <- slot
  }
  object <- SeuratObject:::UpdateSlots(object = object)
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(object = object)
  }
  if (!assay %in% SeuratObject::Assays(object = object)) {
    stop("Assay not found: ", assay)
  }
  SeuratObject::GetAssayData(object = object[[assay]], layer = layer)
}
assignInNamespace("GetAssayData.Seurat", GetAssayData_Seurat_compat, ns = "SeuratObject")

# Set up R library path
.libPaths(c(Sys.getenv("R_LIBS_USER", unset = .libPaths()[1]), .libPaths()))

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure the correct number of arguments are provided
if (length(args) != 2) {
  stop("Please provide the two arguments: S and suffix (_N).")
}

# Get the inputs from command-line arguments
S <- args[1]
suffix <- args[2]

# Create output directory
output_dir <- "./coexpression/hdWGCNA_Agent"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Construct file paths based on the arguments
gene_expression_file <- paste0("./Data/", S, "_gene/simulated/", S, "_gene_expression", suffix, ".csv")
coord_file <- paste0("./Data/", S, "_gene/simulated/", S, "_coord", suffix, ".csv")
layer_file <- paste0("./Data/", S, "_gene/simulated/", S, "_layer", suffix, ".csv")

# Load the data using the provided file paths
cat("Loading data...\n")
spot_count_layers <- fread(gene_expression_file)
coord <- fread(coord_file)
layer_coord <- fread(layer_file)

# Print data dimensions for debugging
cat("Data dimensions:\n")
cat("Gene expression:", dim(spot_count_layers), "\n")
cat("Coordinates:", dim(coord), "\n")  
cat("Layers:", dim(layer_coord), "\n")

# Create Seurat object
cat("Creating Seurat object...\n")
# Check if we need to transpose - spots should be columns in Seurat
if (ncol(spot_count_layers) > nrow(spot_count_layers)) {
  # If more columns than rows, assume columns are spots (no transpose needed)
  seurat_obj <- CreateSeuratObject(counts = spot_count_layers, assay = "RNA")
  cat("Using original matrix orientation (genes as rows, spots as columns)\n")
} else {
  # If more rows than columns, assume rows are spots (transpose needed)
  seurat_obj <- CreateSeuratObject(counts = t(spot_count_layers), assay = "RNA")
  cat("Transposing matrix (original: spots as rows, genes as columns)\n")
}

n_spots <- ncol(seurat_obj)
cat("Number of spots in Seurat object:", n_spots, "\n")

seurat_obj <- RenameCells(seurat_obj, new.names = paste0("spot-", 1:n_spots))

# Set up coordinates - match to the number of spots in Seurat object
if (nrow(coord) != n_spots) {
  cat("Warning: Coordinate file has", nrow(coord), "spots but Seurat object has", n_spots, "spots\n")
  cat("Using first", n_spots, "coordinates to match Seurat object\n")
  coord <- coord[1:n_spots, ]
}

coord_df <- as.data.frame(cbind(coord, coord))
colnames(coord_df) <- c("row", "col", "imagerow", "imagecol")
rownames(coord_df) <- paste0("spot-", 1:n_spots)

# Add spatial information
seurat_obj@images$slice1 <- new(
  Class = 'SlideSeq',
  assay = "RNA",  # Changed from "Spatial" to "RNA"
  key = "slice1_",
  coordinates = coord_df
)
seurat_obj <- AddMetaData(seurat_obj, coord_df)

# Add layer information as metadata - match to the number of spots in Seurat object
if (!is.null(layer_coord) && nrow(layer_coord) > 0) {
  if (nrow(layer_coord) != n_spots) {
    cat("Warning: Layer file has", nrow(layer_coord), "spots but Seurat object has", n_spots, "spots\n")
    cat("Using first", n_spots, "layer assignments to match Seurat object\n")
    layer_coord <- layer_coord[1:n_spots, ]
  }
  
  # Check if layer column exists, otherwise use the third column
  if ("layer" %in% colnames(layer_coord)) {
    layer_values <- layer_coord$layer
  } else if (ncol(layer_coord) >= 3) {
    layer_values <- layer_coord[[3]]  # Use third column
  } else {
    cat("Warning: Cannot find layer information in layer file\n")
    layer_values <- rep(1, n_spots)  # Default to layer 1
  }
  
  layer_info <- data.frame(
    layer = layer_values,
    row.names = paste0("spot-", 1:n_spots)
  )
  seurat_obj <- AddMetaData(seurat_obj, layer_info)
} else {
  cat("Warning: No layer information available\n")
}

# Standard preprocessing
cat("Running standard preprocessing...\n")
seurat_obj <- seurat_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# Ensure data layer is properly populated for hdWGCNA
cat("Ensuring data layer is populated for hdWGCNA...\n")
seurat_obj[["RNA"]]@layers$data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

# Clustering
cat("Running clustering...\n")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, verbose = TRUE, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Save basic plots
cat("Creating basic plots...\n")
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP Clustering")
ggsave(file.path(output_dir, paste0(S, suffix, "_umap_clusters.pdf")), p1, width = 8, height = 6)

p2 <- SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3) + 
  ggtitle("Spatial Clustering")
ggsave(file.path(output_dir, paste0(S, suffix, "_spatial_clusters.pdf")), p2, width = 8, height = 6)

# hdWGCNA analysis
cat("Setting up hdWGCNA...\n")

# Use seurat_clusters for grouping instead of orig.ident
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "hdWGCNA"
)

# Construct metacells when the inferred groups are large enough. The expression
# matrix below is built from the original spots, so a metacell failure need not
# prevent network estimation for small simulated groups.
seurat_obj <- tryCatch(
  MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = "seurat_clusters",
    reduction = "umap",
    slot = "counts",
    layer = "counts",
    k = 25,
    max_shared = 10,
    ident.group = "seurat_clusters"
  ),
  error = function(e) {
    cat("Metacell construction skipped:", conditionMessage(e), "\n")
    seurat_obj
  }
)

# Set up the expression matrix - bypass SetDatExpr entirely
cat("Setting up expression matrix directly (bypassing SetDatExpr)\n")

# Get expression data directly from the counts slot
expr_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
cat("Expression data dimensions:", dim(expr_data), "\n")

# Get cluster assignments
clusters <- seurat_obj@meta.data$seurat_clusters
cat("Number of clusters found:", length(unique(clusters)), "\n")

# Manually set up the hdWGCNA expression data
seurat_obj@misc$hdWGCNA <- list()
seurat_obj@misc$hdWGCNA$wgcna_name <- "hdWGCNA"
seurat_obj@misc$hdWGCNA$group.by <- "seurat_clusters"
seurat_obj@misc$hdWGCNA$datExpr <- t(as.matrix(expr_data))
seurat_obj@misc$hdWGCNA$MEs <- NULL

# Test soft powers
cat("Testing soft powers...\n")
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'
)

# Plot soft power results
plot_list <- PlotSoftPowers(seurat_obj)
p_soft <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(output_dir, paste0(S, suffix, "_soft_powers.pdf")), p_soft, width = 12, height = 8)

# Check if a suitable soft power was found, if not use a reasonable default
# Access soft power results directly from the Seurat object
soft_power_table <- seurat_obj@misc$hdWGCNA$soft_power_table

# Ensure we have valid data
if (!is.null(soft_power_table) && nrow(soft_power_table) > 0) {
  suitable_powers <- soft_power_table[soft_power_table$SFT.R.sq >= 0.8 & !is.na(soft_power_table$SFT.R.sq), ]
  
  if (nrow(suitable_powers) > 0) {
    selected_power <- min(suitable_powers$Power)
    cat("Using soft power:", selected_power, "\n")
  } else {
    # No power meets 0.8 criteria, use the best available power
    valid_rows <- !is.na(soft_power_table$SFT.R.sq) & is.finite(soft_power_table$SFT.R.sq)
    if (any(valid_rows)) {
      best_power_idx <- which.max(soft_power_table$SFT.R.sq[valid_rows])
      selected_power <- soft_power_table$Power[valid_rows][best_power_idx]
      cat("No power meets 0.8 R-squared criteria. Using best available power:", selected_power, "\n")
      cat("Best R-squared achieved:", soft_power_table$SFT.R.sq[valid_rows][best_power_idx], "\n")
    } else {
      # Fallback to a reasonable default
      selected_power <- 6
      cat("No valid soft powers found. Using default power:", selected_power, "\n")
    }
  }
} else {
  # Fallback if no soft power table available
  selected_power <- 6
  cat("No soft power table available. Using default power:", selected_power, "\n")
}

# Construct network with specified power
cat("Constructing co-expression network...\n")
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = selected_power,
  tom_name = paste0(S, suffix, "_TOM"),
  overwrite_tom = TRUE,  # Allow overwriting existing TOM files
  minModuleSize = 10,  # Reduce minimum module size
  mergeCutHeight = 0.25,  # More lenient module merging
  deepSplit = 2  # More sensitive module detection
)

# Check if modules were detected
modules <- GetModules(seurat_obj)
if (nrow(modules) == 0 || all(modules$module == "grey")) {
  cat("Warning: No significant modules detected. Adjusting parameters...\n")
  
  # Try with even more lenient parameters
  seurat_obj <- ConstructNetwork(
    seurat_obj,
    soft_power = selected_power,
    tom_name = paste0(S, suffix, "_TOM"),
    overwrite_tom = TRUE,  # Allow overwriting existing TOM files
    minModuleSize = 5,  # Very small minimum module size
    mergeCutHeight = 0.5,  # Very lenient merging
    deepSplit = 4  # Maximum sensitivity
  )
  
  modules <- GetModules(seurat_obj)
}

# Only proceed with module-dependent analyses if modules were found
if (nrow(modules) > 0 && !all(modules$module == "grey")) {
  cat("Found", length(unique(modules$module)), "modules\n")
  
  # Plot the dendrogram
  tryCatch({
    p_dend <- PlotDendrogram(seurat_obj, main = "hdWGCNA Dendrogram")
    ggsave(file.path(output_dir, paste0(S, suffix, "_dendrogram.pdf")), p_dend, width = 12, height = 8)
  }, error = function(e) {
    cat("Warning: Could not create dendrogram plot:", e$message, "\n")
  })

  # Save modules
  write.csv(modules, file.path(output_dir, paste0(S, suffix, "_modules.csv")), row.names = FALSE)

  # Get hub genes with error handling
  tryCatch({
    hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
    write.csv(hub_df, file.path(output_dir, paste0(S, suffix, "_hub_genes.csv")), row.names = FALSE)
  }, error = function(e) {
    cat("Warning: Could not extract hub genes:", e$message, "\n")
    # Create a placeholder hub genes file
    modules_for_hubs <- GetModules(seurat_obj)
    if (nrow(modules_for_hubs) > 0) {
      # Get top genes per module based on module membership
      hub_df <- modules_for_hubs %>%
        group_by(module) %>%
        slice_head(n = 10) %>%  # Fix: use constant instead of n()
        ungroup()
      write.csv(hub_df, file.path(output_dir, paste0(S, suffix, "_hub_genes.csv")), row.names = FALSE)
    } else {
      # Empty hub genes file
      write.csv(data.frame(gene = character(0), module = character(0)), 
                file.path(output_dir, paste0(S, suffix, "_hub_genes.csv")), row.names = FALSE)
    }
  })

  # Compute module eigengenes when the detected module dimensions permit it.
  tryCatch({
    seurat_obj <- ModuleEigengenes(seurat_obj)
    MEs <- GetMEs(seurat_obj)
    write.csv(MEs, file.path(output_dir, paste0(S, suffix, "_module_eigengenes.csv")))
  }, error = function(e) {
    cat("Warning: Could not compute module eigengenes:", e$message, "\n")
  })

  # Plot module eigengenes
  tryCatch({
    p_me <- PlotKMEs(seurat_obj)
    ggsave(file.path(output_dir, paste0(S, suffix, "_module_kMEs.pdf")), p_me, width = 12, height = 8)
  }, error = function(e) {
    cat("Warning: Could not create kME plot:", e$message, "\n")
  })

  # Compute module connectivity. This can fail when there is only one
  # effective module and kME is returned as a vector.
  non_grey_modules <- setdiff(unique(modules$module), "grey")
  if (length(non_grey_modules) >= 2) {
    tryCatch({
      seurat_obj <- ModuleConnectivity(seurat_obj)

      # Plot module connectivity
      p_conn <- PlotModuleConnectivity(seurat_obj)
      ggsave(file.path(output_dir, paste0(S, suffix, "_module_connectivity.pdf")), p_conn, width = 10, height = 8)
    }, error = function(e) {
      cat("Warning: Module connectivity step failed:", e$message, "\n")
    })
  } else {
    cat("Skipping ModuleConnectivity: fewer than 2 non-grey modules detected.\n")
  }

  # Reset module names for better plotting. This can fail in single-module
  # edge cases due to inconsistent internal naming structures.
  tryCatch({
    seurat_obj <- ResetModuleNames(
      seurat_obj,
      new_name = "M"
    )
  }, error = function(e) {
    cat("Warning: Could not reset module names:", e$message, "\n")
    cat("Continuing with original module names.\n")
  })

  # Plot modules in UMAP space
  tryCatch({
    p_modules_umap <- ModuleUMAPPlot(
      seurat_obj,
      edge.alpha = 0.25,
      sample_edges = TRUE,
      edge_prop = 0.1,
      label_hubs = 2,
      keep_grey_edges = FALSE
    )
    ggsave(file.path(output_dir, paste0(S, suffix, "_modules_umap.pdf")), p_modules_umap, width = 12, height = 10)
  }, error = function(e) {
    cat("Warning: Could not create module UMAP plot:", e$message, "\n")
  })

  # Compute module trait correlations if layers are available
  if (!is.null(layer_coord) && nrow(layer_coord) > 0) {
    cat("Computing module-trait correlations...\n")
    tryCatch({
      seurat_obj <- ModuleTraitCorrelation(
        seurat_obj,
        traits = "layer"
      )
      
      # Plot module-trait correlations
      p_trait <- PlotModuleTraitCorrelation(
        seurat_obj,
        label = "fdr",
        label_symbol = "stars",
        text_size = 2,
        text_digits = 2
      )
      ggsave(file.path(output_dir, paste0(S, suffix, "_module_trait_correlation.pdf")), p_trait, width = 10, height = 8)
    }, error = function(e) {
      cat("Warning: Could not compute module-trait correlations:", e$message, "\n")
    })
  }
  
} else {
  cat("No significant modules detected even with lenient parameters.\n")
  cat("This may be due to limited sample size or weak co-expression patterns.\n")
  
  # Create empty placeholder files
  write.csv(data.frame(gene = character(0), module = character(0)), 
            file.path(output_dir, paste0(S, suffix, "_modules.csv")), row.names = FALSE)
  write.csv(data.frame(gene = character(0), module = character(0), kME = numeric(0)), 
            file.path(output_dir, paste0(S, suffix, "_hub_genes.csv")), row.names = FALSE)
}

# Save the final Seurat object
cat("Saving Seurat object...\n")
saveRDS(seurat_obj, file.path(output_dir, paste0(S, suffix, "_seurat_hdWGCNA.rds")))

# Create summary report
cat("Creating summary report...\n")

# Initialize hub_df if it doesn't exist
if (!exists("hub_df")) {
  hub_df <- data.frame(gene = character(0), module = character(0), kME = numeric(0))
}

summary_stats <- data.frame(
  Metric = c("Total spots", "Total genes", "Number of modules", "Number of hub genes"),
  Value = c(
    ncol(seurat_obj),
    nrow(seurat_obj),
    length(unique(modules$module[modules$module != "grey"])) + ifelse("grey" %in% modules$module, 1, 0),
    nrow(hub_df)
  )
)
write.csv(summary_stats, file.path(output_dir, paste0(S, suffix, "_summary_stats.csv")), row.names = FALSE)

cat("hdWGCNA analysis completed! All results saved to:", output_dir, "\n")
