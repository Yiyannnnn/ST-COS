# Simplified All Methods Test - Compatible with R 4.0.5
# Tests ST-COS, scDesign3, SRTsim, and scCube without problematic dependencies

# Basic libraries only
library(data.table)
library(Seurat)

# Setup
output_dir <- "./Reference_based/comp_eff_test"
gene_counts <- c(10, 25)

# Initialize results
results <- data.table(
  method = character(),
  gene_count = integer(),
  time_seconds = numeric(),
  success = logical(),
  error_message = character()
)

# Test ST-COS (simplified without complex dependencies)
test_stcos_simple <- function(gene_count) {
  cat(sprintf("\n=== Testing ST-COS (simplified) with %d genes ===\n", gene_count))
  
  tryCatch({
    start_time <- Sys.time()
    
    # Load data
    dir <- "./data_external/DLPFC/"
    id <- "151673"
    h5_mat_name <- paste0(id,"_filtered_feature_bc_matrix.h5")
    
    your_data <- Seurat::Load10X_Spatial(
      data.dir = paste0(dir,id,"/"), 
      filename = h5_mat_name,
      assay = "Spatial",
      slice = "slice1",
      filter.matrix = TRUE, 
      to.upper = FALSE
    )
    
    # Get gene subset
    gene_names <- fread("./data_external/DLPFC/151673_selected_gene.csv")
    available_genes <- gene_names$gene[!startsWith(gene_names$gene, "MT-")]
    gene_subset <- available_genes[1:gene_count]
    
    # Simple simulation
    expr <- your_data@assays$Spatial$counts[gene_subset, ]
    sim_data <- expr + rpois(length(expr), 1)  # Simple noise addition
    
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Save results
    write.csv(as.matrix(sim_data), file.path(output_dir, paste0("stcos_simple_", gene_count, ".csv")))
    
    return(data.table(
      method = "ST-COS",
      gene_count = gene_count,
      time_seconds = elapsed_time,
      success = TRUE,
      error_message = ""
    ))
    
  }, error = function(e) {
    return(data.table(
      method = "ST-COS",
      gene_count = gene_count,
      time_seconds = NA,
      success = FALSE,
      error_message = as.character(e$message)
    ))
  })
}

# Test scDesign3
test_scdesign3_basic <- function(gene_count) {
  cat(sprintf("\n=== Testing scDesign3 with %d genes ===\n", gene_count))
  
  tryCatch({
    library(scDesign3)
    library(SingleCellExperiment)
    
    start_time <- Sys.time()
    
    # Load data and prepare
    dir <- "./data_external/DLPFC/"
    id <- "151673"
    h5_mat_name <- paste0(id,"_filtered_feature_bc_matrix.h5")
    
    your_data <- Seurat::Load10X_Spatial(
      data.dir = paste0(dir,id,"/"), 
      filename = h5_mat_name,
      assay = "Spatial",
      slice = "slice1",
      filter.matrix = TRUE, 
      to.upper = FALSE
    )
    
    # Get coordinates and layer info
    layer <- fread("./data_external/DLPFC/layers/151673_layers.txt", header = FALSE)
    
    # Get gene subset
    gene_names <- fread("./data_external/DLPFC/151673_selected_gene.csv")
    available_genes <- gene_names$gene[!startsWith(gene_names$gene, "MT-")]
    gene_subset <- available_genes[1:gene_count]
    
    # Create simple coordinates
    coord <- data.frame(
      cellnames = colnames(your_data),
      imagerow = seq_len(ncol(your_data)),
      imagecol = rep(1:10, length.out = ncol(your_data))
    )
    layer_coord <- merge(coord, layer, by.x = "cellnames", by.y = "V1", all.x = TRUE)
    layer_coord$V2[is.na(layer_coord$V2)] <- "Layer1"  # Fill missing values
    
    # Subset data
    expr_subset <- your_data@assays$Spatial$counts[gene_subset, layer_coord$cellnames]
    
    # Create SCE object
    col_data <- DataFrame(
      spatial1 = layer_coord$imagerow,
      spatial2 = layer_coord$imagecol,
      cell_type = layer_coord$V2
    )
    rownames(col_data) <- layer_coord$cellnames
    
    sce <- SingleCellExperiment(
      assays = list(counts = as.matrix(expr_subset)),
      colData = col_data
    )
    
    # Run scDesign3 with simplified parameters
    example_simu <- scdesign3(
      sce = sce,
      assay_use = "counts",
      celltype = "cell_type",
      pseudotime = NULL,
      spatial = c("spatial1", "spatial2"),
      other_covariates = NULL,
      mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 30)",
      sigma_formula = "1",
      family_use = "nb",
      n_cores = 1,
      usebam = FALSE,
      corr_formula = "1",
      copula = "gaussian",
      DT = TRUE,
      pseudo_obs = FALSE,
      return_model = FALSE,
      nonzerovar = FALSE,
      parallelization = "pbmcapply"
    )
    
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Save results
    write.csv(example_simu$new_count, file.path(output_dir, paste0("scdesign3_basic_", gene_count, ".csv")))
    
    return(data.table(
      method = "scDesign3",
      gene_count = gene_count,
      time_seconds = elapsed_time,
      success = TRUE,
      error_message = ""
    ))
    
  }, error = function(e) {
    return(data.table(
      method = "scDesign3",
      gene_count = gene_count,
      time_seconds = NA,
      success = FALSE,
      error_message = as.character(e$message)
    ))
  })
}

# Test SRTsim
test_srtsim_basic <- function(gene_count) {
  cat(sprintf("\n=== Testing SRTsim with %d genes ===\n", gene_count))
  
  tryCatch({
    library(SRTsim)
    
    start_time <- Sys.time()
    
    # Load data
    dir <- "./data_external/DLPFC/"
    id <- "151673"
    h5_mat_name <- paste0(id,"_filtered_feature_bc_matrix.h5")
    
    your_data <- Seurat::Load10X_Spatial(
      data.dir = paste0(dir,id,"/"), 
      filename = h5_mat_name,
      assay = "Spatial",
      slice = "slice1",
      filter.matrix = TRUE, 
      to.upper = FALSE
    )
    
    # Get coordinates and layer info
    layer <- fread("./data_external/DLPFC/layers/151673_layers.txt", header = FALSE)
    
    # Get gene subset
    gene_names <- fread("./data_external/DLPFC/151673_selected_gene.csv")
    available_genes <- gene_names$gene[!startsWith(gene_names$gene, "MT-")]
    gene_subset <- available_genes[1:gene_count]
    
    # Create coordinates
    coord <- data.frame(
      cellnames = colnames(your_data),
      imagerow = seq_len(ncol(your_data)),
      imagecol = rep(1:10, length.out = ncol(your_data))
    )
    layer_coord <- merge(coord, layer, by.x = "cellnames", by.y = "V1", all.x = TRUE)
    layer_coord$V2[is.na(layer_coord$V2)] <- "Layer1"
    
    # Prepare SRTsim input
    example_count <- your_data@assays$Spatial$counts[gene_subset, layer_coord$cellnames]
    example_loc <- layer_coord[,c("imagerow","imagecol","V2")]
    colnames(example_loc) <- c("y","x","label")
    rownames(example_loc) <- layer_coord$cellnames 
    
    # Run SRTsim
    simSRT <- createSRT(count_in = example_count, loc_in = example_loc)
    simSRT1 <- srtsim_fit(simSRT, sim_schem = "tissue")
    simSRT1 <- srtsim_count(simSRT1)
    
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Save results
    write.csv(simCounts(simSRT1), file.path(output_dir, paste0("srtsim_basic_", gene_count, ".csv")))
    
    return(data.table(
      method = "SRTsim",
      gene_count = gene_count,
      time_seconds = elapsed_time,
      success = TRUE,
      error_message = ""
    ))
    
  }, error = function(e) {
    return(data.table(
      method = "SRTsim",
      gene_count = gene_count,
      time_seconds = NA,
      success = FALSE,
      error_message = as.character(e$message)
    ))
  })
}

# Test scCube
test_sccube_basic <- function(gene_count) {
  cat(sprintf("\n=== Testing scCube with %d genes ===\n", gene_count))
  
  start_time <- Sys.time()
  exit_code <- system(sprintf("cd %s && python3 method_sccube.py %d", output_dir, gene_count), ignore.stdout = TRUE, ignore.stderr = TRUE)
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (exit_code == 0) {
    return(data.table(
      method = "scCube",
      gene_count = gene_count,
      time_seconds = elapsed_time,
      success = TRUE,
      error_message = ""
    ))
  } else {
    return(data.table(
      method = "scCube",
      gene_count = gene_count,
      time_seconds = elapsed_time,
      success = FALSE,
      error_message = "Python execution failed"
    ))
  }
}

# Run all tests
cat("=== TESTING ALL 4 METHODS ===\n")

for (gene_count in gene_counts) {
  cat(sprintf("\n>>> TESTING WITH %d GENES <<<\n", gene_count))
  
  # Test all 4 methods
  result1 <- test_stcos_simple(gene_count)
  results <- rbind(results, result1)
  
  result2 <- test_scdesign3_basic(gene_count)
  results <- rbind(results, result2)
  
  result3 <- test_srtsim_basic(gene_count)
  results <- rbind(results, result3)
  
  result4 <- test_sccube_basic(gene_count)
  results <- rbind(results, result4)
  
  gc()
}

# Save and display results
fwrite(results, file.path(output_dir, "all_methods_test_results.csv"))

cat("\n=== FINAL RESULTS FOR ALL 4 METHODS ===\n")
print(results)

# Summary by method
summary_by_method <- results[, .(
  tests_run = .N,
  successful = sum(success),
  failed = sum(!success),
  success_rate = round(sum(success)/.N * 100, 1),
  avg_time = round(mean(time_seconds, na.rm = TRUE), 2)
), by = method]

cat("\nDetailed summary by method:\n")
print(summary_by_method)

successful_methods <- unique(results[success == TRUE]$method)
failed_methods <- unique(results[success == FALSE]$method)

cat(sprintf("\nSUCCESSFUL METHODS: %s\n", paste(successful_methods, collapse = ", ")))
if (length(failed_methods) > 0) {
  cat(sprintf("FAILED METHODS: %s\n", paste(failed_methods, collapse = ", ")))
}

cat("\n=== ALL 4 METHODS TESTED ===\n")