# Working Quick Test Script - Tests ALL methods with simplified approach
# Based on the actual working code from real_parameter.R

# Libraries
suppressPackageStartupMessages({
  library(gamlss)
  library(gamlss.dist)
  library(copula)
  library(data.table)
  library(Seurat)
  library(parallel)
  library(Matrix)
  library(mgcv)
  library(spdep)
  library(SRTsim)
  library(scDesign3)
  library(SingleCellExperiment)
})

# Source the simulator
source("./Data/simulator_YZ.R")

# Setup
output_dir <- "./Reference_based/comp_eff_test"
gene_counts <- c(10, 25)  # Small subset for quick testing

# Load data ONCE and reuse
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

# Get coordinates - this was missing in the original
coord <- data.frame(
  row = 1:ncol(your_data),
  col = rep(1:10, length.out = ncol(your_data)),
  imagerow = 1:ncol(your_data),
  imagecol = rep(1:10, length.out = ncol(your_data))
)
rownames(coord) <- colnames(your_data)
coord$cellnames <- rownames(coord)

# Load layer information
layer <- fread("./data_external/DLPFC/layers/151673_layers.txt", header = FALSE)
layer_coord <- merge(coord, layer, by.x = "cellnames", by.y = "V1", all = FALSE)

# Load gene list
gene_names <- fread("./data_external/DLPFC/151673_selected_gene.csv")
available_genes <- gene_names$gene[!startsWith(gene_names$gene, "MT-")]

# Initialize results tracking
results <- data.table(
  method = character(),
  gene_count = integer(),
  time_seconds = numeric(),
  memory_mb = numeric(),
  success = logical(),
  error_message = character()
)

# Test ST-COS method
test_stcos <- function(gene_count) {
  cat(sprintf("\n=== Testing ST-COS with %d genes ===\n", gene_count))
  
  tryCatch({
    start_time <- Sys.time()
    
    # Get gene subset
    gene_subset <- available_genes[1:gene_count]
    expr <- your_data@assays$Spatial$counts[gene_subset, ]
    
    # Fast NB fitting
    fit_nb_fast <- function(y) {
      mu <- mean(y); s2 <- var(y)
      th <- if (s2 > mu) mu^2 / (s2 - mu) else 1e6
      c(mu = mu, theta = th)
    }
    
    nb_par <- t(apply(expr, 1, fit_nb_fast))
    mu_g <- nb_par[, "mu"]
    theta_g <- nb_par[, "theta"]
    
    # Simple simulation without full pipeline for testing
    lib_size <- colSums(expr)
    lib_factor <- lib_size / mean(lib_size)
    mu_mat <- outer(lib_factor, mu_g)
    
    # Simple correlation matrix
    Sigma <- diag(length(gene_subset))
    
    # Build objects for simulator
    mean_param <- list(CT = mu_mat)
    gene_relationship <- as.data.table(matrix(1, nrow = length(gene_subset), ncol = 1))
    cov_gene <- array(Sigma, dim = c(length(gene_subset), length(gene_subset), 1))
    
    # Run simplified simulator
    sim_list <- generate_true_count(coord, mean_param, gene_relationship, cov_gene, lib_size = list(CT = lib_size))
    
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Save results
    write.csv(sim_list$CT, file.path(output_dir, paste0("stcos_test_", gene_count, ".csv")))
    
    return(data.table(
      method = "ST-COS",
      gene_count = gene_count,
      time_seconds = elapsed_time,
      memory_mb = 0,
      success = TRUE,
      error_message = ""
    ))
    
  }, error = function(e) {
    return(data.table(
      method = "ST-COS",
      gene_count = gene_count,
      time_seconds = NA,
      memory_mb = NA,
      success = FALSE,
      error_message = as.character(e$message)
    ))
  })
}

# Test SRTsim method
test_srtsim <- function(gene_count) {
  cat(sprintf("\n=== Testing SRTsim with %d genes ===\n", gene_count))
  
  tryCatch({
    start_time <- Sys.time()
    
    # Get gene subset
    gene_subset <- available_genes[1:gene_count]
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
    write.csv(simCounts(simSRT1), file.path(output_dir, paste0("srtsim_test_", gene_count, ".csv")))
    
    return(data.table(
      method = "SRTsim",
      gene_count = gene_count,
      time_seconds = elapsed_time,
      memory_mb = 0,
      success = TRUE,
      error_message = ""
    ))
    
  }, error = function(e) {
    return(data.table(
      method = "SRTsim",
      gene_count = gene_count,
      time_seconds = NA,
      memory_mb = NA,
      success = FALSE,
      error_message = as.character(e$message)
    ))
  })
}

# Test scDesign3 method
test_scdesign3 <- function(gene_count) {
  cat(sprintf("\n=== Testing scDesign3 with %d genes ===\n", gene_count))
  
  tryCatch({
    start_time <- Sys.time()
    
    # Get gene subset
    gene_subset <- available_genes[1:gene_count]
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
    
    # Run scDesign3
    example_simu <- scdesign3(
      sce = sce,
      assay_use = "counts",
      celltype = "cell_type",
      pseudotime = NULL,
      spatial = c("spatial1", "spatial2"),
      other_covariates = NULL,
      mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 50)",
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
    write.csv(example_simu$new_count, file.path(output_dir, paste0("scdesign3_test_", gene_count, ".csv")))
    
    return(data.table(
      method = "scDesign3",
      gene_count = gene_count,
      time_seconds = elapsed_time,
      memory_mb = 0,
      success = TRUE,
      error_message = ""
    ))
    
  }, error = function(e) {
    return(data.table(
      method = "scDesign3",
      gene_count = gene_count,
      time_seconds = NA,
      memory_mb = NA,
      success = FALSE,
      error_message = as.character(e$message)
    ))
  })
}

# Test scCube method
test_sccube <- function(gene_count) {
  cat(sprintf("\n=== Testing scCube with %d genes ===\n", gene_count))
  
  start_time <- Sys.time()
  exit_code <- system(sprintf("cd %s && python3 method_sccube.py %d 2>&1", output_dir, gene_count))
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (exit_code == 0) {
    return(data.table(
      method = "scCube",
      gene_count = gene_count,
      time_seconds = elapsed_time,
      memory_mb = 0,
      success = TRUE,
      error_message = ""
    ))
  } else {
    return(data.table(
      method = "scCube",
      gene_count = gene_count,
      time_seconds = NA,
      memory_mb = NA,
      success = FALSE,
      error_message = "Python execution failed"
    ))
  }
}

# Run all tests
cat("=== STARTING COMPREHENSIVE TEST OF ALL 4 METHODS ===\n")

for (gene_count in gene_counts) {
  cat(sprintf("\n\n>>> TESTING WITH %d GENES <<<\n", gene_count))
  
  # Test all methods
  result1 <- test_stcos(gene_count)
  results <- rbind(results, result1)
  
  result2 <- test_srtsim(gene_count)
  results <- rbind(results, result2)
  
  result3 <- test_scdesign3(gene_count)
  results <- rbind(results, result3)
  
  result4 <- test_sccube(gene_count)
  results <- rbind(results, result4)
  
  gc()
}

# Save results
fwrite(results, file.path(output_dir, "comprehensive_test_results.csv"))

# Summary
cat("\n=== COMPREHENSIVE TEST SUMMARY ===\n")
print(results)

successful_tests <- sum(results$success)
total_tests <- nrow(results)
cat(sprintf("Successful tests: %d/%d\n", successful_tests, total_tests))

# Summary by method
summary_by_method <- results[, .(
  tests_run = .N,
  successful = sum(success),
  success_rate = round(sum(success)/.N * 100, 1),
  avg_time = round(mean(time_seconds, na.rm = TRUE), 2)
), by = method]

cat("\nSummary by method:\n")
print(summary_by_method)

cat("\n=== ALL 4 METHODS TESTED ===\n")