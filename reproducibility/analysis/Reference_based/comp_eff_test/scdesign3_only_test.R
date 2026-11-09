# scDesign3 Only Efficiency Test
# Focused testing script for scDesign3 method only
# Based on comprehensive_efficiency_test.R
# WITH peakRAM PACKAGE FOR ACCURATE MEMORY MONITORING

library(data.table)
library(Seurat)
library(scDesign3)
library(SingleCellExperiment)
library(mgcv)
library(peakRAM)  # For accurate memory monitoring

cat("=== scDesign3 ONLY EFFICIENCY TEST WITH peakRAM MONITORING ===\n")
cat("Testing: scDesign3 method exclusively\n")
cat("Using peakRAM package for accurate memory measurement\n")
cat("COMPLETE workflow - no shortcuts\n\n")

# Setup
output_dir <- "./Reference_based/comp_eff_test"
gene_counts <- c(50, 100, 500, 1000, 3000)

# Initialize results tracking with peakRAM monitoring
results <- data.table(
  method = character(),
  gene_count = integer(),
  time_seconds = numeric(),
  peak_ram_mb = numeric(),
  total_ram_mb = numeric(),
  function_ram_mb = numeric(),
  output_size_mb = numeric(),
  success = logical(),
  error_message = character(),
  status = character()
)

# Load real data once (following real_parameter.R setup)
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

# Load coordinates and layers
coord <- as.data.frame(GetTissueCoordinates(your_data, image = "slice1"))
colnames(coord) <- c("row", "col", "cellnames")

layer <- fread("./data_external/DLPFC/layers/151673_layers.txt", header = FALSE)
layer_coord <- merge(coord, layer, by.x = "cellnames", by.y = "V1", all = FALSE)

# Load gene names
gene_names <- fread("./data_external/DLPFC/151673_selected_gene.csv")
available_genes <- gene_names$gene
cat("Data loaded successfully. Available genes:", length(available_genes), "\n")
cat("Spots with layer info:", nrow(layer_coord), "\n\n")

# Helper function to create standardized result data.table
create_result_dt <- function(method, gene_count, time_seconds = NA, peak_ram_mb = NA, 
                           total_ram_mb = NA, function_ram_mb = NA, output_size_mb = NA, 
                           success = FALSE, error_message = "", status = "FAILED") {
  return(data.table(
    method = method,
    gene_count = as.integer(gene_count),
    time_seconds = as.numeric(time_seconds),
    peak_ram_mb = as.numeric(peak_ram_mb),
    total_ram_mb = as.numeric(total_ram_mb),
    function_ram_mb = as.numeric(function_ram_mb),
    output_size_mb = as.numeric(output_size_mb),
    success = as.logical(success),
    error_message = as.character(error_message),
    status = as.character(status)
  ))
}

# scDesign3 Implementation with COMPLETE workflow and peakRAM monitoring
test_scdesign3 <- function(gene_count) {
  cat(sprintf("\n=== scDesign3 - %d genes (COMPLETE WORKFLOW) ===\n", gene_count))
  
  scdesign3_computation <- function() {
    cat(sprintf("  Starting COMPLETE scDesign3 computation for %d genes\n", gene_count))
    
    # Prepare data
    selected_genes <- available_genes[1:gene_count]
    expr_subset <- your_data@assays$Spatial$counts[selected_genes, layer_coord$cellnames]
    rownames(expr_subset) <- selected_genes
    colnames(expr_subset) <- layer_coord$cellnames
    cat(sprintf("  Data extracted: %d genes x %d spots\n", nrow(expr_subset), ncol(expr_subset)))
    
    # Create colData with spatial coordinates
    col_data <- DataFrame(
      spatial1 = layer_coord$row,
      spatial2 = layer_coord$col,
      cell_type = layer_coord$V2
    )
    rownames(col_data) <- layer_coord$cellnames
    
    # Create SCE object
    cat(sprintf("  Step 1: Creating SingleCellExperiment object\n"))
    sce <- SingleCellExperiment(
      assays = list(counts = as.matrix(expr_subset)),
      colData = col_data
    )
    
    # Use appropriate k parameter that scales with gene count for COMPLETE analysis
    k_param <- as.integer(min(400, max(50, gene_count/3)))  # Ensure integer value
    cat(sprintf("  Step 2: Running COMPLETE scDesign3 with k=%d (no shortcuts)\n", k_param))
    
    # COMPLETE scDesign3 simulation - no shortcuts
    example_simu <- scdesign3(
      sce = sce,
      assay_use = "counts",
      celltype = "cell_type",
      pseudotime = NULL,
      spatial = c("spatial1", "spatial2"),
      other_covariates = NULL,
      mu_formula = paste0("s(spatial1, spatial2, bs = 'gp', k = ", k_param, ")"),
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
    
    sim_counts <- example_simu$new_count
    
    cat(sprintf("  COMPLETE scDesign3 computation finished: %d genes, %d spots\n", 
                nrow(sim_counts), ncol(sim_counts)))
    return(sim_counts)
  }
  
  tryCatch({
    # Use peakRAM to monitor memory usage
    start_time <- Sys.time()
    peak_result <- peakRAM(result <- scdesign3_computation())
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Save results
    write.csv(t(result), 
              file.path(output_dir, paste0("scdesign3_only_", gene_count, ".csv")))
    
    cat(sprintf("  SUCCESS - Peak RAM: %.2f MB, Total RAM: %.2f MB, Function RAM: %.2f MB\n",
                peak_result$Peak_RAM_Used_MiB,
                peak_result$Total_RAM_Used_MiB,
                peak_result$Function_RAM_Used_MiB))
    
    return(create_result_dt(
      method = "scDesign3",
      gene_count = gene_count,
      time_seconds = elapsed_time,
      peak_ram_mb = peak_result$Peak_RAM_Used_MiB,
      total_ram_mb = peak_result$Total_RAM_Used_MiB,
      function_ram_mb = peak_result$Function_RAM_Used_MiB,
      output_size_mb = as.numeric(object.size(result)) / 1024^2,
      success = TRUE,
      error_message = "",
      status = "SUCCESS"
    ))
    
  }, error = function(e) {
    cat(sprintf("  FAILED: %s\n", e$message))
    return(create_result_dt(
      method = "scDesign3",
      gene_count = gene_count,
      success = FALSE,
      error_message = as.character(e$message),
      status = "FAILED"
    ))
  })
}

# Run scDesign3 focused test
cat("\n=== Starting scDesign3 ONLY Efficiency Test with peakRAM Monitoring ===\n")
cat("Gene counts:", paste(gene_counts, collapse = ", "), "\n")
cat("Complete workflow - all computational steps included\n\n")

for (gene_count in gene_counts) {
  cat(sprintf("\n\n>>> TESTING scDesign3 WITH %d GENES <<<\n", gene_count))
  
  # Test scDesign3 with COMPLETE workflow
  result <- test_scdesign3(gene_count)
  results <- rbindlist(list(results, result), fill = TRUE)
  
  # Memory cleanup
  gc(verbose = FALSE)
  
  # Save intermediate results
  fwrite(results, file.path(output_dir, "scdesign3_only_results.csv"))
  
  cat(sprintf("Completed scDesign3 testing with %d genes\n", gene_count))
}

# Final analysis and summary
cat("\n=== scDesign3 EFFICIENCY TEST RESULTS ===\n")
print(results)

# Summary statistics with peakRAM results
if (nrow(results[success == TRUE]) > 0) {
  summary_stats <- results[success == TRUE, .(
    tests_completed = .N,
    avg_time = round(mean(time_seconds, na.rm = TRUE), 2),
    min_time = round(min(time_seconds, na.rm = TRUE), 2),
    max_time = round(max(time_seconds, na.rm = TRUE), 2),
    avg_peak_ram = round(mean(peak_ram_mb, na.rm = TRUE), 2),
    max_peak_ram = round(max(peak_ram_mb, na.rm = TRUE), 2),
    avg_total_ram = round(mean(total_ram_mb, na.rm = TRUE), 2),
    avg_function_ram = round(mean(function_ram_mb, na.rm = TRUE), 2),
    max_function_ram = round(max(function_ram_mb, na.rm = TRUE), 2),
    avg_output_size = round(mean(output_size_mb, na.rm = TRUE), 2),
    total_time = round(sum(time_seconds, na.rm = TRUE), 2)
  )]
  
  cat("\nscDesign3 Summary Statistics (with peakRAM monitoring):\n")
  print(summary_stats)
  
  # Performance analysis by gene count
  performance_by_genes <- results[success == TRUE, .(
    time_seconds = round(time_seconds, 2),
    peak_ram_mb = round(peak_ram_mb, 2),
    function_ram_mb = round(function_ram_mb, 2),
    output_size_mb = round(output_size_mb, 2),
    time_per_gene = round(time_seconds / gene_count, 6),
    ram_per_gene = round(function_ram_mb / gene_count, 6)
  ), by = gene_count]
  
  cat("\nscDesign3 Performance by Gene Count:\n")
  print(performance_by_genes)
  
  # Scalability analysis
  if (nrow(results[success == TRUE]) > 1) {
    cat("\n=== scDesign3 SCALABILITY ANALYSIS ===\n")
    
    # Time complexity analysis
    time_model <- lm(time_seconds ~ gene_count, data = results[success == TRUE])
    time_slope <- round(coef(time_model)[2], 6)
    time_r_squared <- round(summary(time_model)$r.squared, 4)
    
    cat(sprintf("Time Scaling: %.6f seconds per additional gene (R² = %.4f)\n", 
                time_slope, time_r_squared))
    
    # RAM scaling analysis
    if (!all(is.na(results[success == TRUE]$function_ram_mb))) {
      ram_model <- lm(function_ram_mb ~ gene_count, data = results[success == TRUE])
      ram_slope <- round(coef(ram_model)[2], 6)
      ram_r_squared <- round(summary(ram_model)$r.squared, 4)
      
      cat(sprintf("RAM Scaling: %.6f MB per additional gene (R² = %.4f)\n", 
                  ram_slope, ram_r_squared))
    }
    
    # Efficiency score calculation
    efficiency_analysis <- results[success == TRUE, .(
      gene_count = gene_count,
      time_per_gene = round(time_seconds / gene_count, 6),
      ram_per_gene = round(function_ram_mb / gene_count, 6),
      efficiency_score = round(gene_count^2 / (time_seconds * function_ram_mb + 1), 4)
    )]
    
    cat("\nscDesign3 Efficiency Analysis:\n")
    print(efficiency_analysis)
    
    # Save scalability results
    fwrite(performance_by_genes, file.path(output_dir, "scdesign3_only_performance.csv"))
    fwrite(efficiency_analysis, file.path(output_dir, "scdesign3_only_efficiency.csv"))
  }
  
  # Save summary statistics
  fwrite(summary_stats, file.path(output_dir, "scdesign3_only_summary.csv"))
}

# Save final comprehensive results
fwrite(results, file.path(output_dir, "scdesign3_only_complete_results.csv"))

cat("\n=== scDesign3 EFFICIENCY TEST COMPLETED ===\n")
cat("Results saved to:", output_dir, "\n")
cat("Files generated:\n")
cat("- scdesign3_only_complete_results.csv (main results)\n")
cat("- scdesign3_only_summary.csv (summary statistics)\n")
if (exists("performance_by_genes")) {
  cat("- scdesign3_only_performance.csv (performance by gene count)\n")
}
if (exists("efficiency_analysis")) {
  cat("- scdesign3_only_efficiency.csv (efficiency analysis)\n")
}
cat("- scdesign3_only_[gene_count].csv (individual simulation results)\n")

# Print final summary
if (nrow(results[success == TRUE]) > 0) {
  cat("\nFinal Summary:\n")
  cat(sprintf("Successful tests: %d/%d\n", 
              nrow(results[success == TRUE]), nrow(results)))
  cat(sprintf("Total execution time: %.2f seconds\n", 
              sum(results[success == TRUE]$time_seconds, na.rm = TRUE)))
  cat(sprintf("Average time per test: %.2f seconds\n", 
              mean(results[success == TRUE]$time_seconds, na.rm = TRUE)))
  cat(sprintf("Peak RAM usage: %.2f MB\n", 
              max(results[success == TRUE]$peak_ram_mb, na.rm = TRUE)))
} else {
  cat("\nNo successful tests completed. Check error messages above.\n")
}

cat("\nscDesign3 focused efficiency testing completed!\n")