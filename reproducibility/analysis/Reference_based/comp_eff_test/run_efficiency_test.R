# Comprehensive Computational Efficiency Testing Script
# Tests 3 methods: ST-COS, scDesign3, and SRTsim (scCube disabled)

# Libraries
library(data.table)
library(Seurat)
library(microbenchmark)

# Setup
output_dir <- "./Reference_based/comp_eff_test"
gene_counts <- c(50, 100, 500, 1000, 3000)

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

# Function to measure memory usage
measure_memory <- function(expr) {
  gc_before <- gc(reset = TRUE)
  result <- eval(expr)
  gc_after <- gc()
  memory_used <- sum(gc_after[,2] - gc_before[,2]) * 1024 * 1024 / 1024^2  # Convert to MB
  return(list(result = result, memory_mb = memory_used))
}

# Function to run efficiency test for R methods
run_r_method <- function(method_name, gene_count, gene_subset) {
  cat(sprintf("\n=== Testing %s with %d genes ===\n", method_name, gene_count))
  
  tryCatch({
    start_time <- Sys.time()
    
    if (method_name == "ST-COS") {
      source(file.path(output_dir, "method_stcos.R"))
      mem_result <- measure_memory(run_stcos(gene_subset, output_dir, gene_count))
    } else if (method_name == "scDesign3") {
      source(file.path(output_dir, "method_scdesign3.R"))
      mem_result <- measure_memory(run_scdesign3(gene_subset, output_dir, gene_count))
    } else if (method_name == "SRTsim") {
      source(file.path(output_dir, "method_srtsim.R"))
      mem_result <- measure_memory(run_srtsim(gene_subset, output_dir, gene_count))
    }
    
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    new_result <- data.table(
      method = method_name,
      gene_count = gene_count,
      time_seconds = elapsed_time,
      memory_mb = mem_result$memory_mb,
      success = TRUE,
      error_message = ""
    )
    
    cat(sprintf("✓ %s completed in %.2f seconds, using %.2f MB\n", 
                method_name, elapsed_time, mem_result$memory_mb))
    
    return(new_result)
    
  }, error = function(e) {
    cat(sprintf("✗ %s failed: %s\n", method_name, e$message))
    return(data.table(
      method = method_name,
      gene_count = gene_count,
      time_seconds = NA,
      memory_mb = NA,
      success = FALSE,
      error_message = as.character(e$message)
    ))
  })
}

# Main testing loop - 3 METHODS
for (gene_count in gene_counts) {
  cat(sprintf("\n\n>>> TESTING WITH %d GENES <<<\n", gene_count))
  
  # Select gene subset
  gene_subset <- available_genes[1:min(gene_count, length(available_genes))]
  
  # Test ALL R methods - MANDATORY
  r_methods <- c("ST-COS", "scDesign3", "SRTsim")
  
  for (method in r_methods) {
    result <- run_r_method(method, gene_count, gene_subset)
    results <- rbind(results, result)
    
    # Save intermediate results
    fwrite(results, file.path(output_dir, "efficiency_results_intermediate.csv"))
    
    # Clean up memory between tests
    gc()
  }
  
  # Save intermediate results
  fwrite(results, file.path(output_dir, "efficiency_results_intermediate.csv"))
  
  # Clean up
  gc()
}

# Save final results
fwrite(results, file.path(output_dir, "efficiency_results_final.csv"))

# Generate summary report
source(file.path(output_dir, "generate_report.R"))
generate_efficiency_report(results, output_dir)

cat("\n=== Efficiency testing completed ===\n")
cat("Results saved to:", file.path(output_dir, "efficiency_results_final.csv"), "\n")

# Print final summary
successful_methods <- unique(results[success == TRUE]$method)
failed_methods <- unique(results[success == FALSE]$method)

cat("\nFINAL SUMMARY:\n")
cat("Successful methods:", paste(successful_methods, collapse = ", "), "\n")
if (length(failed_methods) > 0) {
  cat("Failed methods:", paste(failed_methods, collapse = ", "), "\n")
}

# Show detailed results by method
summary_by_method <- results[, .(
  tests_run = .N,
  successful = sum(success),
  failed = sum(!success),
  success_rate = round(sum(success)/.N * 100, 1),
  avg_time = round(mean(time_seconds, na.rm = TRUE), 2),
  min_time = round(min(time_seconds, na.rm = TRUE), 2),
  max_time = round(max(time_seconds, na.rm = TRUE), 2)
), by = method]

cat("\nDetailed results by method:\n")
print(summary_by_method)
