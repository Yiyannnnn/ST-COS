# Quick Test Script - Minimal gene counts
# Tests ST-COS, scDesign3, and SRTsim (no scCube)

# Libraries
library(data.table)
library(Seurat)

# Setup
output_dir <- "./Reference_based/comp_eff_test"
gene_counts <- c(10, 25)  # Small subset for quick testing

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

# Quick test function
quick_test_method <- function(method_name, gene_count, gene_subset) {
  cat(sprintf("\n=== Quick Testing %s with %d genes ===\n", method_name, gene_count))
  
  tryCatch({
    start_time <- Sys.time()
    
    if (method_name == "ST-COS") {
      source(file.path(output_dir, "method_stcos.R"))
      result <- run_stcos(gene_subset, output_dir, gene_count)
    } else if (method_name == "scDesign3") {
      source(file.path(output_dir, "method_scdesign3.R"))
      result <- run_scdesign3(gene_subset, output_dir, gene_count)
    } else if (method_name == "SRTsim") {
      source(file.path(output_dir, "method_srtsim.R"))
      result <- run_srtsim(gene_subset, output_dir, gene_count)
    }
    
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    new_result <- data.table(
      method = method_name,
      gene_count = gene_count,
      time_seconds = elapsed_time,
      memory_mb = NA,  # Simplified for quick test
      success = TRUE,
      error_message = ""
    )
    
    cat(sprintf("✓ %s completed in %.2f seconds\n", method_name, elapsed_time))
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

# Quick test loop - TEST R METHODS ONLY
for (gene_count in gene_counts) {
  cat(sprintf("\n\n>>> QUICK TESTING WITH %d GENES <<<\n", gene_count))
  
  # Select gene subset
  gene_subset <- available_genes[1:min(gene_count, length(available_genes))]
  
  # Test ALL R methods
  methods <- c("ST-COS", "scDesign3", "SRTsim")
  
  for (method in methods) {
    result <- quick_test_method(method, gene_count, gene_subset)
    results <- rbind(results, result)
    gc()
  }
}

# Save quick test results
fwrite(results, file.path(output_dir, "quick_test_results.csv"))

# Summary
cat("\n=== QUICK TEST SUMMARY (R METHODS ONLY) ===\n")
print(results)

successful_tests <- sum(results$success)
total_tests <- nrow(results)
cat(sprintf("Successful tests: %d/%d\n", successful_tests, total_tests))

# Show results by method
summary_by_method <- results[, .(
  tests_run = .N,
  successful = sum(success),
  avg_time = mean(time_seconds, na.rm = TRUE)
), by = method]

cat("\nSummary by method:\n")
print(summary_by_method)

if (successful_tests == total_tests) {
  cat("\n✓ All quick tests passed! Ready to run full efficiency test.\n")
  cat("Run: GENE_COUNTS=50,100,500,1000,3000 Rscript comprehensive_efficiency_test_fixed.R\n")
} else {
  cat("\n⚠ Some tests failed. Check error messages above.\n")
  failed_methods <- unique(results[success == FALSE]$method)
  cat("Failed methods:", paste(failed_methods, collapse = ", "), "\n")
}
