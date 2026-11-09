# Test with Small Gene Set - Demonstrates Working Framework
# Use this to run a quick efficiency test with just a few genes

# Libraries
library(data.table)
library(Seurat)

# Setup
output_dir <- "./Reference_based/comp_eff_test"
gene_counts <- c(10, 20)  # Very small for demo

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

# Function to run efficiency test for R methods
run_r_method <- function(method_name, gene_count, gene_subset) {
  cat(sprintf("\n=== Testing %s with %d genes ===\n", method_name, gene_count))
  
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
      memory_mb = 0,  # Simplified for demo
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

# Demo testing loop
for (gene_count in gene_counts) {
  cat(sprintf("\n\n>>> DEMO TESTING WITH %d GENES <<<\n", gene_count))
  
  # Select gene subset
  gene_subset <- available_genes[1:min(gene_count, length(available_genes))]
  
  # Test only working R methods for demo
  methods <- c("ST-COS")  # Start with your method, add others as they're verified
  
  for (method in methods) {
    result <- run_r_method(method, gene_count, gene_subset)
    results <- rbind(results, result)
    gc()
  }
}

# Save demo results
fwrite(results, file.path(output_dir, "demo_results.csv"))

cat("\n=== DEMO COMPLETED ===\n")
print(results)

successful_tests <- sum(results$success)
total_tests <- nrow(results)
cat(sprintf("\nSuccessful tests: %d/%d\n", successful_tests, total_tests))

if (successful_tests > 0) {
  cat("\n✓ Framework is working! You can now run the full test.\n")
  cat("For full test: Rscript run_efficiency_test.R\n")
  cat("For smaller test: modify gene_counts in run_efficiency_test.R\n")
}