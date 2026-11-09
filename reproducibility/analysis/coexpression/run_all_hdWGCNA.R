#!/usr/bin/env Rscript

# Script to run hdWGCNA analysis for all scenarios
# This script will execute the hdWGCNA.R script for all S1-S6 scenarios with replicates 1-5

library(parallel)

# Define scenarios and replicates
scenarios <- paste0("S", 1:6)
replicates <- paste0("_", 1:5)

# Create all combinations
combinations <- expand.grid(scenarios = scenarios, replicates = replicates, stringsAsFactors = FALSE)

# Function to run hdWGCNA for a single combination
run_hdwgcna <- function(scenario, replicate) {
  cat("Running hdWGCNA for", scenario, "replicate", replicate, "\n")
  
  # Construct the command
  cmd <- sprintf("Rscript ./coexpression/hdWGCNA.R %s %s", 
                 scenario, replicate)
  
  # Execute the command
  result <- tryCatch({
    system(cmd, wait = TRUE)
    cat("Completed:", scenario, replicate, "\n")
    return(paste(scenario, replicate, "SUCCESS"))
  }, error = function(e) {
    cat("Error in:", scenario, replicate, "-", e$message, "\n")
    return(paste(scenario, replicate, "ERROR:", e$message))
  })
  
  return(result)
}

# Function to run all combinations sequentially
run_all_sequential <- function() {
  cat("Starting hdWGCNA analysis for all scenarios...\n")
  cat("Total combinations to process:", nrow(combinations), "\n\n")
  
  results <- list()
  
  for (i in 1:nrow(combinations)) {
    scenario <- combinations$scenarios[i]
    replicate <- combinations$replicates[i]
    
    cat(sprintf("Processing %d/%d: %s%s\n", i, nrow(combinations), scenario, replicate))
    
    result <- run_hdwgcna(scenario, replicate)
    results[[i]] <- result
    
    cat("Progress:", i, "/", nrow(combinations), "completed\n\n")
  }
  
  return(results)
}

# Function to run combinations in parallel (optional)
run_all_parallel <- function(n_cores = 4) {
  cat("Starting hdWGCNA analysis for all scenarios in parallel...\n")
  cat("Using", n_cores, "cores\n")
  cat("Total combinations to process:", nrow(combinations), "\n\n")
  
  # Create cluster
  cl <- makeCluster(n_cores)
  
  # Export necessary objects to cluster
  clusterExport(cl, c("run_hdwgcna"))
  
  # Run in parallel
  results <- clusterMap(cl, 
                       run_hdwgcna, 
                       combinations$scenarios, 
                       combinations$replicates,
                       SIMPLIFY = FALSE)
  
  # Stop cluster
  stopCluster(cl)
  
  return(results)
}

# Main execution
main <- function() {
  # Check if hdWGCNA.R exists
  hdwgcna_script <- "./coexpression/hdWGCNA.R"
  if (!file.exists(hdwgcna_script)) {
    stop("hdWGCNA.R script not found at: ", hdwgcna_script)
  }
  
  # Check if output directory exists, create if not
  output_dir <- "./coexpression/hdWGCNA_Agent"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Get command line arguments to decide execution mode
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) > 0 && args[1] == "parallel") {
    n_cores <- ifelse(length(args) > 1, as.numeric(args[2]), 4)
    cat("Running in parallel mode with", n_cores, "cores\n")
    results <- run_all_parallel(n_cores)
  } else {
    cat("Running in sequential mode\n")
    cat("To run in parallel, use: Rscript run_all_hdWGCNA.R parallel [n_cores]\n\n")
    results <- run_all_sequential()
  }
  
  # Save results summary
  results_summary <- data.frame(
    scenario = combinations$scenarios,
    replicate = combinations$replicates,
    result = unlist(results),
    timestamp = Sys.time()
  )
  
  # Write results summary
  summary_file <- file.path(output_dir, paste0("hdWGCNA_run_summary_", 
                                              format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
  write.csv(results_summary, summary_file, row.names = FALSE)
  
  cat("\n=== SUMMARY ===\n")
  cat("Total runs:", nrow(results_summary), "\n")
  cat("Successful runs:", sum(grepl("SUCCESS", results_summary$result)), "\n")
  cat("Failed runs:", sum(grepl("ERROR", results_summary$result)), "\n")
  cat("Results summary saved to:", summary_file, "\n")
  
  # Print any errors
  errors <- results_summary[grepl("ERROR", results_summary$result), ]
  if (nrow(errors) > 0) {
    cat("\nErrors encountered:\n")
    for (i in 1:nrow(errors)) {
      cat("-", errors$scenario[i], errors$replicate[i], ":", errors$result[i], "\n")
    }
  }
  
  cat("\nAll hdWGCNA analyses completed!\n")
  cat("Results are saved in:", output_dir, "\n")
}

# Execute main function
if (!interactive()) {
  main()
}