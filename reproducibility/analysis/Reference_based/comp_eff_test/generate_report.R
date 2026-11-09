# Report Generation Script for Computational Efficiency Testing

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

generate_efficiency_report <- function(results, output_dir) {
  
  cat("Generating computational efficiency report...\n")
  
  # Load and combine scCube timing results if they exist
  sccube_files <- list.files(output_dir, pattern = "sccube_timing_.*\\.csv", full.names = TRUE)
  if (length(sccube_files) > 0) {
    sccube_results <- rbindlist(lapply(sccube_files, fread))
    sccube_results[, memory_mb := NA]  # Memory not tracked for Python script
    sccube_results[, error_message := ""]
    
    # Combine with R results
    results <- rbind(results, sccube_results[, .(method, gene_count, time_seconds, memory_mb, success, error_message)])
  }
  
  # Filter successful runs
  successful_results <- results[success == TRUE]
  
  # Summary statistics
  summary_stats <- successful_results[, .(
    mean_time = mean(time_seconds, na.rm = TRUE),
    median_time = median(time_seconds, na.rm = TRUE),
    min_time = min(time_seconds, na.rm = TRUE),
    max_time = max(time_seconds, na.rm = TRUE),
    mean_memory = mean(memory_mb, na.rm = TRUE),
    median_memory = median(memory_mb, na.rm = TRUE),
    runs_completed = .N
  ), by = .(method, gene_count)]
  
  # Save summary
  fwrite(summary_stats, file.path(output_dir, "efficiency_summary.csv"))
  
  # Create visualizations
  pdf(file.path(output_dir, "efficiency_comparison.pdf"), width = 16, height = 12)
  
  # Time comparison plot
  if (nrow(successful_results) > 0) {
    p1 <- ggplot(successful_results, aes(x = gene_count, y = time_seconds, color = method)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      scale_x_log10() +
      scale_y_log10() +
      labs(title = "Computation Time vs Number of Genes",
           x = "Number of Genes (log scale)",
           y = "Time (seconds, log scale)",
           color = "Method") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    print(p1)
    
    # Memory comparison plot (excluding scCube since memory not tracked)
    memory_results <- successful_results[!is.na(memory_mb)]
    if (nrow(memory_results) > 0) {
      p2 <- ggplot(memory_results, aes(x = gene_count, y = memory_mb, color = method)) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        scale_x_log10() +
        scale_y_log10() +
        labs(title = "Memory Usage vs Number of Genes",
             x = "Number of Genes (log scale)",
             y = "Memory (MB, log scale)",
             color = "Method") +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      print(p2)
    }
    
    # Bar plot for latest run
    latest_results <- successful_results[gene_count == max(gene_count)]
    if (nrow(latest_results) > 0) {
      p3 <- ggplot(latest_results, aes(x = method, y = time_seconds, fill = method)) +
        geom_bar(stat = "identity") +
        labs(title = paste("Computation Time for", max(latest_results$gene_count), "Genes"),
             x = "Method",
             y = "Time (seconds)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
      
      print(p3)
    }
    
    # Scaling efficiency plot
    scaling_data <- successful_results[, .(
      scaling_factor = time_seconds / time_seconds[gene_count == min(gene_count)],
      gene_ratio = gene_count / min(gene_count)
    ), by = method]
    
    if (nrow(scaling_data) > 0) {
      p4 <- ggplot(scaling_data, aes(x = gene_ratio, y = scaling_factor, color = method)) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
        labs(title = "Scaling Efficiency (compared to linear scaling)",
             x = "Gene Count Ratio",
             y = "Time Scaling Factor",
             color = "Method") +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      print(p4)
    }
  }
  
  dev.off()
  
  # Generate text report
  report_text <- paste0(
    "=== COMPUTATIONAL EFFICIENCY TESTING REPORT ===\n\n",
    "Test Configuration:\n",
    "- Dataset: 151673 (DLPFC)\n",
    "- Gene counts tested: ", paste(unique(results$gene_count), collapse = ", "), "\n",
    "- Methods tested: ", paste(unique(results$method), collapse = ", "), "\n\n",
    
    "=== SUMMARY STATISTICS ===\n\n"
  )
  
  for (method in unique(summary_stats$method)) {
    method_stats <- summary_stats[method == method]
    report_text <- paste0(report_text, 
      method, ":\n",
      "  Successful runs: ", sum(method_stats$runs_completed), "\n"
    )
    
    if (nrow(method_stats) > 0) {
      report_text <- paste0(report_text,
        "  Time range: ", round(min(method_stats$min_time), 2), " - ", 
        round(max(method_stats$max_time), 2), " seconds\n"
      )
      
      if (!all(is.na(method_stats$mean_memory))) {
        report_text <- paste0(report_text,
          "  Memory range: ", round(min(method_stats$mean_memory, na.rm = TRUE), 2), " - ", 
          round(max(method_stats$mean_memory, na.rm = TRUE), 2), " MB\n"
        )
      }
    }
    report_text <- paste0(report_text, "\n")
  }
  
  # Failed runs
  failed_runs <- results[success == FALSE]
  if (nrow(failed_runs) > 0) {
    report_text <- paste0(report_text, 
      "=== FAILED RUNS ===\n\n"
    )
    
    for (i in 1:nrow(failed_runs)) {
      report_text <- paste0(report_text,
        failed_runs$method[i], " (", failed_runs$gene_count[i], " genes): ",
        failed_runs$error_message[i], "\n"
      )
    }
  }
  
  # Write report
  writeLines(report_text, file.path(output_dir, "efficiency_report.txt"))
  
  cat("Report generated successfully!\n")
  cat("Files created:\n")
  cat("- efficiency_summary.csv\n")
  cat("- efficiency_comparison.pdf\n") 
  cat("- efficiency_report.txt\n")
  
  return(summary_stats)
}