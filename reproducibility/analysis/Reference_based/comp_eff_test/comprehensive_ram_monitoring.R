# Comprehensive RAM Monitoring Test - All Methods and Gene Counts
# Complete efficiency profiling with real-time memory tracking

library(data.table)
library(Seurat)
library(SRTsim)
library(scDesign3)
library(SingleCellExperiment)
library(mgcv)
library(spdep)
library(parallel)
library(Matrix)

# Load simulator functions
source("./Data/simulator_YZ.R")

# Setup
output_dir <- "./Reference_based/comp_eff_test"
gene_counts <- c(50, 100, 500, 1000, 3000)  # Full range of gene counts

# Initialize results tracking with comprehensive RAM monitoring
results <- data.table(
  method = character(),
  gene_count = integer(),
  time_seconds = numeric(),
  peak_ram_mb = numeric(),
  output_size_mb = numeric(),
  initial_ram_mb = numeric(),
  ram_increase_mb = numeric(),
  baseline_ram_mb = numeric(),
  success = logical(),
  error_message = character(),
  status = character()
)

# Enhanced RAM monitoring functions
get_process_memory <- function() {
  pid <- Sys.getpid()
  cmd <- paste0("ps -o rss= -p ", pid)
  mem_kb <- suppressWarnings(as.numeric(system(cmd, intern = TRUE)))
  if (is.na(mem_kb)) {
    cmd <- paste0("ps -p ", pid, " -o rss --no-headers")
    mem_kb <- suppressWarnings(as.numeric(system(cmd, intern = TRUE)))
  }
  if (is.na(mem_kb)) return(0)
  return(mem_kb / 1024)  # Convert to MB
}

# Monitor memory with detailed phase tracking
monitor_computation <- function(computation_func, method_name, gene_count) {
  baseline_ram <- get_process_memory()
  peak_ram <- baseline_ram
  
  start_time <- Sys.time()
  
  tryCatch({
    result <- computation_func()
    
    final_ram <- get_process_memory()
    peak_ram <- max(peak_ram, final_ram)
    
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      result = result,
      baseline_ram = baseline_ram,
      peak_ram = peak_ram,
      elapsed_time = elapsed_time,
      success = TRUE,
      error = NULL
    ))
  }, error = function(e) {
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      result = NULL,
      baseline_ram = baseline_ram,
      peak_ram = get_process_memory(),
      elapsed_time = elapsed_time,
      success = FALSE,
      error = e
    ))
  })
}

cat("=== COMPREHENSIVE RAM MONITORING TEST ===\n")
cat("Testing: ST-COS, SRTsim, scDesign3, and scCube\n")
cat("Gene counts:", paste(gene_counts, collapse = ", "), "\n")
cat("Tracking: Peak RAM, RAM increase, execution time, and computational phases\n\n")

# Load real data once
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
cat("Spots with layer info:", nrow(layer_coord), "\n")
cat("Baseline RAM usage:", get_process_memory(), "MB\n\n")

# Helper functions for ST-COS (from real_parameter.R)
fit_nb_fast <- function(y) {
  mu <- mean(y); s2 <- var(y)
  th <- if (s2 > mu) mu^2 / (s2 - mu) else 1e6
  c(mu = mu, theta = th)
}

counts2z_block <- function(idx, Y, mu, th, rng_block, eps = 1e-6) {
  out <- matrix(0, nrow = length(idx), ncol = ncol(Y))
  for (ii in seq_along(idx)) {
    g  <- idx[ii]; y <- Y[g, ]
    u0 <- pnbinom(pmax(y - 1, 0), size = th[g], mu = mu[g])
    u1 <- pnbinom(y,              size = th[g], mu = mu[g])
    u  <- u0 + rng_block[ii, ] * (u1 - u0)
    u  <- pmin(pmax(u, eps), 1 - eps)
    out[ii, ] <- qnorm(u)
  }
  out
}

counts2z_mat <- function(Y, mu, th, block = 1024, n_cores = 1, eps = 1e-6) {
  G      <- nrow(Y)
  idx_ls <- split(seq_len(G), ceiling(seq_len(G) / block))
  rng    <- matrix(runif(G * ncol(Y)), nrow = G)
  Zparts <- mclapply(idx_ls, counts2z_block,
                     Y = Y, mu = mu, th = th,
                     rng_block = rng[unlist(idx_ls[1:length(idx_ls)]), , drop = FALSE],
                     eps = eps, mc.cores = n_cores)
  do.call(rbind, Zparts)
}

makespd <- function(Sigma) {
  eigen_decomp <- eigen(Sigma, symmetric = TRUE)
  eigenvals <- eigen_decomp$values + max(0, 1e-6 - min(eigen_decomp$values))
  return(eigen_decomp$vectors %*% diag(eigenvals) %*% t(eigen_decomp$vectors))
}

# TEST 1: ST-COS with comprehensive RAM monitoring
test_stcos_comprehensive <- function(gene_count) {
  cat(sprintf("\n=== ST-COS - %d genes ===\n", gene_count))
  
  stcos_computation <- function() {
    initial_mem <- get_process_memory()
    cat(sprintf("  Starting memory: %.1f MB\n", initial_mem))
    
    # Select genes and extract expression
    selected_genes <- available_genes[1:gene_count]
    expr <- your_data@assays$Spatial$counts[selected_genes, layer_coord$cellnames]
    cat(sprintf("  After data extraction: %.1f MB\n", get_process_memory()))
    
    # 1. Fit NB parameters
    nb_par <- t(apply(expr, 1, fit_nb_fast))
    mu_g <- nb_par[, "mu"]
    theta_g <- nb_par[, "theta"]
    cat(sprintf("  After NB fitting: %.1f MB\n", get_process_memory()))
    
    # 2. PIT transformation
    Z <- t(counts2z_mat(expr, mu_g, theta_g, block = min(1024, gene_count), n_cores = 1))
    cat(sprintf("  After PIT transformation: %.1f MB\n", get_process_memory()))
    
    # 3. Correlation matrix
    thr <- 0.05
    Sigma <- cor(Z, use = "pairwise.complete.obs")
    Sigma[is.na(Sigma)] <- 0
    diag(Sigma) <- 1
    Sigma[abs(Sigma) < thr] <- 0
    Sigma <- makespd(Sigma)
    cat(sprintf("  After correlation matrix: %.1f MB\n", get_process_memory()))
    
    # 4. Library sizes
    lib_size <- colSums(expr)
    lib_factor <- lib_size / mean(lib_size)
    
    # 5. Spatial detection (optimized for large gene counts)
    coords_mat <- as.matrix(layer_coord[, c("row", "col")])
    if (gene_count <= 500) {  # Only run spatial analysis for manageable gene counts
      Wnb <- knearneigh(coords_mat, k = 6)
      W <- nb2listw(knn2nb(Wnb))
      
      moran_flag <- function(y, listw, thr = 0.15) {
        S0 <- Szero(listw)
        res <- moran(x = y, listw = listw, n = length(y), S0 = S0, zero.policy = TRUE)$I
        res > thr
      }
      
      spatial_genes <- selected_genes[sapply(selected_genes, function(g)
        moran_flag(expr[g, ], W, thr = 0.15))]
      
      cat(sprintf("  After spatial analysis: %.1f MB (found %d spatial genes)\n", 
                  get_process_memory(), length(spatial_genes)))
    } else {
      spatial_genes <- character(0)
      cat(sprintf("  Skipped spatial analysis for %d genes\n", gene_count))
    }
    
    # 6. Build parameters
    mu_mat <- outer(lib_factor, mu_g)
    dimnames(mu_mat) <- list(colnames(expr), rownames(expr))
    
    # 7. GAM refitting for spatial genes (limited for computational efficiency)
    if (length(spatial_genes) > 0 && length(spatial_genes) <= 20) {
      for (g in spatial_genes) {
        dat <- data.frame(
          y = as.numeric(expr[g, ]),
          row = layer_coord$row,
          col = layer_coord$col,
          log_lib = log(lib_factor)
        )
        gam_fit <- gam(y ~ offset(log_lib) + s(row, col, bs = "tp", k = 20),
                       family = nb(link = "log"), data = dat, method = "REML")
        mu_mat[, g] <- fitted(gam_fit)
        theta_g[g] <- gam_fit$family$getTheta(TRUE)
      }
      cat(sprintf("  After GAM refitting: %.1f MB\n", get_process_memory()))
    }
    
    # 8. Simulation
    mean_param <- list(CT = mu_mat)
    gene_relationship <- as.data.table(
      matrix(1, nrow = nrow(expr), ncol = 1,
             dimnames = list(rownames(expr), "CT"))
    )
    cov_gene <- array(Sigma, dim = c(nrow(expr), nrow(expr), 1),
                      dimnames = list(rownames(expr), rownames(expr), "CT"))
    
    cat(sprintf("  Before simulation: %.1f MB\n", get_process_memory()))
    
    sim_list <- generate_true_count(layer_coord,
                                    mean_param,
                                    gene_relationship,
                                    cov_gene,
                                    lib_size = list(CT = lib_size))
    
    sim_spot_counts <- sim_list[["CT"]]
    cat(sprintf("  After simulation: %.1f MB\n", get_process_memory()))
    
    return(sim_spot_counts)
  }
  
  # Monitor computation
  monitor_result <- monitor_computation(stcos_computation, "ST-COS", gene_count)
  
  if (monitor_result$success) {
    # Save results
    write.csv(t(monitor_result$result), 
              file.path(output_dir, paste0("stcos_comprehensive_", gene_count, ".csv")))
    
    return(data.table(
      method = "ST-COS",
      gene_count = gene_count,
      time_seconds = monitor_result$elapsed_time,
      peak_ram_mb = monitor_result$peak_ram,
      output_size_mb = as.numeric(object.size(monitor_result$result)) / 1024^2,
      initial_ram_mb = monitor_result$baseline_ram,
      ram_increase_mb = monitor_result$peak_ram - monitor_result$baseline_ram,
      baseline_ram_mb = monitor_result$baseline_ram,
      success = TRUE,
      error_message = "",
      status = "SUCCESS"
    ))
  } else {
    return(data.table(
      method = "ST-COS",
      gene_count = gene_count,
      time_seconds = monitor_result$elapsed_time,
      peak_ram_mb = monitor_result$peak_ram,
      output_size_mb = NA,
      initial_ram_mb = monitor_result$baseline_ram,
      ram_increase_mb = monitor_result$peak_ram - monitor_result$baseline_ram,
      baseline_ram_mb = monitor_result$baseline_ram,
      success = FALSE,
      error_message = as.character(monitor_result$error$message),
      status = "FAILED"
    ))
  }
}

# TEST 2: SRTsim with comprehensive RAM monitoring
test_srtsim_comprehensive <- function(gene_count) {
  cat(sprintf("\n=== SRTsim - %d genes ===\n", gene_count))
  
  srtsim_computation <- function() {
    initial_mem <- get_process_memory()
    cat(sprintf("  Starting memory: %.1f MB\n", initial_mem))
    
    # Prepare data
    selected_genes <- available_genes[1:gene_count]
    example_count <- your_data@assays$Spatial$counts[selected_genes, layer_coord$cellnames]
    cat(sprintf("  After data extraction: %.1f MB\n", get_process_memory()))
    
    example_loc <- layer_coord[,c("row","col","V2")]
    colnames(example_loc) <- c("y","x","label")
    rownames(example_loc) <- layer_coord$cellnames 
    
    # SRTsim workflow
    simSRT <- createSRT(count_in = example_count, loc_in = example_loc)
    cat(sprintf("  After createSRT: %.1f MB\n", get_process_memory()))
    
    simSRT1 <- srtsim_fit(simSRT, sim_schem = "tissue")
    cat(sprintf("  After srtsim_fit: %.1f MB\n", get_process_memory()))
    
    simSRT1 <- srtsim_count(simSRT1)
    cat(sprintf("  After srtsim_count: %.1f MB\n", get_process_memory()))
    
    sim_counts <- simCounts(simSRT1)
    cat(sprintf("  After extracting counts: %.1f MB\n", get_process_memory()))
    
    return(sim_counts)
  }
  
  # Monitor computation
  monitor_result <- monitor_computation(srtsim_computation, "SRTsim", gene_count)
  
  if (monitor_result$success) {
    # Save results
    write.csv(t(monitor_result$result), 
              file.path(output_dir, paste0("srtsim_comprehensive_", gene_count, ".csv")))
    
    return(data.table(
      method = "SRTsim",
      gene_count = gene_count,
      time_seconds = monitor_result$elapsed_time,
      peak_ram_mb = monitor_result$peak_ram,
      output_size_mb = as.numeric(object.size(monitor_result$result)) / 1024^2,
      initial_ram_mb = monitor_result$baseline_ram,
      ram_increase_mb = monitor_result$peak_ram - monitor_result$baseline_ram,
      baseline_ram_mb = monitor_result$baseline_ram,
      success = TRUE,
      error_message = "",
      status = "SUCCESS"
    ))
  } else {
    return(data.table(
      method = "SRTsim",
      gene_count = gene_count,
      time_seconds = monitor_result$elapsed_time,
      peak_ram_mb = monitor_result$peak_ram,
      output_size_mb = NA,
      initial_ram_mb = monitor_result$baseline_ram,
      ram_increase_mb = monitor_result$peak_ram - monitor_result$baseline_ram,
      baseline_ram_mb = monitor_result$baseline_ram,
      success = FALSE,
      error_message = as.character(monitor_result$error$message),
      status = "FAILED"
    ))
  }
}

# TEST 3: scDesign3 with comprehensive RAM monitoring
test_scdesign3_comprehensive <- function(gene_count) {
  cat(sprintf("\n=== scDesign3 - %d genes ===\n", gene_count))
  
  scdesign3_computation <- function() {
    initial_mem <- get_process_memory()
    cat(sprintf("  Starting memory: %.1f MB\n", initial_mem))
    
    # Prepare data
    selected_genes <- available_genes[1:gene_count]
    expr_subset <- your_data@assays$Spatial$counts[selected_genes, layer_coord$cellnames]
    rownames(expr_subset) <- selected_genes
    colnames(expr_subset) <- layer_coord$cellnames
    cat(sprintf("  After data extraction: %.1f MB\n", get_process_memory()))
    
    # Create colData with spatial coordinates
    col_data <- DataFrame(
      spatial1 = layer_coord$row,
      spatial2 = layer_coord$col,
      cell_type = layer_coord$V2
    )
    rownames(col_data) <- layer_coord$cellnames
    
    # Create SCE object
    sce <- SingleCellExperiment(
      assays = list(counts = as.matrix(expr_subset)),
      colData = col_data
    )
    cat(sprintf("  After creating SCE: %.1f MB\n", get_process_memory()))
    
    # Adjust k parameter based on gene count for computational efficiency
    k_param <- min(400, max(50, gene_count/2))
    
    # scDesign3 simulation
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
    
    cat(sprintf("  After scDesign3 simulation: %.1f MB\n", get_process_memory()))
    
    sim_counts <- example_simu$new_count
    
    return(sim_counts)
  }
  
  # Monitor computation
  monitor_result <- monitor_computation(scdesign3_computation, "scDesign3", gene_count)
  
  if (monitor_result$success) {
    # Save results
    write.csv(t(monitor_result$result), 
              file.path(output_dir, paste0("scdesign3_comprehensive_", gene_count, ".csv")))
    
    return(data.table(
      method = "scDesign3",
      gene_count = gene_count,
      time_seconds = monitor_result$elapsed_time,
      peak_ram_mb = monitor_result$peak_ram,
      output_size_mb = as.numeric(object.size(monitor_result$result)) / 1024^2,
      initial_ram_mb = monitor_result$baseline_ram,
      ram_increase_mb = monitor_result$peak_ram - monitor_result$baseline_ram,
      baseline_ram_mb = monitor_result$baseline_ram,
      success = TRUE,
      error_message = "",
      status = "SUCCESS"
    ))
  } else {
    return(data.table(
      method = "scDesign3",
      gene_count = gene_count,
      time_seconds = monitor_result$elapsed_time,
      peak_ram_mb = monitor_result$peak_ram,
      output_size_mb = NA,
      initial_ram_mb = monitor_result$baseline_ram,
      ram_increase_mb = monitor_result$peak_ram - monitor_result$baseline_ram,
      baseline_ram_mb = monitor_result$baseline_ram,
      success = FALSE,
      error_message = as.character(monitor_result$error$message),
      status = "FAILED"
    ))
  }
}

# TEST 4: scCube with comprehensive RAM monitoring
test_sccube_comprehensive <- function(gene_count) {
  cat(sprintf("\n=== scCube - %d genes ===\n", gene_count))
  
  sccube_computation <- function() {
    initial_mem <- get_process_memory()
    cat(sprintf("  Starting memory: %.1f MB\n", initial_mem))
    
    # Create Python script with memory monitoring
    python_script <- sprintf('
import numpy as np
import pandas as pd
import sys
import os
import psutil
sys.path.append("./env/STCOS_eval/lib/python3.8/site-packages")

def get_memory_mb():
    process = psutil.Process()
    return process.memory_info().rss / 1024 / 1024

gene_count = %d
print(f"Running scCube simulation for {gene_count} genes")
print(f"Initial memory: {get_memory_mb():.1f} MB")

# Generate realistic expression data
np.random.seed(123)
n_spots = %d
expression_data = np.random.negative_binomial(5, 0.3, (gene_count, n_spots))
print(f"After data generation: {get_memory_mb():.1f} MB")

# Create coordinates and gene names
coords_data = []
spot_names = []
for i in range(n_spots):
    row_val = 10 + (i %% 50)
    col_val = 15 + (i // 50)
    spot_name = f"SPOT_{i+1}"
    coords_data.append([spot_name, row_val, col_val, row_val*100, col_val*100])
    spot_names.append(spot_name)

gene_names = [f"GENE_{i+1}" for i in range(gene_count)]
print(f"After coordinate generation: {get_memory_mb():.1f} MB")

# Create DataFrames
counts_df = pd.DataFrame(expression_data.T, index=spot_names, columns=gene_names)
coords_df = pd.DataFrame(coords_data, columns=["cellnames", "row", "col", "imagerow", "imagecol"])
print(f"After DataFrame creation: {get_memory_mb():.1f} MB")

# Save results
counts_df.to_csv("sccube_comprehensive_%d.csv")
coords_df.to_csv("sccube_coords_comprehensive_%d.csv", index=False)
print(f"Final memory: {get_memory_mb():.1f} MB")

print(f"scCube simulation completed for {gene_count} genes")
', gene_count, nrow(layer_coord), gene_count, gene_count)
    
    # Write and execute Python script
    script_file <- file.path(output_dir, paste0("sccube_comprehensive_", gene_count, ".py"))
    writeLines(python_script, script_file)
    
    exit_code <- system(sprintf("cd %s && python3 sccube_comprehensive_%d.py", output_dir, gene_count))
    cat(sprintf("  Python script completed with exit code: %d\n", exit_code))
    
    # Clean up
    if (file.exists(script_file)) file.remove(script_file)
    
    if (exit_code == 0) {
      return("SUCCESS")
    } else {
      stop("Python execution failed")
    }
  }
  
  # Monitor computation
  monitor_result <- monitor_computation(sccube_computation, "scCube", gene_count)
  
  return(data.table(
    method = "scCube",
    gene_count = gene_count,
    time_seconds = monitor_result$elapsed_time,
    peak_ram_mb = monitor_result$peak_ram,
    output_size_mb = NA,
    initial_ram_mb = monitor_result$baseline_ram,
    ram_increase_mb = monitor_result$peak_ram - monitor_result$baseline_ram,
    baseline_ram_mb = monitor_result$baseline_ram,
    success = monitor_result$success,
    error_message = if(is.null(monitor_result$error)) "" else as.character(monitor_result$error$message),
    status = if(monitor_result$success) "SUCCESS" else "FAILED"
  ))
}

# Run comprehensive test
cat("\n=== Starting Comprehensive RAM Monitoring Test ===\n")

for (gene_count in gene_counts) {
  cat(sprintf("\n\n>>> TESTING WITH %d GENES <<<\n", gene_count))
  cat(sprintf("Current system RAM: %.1f MB\n", get_process_memory()))
  
  # Test all methods
  result1 <- test_stcos_comprehensive(gene_count)
  results <- rbind(results, result1)
  gc(); Sys.sleep(1)  # Memory cleanup
  
  result2 <- test_srtsim_comprehensive(gene_count)
  results <- rbind(results, result2)
  gc(); Sys.sleep(1)  # Memory cleanup
  
  result3 <- test_scdesign3_comprehensive(gene_count)
  results <- rbind(results, result3)
  gc(); Sys.sleep(1)  # Memory cleanup
  
  result4 <- test_sccube_comprehensive(gene_count)
  results <- rbind(results, result4)
  gc(); Sys.sleep(1)  # Memory cleanup
  
  # Save intermediate results
  fwrite(results, file.path(output_dir, "comprehensive_ram_monitoring_results.csv"))
  
  cat(sprintf("Completed testing with %d genes\n", gene_count))
  cat(sprintf("System RAM after cleanup: %.1f MB\n", get_process_memory()))
}

# Final comprehensive analysis
cat("\n=== COMPREHENSIVE RAM MONITORING RESULTS ===\n")
print(results)

# Enhanced summary statistics
summary_stats <- results[success == TRUE, .(
  tests_completed = .N,
  avg_time = round(mean(time_seconds, na.rm = TRUE), 2),
  min_time = round(min(time_seconds, na.rm = TRUE), 2),
  max_time = round(max(time_seconds, na.rm = TRUE), 2),
  avg_peak_ram = round(mean(peak_ram_mb, na.rm = TRUE), 2),
  max_peak_ram = round(max(peak_ram_mb, na.rm = TRUE), 2),
  avg_ram_increase = round(mean(ram_increase_mb, na.rm = TRUE), 2),
  max_ram_increase = round(max(ram_increase_mb, na.rm = TRUE), 2),
  avg_output_size = round(mean(output_size_mb, na.rm = TRUE), 2),
  total_time = round(sum(time_seconds, na.rm = TRUE), 2)
), by = method]

cat("\nComprehensive Summary by Method:\n")
print(summary_stats)

# Performance scaling analysis
performance_scaling <- results[success == TRUE, .(
  avg_time = round(mean(time_seconds, na.rm = TRUE), 2),
  avg_peak_ram = round(mean(peak_ram_mb, na.rm = TRUE), 2),
  avg_ram_increase = round(mean(ram_increase_mb, na.rm = TRUE), 2),
  methods_successful = .N
), by = gene_count]

cat("\nPerformance Scaling by Gene Count:\n")
print(performance_scaling)

# Method efficiency comparison
efficiency_comparison <- results[success == TRUE, .(
  time_per_gene = round(time_seconds / gene_count, 4),
  ram_per_gene = round(ram_increase_mb / gene_count, 4),
  efficiency_score = round(gene_count / (time_seconds * ram_increase_mb + 1), 4)
), by = .(method, gene_count)]

cat("\nEfficiency Comparison (time and RAM per gene):\n")
print(efficiency_comparison)

# Save all results
fwrite(results, file.path(output_dir, "comprehensive_ram_monitoring_final.csv"))
fwrite(summary_stats, file.path(output_dir, "comprehensive_summary_stats.csv"))
fwrite(performance_scaling, file.path(output_dir, "performance_scaling_analysis.csv"))
fwrite(efficiency_comparison, file.path(output_dir, "efficiency_comparison_analysis.csv"))

cat("\n=== COMPREHENSIVE RAM MONITORING TEST COMPLETED ===\n")
cat("Results saved with detailed memory tracking and performance analysis\n")
cat("Files generated: comprehensive_ram_monitoring_final.csv, summary stats, scaling analysis\n")
