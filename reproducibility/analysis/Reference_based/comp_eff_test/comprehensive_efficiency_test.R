# Comprehensive Computational Efficiency Test - All Methods Working
# Based on successful real_parameter.R implementations
# NOW WITH peakRAM PACKAGE FOR ACCURATE MEMORY MONITORING
# NO SHORTCUTS - ALL COMPUTATIONAL STEPS INCLUDED

library(data.table)
library(Seurat)
library(SRTsim)
library(scDesign3)
library(SingleCellExperiment)
library(mgcv)
library(spdep)
library(parallel)
library(Matrix)
library(peakRAM)  # For accurate memory monitoring

# Load simulator functions
source("./Data/simulator_YZ.R")

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

cat("=== COMPUTATIONAL EFFICIENCY TEST WITH peakRAM MONITORING ===\n")
cat("Testing: ST-COS, SRTsim, and scDesign3 (scCube removed)\n")
cat("Gene counts:", paste(gene_counts, collapse = ", "), "\n")
cat("Using peakRAM package for accurate memory measurement\n")
cat("NO SHORTCUTS - ALL COMPUTATIONAL STEPS INCLUDED\n\n")

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

# TEST 1: ST-COS Implementation with COMPLETE workflow and peakRAM monitoring
test_stcos <- function(gene_count) {
  cat(sprintf("\n=== ST-COS - %d genes (COMPLETE WORKFLOW) ===\n", gene_count))
  
  stcos_computation <- function() {
    cat(sprintf("  Starting COMPLETE ST-COS computation for %d genes\n", gene_count))
    
    # Select genes and extract expression
    selected_genes <- available_genes[1:gene_count]
    expr <- your_data@assays$Spatial$counts[selected_genes, layer_coord$cellnames]
    cat(sprintf("  Data extracted: %d genes x %d spots\n", nrow(expr), ncol(expr)))
    
    # 1. Fit NB parameters for ALL genes
    cat(sprintf("  Step 1: Fitting NB parameters for %d genes\n", gene_count))
    nb_par <- t(apply(expr, 1, fit_nb_fast))
    mu_g <- nb_par[, "mu"]
    theta_g <- nb_par[, "theta"]
    
    # 2. PIT transformation for ALL genes
    cat(sprintf("  Step 2: PIT transformation for %d genes\n", gene_count))
    Z <- t(counts2z_mat(expr, mu_g, theta_g, block = min(1024, gene_count), n_cores = 1))
    
    # 3. Correlation matrix for ALL genes
    cat(sprintf("  Step 3: Computing correlation matrix (%d x %d)\n", gene_count, gene_count))
    thr <- 0.05
    Sigma <- cor(Z, use = "pairwise.complete.obs")
    Sigma[is.na(Sigma)] <- 0
    diag(Sigma) <- 1
    Sigma[abs(Sigma) < thr] <- 0
    Sigma <- makespd(Sigma)
    
    # 4. Library sizes
    lib_size <- colSums(expr)
    lib_factor <- lib_size / mean(lib_size)
    
    # 5. Spatial detection for ALL genes (NO SHORTCUTS)
    cat(sprintf("  Step 5: Spatial analysis for ALL %d genes\n", gene_count))
    coords_mat <- as.matrix(layer_coord[, c("row", "col")])
    Wnb <- knearneigh(coords_mat, k = 6)
    W <- nb2listw(knn2nb(Wnb))
    
    moran_flag <- function(y, listw, thr = 0.15) {
      S0 <- Szero(listw)
      res <- moran(x = y, listw = listw, n = length(y), S0 = S0, zero.policy = TRUE)$I
      res > thr
    }
    
    # Check ALL genes for spatial patterns (no shortcuts)
    spatial_genes <- selected_genes[sapply(selected_genes, function(g)
      moran_flag(expr[g, ], W, thr = 0.15))]
    
    cat(sprintf("  Found %d spatial genes out of %d total genes\n", length(spatial_genes), gene_count))
    
    # 6. Build initial parameters
    mu_mat <- outer(lib_factor, mu_g)
    dimnames(mu_mat) <- list(colnames(expr), rownames(expr))
    
    # 7. GAM refitting for ALL spatial genes (NO SHORTCUTS)
    if (length(spatial_genes) > 0) {
      cat(sprintf("  Step 7: GAM refitting for ALL %d spatial genes (no shortcuts)\n", length(spatial_genes)))
      for (i in seq_along(spatial_genes)) {
        g <- spatial_genes[i]
        if (i %% 50 == 0) {
          cat(sprintf("    Processing spatial gene %d/%d (%s)\n", i, length(spatial_genes), g))
        }
        
        dat <- data.frame(
          y = as.numeric(expr[g, ]),
          row = layer_coord$row,
          col = layer_coord$col,
          log_lib = log(lib_factor)
        )
        
        # Fit GAM for this spatial gene
        gam_fit <- gam(y ~ offset(log_lib) + s(row, col, bs = "tp", k = 20),
                       family = nb(link = "log"), data = dat, method = "REML")
        mu_mat[, g] <- fitted(gam_fit)
        theta_g[g] <- gam_fit$family$getTheta(TRUE)
      }
      cat(sprintf("  Completed GAM refitting for all %d spatial genes\n", length(spatial_genes)))
    } else {
      cat(sprintf("  No spatial genes found - skipping GAM refitting\n"))
    }
    
    # 8. Final simulation with ALL computed parameters
    cat(sprintf("  Step 8: Final simulation with all parameters\n"))
    mean_param <- list(CT = mu_mat)
    gene_relationship <- as.data.table(
      matrix(1, nrow = nrow(expr), ncol = 1,
             dimnames = list(rownames(expr), "CT"))
    )
    cov_gene <- array(Sigma, dim = c(nrow(expr), nrow(expr), 1),
                      dimnames = list(rownames(expr), rownames(expr), "CT"))
    
    sim_list <- generate_true_count(layer_coord,
                                    mean_param,
                                    gene_relationship,
                                    cov_gene,
                                    lib_size = list(CT = lib_size))
    
    sim_spot_counts <- sim_list[["CT"]]
    
    cat(sprintf("  COMPLETE ST-COS computation finished: %d genes, %d spots\n", 
                nrow(sim_spot_counts), ncol(sim_spot_counts)))
    return(sim_spot_counts)
  }
  
  tryCatch({
    # Use peakRAM to monitor memory usage
    start_time <- Sys.time()
    peak_result <- peakRAM(result <- stcos_computation())
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Save results
    write.csv(t(result), 
              file.path(output_dir, paste0("stcos_complete_", gene_count, ".csv")))
    
    cat(sprintf("  SUCCESS - Peak RAM: %.2f MB, Total RAM: %.2f MB, Function RAM: %.2f MB\n",
                peak_result$Peak_RAM_Used_MiB,
                peak_result$Total_RAM_Used_MiB,
                peak_result$Function_RAM_Used_MiB))
    
    return(create_result_dt(
      method = "ST-COS",
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
      method = "ST-COS",
      gene_count = gene_count,
      success = FALSE,
      error_message = as.character(e$message),
      status = "FAILED"
    ))
  })
}

# TEST 2: SRTsim Implementation with COMPLETE workflow and peakRAM monitoring
test_srtsim <- function(gene_count) {
  cat(sprintf("\n=== SRTsim - %d genes (COMPLETE WORKFLOW) ===\n", gene_count))
  
  srtsim_computation <- function() {
    cat(sprintf("  Starting COMPLETE SRTsim computation for %d genes\n", gene_count))
    
    # Prepare data
    selected_genes <- available_genes[1:gene_count]
    example_count <- your_data@assays$Spatial$counts[selected_genes, layer_coord$cellnames]
    cat(sprintf("  Data extracted: %d genes x %d spots\n", nrow(example_count), ncol(example_count)))
    
    example_loc <- layer_coord[,c("row","col","V2")]
    colnames(example_loc) <- c("y","x","label")
    rownames(example_loc) <- layer_coord$cellnames 
    
    # COMPLETE SRTsim workflow - no shortcuts
    cat(sprintf("  Step 1: Creating SRT object\n"))
    simSRT <- createSRT(count_in = example_count, loc_in = example_loc)
    
    cat(sprintf("  Step 2: SRTsim fitting (complete analysis)\n"))
    simSRT1 <- srtsim_fit(simSRT, sim_schem = "tissue")
    
    cat(sprintf("  Step 3: SRTsim count simulation\n"))
    simSRT1 <- srtsim_count(simSRT1)
    
    cat(sprintf("  Step 4: Extracting simulated counts\n"))
    sim_counts <- simCounts(simSRT1)
    
    cat(sprintf("  COMPLETE SRTsim computation finished: %d genes, %d spots\n", 
                nrow(sim_counts), ncol(sim_counts)))
    return(sim_counts)
  }
  
  tryCatch({
    # Use peakRAM to monitor memory usage
    start_time <- Sys.time()
    peak_result <- peakRAM(result <- srtsim_computation())
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Save results
    write.csv(t(result), 
              file.path(output_dir, paste0("srtsim_complete_", gene_count, ".csv")))
    
    cat(sprintf("  SUCCESS - Peak RAM: %.2f MB, Total RAM: %.2f MB, Function RAM: %.2f MB\n",
                peak_result$Peak_RAM_Used_MiB,
                peak_result$Total_RAM_Used_MiB,
                peak_result$Function_RAM_Used_MiB))
    
    return(create_result_dt(
      method = "SRTsim",
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
      method = "SRTsim",
      gene_count = gene_count,
      success = FALSE,
      error_message = as.character(e$message),
      status = "FAILED"
    ))
  })
}

# TEST 3: scDesign3 Implementation with COMPLETE workflow and peakRAM monitoring
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
              file.path(output_dir, paste0("scdesign3_complete_", gene_count, ".csv")))
    
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

# Run COMPLETE comprehensive test
cat("\n=== Starting COMPLETE Efficiency Test with peakRAM Monitoring ===\n")
cat("NO shortcuts or optimizations - all computational steps included\n\n")

for (gene_count in gene_counts) {
  cat(sprintf("\n\n>>> TESTING WITH %d GENES (COMPLETE WORKFLOWS) <<<\n", gene_count))
  
  # Test all methods with COMPLETE workflows using rbindlist with fill=TRUE
  result1 <- test_stcos(gene_count)
  results <- rbindlist(list(results, result1), fill = TRUE)
  gc(verbose = FALSE)  # Memory cleanup
  
  result2 <- test_srtsim(gene_count)
  results <- rbindlist(list(results, result2), fill = TRUE)
  gc(verbose = FALSE)  # Memory cleanup
  
  result3 <- test_scdesign3(gene_count)
  results <- rbindlist(list(results, result3), fill = TRUE)
  gc(verbose = FALSE)  # Memory cleanup
  
  # scCube test removed for this run
  
  # Save intermediate results
  fwrite(results, file.path(output_dir, "complete_efficiency_results.csv"))
  
  cat(sprintf("Completed COMPLETE testing with %d genes\n", gene_count))
}

# Final comprehensive analysis
cat("\n=== COMPLETE COMPUTATIONAL EFFICIENCY RESULTS ===\n")
print(results)

# Enhanced summary statistics with peakRAM results
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
), by = method]

cat("\nCOMPLETE Summary by Method (with peakRAM monitoring):\n")
print(summary_stats)

# Performance analysis by gene count
performance_by_genes <- results[success == TRUE, .(
  avg_time = round(mean(time_seconds, na.rm = TRUE), 2),
  avg_peak_ram = round(mean(peak_ram_mb, na.rm = TRUE), 2),
  avg_function_ram = round(mean(function_ram_mb, na.rm = TRUE), 2),
  methods_successful = .N
), by = gene_count]

cat("\nCOMPLETE Performance by Gene Count:\n")
print(performance_by_genes)

# Comprehensive scalability analysis
if (nrow(results[success == TRUE]) > 0) {
  cat("\n=== COMPLETE SCALABILITY ANALYSIS ===\n")
  
  # Time complexity analysis
  time_scaling <- results[success == TRUE, .(
    time_slope = round(coef(lm(time_seconds ~ gene_count))[2], 6),
    time_r_squared = round(summary(lm(time_seconds ~ gene_count))$r.squared, 4)
  ), by = method]
  
  cat("Time Scaling (seconds per additional gene):\n")
  print(time_scaling)
  
  # RAM scaling analysis
  ram_scaling <- results[success == TRUE & !is.na(function_ram_mb), .(
    ram_slope = round(coef(lm(function_ram_mb ~ gene_count))[2], 6),
    ram_r_squared = round(summary(lm(function_ram_mb ~ gene_count))$r.squared, 4)
  ), by = method]
  
  cat("\nRAM Scaling (MB per additional gene):\n")
  print(ram_scaling)
  
  # Comprehensive efficiency analysis
  efficiency_analysis <- results[success == TRUE, .(
    time_per_gene = round(time_seconds / gene_count, 6),
    ram_per_gene = round(function_ram_mb / gene_count, 6),
    efficiency_score = round(gene_count^2 / (time_seconds * function_ram_mb + 1), 4)
  ), by = .(method, gene_count)]
  
  cat("\nCOMPLETE Efficiency Analysis (per gene):\n")
  print(efficiency_analysis)
}

# Save COMPLETE comprehensive results
fwrite(results, file.path(output_dir, "complete_peakram_efficiency_results.csv"))
fwrite(summary_stats, file.path(output_dir, "complete_method_summary_stats.csv"))
fwrite(performance_by_genes, file.path(output_dir, "complete_performance_by_genes.csv"))

if (exists("time_scaling")) {
  fwrite(time_scaling, file.path(output_dir, "complete_time_scaling_analysis.csv"))
}
if (exists("ram_scaling")) {
  fwrite(ram_scaling, file.path(output_dir, "complete_ram_scaling_analysis.csv"))
}
if (exists("efficiency_analysis")) {
  fwrite(efficiency_analysis, file.path(output_dir, "complete_efficiency_analysis.csv"))
}

cat("\n=== COMPLETE COMPUTATIONAL EFFICIENCY TEST COMPLETED ===\n")
cat("Results saved to:", output_dir, "\n")
cat("COMPLETE efficiency analysis with peakRAM monitoring finished\n")
cat("NO shortcuts taken - all computational steps included for all methods\n")
cat("Files generated:\n")
cat("- complete_peakram_efficiency_results.csv\n")
cat("- complete_method_summary_stats.csv\n")
cat("- complete_performance_by_genes.csv\n")
cat("- complete_time_scaling_analysis.csv\n")
cat("- complete_ram_scaling_analysis.csv\n")
cat("- complete_efficiency_analysis.csv\n")
