
library(SPARK)
library(data.table)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure the correct number of arguments are provided
if (length(args) != 2) {
  stop("Please provide the two arguments: S and suffix (_N).")
}

# Get the inputs from command-line arguments
S <- args[1]
suffix <- args[2]

# Construct file paths based on the arguments
gene_expression_file <- paste0("./Data/", S, "_gene/simulated/", S, "_gene_expression", suffix, ".csv")
coord_file <- paste0("./Data/", S, "_gene/simulated/", S, "_coord", suffix, ".csv")
spark_output_file <- paste0("./Data/", S, "_gene/SPARK", suffix, ".csv")
sparkx_output_file <- paste0("./Data/", S, "_gene/SPARKX", suffix, ".csv")

# Load the data using the provided file paths
rawcount <- t(fread(gene_expression_file))
info <- fread(coord_file)
colnames(info) <- c("x", "y")
colnames(rawcount) <- 1:ncol(rawcount)

# SPARK analysis
spark <- CreateSPARKObject(counts = rawcount, 
                           location = info[, 1:2], 
                           percentage = 0.1, 
                           min_total_counts = 10)

# Calculate total counts for each cell/spot
spark@lib_size <- apply(spark@counts, 2, sum)

# Variance component estimation
spark <- spark.vc(spark, 
                  covariates = NULL, 
                  lib_size = spark@lib_size, 
                  num_core = 1, 
                  verbose = FALSE)

# Hypothesis testing
spark <- spark.test(spark, 
                    check_positive = TRUE, 
                    verbose = FALSE)

# Save SPARK results
spark_output <- cbind(gene = rownames(spark@res_mtest), spark@res_mtest[, c("combined_pvalue", "adjusted_pvalue")])
fwrite(spark_output, spark_output_file)

# SPARKX analysis
sp_count <- t(fread(gene_expression_file))
pos <- fread(coord_file)
colnames(pos) <- c("x", "y")
sparkX <- sparkx(sp_count, as.matrix(pos), numCores = 1, option = "mixture")

# Save SPARKX results
sparkx_output <- cbind(gene = rownames(sparkX$res_mtest), sparkX$res_mtest)
fwrite(sparkx_output, sparkx_output_file)


