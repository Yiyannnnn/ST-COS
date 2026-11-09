#devtools::install_github('https://github.com/cz0316/HEARTSVG',force = T)
library(HEARTSVG)
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
output_file <- paste0("./Data/", S, "_gene/HEARTSVG", suffix, ".csv")

# Load the data using the provided file paths
spot_count <- fread(gene_expression_file)
coord <- fread(coord_file)
colnames(coord) <- c("row", "col")

# Combine the coordinate and spot count data
input <- cbind(coord, spot_count)

# Perform the HEARTSVG analysis
result <- heartsvg(input, scale = FALSE)

# Save the results to the specified file
fwrite(result, output_file)
