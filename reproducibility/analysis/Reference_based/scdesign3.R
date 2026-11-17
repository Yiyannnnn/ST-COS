#!/usr/bin/env Rscript
# scDesign3 reference-based simulation (NB marginals, Gaussian-process spatial
# mean, Gaussian copula).
#
# Usage: Rscript scdesign3.R counts.csv coords.csv output_dir [seed]
#
# Input and output follow the shared reference-based format described in
# reproducibility/README.md. The optional `label` column of coords.csv is
# passed to scDesign3 as the cell-type covariate.
# Writes output_dir/simulated_counts.csv in the same layout as counts.csv.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(scDesign3)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript scdesign3.R counts.csv coords.csv output_dir [seed]")
}
seed <- if (length(args) >= 4L) as.integer(args[4L]) else 20260922L

counts_dt <- fread(args[1L], check.names = FALSE)
counts <- as.matrix(counts_dt[, -1L, with = FALSE])
storage.mode(counts) <- "double"
rownames(counts) <- as.character(counts_dt[[1L]])
coords <- fread(args[2L], check.names = FALSE)
setnames(coords, 1L, "id")
coords <- coords[match(colnames(counts), as.character(coords$id))]
if (anyNA(coords$id)) stop("Every observation in counts.csv needs a row in coords.csv.")
label <- if ("label" %in% names(coords)) as.character(coords$label) else rep("CT", nrow(coords))

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = DataFrame(spatial1 = coords$row, spatial2 = coords$col,
                      cell_type = label, row.names = colnames(counts))
)

set.seed(seed)
fit <- scdesign3(
  sce = sce, assay_use = "counts", celltype = "cell_type", pseudotime = NULL,
  spatial = c("spatial1", "spatial2"), other_covariates = NULL,
  mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 200)", sigma_formula = "1",
  family_use = "nb", n_cores = 1, usebam = FALSE, corr_formula = "1",
  copula = "gaussian", DT = TRUE, pseudo_obs = FALSE, return_model = FALSE,
  nonzerovar = FALSE, parallelization = "pbmcapply"
)
simulated <- as.matrix(fit$new_count)
if (!all(rownames(counts) %in% rownames(simulated))) simulated <- t(simulated)
simulated <- simulated[rownames(counts), colnames(counts), drop = FALSE]

dir.create(args[3L], recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(gene = rownames(simulated), as.data.table(simulated)),
       file.path(args[3L], "simulated_counts.csv"))
writeLines(capture.output(sessionInfo()), file.path(args[3L], "session_info.txt"))
