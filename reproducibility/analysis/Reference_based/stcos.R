#!/usr/bin/env Rscript
# ST-COS reference-based simulation.
#
# Usage: Rscript stcos.R counts.csv coords.csv output_dir [seed]
#
# Input and output follow the shared reference-based format described in
# reproducibility/README.md: counts.csv is genes x observations (first column
# gene IDs), coords.csv has the observation ID first plus `row` and `col`.
# Writes output_dir/simulated_counts.csv in the same layout as counts.csv.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript stcos.R counts.csv coords.csv output_dir [seed]")
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
locations <- data.frame(row = coords$row, col = coords$col, row.names = colnames(counts))

set.seed(seed)
fit <- STCOS::get_real_param(counts, locations, rownames(counts))
simulated <- STCOS::generate_true_count(
  coord = locations,
  mean_param = fit$mean_param,
  gene_relationship = fit$gene_relationship,
  cov_gene = fit$cov_gene,
  theta_g = fit$theta_g,
  library_size = fit$library_size
)$CT

dir.create(args[3L], recursive = TRUE, showWarnings = FALSE)
out <- data.table(gene = rownames(counts), t(simulated[colnames(counts), rownames(counts)]))
setnames(out, c("gene", colnames(counts)))
fwrite(out, file.path(args[3L], "simulated_counts.csv"))
writeLines(capture.output(sessionInfo()), file.path(args[3L], "session_info.txt"))
