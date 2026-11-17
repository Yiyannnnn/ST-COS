#!/usr/bin/env Rscript
# SRTsim reference-based simulation (tissue-based scheme).
#
# Usage: Rscript srtsim.R counts.csv coords.csv output_dir [seed]
#
# Input and output follow the shared reference-based format described in
# reproducibility/README.md. The optional `label` column of coords.csv is
# passed to SRTsim as the domain label.
# Writes output_dir/simulated_counts.csv in the same layout as counts.csv.

suppressPackageStartupMessages({
  library(data.table)
  library(SRTsim)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript srtsim.R counts.csv coords.csv output_dir [seed]")
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
locations <- data.frame(x = coords$col, y = coords$row, label = label,
                        row.names = colnames(counts))

set.seed(seed)
fit <- createSRT(count_in = counts, loc_in = locations)
fit <- srtsim_fit(fit, sim_schem = "tissue")
fit <- srtsim_count(fit)
simulated <- as.matrix(simCounts(fit))
if (!all(rownames(counts) %in% rownames(simulated))) simulated <- t(simulated)
simulated <- simulated[rownames(counts), colnames(counts), drop = FALSE]

dir.create(args[3L], recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(gene = rownames(simulated), as.data.table(simulated)),
       file.path(args[3L], "simulated_counts.csv"))
writeLines(capture.output(sessionInfo()), file.path(args[3L], "session_info.txt"))
