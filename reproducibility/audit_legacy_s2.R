#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else "."
repro_dir <- normalizePath(dirname(script_file), mustWork = TRUE)
source(file.path(repro_dir, "historical", "simulator_YZ_manuscript.R"))

legacy_batch_count <- function(true_count, p0, celltype_proportion, seq_depth) {
  count <- lapply(true_count, function(ct) {
    apply(ct, c(1, 2), function(n) rbinom(1, n, 1 - p0))
  })
  count_failed <- lapply(count, function(ct) {
    vapply(seq_len(nrow(ct)), function(n) {
      sum(ct[n, ]) > seq_depth * celltype_proportion[n]
    }, logical(1))
  })
  count_obs <- count
  for (ct in names(true_count)) {
    failed <- which(count_failed[[ct]])
    if (length(failed)) {
      count_obs[[ct]][failed, ] <- round(
        proportions(count[[ct]][failed, ], 1) *
          seq_depth * celltype_proportion[failed, ct]
      )
    }
  }
  count_obs
}

nGenes <- 150L
gene_pattern <- data.frame(
  ct1 = rep("no_pattern", nGenes),
  ct2 = rep("no_pattern", nGenes),
  ct3 = rep("no_pattern", nGenes)
)
rownames(gene_pattern) <- paste0("gene", seq_len(nGenes))
layer_param <- data.frame(
  ct1 = c(1 / 2, 0, 0),
  ct2 = c(1 / 4, 3 / 4, 0),
  ct3 = c(1 / 3, 0, 1)
)

coord <- generate_coordinates(100, 100, "hex", spot_distance = 5, seed = 123)
mean_param <- generate_gene_pattern(
  coord = coord, gene_pattern = gene_pattern, seed = 123
)
true_count <- generate_true_count(coord, mean_param, gene_relationship = NULL)
layer_coord <- generate_layer(
  coord,
  xmin = c(0, 33, 66), xmax = c(33, 66, 100),
  ymin = c(0, 0, 0), ymax = c(100, 100, 100)
)
celltype_proportion <- generate_cell_prop(layer_coord, layer_param)
observed_count <- legacy_batch_count(
  true_count, p0 = 0.05, celltype_proportion = celltype_proportion,
  seq_depth = 1e6
)
spot_count <- generate_spot_count_legacy(
  observed_count, nGenes = nGenes, cell_number = 10,
  celltype_proportion = celltype_proportion
)
rownames(spot_count) <- paste0("spot_", seq_len(nrow(spot_count)))
colnames(spot_count) <- paste0("gene_", seq_len(ncol(spot_count)))

out_dir <- file.path(repro_dir, "results", "legacy_provenance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "S2_gene_expression_seed123.csv")
fwrite(spot_count, out_file)

sha_output <- system2("sha256sum", out_file, stdout = TRUE)
observed_sha <- strsplit(sha_output, "[[:space:]]+")[[1L]][1L]
expected_sha <- "7cfd5760c886f20f124d8d48f9524c13fc27241bea270ec96253e2c665e57c4f"
summary <- data.table(
  scenario = "S2",
  seed = 123L,
  thinning_probability = 0.05,
  aggregation = "abundance-weighted spot-level mixture",
  expected_sha256 = expected_sha,
  observed_sha256 = observed_sha,
  exact_match = identical(observed_sha, expected_sha)
)
fwrite(summary, file.path(out_dir, "legacy_s2_hash_audit.csv"))
print(summary)
if (!summary$exact_match) quit(status = 1L)
