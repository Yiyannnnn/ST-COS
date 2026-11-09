#!/usr/bin/env Rscript

repo_dir <- normalizePath(".", mustWork = TRUE)
source_s7_dir <- file.path(repo_dir, "Data", "S7_gene", "simulated")
source_srtsim_dir <- file.path(repo_dir, "Independent_simulator_anchor", "SRTsim_S7_celltype_aggregate")
target_dir <- file.path(repo_dir, "Data", "S7SRT_gene")
target_sim_dir <- file.path(target_dir, "simulated")

dir.create(target_sim_dir, recursive = TRUE, showWarnings = FALSE)

copy_checked <- function(from, to) {
  if (!file.exists(from)) {
    stop("Missing input: ", from)
  }
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok) {
    stop("Could not copy ", from, " to ", to)
  }
  data.frame(source = from, target = to, stringsAsFactors = FALSE)
}

manifest <- list()

for (rep_id in 1:5) {
  manifest[[length(manifest) + 1]] <- copy_checked(
    file.path(source_srtsim_dir, sprintf("S7_CTaggregate_SRTsim_counts_rep%d.csv", rep_id)),
    file.path(target_sim_dir, sprintf("S7SRT_gene_expression_%d.csv", rep_id))
  )
  for (kind in c("coord", "layer", "region", "scref")) {
    manifest[[length(manifest) + 1]] <- copy_checked(
      file.path(source_s7_dir, sprintf("S7_%s_%d.csv", kind, rep_id)),
      file.path(target_sim_dir, sprintf("S7SRT_%s_%d.csv", kind, rep_id))
    )
  }
}

for (kind in c("marker_gene_truth", "coexpression_truth")) {
  manifest[[length(manifest) + 1]] <- copy_checked(
    file.path(source_s7_dir, sprintf("S7_%s.csv", kind)),
    file.path(target_sim_dir, sprintf("S7SRT_%s.csv", kind))
  )
}

manifest <- do.call(rbind, manifest)
manifest$created_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
manifest$note <- "S7SRT uses SRTsim cell-type aggregate spatial expression with S7 coordinates, labels, proportions, marker truth, co-expression truth, and scRNA references."

write.csv(manifest, file.path(target_dir, "S7SRT_manifest.csv"), row.names = FALSE)

message("Prepared S7SRT scenario at: ", target_sim_dir)
