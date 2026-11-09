#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else "."
repro_dir <- normalizePath(dirname(script_file), mustWork = TRUE)

core_packages <- c("data.table", "ggplot2", "Matrix", "mgcv", "spdep")
reference_packages <- c(
  "Seurat", "SRTsim", "scDesign3", "SingleCellExperiment", "cowplot"
)
downstream_packages <- c(
  "CARD", "spacexr", "STdeconvolve", "BayesSpace", "BASS", "SPARK",
  "aricode", "dplyr", "tidyr"
)

groups <- list(
  core = core_packages,
  reference_based = reference_packages,
  downstream = downstream_packages
)
rows <- do.call(rbind, lapply(names(groups), function(group) {
  packages <- groups[[group]]
  data.frame(
    group = group,
    package = packages,
    installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    version = vapply(packages, function(package) {
      if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
      as.character(utils::packageVersion(package))
    }, character(1)),
    stringsAsFactors = FALSE
  )
}))

out_dir <- file.path(repro_dir, "results", "session")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(rows, file.path(out_dir, "package_availability.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(rows, row.names = FALSE)
cat("\nPython:", Sys.which("python3"), "\n")
cat("External data directory:", file.path(repro_dir, "data_external"), "\n")
