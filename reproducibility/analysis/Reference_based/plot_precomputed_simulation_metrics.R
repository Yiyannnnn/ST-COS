#!/usr/bin/env Rscript

## This script loads the precomputed simulation matrices (ST-COS, SRTsim, scDesign3,
## + scCube pretrained) and reuses the metric/plot workflow from real_parameter.R.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(Seurat)
})

gene_list <- fread("./data_external/DLPFC/151673_selected_gene.csv")
eval_genes <- gene_list$gene[1:3000]

read_matrix <- function(path, transpose = FALSE) {
  mat <- fread(path)
  rownames(mat) <- mat[[1]]
  mat <- as.matrix(mat[, -1, with = FALSE])
  if (transpose) {
    mat <- t(mat)
  }
  mat
}

subset_genes <- function(mat, genes) {
  rn <- rownames(mat)
  keep <- intersect(genes, rn)
  if (length(keep) == 0) {
    keep <- rn[seq_len(min(length(genes), length(rn)))]
  }
  mat[keep, , drop = FALSE]
}

align_spots <- function(mat, spot_ids) {
  keep <- intersect(spot_ids, colnames(mat))
  if (length(keep) == 0) return(mat[, integer(0), drop = FALSE])
  mat[, keep, drop = FALSE]
}

stcos_mat <- read_matrix(
  "./Reference_based/comp_eff_test/stcos_complete_3000_expr_gene_x_spot.csv",
  transpose = FALSE
)
srt_mat <- read_matrix(
  "./Reference_based/comp_eff_test/srtsim_counts_3000_expr_gene_x_spot.csv",
  transpose = FALSE
)
scd3_mat <- read_matrix(
  "./Reference_based/comp_eff_test/scdesign3_complete_3000.csv",
  transpose = TRUE
)
scCube_mat <- read_matrix("./data_external/DLPFC/scCube_151673_expr_cell_x_gene.csv", transpose = TRUE)
rownames(stcos_mat) <- eval_genes
rownames(srt_mat) <- eval_genes
rownames(scd3_mat) <- eval_genes
colnames(stcos_mat) <- colnames(srt_mat)
colnames(scd3_mat) <- colnames(srt_mat)
colnames(scCube_mat) <- colnames(stcos_mat)
spot_ids <- colnames(stcos_mat)

stcos_mat <- subset_genes(stcos_mat, eval_genes)
srt_mat <- subset_genes(srt_mat, eval_genes)
scd3_mat <- subset_genes(scd3_mat, eval_genes)
scCube_mat <- subset_genes(scCube_mat, eval_genes)

stcos_mat <- align_spots(stcos_mat, spot_ids)
srt_mat <- align_spots(srt_mat, spot_ids)
scd3_mat <- align_spots(scd3_mat, spot_ids)
scCube_mat <- align_spots(scCube_mat, spot_ids)

dir <- "./data_external/DLPFC/"
id <- "151673"
h5_mat_name <- paste0(id, "_filtered_feature_bc_matrix.h5")

your_data <- Load10X_Spatial(
  data.dir = paste0(dir, id, "/"),
  filename = h5_mat_name,
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE
)
expr_true <- as.matrix(your_data@assays$Spatial$counts)
truth_mat <- subset_genes(expr_true, eval_genes)
truth_mat <- align_spots(truth_mat, spot_ids)

combine_long <- function(mat, method) {
  dt <- as.data.table(mat)
  if (is.null(rownames(mat))) {
    stop("Matrix must have row names for genes")
  }
  dt[, gene := rownames(mat)]
  dt_long <- melt(dt, id.vars = "gene", variable.name = "variable", value.name = "count")
  dt_long[, source := method]
  setnames(dt_long, "variable", "cellnames")
  dt_long
}

combined <- rbindlist(list(
  combine_long(stcos_mat, "ST-COS"),
  combine_long(srt_mat, "SRTsim"),
  combine_long(scd3_mat, "scDesign3"),
  combine_long(scCube_mat, "scCube")
), use.names = TRUE)

if (nrow(truth_mat) > 0 && ncol(truth_mat) > 0) {
  truth_dt <- combine_long(truth_mat, "Truth")
  combined <- rbindlist(list(combined, truth_dt), use.names = TRUE)
}

combined$method <- factor(combined$source,
                          levels = c("Truth", "ST-COS", "SRTsim", "scDesign3", "scCube"))

gene_metrics <- combined %>%
  group_by(gene, method) %>%
  summarise(
    GeneMean = mean(count, na.rm = TRUE),
    GeneVar = var(count, na.rm = TRUE),
    GeneCV = sd(count, na.rm = TRUE) / mean(count, na.rm = TRUE),
    GeneZeroProp = mean(count == 0, na.rm = TRUE)
  ) %>%
  ungroup()

location_metrics <- combined %>%
  group_by(cellnames, method) %>%
  summarise(
    LocZeroProp = mean(count == 0, na.rm = TRUE),
    LocLibSize = sum(count, na.rm = TRUE)
  ) %>%
  ungroup()

gene_summary <- gene_metrics %>%
  group_by(method) %>%
  summarise(
    Mean_GeneMean = mean(GeneMean, na.rm = TRUE),
    Mean_GeneVar = mean(GeneVar, na.rm = TRUE),
    Mean_GeneCV = mean(GeneCV, na.rm = TRUE),
    Mean_GeneZeroProp = mean(GeneZeroProp, na.rm = TRUE)
  )

location_summary <- location_metrics %>%
  group_by(method) %>%
  summarise(
    Mean_LocZeroProp = mean(LocZeroProp, na.rm = TRUE),
    Mean_LocLibSize = mean(LocLibSize, na.rm = TRUE)
  )

cat("Gene-level metrics summary:\n")
print(gene_summary)
cat("\nLocation-level metrics summary:\n")
print(location_summary)

write.csv(gene_metrics, file = "precomputed_gene_metrics.csv", row.names = FALSE)
write.csv(location_metrics, file = "precomputed_location_metrics.csv", row.names = FALSE)

plot_metrics <- function(dt, y_label, title) {
  dt %>%
    ggplot(aes(x = method, y = value, fill = method)) +
    geom_violin(trim = FALSE, alpha = 0.9) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.6, outlier.shape = NA) +
    facet_wrap(~metric, scales = "free_y") +
    theme_minimal(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 20),
          legend.position = "none",
          strip.text = element_text(size = 16)) +
    labs(title = title, x = "Method", y = y_label)
}

gene_metrics_plot <- gene_metrics %>%
  pivot_longer(cols = c(GeneMean, GeneVar, GeneCV, GeneZeroProp),
               names_to = "metric", values_to = "value") %>%
  plot_metrics("Metric value", "Gene-level Metrics Comparison (precomputed)")

location_metrics_plot <- location_metrics %>%
  pivot_longer(cols = c(LocZeroProp, LocLibSize),
               names_to = "metric", values_to = "value") %>%
  plot_metrics("Metric value", "Location-level Metrics Comparison (precomputed)")

pdf("precomputed_metrics_comparison.pdf", width = 14, height = 10)
print(gene_metrics_plot)
print(location_metrics_plot)
dev.off()
