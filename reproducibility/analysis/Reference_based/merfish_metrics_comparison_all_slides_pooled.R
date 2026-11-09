#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Metrics comparison pooled across slides
## Compares: True vs ST-COS, SRTsim, scDesign3, scCube
## Usage: Rscript merfish_metrics_comparison_all_slides_pooled.R 04 09 14 19 24
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  slides <- sprintf("%02d", as.integer(args))
} else {
  slides <- c("04", "09", "14", "19", "24")
}

base_count_dir <- "./data_external/MERFISH/count"
base_mapping_dir <- "./data_external/MERFISH/mapping"
output_dir <- "./Reference_based"

map_scube_rows <- function(sc_expr, coords, spot_ids, meta_dt) {
  meta_dt$key_xy <- paste0(round(meta_dt$x), "x", round(meta_dt$y))
  row_ids <- rownames(sc_expr)
  spot_key_xy <- paste0(round(coords[, 1]), "x", round(coords[, 2]))
  spot_key_yx <- paste0(round(coords[, 2]), "x", round(coords[, 1]))
  spot_map_xy <- setNames(spot_ids, spot_key_xy)
  spot_map_yx <- setNames(spot_ids, spot_key_yx)
  mapped_xy <- spot_map_xy[row_ids]
  mapped_yx <- spot_map_yx[row_ids]
  n_xy <- sum(!is.na(mapped_xy))
  n_yx <- sum(!is.na(mapped_yx))
  mapped <- if (n_yx > n_xy) mapped_yx else mapped_xy
  rownames(sc_expr) <- ifelse(is.na(mapped), rownames(sc_expr), mapped)
  sc_expr
}

orient_by_names <- function(mat, gene_ids, spot_ids) {
  gene_in_rows <- length(intersect(rownames(mat), gene_ids))
  gene_in_cols <- length(intersect(colnames(mat), gene_ids))
  spot_in_rows <- length(intersect(rownames(mat), spot_ids))
  spot_in_cols <- length(intersect(colnames(mat), spot_ids))
  if (gene_in_cols > gene_in_rows && spot_in_rows > spot_in_cols) {
    mat <- t(mat)
  }
  mat
}

as_gene_by_spot <- function(mat, gene_ids, spot_ids) {
  mat <- as.matrix(mat)
  if (length(intersect(gene_ids, colnames(mat))) > 0 &&
      length(intersect(spot_ids, rownames(mat))) > 0) {
    mat <- t(mat)
  }
  gene_keep <- intersect(gene_ids, rownames(mat))
  spot_keep <- intersect(spot_ids, colnames(mat))
  if (length(gene_keep) == 0 || length(spot_keep) == 0) {
    stop("No overlapping genes or spots found for metrics comparison.")
  }
  mat <- mat[gene_keep, spot_keep, drop = FALSE]
  mat
}

get_gene_ids <- function(mat, gene_names) {
  row_hits <- length(intersect(rownames(mat), gene_names))
  col_hits <- length(intersect(colnames(mat), gene_names))
  if (col_hits > row_hits) colnames(mat) else rownames(mat)
}

load_slide <- function(slide_id) {
  count_path <- file.path(base_count_dir, paste0("X0_", slide_id, "_count.csv"))
  mapping_path <- file.path(base_mapping_dir, paste0("X0_", slide_id, "_mapping.csv"))
  output_prefix <- file.path(output_dir, paste0("MERFISH_X0_", slide_id))

  stcos_path <- paste0(output_prefix, "_simulated_counts.csv")
  srt_path <- paste0(output_prefix, "_SRTsim_counts.csv")
  scd3_path <- paste0(output_prefix, "_scDesign3_counts.csv")

  scCube_path <- file.path(output_dir, paste0("scCube_MERFISH_", slide_id, ".csv"))
  scCube_meta_path <- file.path(output_dir, paste0("scCube_MERFISH_", slide_id, "_meta.csv"))

  message(sprintf("[INFO] Slide %s: loading counts", slide_id))
  counts_dt <- fread(count_path)
  spot_ids <- counts_dt[[1]]
  gene_names <- colnames(counts_dt)[-1]
  expr_mat_spot_gene <- as.matrix(counts_dt[, -1, with = FALSE])
  expr_mat <- t(expr_mat_spot_gene)
  rownames(expr_mat) <- gene_names
  colnames(expr_mat) <- spot_ids

  message(sprintf("[INFO] Slide %s: loading coordinates", slide_id))
  mapping_dt <- fread(mapping_path)
  coords_mat <- as.matrix(mapping_dt[match(spot_ids, mapping_dt$row_name), c("x", "y"), with = FALSE])
  if (any(is.na(coords_mat))) {
    stop("Missing coordinates for some spots after matching mapping table.")
  }
  rownames(coords_mat) <- spot_ids

  if (!file.exists(stcos_path)) stop(sprintf("Missing ST-COS output: %s", stcos_path))
  if (!file.exists(srt_path)) stop(sprintf("Missing SRTsim output: %s", srt_path))
  if (!file.exists(scd3_path)) stop(sprintf("Missing scDesign3 output: %s", scd3_path))
  if (!file.exists(scCube_path)) stop(sprintf("Missing scCube output: %s", scCube_path))

  sim_spot_counts <- as.matrix(read.csv(stcos_path, row.names = 1, check.names = FALSE))
  srt_sim_counts <- as.matrix(read.csv(srt_path, row.names = 1, check.names = FALSE))
  scd3_counts <- as.matrix(read.csv(scd3_path, row.names = 1, check.names = FALSE))

  scCube_expr <- read.csv(scCube_path, row.names = 1, check.names = FALSE)
  scCube_expr <- as.matrix(scCube_expr)
  if (file.exists(scCube_meta_path)) {
    scCube_meta <- read.csv(scCube_meta_path, check.names = FALSE)
    if (!("key_xy" %in% colnames(scCube_meta))) {
      scCube_meta$key_xy <- paste0(round(scCube_meta$x), "x", round(scCube_meta$y))
    }
    if (length(intersect(rownames(scCube_expr), scCube_meta$key_xy)) > 0) {
      scCube_expr <- map_scube_rows(scCube_expr, coords_mat, spot_ids, scCube_meta)
    } else if (length(intersect(rownames(scCube_expr), scCube_meta$Cell)) > 0) {
      spot_key <- paste0(round(coords_mat[, 1]), "x", round(coords_mat[, 2]))
      spot_key_map <- setNames(spot_ids, spot_key)
      cell_keys <- scCube_meta$key_xy[match(rownames(scCube_expr), scCube_meta$Cell)]
      mapped <- spot_key_map[cell_keys]
      rownames(scCube_expr) <- ifelse(is.na(mapped), rownames(scCube_expr), mapped)
    }
  } else {
    message(sprintf("[WARN] scCube meta not found at %s; skipping coordinate mapping.", scCube_meta_path))
  }
  scCube_expr <- orient_by_names(scCube_expr, gene_names, spot_ids)

  common_genes <- Reduce(intersect, list(
    gene_names,
    get_gene_ids(sim_spot_counts, gene_names),
    get_gene_ids(srt_sim_counts, gene_names),
    get_gene_ids(scd3_counts, gene_names),
    get_gene_ids(scCube_expr, gene_names)
  ))

  list(
    slide_id = slide_id,
    spot_ids = spot_ids,
    gene_names = gene_names,
    expr_mat = expr_mat,
    sim_spot_counts = sim_spot_counts,
    srt_sim_counts = srt_sim_counts,
    scd3_counts = scd3_counts,
    scCube_expr = scCube_expr,
    common_genes = common_genes
  )
}

slide_data <- list()
common_genes_list <- list()
for (slide_id in slides) {
  data <- load_slide(slide_id)
  slide_data[[slide_id]] <- data
  common_genes_list[[slide_id]] <- data$common_genes
}

global_genes <- Reduce(intersect, common_genes_list)
if (length(global_genes) == 0) {
  stop("No common genes found across slides and methods.")
}
message(sprintf("[INFO] Using %d common genes across slides", length(global_genes)))

make_long <- function(expr, label) {
  dt <- as.data.table(t(expr), keep.rownames = "spot")
  dt <- melt(
    dt,
    id.vars = "spot",
    variable.name = "gene",
    value.name = "count"
  )
  dt <- as.data.table(dt)
  dt[, source := label]
  dt
}

plot_dt_list <- list()
for (slide_id in slides) {
  data <- slide_data[[slide_id]]
  spot_ids <- data$spot_ids

  expr_true <- as_gene_by_spot(data$expr_mat, global_genes, spot_ids)
  expr_stcos <- as_gene_by_spot(data$sim_spot_counts, global_genes, spot_ids)
  expr_srt <- as_gene_by_spot(data$srt_sim_counts, global_genes, spot_ids)
  expr_scd3 <- as_gene_by_spot(data$scd3_counts, global_genes, spot_ids)
  expr_scube <- as_gene_by_spot(data$scCube_expr, global_genes, spot_ids)

  expr_true <- expr_true[global_genes, , drop = FALSE]
  expr_stcos <- expr_stcos[global_genes, , drop = FALSE]
  expr_srt <- expr_srt[global_genes, , drop = FALSE]
  expr_scd3 <- expr_scd3[global_genes, , drop = FALSE]
  expr_scube <- expr_scube[global_genes, , drop = FALSE]

  prefix_spots <- paste0("X0_", slide_id, "_", colnames(expr_true))
  colnames(expr_true) <- prefix_spots
  colnames(expr_stcos) <- prefix_spots
  colnames(expr_srt) <- prefix_spots
  colnames(expr_scd3) <- prefix_spots
  colnames(expr_scube) <- prefix_spots

  plot_dt_list[[length(plot_dt_list) + 1]] <- make_long(expr_true, "True")
  plot_dt_list[[length(plot_dt_list) + 1]] <- make_long(expr_stcos, "ST-COS")
  plot_dt_list[[length(plot_dt_list) + 1]] <- make_long(expr_srt, "SRTsim")
  plot_dt_list[[length(plot_dt_list) + 1]] <- make_long(expr_scd3, "scDesign3")
  plot_dt_list[[length(plot_dt_list) + 1]] <- make_long(expr_scube, "scCube")
}

plot_dt <- rbindlist(plot_dt_list)
plot_dt$source <- factor(plot_dt$source,
                         levels = c("True", "ST-COS", "SRTsim", "scDesign3", "scCube"))

message("[INFO] Computing pooled metrics ...")
gene_metrics <- plot_dt[, .(
  GeneMean = mean(count, na.rm = TRUE),
  GeneVar = var(count, na.rm = TRUE),
  GeneCV = sd(count, na.rm = TRUE) / mean(count, na.rm = TRUE),
  GeneZeroProp = mean(count == 0, na.rm = TRUE)
), by = .(gene, source)]

location_metrics <- plot_dt[, .(
  LocZeroProp = mean(count == 0, na.rm = TRUE),
  LocLibSize = sum(count, na.rm = TRUE)
), by = .(spot, source)]

gene_metrics_plot <- gene_metrics %>%
  pivot_longer(cols = c(GeneMean, GeneVar, GeneCV, GeneZeroProp),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = source, y = value, fill = source)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene-level Metrics Comparison (Pooled Slides)")

location_metrics_plot <- location_metrics %>%
  pivot_longer(cols = c(LocZeroProp, LocLibSize),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = source, y = value, fill = source)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Location-level Metrics Comparison (Pooled Slides)")

metrics_out <- file.path(output_dir, "MERFISH_all_slides_metrics_comparison_pooled.pdf")
pdf(metrics_out, width = 12, height = 8)
print(gene_metrics_plot)
print(location_metrics_plot)
dev.off()

message(sprintf("[INFO] Pooled metrics comparison saved to %s", metrics_out))
