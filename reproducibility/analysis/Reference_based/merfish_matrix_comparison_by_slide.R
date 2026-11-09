#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Matrix comparison plots for a given MERFISH slide
## Compares: True vs ST-COS, SRTsim, scDesign3, scCube, scCube_train
## Usage: Rscript merfish_matrix_comparison_by_slide.R 04
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(corrplot)
  library(ggplot2)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide slide number, e.g. 04")
}
slide_id <- sprintf("%02d", as.integer(args[1]))

count_path <- file.path("./data_external/MERFISH/count",
                        paste0("X0_", slide_id, "_count.csv"))
mapping_path <- file.path("./data_external/MERFISH/mapping",
                          paste0("X0_", slide_id, "_mapping.csv"))
output_prefix <- file.path("./Reference_based",
                           paste0("MERFISH_X0_", slide_id))

stcos_path <- paste0(output_prefix, "_simulated_counts.csv")
srt_path <- paste0(output_prefix, "_SRTsim_counts.csv")
scd3_path <- paste0(output_prefix, "_scDesign3_counts.csv")

scCube_path <- file.path("./Reference_based",
                         paste0("scCube_MERFISH_", slide_id, ".csv"))
scCube_meta_path <- file.path("./Reference_based",
                              paste0("scCube_MERFISH_", slide_id, "_meta.csv"))
scCube_train_dir <- "./Reference_based/scCube_train_MERFISH_X0_04"
scCube_train_path <- file.path(scCube_train_dir,
                               paste0("MERFISH_X0_", slide_id, "_sim_data.csv.gz"))
scCube_train_meta_path <- file.path(scCube_train_dir,
                                    paste0("MERFISH_X0_", slide_id, "_sim_meta.csv.gz"))

message(sprintf("[INFO] Slide %s: loading true counts", slide_id))
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

map_scube_cols <- function(sc_expr, coords, spot_ids, meta_dt) {
  meta_dt$key_xy <- paste0(round(meta_dt$x), "x", round(meta_dt$y))
  col_ids <- colnames(sc_expr)
  spot_key_xy <- paste0(round(coords[, 1]), "x", round(coords[, 2]))
  spot_key_yx <- paste0(round(coords[, 2]), "x", round(coords[, 1]))
  spot_map_xy <- setNames(spot_ids, spot_key_xy)
  spot_map_yx <- setNames(spot_ids, spot_key_yx)
  cell_keys <- meta_dt$key_xy[match(col_ids, meta_dt$Cell)]
  mapped_xy <- spot_map_xy[cell_keys]
  mapped_yx <- spot_map_yx[cell_keys]
  n_xy <- sum(!is.na(mapped_xy))
  n_yx <- sum(!is.na(mapped_yx))
  mapped <- if (n_yx > n_xy) mapped_yx else mapped_xy
  colnames(sc_expr) <- ifelse(is.na(mapped), colnames(sc_expr), mapped)
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

subset_genes_matrix <- function(mat, genes) {
  if (length(intersect(rownames(mat), genes)) > 0) {
    return(mat[intersect(rownames(mat), genes), , drop = FALSE])
  }
  if (length(intersect(colnames(mat), genes)) > 0) {
    return(mat[, intersect(colnames(mat), genes), drop = FALSE])
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

message("[INFO] Loading method outputs ...")
if (!file.exists(stcos_path)) stop(sprintf("Missing ST-COS output: %s", stcos_path))
if (!file.exists(srt_path)) stop(sprintf("Missing SRTsim output: %s", srt_path))
if (!file.exists(scd3_path)) stop(sprintf("Missing scDesign3 output: %s", scd3_path))

sim_spot_counts <- as.matrix(read.csv(stcos_path, row.names = 1, check.names = FALSE))
srt_sim_counts <- as.matrix(read.csv(srt_path, row.names = 1, check.names = FALSE))
scd3_counts <- as.matrix(read.csv(scd3_path, row.names = 1, check.names = FALSE))

scCube_expr <- NULL
if (file.exists(scCube_path)) {
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
} else {
  message(sprintf("[WARN] scCube results not found at %s; skipping scCube plots.", scCube_path))
}

scCube_train_expr <- NULL
if (file.exists(scCube_train_path)) {
  scCube_train_dt <- fread(cmd = paste("zcat", shQuote(scCube_train_path)))
  gene_names_scCube_train <- scCube_train_dt[[1]]
  rownames(scCube_train_dt) <- gene_names_scCube_train
  scCube_train_expr <- as.matrix(scCube_train_dt[, -1, with = FALSE])
  if (file.exists(scCube_train_meta_path)) {
    scCube_train_meta <- fread(cmd = paste("zcat", shQuote(scCube_train_meta_path)))
    if (length(intersect(colnames(scCube_train_expr), scCube_train_meta$Cell)) > 0) {
      scCube_train_expr <- map_scube_cols(scCube_train_expr, coords_mat, spot_ids, scCube_train_meta)
    } else if (length(intersect(rownames(scCube_train_expr), scCube_train_meta$Cell)) > 0) {
      scCube_train_expr <- map_scube_rows(scCube_train_expr, coords_mat, spot_ids, scCube_train_meta)
    }
  } else {
    message(sprintf("[WARN] scCube self-trained meta not found at %s; skipping coordinate mapping.", scCube_train_meta_path))
  }
  scCube_train_expr <- orient_by_names(scCube_train_expr, gene_names, spot_ids)
} else {
  message(sprintf("[WARN] scCube self-trained results not found at %s; skipping scCube_train plots.", scCube_train_path))
}

common_genes <- Reduce(intersect, list(
  rownames(expr_mat),
  colnames(sim_spot_counts),
  rownames(srt_sim_counts),
  colnames(scd3_counts),
  if (!is.null(scCube_expr)) rownames(scCube_expr) else rownames(expr_mat),
  if (!is.null(scCube_train_expr)) gene_names_scCube_train else rownames(expr_mat)
))
rownames(scCube_train_expr) <- gene_names_scCube_train
expr_mat <- expr_mat[common_genes, , drop = FALSE]
sim_spot_counts <- subset_genes_matrix(sim_spot_counts, common_genes)
srt_sim_counts <- subset_genes_matrix(srt_sim_counts, common_genes)
scd3_counts <- subset_genes_matrix(scd3_counts, common_genes)
if (!is.null(scCube_expr)) {
  scCube_expr <- subset_genes_matrix(scCube_expr, common_genes)
}
if (!is.null(scCube_train_expr)) {
  scCube_train_expr <- subset_genes_matrix(scCube_train_expr, common_genes)
}

cor_pal <- colorRampPalette(c("#3B4CC0", "#BFD3E6", "#FFFFFF", "#F4B6C2", "#B40426"))(200)

plot_corr_dual <- function(base_mat, method_mat, title, ord = NULL) {
  if (!is.null(ord)) {
    base_mat <- base_mat[ord, ord, drop = FALSE]
    method_mat <- method_mat[ord, ord, drop = FALSE]
  }
  corrplot::corrplot(
    base_mat,
    method = "color",
    col = cor_pal,
    type = "lower",
    order = "original",
    addrect = 4,
    tl.pos = "n",
    cl.pos = "n",
    diag = FALSE,
    mar = c(0, 0, 2, 0),
    title = title
  )
  corrplot::corrplot(
    method_mat,
    method = "color",
    col = cor_pal,
    type = "upper",
    order = "original",
    add = TRUE,
    tl.pos = "n",
    cl.cex = 0.7,
    diag = FALSE
  )
}

base_corr <- cor(t(expr_mat))
base_ord <- corrplot::corrMatOrder(base_corr, order = "hclust")

out_pdf <- paste0(output_prefix, "_matrix_comparison.pdf")
pdf(out_pdf, width = 8, height = 8)
par(mfrow = c(2, 2))
plot_corr_dual(base_corr, cor(as.matrix(sim_spot_counts)), "ST-COS vs True", ord = base_ord)
plot_corr_dual(base_corr, cor(as.matrix(t(srt_sim_counts))), "SRTsim vs True", ord = base_ord)
plot_corr_dual(base_corr, cor(as.matrix(scd3_counts)), "scDesign3 vs True", ord = base_ord)
if (!is.null(scCube_expr)) {
  plot_corr_dual(base_corr, cor(as.matrix(scCube_expr)), "scCube vs True", ord = base_ord)
}
if (!is.null(scCube_train_expr)) {
  plot_corr_dual(base_corr, cor(as.matrix(scCube_train_expr)), "scCube_train vs True", ord = base_ord)
}
dev.off()

message(sprintf("[INFO] Matrix comparison saved to %s", out_pdf))

message("[INFO] Building metrics comparison plots ...")
expr_true <- as_gene_by_spot(expr_mat, common_genes, spot_ids)
expr_stcos <- as_gene_by_spot(sim_spot_counts, common_genes, spot_ids)
expr_srt <- as_gene_by_spot(srt_sim_counts, common_genes, spot_ids)
expr_scd3 <- as_gene_by_spot(scd3_counts, common_genes, spot_ids)
expr_scube <- NULL
if (!is.null(scCube_expr)) {
  expr_scube <- as_gene_by_spot(scCube_expr, common_genes, spot_ids)
}
expr_scube_train <- NULL
if (!is.null(scCube_train_expr)) {
  expr_scube_train <- as_gene_by_spot(scCube_train_expr, common_genes, spot_ids)
}

make_long <- function(expr, label, spot_ids) {
  dt <- as.data.table(t(expr))
  dt[, spot := spot_ids]
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

plot_dt <- rbind(
  make_long(expr_true, "True", spot_ids),
  make_long(expr_stcos, "ST-COS", spot_ids),
  make_long(expr_srt, "SRTsim", spot_ids),
  make_long(expr_scd3, "scDesign3", spot_ids)
)
if (!is.null(expr_scube)) {
  plot_dt <- rbind(plot_dt, make_long(expr_scube, "scCube", spot_ids))
}
if (!is.null(expr_scube_train)) {
  plot_dt <- rbind(plot_dt, make_long(expr_scube_train, "scCube_train", spot_ids))
}

plot_dt$source <- factor(plot_dt$source,
                         levels = c("True", "ST-COS", "SRTsim", "scDesign3", "scCube", "scCube_train"))

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
  labs(title = "Gene-level Metrics Comparison")

location_metrics_plot <- location_metrics %>%
  pivot_longer(cols = c(LocZeroProp, LocLibSize),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = source, y = value, fill = source)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Location-level Metrics Comparison")

metrics_out <- paste0(output_prefix, "_metrics_comparison.pdf")
pdf(metrics_out, width = 12, height = 8)
print(gene_metrics_plot)
print(location_metrics_plot)
dev.off()

message(sprintf("[INFO] Metrics comparison saved to %s", metrics_out))
