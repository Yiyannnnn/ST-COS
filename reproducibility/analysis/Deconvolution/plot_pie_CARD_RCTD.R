library(data.table)
library(plotrix)

args <- commandArgs(trailingOnly = TRUE)
scenarios <- if (length(args) > 0) args else c("S7", "S8")

get_reps <- function(S) {
  if (S %in% c("S7", "S8")) return(1:3)
  1:5
}

plot_pies <- function(prop_mat, coords, out_pdf, colors = NULL, radius = 1) {
  common <- intersect(rownames(prop_mat), rownames(coords))
  if (length(common) == 0) {
    message(sprintf("[WARN] No overlapping spots for %s; skipping.", out_pdf))
    return(invisible(NULL))
  }
  prop_mat <- prop_mat[common, , drop = FALSE]
  coords <- coords[common, , drop = FALSE]
  if (is.null(colors)) {
    colors <- grDevices::rainbow(ncol(prop_mat))
  }
  pdf(out_pdf)
  plot(
    coords$x, coords$y,
    type = "n",
    xlab = "", ylab = "",
    axes = FALSE,
    asp = 1
  )
  for (idx in seq_len(nrow(coords))) {
    vals <- as.numeric(prop_mat[idx, ])
    if (all(is.na(vals)) || sum(vals, na.rm = TRUE) <= 0) next
    plotrix::floating.pie(
      xpos = coords$x[idx],
      ypos = coords$y[idx],
      x = vals,
      radius = radius,
      col = colors
    )
  }
  dev.off()
}

for (S in scenarios) {
  for (i in get_reps(S)) {
    coord_file <- sprintf("./Data/%s_gene/simulated/%s_coord_%d.csv", S, S, i)
    coords <- fread(coord_file)
    coords <- as.data.frame(coords)
    rownames(coords) <- paste0("spot_", seq_len(nrow(coords)))
    coords <- coords[, c("x_coord", "y_coord")]
    colnames(coords) <- c("x", "y")

    card_file <- sprintf("./Data/%s_gene/CARD_%d.csv", S, i)
    if (file.exists(card_file)) {
      card_dt <- as.data.frame(fread(card_file))
      if ("spot" %in% names(card_dt)) {
        rownames(card_dt) <- card_dt$spot
        card_dt$spot <- NULL
      }
      card_mat <- as.matrix(card_dt)
      out_pdf <- sprintf("./Deconvolution/%s_CARD_pie_%d.pdf", S, i)
      plot_pies(card_mat, coords, out_pdf)
    }

    rctd_file <- sprintf("./Data/%s_gene/RCTD_%d.csv", S, i)
    if (file.exists(rctd_file)) {
      rctd_dt <- as.data.frame(fread(rctd_file))
      if ("barcodes" %in% names(rctd_dt)) {
        rownames(rctd_dt) <- rctd_dt$barcodes
        rctd_dt$barcodes <- NULL
      }
      rctd_mat <- as.matrix(rctd_dt)
      out_pdf <- sprintf("./Deconvolution/%s_RCTD_pie_%d.pdf", S, i)
      plot_pies(rctd_mat, coords, out_pdf)
    }
  }
}
