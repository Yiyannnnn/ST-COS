#!/usr/bin/env bash
set -euo pipefail

# Wrapper to run the image-driven simulator on Longleaf

IMG_PATH="./data_external/HER2/9769_C1_HE_small.jpg"
OUT_PREFIX="./Image_feature/simulator_image"
K_DOMAINS=3
SPOT_DISTANCE=6
NGENES=150
N_CELLTYPES=3
SVG_PER_CT=5
SEQ_DEPTH=1000000
XY_WEIGHT=0.15
SEED=123
SAMPLE="A1"
SEG_METHOD="slic"            # options: kmeans, gmm, slic
SLIC_SUPERPIXELS=1500
SLIC_COMPACTNESS=25
export SEG_METHOD SLIC_SUPERPIXELS SLIC_COMPACTNESS

mkdir -p ./temp
export TMPDIR=./temp
export TMP=$TMPDIR
export TEMP=$TMPDIR
export R_SESSION_TMPDIR=$TMPDIR

module --ignore_cache load r/4.5.0 || module --ignore_cache load r/4.4.0

Rscript - <<'RS'
source("STSimulator/Image_feature/simulator_image.R")
simulator_image(
  image_path      = "./data_external/HER2/9769_C1_HE_small.jpg",
  out_prefix      = "./Image_feature/simulator_image",
  K_domains       = 3,
  spot_distance   = 6,
  nGenes          = 150,
  n_celltypes     = 3,
  svg_per_ct      = 5,
  seq_depth       = 1e6,
  xy_weight       = 0.08,
  seed            = 123,
  mask_background = TRUE,
  white_thr       = 0.94,
  use_mapping     = TRUE,
  sample          = "A1",
  mapping_x_col   = "pixel_x",
  mapping_y_col   = "pixel_y",
  spots_only      = TRUE,
  spot_radius     = 8,
  segment_only    = TRUE,
  segmentation_method = Sys.getenv("SEG_METHOD", unset = "kmeans"),
  slic_superpixel_count = as.numeric(Sys.getenv("SLIC_SUPERPIXELS", unset = "1200")),
  slic_compactness = as.numeric(Sys.getenv("SLIC_COMPACTNESS", unset = "20"))
)
RS

echo "Done. Outputs at ${OUT_PREFIX}_*.csv and ${OUT_PREFIX}_segmentation.png"
