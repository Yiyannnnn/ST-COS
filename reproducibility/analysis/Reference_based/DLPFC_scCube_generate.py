#!/usr/bin/env python
"""
Generate scCube simulated expression for one DLPFC slide using a pretrained VAE.

Usage:
    python DLPFC_scCube_generate.py <slide_id>
Example:
    python DLPFC_scCube_generate.py 151673

Outputs (in out_dir):
    {slide}_scCube_lognorm.csv.gz   -- genes x spots, log-normalized (as generated)
    {slide}_scCube_meta.csv         -- spots metadata: Cell, Cell_type, x, y
    {slide}_scCube_timing.csv       -- elapsed seconds
"""

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

warnings.filterwarnings("ignore")

# ---- args ----------------------------------------------------------------
if len(sys.argv) < 2:
    sys.exit("Usage: python DLPFC_scCube_generate.py <slide_id>")
slide_id = sys.argv[1]
print(f"[{slide_id}] scCube generation starting ...")

# ---- paths ---------------------------------------------------------------
DATA_DIR    = Path("./data_external/DLPFC")
LAYER_DIR   = DATA_DIR / "layers"
MODEL_PATH  = f"./data_external/DLPFC/models/slice{slide_id}_epoch100000.pth"
GENE_LIST   = "./Reference_based/sccube_18094_genes.csv"
OUT_DIR     = Path("./Reference_based/DLPFC")

# ---- load gene whitelist (18094 training genes) --------------------------
train_genes = pd.read_csv(GENE_LIST, header=None)[0].tolist()
print(f"[{slide_id}] Training gene set: {len(train_genes)} genes")

# ---- load 10X h5 ---------------------------------------------------------
h5_path = DATA_DIR / slide_id / f"{slide_id}_filtered_feature_bc_matrix.h5"
print(f"[{slide_id}] Loading {h5_path} ...")
adata = sc.read_10x_h5(str(h5_path))
adata.var_names_make_unique()

# ---- spatial coordinates (array_row / array_col) -------------------------
spatial_dir = DATA_DIR / slide_id / "spatial"
for fname in ["tissue_positions.txt", "tissue_positions_list.txt"]:
    tp = spatial_dir / fname
    if tp.exists():
        try:
            df = pd.read_csv(tp)
            if "barcode" not in df.columns:
                raise ValueError
        except Exception:
            df = pd.read_csv(tp, header=None)
            df.columns = ["barcode","in_tissue","array_row","array_col",
                          "pxl_row_in_fullres","pxl_col_in_fullres"]
        break
coords = df[df["in_tissue"] == 1].set_index("barcode")
common_bc = adata.obs_names.intersection(coords.index)
adata = adata[common_bc].copy()
adata.obs["x"] = coords.loc[adata.obs_names, "array_col"].values
adata.obs["y"] = coords.loc[adata.obs_names, "array_row"].values

# ---- layer labels --------------------------------------------------------
layer_path = LAYER_DIR / f"{slide_id}_layers.txt"
if layer_path.exists():
    layers = pd.read_csv(layer_path, header=None, sep=None, engine="python")
    layer_map = layers.set_index(layers.columns[0]).iloc[:, 0]
    layer_map.index = layer_map.index.astype(str)
    keep = adata.obs_names.intersection(layer_map.index)
    adata = adata[keep].copy()
    adata.obs["Cell_type"] = layer_map.loc[adata.obs_names].astype(str).values
    print(f"[{slide_id}] {len(keep)} spots with layer labels")
else:
    adata.obs["Cell_type"] = "CT"
    print(f"[{slide_id}] No layer file found — using single label CT")

adata.obs["Cell"] = adata.obs_names.astype(str)

# ---- align to training genes (zero-pad missing) -------------------------
available = set(adata.var_names)
missing   = [g for g in train_genes if g not in available]
present   = [g for g in train_genes if g in available]

if missing:
    print(f"[{slide_id}] Zero-padding {len(missing)} missing genes (of {len(train_genes)})")
    import scipy.sparse as ssp
    n_spots  = adata.n_obs
    zero_X   = ssp.csr_matrix((n_spots, len(missing)), dtype=np.float32)
    zero_adata = sc.AnnData(
        X   = zero_X,
        obs = adata.obs.copy(),
        var = pd.DataFrame(index=missing),
    )
    adata_present = adata[:, present].copy()
    adata = sc.concat([adata_present, zero_adata], axis=1, merge="same")
    adata = adata[:, train_genes].copy()
else:
    adata = adata[:, train_genes].copy()

print(f"[{slide_id}] Gene matrix: {adata.n_vars} genes x {adata.n_obs} spots")

# ---- normalize: total 1e4 + log1p  (same as scCube pre_process default) --
if sp.issparse(adata.X):
    adata.X = adata.X.toarray()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ---- build sc_data (genes x cells) and sc_meta ---------------------------
sc_data = pd.DataFrame(adata.X.T, index=adata.var_names, columns=adata.obs_names)
sc_meta = adata.obs[["Cell", "Cell_type", "x", "y"]].copy().reset_index(drop=True)
sc_meta["Cell"] = sc_meta["Cell"].astype(str)
sc_data.columns = sc_data.columns.astype(str)

# ---- scCube model --------------------------------------------------------
from scCube import scCube
model = scCube()

t0 = time.time()
print(f"[{slide_id}] load_vae_and_generate_cell ...")
generate_sc_meta, generate_sc_data = model.load_vae_and_generate_cell(
    sc_adata    = adata,
    celltype_key = "Cell_type",
    cell_key     = "Cell",
    target_num   = None,
    hidden_size  = 128,
    load_path    = MODEL_PATH,
    used_device  = "cpu",
)
t1 = time.time()
print(f"[{slide_id}] generate_sc_meta shape: {generate_sc_meta.shape}, "
      f"generate_sc_data shape: {generate_sc_data.shape}")

print(f"[{slide_id}] generate_pattern_reference ...")
generate_sc_data, generate_sc_meta = model.generate_pattern_reference(
    sc_adata          = adata,
    generate_sc_data  = generate_sc_data,
    generate_sc_meta  = generate_sc_meta,
    celltype_key      = "Cell_type",
    spatial_key       = ["x", "y"],
    cost_metric       = "sqeuclidean",
)
t2 = time.time()

total_sec = t2 - t0
print(f"[{slide_id}] Done in {total_sec:.1f}s")

# ---- save ----------------------------------------------------------------
lognorm_out = OUT_DIR / f"{slide_id}_scCube_lognorm.csv.gz"
meta_out    = OUT_DIR / f"{slide_id}_scCube_meta.csv"
timing_out  = OUT_DIR / f"{slide_id}_scCube_timing.csv"

# generate_sc_data: genes x cells  — save as-is (log-normalized)
generate_sc_data.to_csv(lognorm_out, compression="gzip")
generate_sc_meta.to_csv(meta_out, index=False)

pd.DataFrame({"slide": [slide_id], "n_genes": [adata.n_vars],
              "n_spots": [adata.n_obs], "elapsed_sec": [total_sec]}
             ).to_csv(timing_out, index=False)

print(f"[{slide_id}] Saved: {lognorm_out}")
print(f"[{slide_id}] Saved: {meta_out}")
print(f"[{slide_id}] Saved: {timing_out}")
