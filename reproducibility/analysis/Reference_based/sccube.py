#!/usr/bin/env python3
"""scCube reference-based generation from a pretrained model.

Usage:
    python sccube.py counts.csv coords.csv output_dir model.pth \
        [--training-genes genes.txt] [--seed 20260922]

Input follows the shared reference-based format described in
reproducibility/README.md, except for the scCube-specific DLPFC input in
section 4 of reproducibility/EXTERNAL_INPUTS.txt (all genes, labelled spots
only). The optional ``label`` column of coords.csv is used as the cell-type
key. scCube generates library-size-normalized log1p expression, so the output
is output_dir/simulated_lognorm.csv (genes x generated observations, first
column gene IDs) rather than counts. ``--training-genes`` gives the model's
ordered training-gene list (one gene per line, no header). Training genes
absent from counts.csv are filled with zeros.
"""

import argparse
import random
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch
from scCube import scCube

parser = argparse.ArgumentParser()
parser.add_argument("counts")
parser.add_argument("coords")
parser.add_argument("output_dir")
parser.add_argument("model")
parser.add_argument("--training-genes", default=None)
parser.add_argument("--seed", type=int, default=20260922)
args = parser.parse_args()

random.seed(args.seed)
np.random.seed(args.seed)
torch.manual_seed(args.seed)

counts = pd.read_csv(args.counts, index_col=0)
counts.index = counts.index.astype(str)
counts.columns = counts.columns.astype(str)
coords = pd.read_csv(args.coords)
coords = coords.rename(columns={coords.columns[0]: "id"})
coords["id"] = coords["id"].astype(str)
coords = coords.set_index("id").reindex(counts.columns)
if coords[["row", "col"]].isna().any().any():
    raise ValueError("Every observation in counts.csv needs a row in coords.csv.")
labels = coords["label"].astype(str) if "label" in coords.columns else "CT"

if args.training_genes is not None:
    genes = pd.read_csv(args.training_genes, header=None).iloc[:, 0].astype(str)
    counts = counts.reindex(genes, fill_value=0)

adata = sc.AnnData(
    X=sp.csr_matrix(counts.T.to_numpy(dtype=np.float32)),
    obs=pd.DataFrame(index=counts.columns),
    var=pd.DataFrame(index=counts.index),
)
adata.obs["Cell"] = adata.obs_names
adata.obs["Cell_type"] = labels
adata.obs["x"] = coords["col"].to_numpy()
adata.obs["y"] = coords["row"].to_numpy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

model = scCube()
meta, data = model.load_vae_and_generate_cell(
    sc_adata=adata, celltype_key="Cell_type", cell_key="Cell", target_num=None,
    hidden_size=128, load_path=args.model, used_device="cpu",
)
data, meta = model.generate_pattern_reference(
    sc_adata=adata, generate_sc_data=data, generate_sc_meta=meta,
    celltype_key="Cell_type", spatial_key=["x", "y"], cost_metric="sqeuclidean",
)

expression = pd.DataFrame(data)
expression.index = expression.index.astype(str)
if not set(counts.index).issubset(expression.index):
    expression = expression.T
    expression.index = expression.index.astype(str)
expression = expression.loc[counts.index]
meta = pd.DataFrame(meta)
if "Cell" in meta.columns and len(meta) == expression.shape[1]:
    expression.columns = meta["Cell"].astype(str)

output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)
expression.to_csv(output_dir / "simulated_lognorm.csv", index_label="gene")
