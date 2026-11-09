import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import math
import sys
import os
import sys
from importlib import import_module
spg = import_module('SpaGCN')

print("SpaGCN is imported from:", spg.__file__)

import random, torch
import warnings
warnings.filterwarnings("ignore")
# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run SpaGCN analysis with dynamic file paths.")
parser.add_argument('S', type=str, help="The dataset identifier (e.g., S1).")
parser.add_argument('suffix', type=str, help="The suffix for file naming (e.g., _1).")
args = parser.parse_args()

# Dynamically create file paths
data_f = f'./Data/{args.S}_gene/simulated/'
gene_f = data_f + f'{args.S}_gene_expression{args.suffix}.csv'
coor_f = data_f + f'{args.S}_coord{args.suffix}.csv'
true_f = data_f + f'{args.S}_layer{args.suffix}.csv'

# Load data
expr = pd.read_csv(gene_f)  # genes as columns
coords = pd.read_csv(coor_f)
layers = pd.read_csv(true_f)

# Create the AnnData object
adata = sc.AnnData(X=expr, dtype=np.float32)

# Add the spatial coordinates
adata.obsm["spatial"] = coords.values

# Add the layer annotations (clusters)
adata.obs["Ground Truth"] = layers.layer.values
adata.obs["x_pixel"] = coords["x_coord"].values
adata.obs["y_pixel"] = coords["y_coord"].values

print(adata)

# Convert the expression matrix to float
adata.X = adata.X.astype(float)

# Select captured samples
adata.var_names = [i.upper() for i in list(adata.var_names)]
adata.var["genename"] = adata.var.index.astype("str")

# Set coordinates
x_pixel = adata.obs["x_pixel"].astype(int).tolist()
y_pixel = adata.obs["y_pixel"].astype(int).tolist()

# Calculate adjacency matrix
adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, histology=False)

# Ensure unique variable names
adata.var_names_make_unique()

# Pre-filter genes
spg.prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)

# Normalize and log transform
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

# Set parameters
p = 0.5
l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

# Methods requiring a fixed cluster count receive the simulated truth count.
n_clusters = int(adata.obs["Ground Truth"].nunique(dropna=True))
if n_clusters < 2:
    raise ValueError(f"Expected at least two spatial domains, found {n_clusters}.")
print(f"Using {n_clusters} spatial domains from the simulation truth.")

# Set seeds
r_seed = t_seed = n_seed = 100

# Search for suitable resolution
res = spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20,
                     r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

# Initialize and train the SpaGCN model
clf = spg.SpaGCN()
clf.set_l(l)

# Set random seed
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)

# Train the model
clf.train(adata, adj, init_spa=True, init="louvain", res=res, tol=5e-3, lr=0.05)
y_pred, prob = clf.predict()

# Store predictions
adata.obs["pred"] = pd.Categorical(y_pred)

# Calculate 2D adjacency matrix for refinement
adj_2d = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, histology=False)

# Refine the predictions
refined_pred = spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d,
                          shape="square")
adata.obs["refined_pred"] = pd.Categorical(refined_pred)


# Save the results
results = adata.obs[['refined_pred', 'Ground Truth']]
results['x'] = coords['x_coord'].values
results['y'] = coords['y_coord'].values

output_file = f'./Data/{args.S}_gene/SpaGCN{args.suffix}.csv'
results.to_csv(output_file, index_label='Row_names')

#adjusted_rand_score(adata.obs['refined_pred'], adata.obs['Ground Truth'])
