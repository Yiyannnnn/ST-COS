#!/usr/bin/env conda run -n STCOS_eval python
import torch
import argparse
import scanpy as sc
import sys
sys.path.append('./Segmentation/GraphST')
from GraphST.utils import clustering
from GraphST import GraphST
import pandas as pd
import numpy as np 



# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run GraphST analysis with dynamic file paths.")
parser.add_argument('S', type=str, help="The dataset identifier (e.g., S1).")
parser.add_argument('suffix', type=str, help="The suffix for file naming (e.g., _1).")
parser.add_argument('--radius', type=int, default=50, help="Neighborhood radius used in GraphST refinement.")
args = parser.parse_args()

# Dynamically create file paths
data_f = f'./Data/{args.S}_gene/simulated/' 
gene_f = data_f + f'{args.S}_gene_expression{args.suffix}.csv'
coor_f = data_f + f'{args.S}_coord{args.suffix}.csv'
true_f = data_f + f'{args.S}_layer{args.suffix}.csv'

# Check if a GPU is available and set the device accordingly
use_gpu = torch.cuda.is_available()
device = torch.device('cuda:0' if use_gpu else 'cpu')

# Load data
expr = pd.read_csv(gene_f)  # genes as columns
coords = pd.read_csv(coor_f)
layers = pd.read_csv(true_f)

# Create the AnnData object
adata = sc.AnnData(X=expr.round(), dtype=np.float32)

# Add the spatial coordinates
adata.obsm["spatial"] = coords.values

# Add the layer annotations (clusters)
adata.obs["Ground Truth"] = layers.layer.values

print(adata)

# Ensure unique variable names
adata.var_names_make_unique()

# Initialize the GraphST model
model = GraphST.GraphST(adata, device=device)

# Train the model
adata = model.train()

# Methods requiring a fixed cluster count receive the simulated truth count.
cluster_num = int(adata.obs["Ground Truth"].nunique(dropna=True))
if cluster_num < 2:
    raise ValueError(f"Expected at least two spatial domains, found {cluster_num}.")
print(f"Using {cluster_num} spatial domains from the simulation truth.")
tool = 'mclust'

# Perform clustering
if tool == 'mclust':
    clustering(adata, cluster_num, method=tool, radius=args.radius, refinement=True)
elif tool in ['leiden', 'louvain']:
    clustering(adata, cluster_num, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

# Save results
results = adata.obs[['domain', 'Ground Truth']]
results['x'] = coords['x_coord'].values
results['y'] = coords['y_coord'].values
output_file = f'./Data/{args.S}_gene/GraphST{args.suffix}.csv'
results.to_csv(output_file, index_label='Row_names')

# from sklearn.metrics.cluster import adjusted_rand_score

# adjusted_rand_score(adata.obs['domain'], adata.obs['Ground Truth'])
