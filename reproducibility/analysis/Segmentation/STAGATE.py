import argparse

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys

from sklearn.metrics.cluster import adjusted_rand_score

import STAGATE_pyG


# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run STAGATE analysis with dynamic file paths.")
parser.add_argument('S', type=str, help="The dataset identifier (e.g., S1).")
parser.add_argument('suffix', type=str, help="The suffix for file naming (e.g., _1).")
parser.add_argument(
    "--graph_model",
    type=str,
    default="KNN",
    choices=["KNN", "Radius"],
    help="Spatial graph construction model."
)
parser.add_argument(
    "--k_cutoff",
    type=int,
    default=8,
    help="Number of nearest neighbors when graph_model=KNN."
)
parser.add_argument(
    "--rad_cutoff",
    type=float,
    default=3.0,
    help="Radius cutoff when graph_model=Radius."
)
args = parser.parse_args()

# Dynamically create file paths
data_f = f'./Data/{args.S}_gene/simulated/{args.S}_'
gene_f = data_f + f'gene_expression{args.suffix}.csv'
coor_f = data_f + f'coord{args.suffix}.csv'
true_f = data_f + f'layer{args.suffix}.csv'

# Read in data
expr = pd.read_csv(gene_f)  # genes as columns
coords = pd.read_csv(coor_f)
layers = pd.read_csv(true_f)

# Create the AnnData object
adata = sc.AnnData(X=expr)

# Add the spatial coordinates
adata.obsm["spatial"] = coords.values

# Add the layer annotations (clusters)
adata.obs["Ground Truth"] = layers.layer.values

print(adata)

# Preprocessing steps
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Calculate spatial network and train the STAGATE model
if args.graph_model == "KNN":
    if args.k_cutoff <= 0:
        raise ValueError("--k_cutoff must be > 0 when --graph_model KNN.")
    print(f"Building spatial graph with KNN, k_cutoff={args.k_cutoff}")
    STAGATE_pyG.Cal_Spatial_Net(adata, model="KNN", k_cutoff=args.k_cutoff)
else:
    if args.rad_cutoff <= 0:
        raise ValueError("--rad_cutoff must be > 0 when --graph_model Radius.")
    print(f"Building spatial graph with Radius, rad_cutoff={args.rad_cutoff}")
    STAGATE_pyG.Cal_Spatial_Net(adata, model="Radius", rad_cutoff=args.rad_cutoff)

adata = STAGATE_pyG.train_STAGATE(adata)

# Perform clustering
sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)

intrim_file = f'./Data/{args.S}_gene/{args.S}_STAGATE{args.suffix}.csv'
df = pd.DataFrame(adata.obsm["STAGATE"])
df.to_csv(intrim_file) 

# Methods requiring a fixed cluster count receive the simulated truth count.
n_clusters = int(adata.obs["Ground Truth"].nunique(dropna=True))
if n_clusters < 2:
    raise ValueError(f"Expected at least two spatial domains, found {n_clusters}.")
print(f"Using {n_clusters} spatial domains from the simulation truth.")
adata = STAGATE_pyG.mclust_R(
    adata, used_obsm='STAGATE', num_cluster=n_clusters
)

# Save results
results = adata.obs[['mclust', 'Ground Truth']]
results['x'] = coords['x_coord'].values
results['y'] = coords['y_coord'].values

# Compute the adjusted Rand index (ARI) between predicted and ground truth clusters
# ari = adjusted_rand_score(results['mclust'], results['Ground Truth'])
# print(f'Adjusted Rand Index (ARI): {ari}')

# Save the results to a CSV file
output_file = f'./Data/{args.S}_gene/STAGATE{args.suffix}.csv'
results.to_csv(output_file, index_label='Row_names')

# data_f = './Data/S2_gene/simulated/S2_' 
# gene_f = data_f + 'gene_expression_3' + '.csv'
# coor_f = data_f + 'coord_3'  + '.csv'
# true_f = data_f + 'layer_3'  + '.csv'

# expr = pd.read_csv(gene_f)  # genes as columns
# coords = pd.read_csv(coor_f)
# layers = pd.read_csv(true_f)

# # Create the AnnData object
# adata = sc.AnnData(X=expr)

# # Add the spatial coordinates
# adata.obsm["spatial"] = coords.values

# # Add the layer annotations (clusters)
# adata.obs["Ground Truth"] = layers.layer.values

# print(adata)

# sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

# STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=150)
# adata = STAGATE_pyG.train_STAGATE(adata)

# sc.pp.neighbors(adata, use_rep='STAGATE')
# sc.tl.umap(adata)

# adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=3)
# results = adata.obs[['mclust','Ground Truth']]
# results['x'] = coords['x_coord'].values
# results['y'] = coords['y_coord'].values
# adjusted_rand_score(results['mclust'], results['Ground Truth'])
# results.to_csv('./Data/S2_gene/STAGATE_3.csv', index_label='Row_names')


# library(data.table)
# library(mclust)
# test <- fread("S4_STAGATE_5.csv")[-1,]
# res <- Mclust(test,G = 4, moduleName = "EEE")
# layer <- fread("./Data/S4_gene/simulated/S4_layer_5.csv")
# res_df <- cbind(mclust = res$classification,`Ground Truth` = layer$layer, x = layer[,1], y = layer[,2])
# fwrite(res_df, "./Data/S4_gene/STAGATE_5.csv")
