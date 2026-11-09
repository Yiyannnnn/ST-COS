"""Run SpatialDE"""
#MSC
import argparse
import pandas as pd
import numpy as np
import squidpy as sq
from anndata import AnnData

import sys
sys.path.append('./env/STCOS_eval/lib/python3.8/site-packages/SpatialDE/')

#import NaiveDE
import SpatialDE
from SpatialDE.base import *

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run SpatialDE analysis with dynamic file paths.")
parser.add_argument('S', type=str, help="The dataset identifier (e.g., S1).")
parser.add_argument('suffix', type=str, help="The suffix for file naming (e.g., _1).")
args = parser.parse_args()

# Dynamically create file paths
data_f = f'./Data/{args.S}_gene/simulated/'
gene_f = data_f + f'{args.S}_gene_expression{args.suffix}.csv'
coord_f = data_f + f'{args.S}_coord{args.suffix}.csv'
output_f =  f'./Data/{args.S}_gene/spatialDE{args.suffix}.csv'
moran_output_f = './Data/' + args.S + '_gene/moranI' + args.suffix + '.csv'

# Read in data
st = pd.read_csv(gene_f)
spatial_coord = pd.read_csv(coord_f)

# Create AnnData object for Moran's I
adata_st = AnnData(st)
adata_st.obsm['spatial'] = spatial_coord.values

# Run Moran's I analysis
sq.gr.spatial_neighbors(adata_st)
sq.gr.spatial_autocorr(
    adata_st,
    mode="moran",
    genes=adata_st.var_names,
    n_perms=100,
    n_jobs=1,
)

# Save Moran's I results
adata_st.uns["moranI"].to_csv(moran_output_f)

# Run SpatialDE
df = run(spatial_coord, st)

# Process and save the results
df.set_index("g", inplace=True)
df = df.loc[~df.index.duplicated(keep='first')]
df.to_csv(output_f)



# """Moran's I """
# import squidpy as sq

# st = pd.read_csv("./Data/S1_gene/simulated/S1_gene_expression_1.csv")
# adata_st = AnnData(st)

# spatial_coord = pd.read_csv("./Data/S1_gene/simulated/S1_coord_1.csv")
# adata_st.obsm['spatial'] = spatial_coord.values

# sq.gr.spatial_neighbors(adata_st)
# sq.gr.spatial_autocorr(
#     adata_st,
#     mode="moran",
#     genes=adata_st.var_names,
#     n_perms=100,
#     n_jobs=1,
# )

# adata_st.uns["moranI"].to_csv("./Data/S1_gene/moranI_1.csv")

# """Run SpatialDE"""

# st = pd.read_csv("./Data/S1_gene/simulated/S1_gene_expression_1.csv")
# spatial_coord = pd.read_csv("./Data/S1_gene/simulated/S1_coord_1.csv")

# df = SpatialDE.run(spatial_coord, st)
    
# df.set_index("g", inplace=True)
    
# df = df.loc[~df.index.duplicated(keep='first')]

# df.to_csv("./Data/S1_gene/spatialDE_1.csv")
