import argparse
import numpy as np
import pandas as pd
from anndata import AnnData
import time

# import tangram for spatial deconvolution
import tangram as tg


# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run Tangram analysis with dynamic file paths.")
parser.add_argument('S', type=str, help="The dataset identifier (e.g., S1).")
parser.add_argument('suffix', type=str, help="The suffix for file naming (e.g., _1).")
args = parser.parse_args()

# Define file paths
data_f = f'./Data/{args.S}_gene/simulated/'
gene_f = data_f + f'{args.S}_gene_expression{args.suffix}.csv'
coord_f = data_f + f'{args.S}_coord{args.suffix}.csv'
layer_f = data_f + f'{args.S}_layer{args.suffix}.csv'
sc_ref_f = data_f + f'{args.S}_scref{args.suffix}.csv'
output_f = f'./Data/{args.S}_gene/Tangram{args.suffix}.csv'

# Read in ST data
st = pd.read_csv(gene_f)
adata_st = AnnData(st)

# Append cell count for each spot
adata_st.obs['cell_count'] = 10

# Create spatial coordinates information
spatial_coord = pd.read_csv(coord_f)
adata_st.obsm['spatial'] = spatial_coord.values

# Read in layer information
layer = pd.read_csv(layer_f)
image_features = pd.DataFrame(index=adata_st.obs.index)
image_features['segmentation_label'] = layer['layer'].values
image_features['segmentation_centroid_x'] = layer['x_coord'].values
image_features['segmentation_centroid_y'] = layer['y_coord'].values

# Assign image features to adata_st
adata_st.obsm['image_features'] = image_features

# Read in scRNA-seq data
scdat = pd.read_csv(sc_ref_f)
sc_meta = scdat.iloc[:, 0].to_frame()
sc_meta.iloc[:, 0] = pd.Categorical(sc_meta.iloc[:, 0])

adata_sc = AnnData(scdat.iloc[:, 1:])
adata_sc.obs_names = sc_meta.index
adata_sc.obs = adata_sc.obs.merge(sc_meta, how='left', left_index=True, right_index=True)

# Start time tracking
start_time = time.time()

# Preprocess: find common genes between sc and st
tg.pp_adatas(adata_sc, adata_st, genes=None)

# Deconvolution
ad_map = tg.map_cells_to_space(
    adata_sc,
    adata_st,
    mode="constrained",
    target_count=adata_st.obs.cell_count.sum(),
    density_prior=np.array(adata_st.obs.cell_count) / adata_st.obs.cell_count.sum(),
    num_epochs=1000,
    device='cpu'
)

# Map the cell type information to the ST AnnData object
tg.project_cell_annotations(ad_map, adata_st, annotation="Cell_type")

# End time tracking
end_time = time.time()
print("--- %.2f seconds ---" % (end_time - start_time))

# Normalize the probability matrix and save as CSV
prob_mat = adata_st.obsm["tangram_ct_pred"]
prob_mat = prob_mat.div(prob_mat.sum(axis=1), axis=0)
prob_mat['spot'] = [f"spot_{i}" for i in range(1, len(prob_mat) + 1)]
prob_mat.to_csv(output_f)


# # read in ST data
# st = pd.read_csv("./Data/S1_gene/S1_gene_expression.csv")
# adata_st = AnnData(st)
# adata_st

# ## append cell count for each spot
# adata_st.obs['cell_count'] = 10

# ## create spatial coordinates information
# spatial_coord = pd.read_csv("./Data/S1_gene/S1_coord.csv")
# adata_st.obsm['spatial'] = spatial_coord.values


# layer = pd.read_csv("./Data/S1_gene/S1_layer.csv")

# # Assign the centroid information to adata_st
# image_features = pd.DataFrame(index=adata_st.obs.index)
# image_features['segmentation_label'] = layer['layer'].values
# image_features['segmentation_centroid_x'] = layer['x_coord'].values
# image_features['segmentation_centroid_y'] = layer['y_coord'].values

# # Assign the image_features to adata_st.obsm
# adata_st.obsm['image_features'] = image_features

# # Read in scRNA-seq data
# scdat = pd.read_csv("./Data/S1_gene/S1_scref.csv")

# sc_meta = scdat.iloc[:, 0].to_frame()
# sc_meta.iloc[:,0] = pd.Categorical(sc_meta.iloc[:,0])

# adata_sc = AnnData(scdat.iloc[:, 1:])
# adata_sc.obs_names = sc_meta.index
# adata_sc.obs = adata_sc.obs.merge(sc_meta, how = 'left', left_index=True, right_index=True)
# adata_sc

# start_time = time.time()
# # preprocessing: find the common genes between sc and st
# tg.pp_adatas(adata_sc, adata_st, genes=None)


# # Deconvolution
# ad_map = tg.map_cells_to_space(
#     adata_sc,
#     adata_st,
#     mode="constrained",
#     target_count=adata_st.obs.cell_count.sum(),
#     density_prior=np.array(adata_st.obs.cell_count) / adata_st.obs.cell_count.sum(),
#     num_epochs=1000,
#     #device="cuda:0",
#     device='cpu'
# )

# # Gather deconvolution results
# ## map the cell type information to the st AnnData object
# ## The output created is the unnormalized probability matrix
# tg.project_cell_annotations(ad_map, adata_st, annotation="Cell_type")
# end_time = time.time()
# print("--- %.2f seconds ---" % (time.time() - start_time))

# ## normalize the probability matrix and save as csv
# prob_mat = adata_st.obsm["tangram_ct_pred"]
# prob_mat = prob_mat.div(prob_mat.sum(axis=1), axis=0)
# prob_mat['spot'] = [f"spot_{i}" for i in range(1, len(prob_mat) + 1)]
# prob_mat.to_csv("./Data/S1_gene/Tangram_1.csv")

