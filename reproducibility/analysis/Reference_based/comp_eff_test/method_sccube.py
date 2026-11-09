#!/usr/bin/env python3
# scCube Method Implementation for Efficiency Testing

import sys
import os
import time
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

def run_sccube(gene_count, output_dir):
    """Run scCube method for efficiency testing"""
    
    # Import required libraries
    import torch
    import scanpy as sc
    
    # Add scCube to path
    sys.path.append('./temp/scCube')
    
    import scCube
    from scCube import scCube as scCubeModel
    
    # Load gene names
    gene_names = pd.read_csv("./data_external/DLPFC/151673_selected_gene.csv")
    gene_names = gene_names.tolist()
    genes = gene_names[:gene_count]  # Select subset of genes
    
    # Load spatial transcriptomics data
    sc_adata = sc.read_10x_h5('./data_external/DLPFC/151673/151673_filtered_feature_bc_matrix.h5')
    
    # Load metadata with cell types and coordinates
    sc_adata_paper = sc.read_h5ad('./data_external/DLPFC/DLPFC_151673_adata.h5ad')
    sc_adata = sc_adata[sc_adata_paper.obs_names]
    sc_adata.obs['Cell_type'] = sc_adata_paper.obs['Cell_type']
    sc_adata.obs['Cell'] = sc_adata_paper.obs['Cell']
    
    # Ensure var_names are unique
    sc_adata_paper.var_names_make_unique()
    
    # Subset to selected genes
    sc_adata = sc_adata[:, genes]
    
    # Create log-transformed layer
    sc_adata.layers["log_transformed"] = sc_adata.X
    sc_data = sc_adata.to_df(layer="log_transformed").T
    sc_meta = sc_adata.obs
    
    # Initialize scCube model
    model = scCubeModel()
    
    # Generate synthetic data
    generate_sc_meta, generate_sc_data = model.load_vae_and_generate_cell(
        sc_adata=sc_adata,
        celltype_key='Cell_type',
        cell_key='Cell',
        target_num=None,
        hidden_size=128,
        load_path='./data_external/DLPFC/models/slice151673_epoch100000.pth',
        used_device='cpu'
    )
    
    # Generate pattern reference
    generate_sc_data, generate_sc_meta = model.generate_pattern_reference(
        sc_adata=sc_adata,
        generate_sc_data=generate_sc_data,
        generate_sc_meta=generate_sc_meta,
        celltype_key='Cell_type',
        spatial_key=['x', 'y'],
        cost_metric='sqeuclidean'
    )
    
    # Save results
    generate_sc_data.T.to_csv(
        os.path.join(output_dir, f"sccube_simulated_counts_{gene_count}.csv")
    )
    
    generate_sc_meta.to_csv(
        os.path.join(output_dir, f"sccube_coordinates_{gene_count}.csv")
    )
    
    return True

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python method_sccube.py <gene_count>")
        sys.exit(1)
    
    gene_count = int(sys.argv[1])
    output_dir = "./Reference_based/comp_eff_test"
    
    try:
        start_time = time.time()
        result = run_sccube(gene_count, output_dir)
        end_time = time.time()
        
        # Save timing information
        timing_info = pd.DataFrame({
            'method': ['scCube'],
            'gene_count': [gene_count],
            'time_seconds': [end_time - start_time],
            'success': [result]
        })
        timing_info.to_csv(
            os.path.join(output_dir, f"sccube_timing_{gene_count}.csv"),
            index=False
        )
        
        print(f"scCube completed for {gene_count} genes in {end_time - start_time:.2f} seconds")
        
    except Exception as e:
        print(f"scCube failed for {gene_count} genes: {str(e)}")
        # Save error information
        error_info = pd.DataFrame({
            'method': ['scCube'],
            'gene_count': [gene_count],
            'time_seconds': [np.nan],
            'success': [False],
            'error': [str(e)]
        })
        error_info.to_csv(
            os.path.join(output_dir, f"sccube_error_{gene_count}.csv"),
            index=False
        )
        sys.exit(1)