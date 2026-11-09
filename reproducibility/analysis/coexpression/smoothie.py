#!/usr/bin/env python3
"""
Simplified Smoothie Spatial Co-expression Analysis Script
This implementation uses only numpy and basic Python libraries
"""

import os
import sys

# Avoid accidentally importing the vendored numpy source tree under STSimulator/env
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
bad_env_path = os.path.join(repo_root, "env")
if bad_env_path in sys.path:
    sys.path.remove(bad_env_path)

import numpy as np
import warnings
warnings.filterwarnings('ignore')

def load_data(gene_expression_file, coord_file, layer_file):
    """
    Load spatial transcriptomics data using basic file reading
    """
    print("Loading data...")
    
    # Load gene expression data (spots x genes)
    gene_expr = []
    gene_names = []
    spot_names = []
    
    with open(gene_expression_file, 'r') as f:
        lines = f.readlines()
        header = lines[0].strip().split(',')[1:]  # Skip first column (spot names)
        gene_names = header
        
        for i, line in enumerate(lines[1:]):
            parts = line.strip().split(',')
            spot_names.append(parts[0])
            values = [float(x) for x in parts[1:]]
            gene_expr.append(values)
    
    gene_expr = np.array(gene_expr)
    print(f"Gene expression shape: {gene_expr.shape}")
    
    # Load coordinates - format: x_coord,y_coord (no index column)
    coords = []
    with open(coord_file, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:  # Skip header
            parts = line.strip().split(',')
            coords.append([float(parts[0]), float(parts[1])])  # x, y coordinates
    
    coords = np.array(coords)
    print(f"Coordinates shape: {coords.shape}")
    
    # Load layer information
    layers = []
    with open(layer_file, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:  # Skip header
            parts = line.strip().split(',')
            if len(parts) >= 2:
                layers.append(int(float(parts[1])))  # layer assignment from second column
            else:
                layers.append(1)  # default layer if parsing fails
    
    layers = np.array(layers)
    print(f"Layer info shape: {layers.shape}")
    
    # Ensure consistent number of spots
    n_spots = min(gene_expr.shape[0], coords.shape[0], len(layers))
    gene_expr = gene_expr[:n_spots, :]
    coords = coords[:n_spots, :]
    layers = layers[:n_spots]
    
    print(f"Using {n_spots} spots and {gene_expr.shape[1]} genes")
    
    return gene_expr, coords, layers, gene_names, spot_names[:n_spots]

def gaussian_smoothing(gene_expr, coords, sigma=46.37, min_spots=25):
    """
    Apply Gaussian smoothing to gene expression data
    """
    print("Running Gaussian smoothing...")
    
    n_spots, n_genes = gene_expr.shape
    smoothed_expr = np.zeros_like(gene_expr)
    
    for i in range(n_spots):
        # Calculate distances to all other spots
        distances = np.sqrt(np.sum((coords - coords[i])**2, axis=1))
        
        # Calculate Gaussian weights
        weights = np.exp(-(distances**2) / (2 * sigma**2))
        
        # Only use spots within reasonable distance (weight > 0.01)
        valid_spots = weights > 0.01
        
        if np.sum(valid_spots) >= min_spots:
            # Normalize weights
            weights = weights[valid_spots]
            weights = weights / np.sum(weights)
            
            # Apply weighted average
            smoothed_expr[i, :] = np.average(gene_expr[valid_spots, :], weights=weights, axis=0)
        else:
            # Not enough neighboring spots, use original expression
            smoothed_expr[i, :] = gene_expr[i, :]
    
    print(f"Smoothed data shape: {smoothed_expr.shape}")
    return smoothed_expr

def compute_correlation_matrix(expression_data):
    """
    Compute pairwise gene correlation matrix
    """
    print("Calculating pairwise gene correlation matrices...")
    
    # Transpose to have genes as rows, spots as columns
    expr_t = expression_data.T
    n_genes = expr_t.shape[0]
    
    # Calculate correlation matrix
    correlation_matrix = np.corrcoef(expr_t)
    
    # Handle NaN values (set to 0)
    correlation_matrix = np.nan_to_num(correlation_matrix)
    
    return correlation_matrix

def create_shuffled_data(gene_expr, coords):
    """
    Create shuffled control dataset
    """
    print("Creating shuffled dataset for null hypothesis testing...")
    
    shuffled_coords = coords.copy()
    np.random.seed(0)
    np.random.shuffle(shuffled_coords)
    
    return shuffled_coords

def simple_clustering(correlation_matrix, gene_names, pcc_cutoff=0.4):
    """
    Simple clustering based on correlation thresholding
    """
    print("Building spatial gene correlation network...")
    
    n_genes = correlation_matrix.shape[0]
    edges = []
    
    # Extract edges above threshold
    for i in range(n_genes):
        for j in range(i+1, n_genes):
            corr = correlation_matrix[i, j]
            if abs(corr) >= pcc_cutoff:
                edges.append([gene_names[i], gene_names[j], corr, corr])
    
    print(f"Network created with {len(edges)} edges")
    
    # Simple community detection using connected components
    if len(edges) > 0:
        # Create adjacency matrix
        adj_matrix = np.abs(correlation_matrix) >= pcc_cutoff
        
        # Find connected components (communities)
        visited = [False] * n_genes
        communities = []
        community_id = 1
        
        def dfs(node, current_community):
            visited[node] = True
            current_community.append(node)
            for neighbor in range(n_genes):
                if adj_matrix[node, neighbor] and not visited[neighbor]:
                    dfs(neighbor, current_community)
        
        for i in range(n_genes):
            if not visited[i]:
                current_community = []
                dfs(i, current_community)
                communities.append((community_id, current_community))
                community_id += 1
        
        # Create node labels
        node_labels = []
        gene_to_community = {}
        for comm_id, genes_in_comm in communities:
            for gene_idx in genes_in_comm:
                gene_to_community[gene_idx] = comm_id
        
        for i, gene in enumerate(gene_names):
            community = gene_to_community.get(i, 1)
            node_labels.append([gene, community])
    else:
        # No significant edges, assign all genes to one cluster
        node_labels = [[gene, 1] for gene in gene_names]
    
    return edges, node_labels

def save_results(output_dir, S, suffix, correlation_matrix, shuffled_correlation_matrix, 
                edges, node_labels, gene_expr, smoothed_expr, coords, layers, gene_names):
    """
    Save all analysis results
    """
    print("Saving results...")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Save correlation matrices
    np.save(os.path.join(output_dir, f"{S}{suffix}_correlation_matrix.npy"), correlation_matrix)
    np.save(os.path.join(output_dir, f"{S}{suffix}_shuffled_correlation_matrix.npy"), shuffled_correlation_matrix)
    
    # Save edge list
    edge_file = os.path.join(output_dir, f"{S}{suffix}_edge_list.csv")
    with open(edge_file, 'w') as f:
        f.write("gene1,gene2,PCC,Rescaled_PCC\n")
        for edge in edges:
            f.write(f"{edge[0]},{edge[1]},{edge[2]},{edge[3]}\n")
    
    # Save node labels
    node_file = os.path.join(output_dir, f"{S}{suffix}_node_labels.csv")
    with open(node_file, 'w') as f:
        f.write("gene_name,community\n")
        for node in node_labels:
            f.write(f"{node[0]},{node[1]}\n")
    
    # Calculate and save statistics
    lower_tri_indices = np.tril_indices(correlation_matrix.shape[0], -1)
    true_values = correlation_matrix[lower_tri_indices]
    shuffled_values = shuffled_correlation_matrix[lower_tri_indices]
    
    stats = {
        'True_95th': np.percentile(true_values, 95),
        'True_99th': np.percentile(true_values, 99),
        'True_99.9th': np.percentile(true_values, 99.9),
        'Shuffled_95th': np.percentile(shuffled_values, 95),
        'Shuffled_99th': np.percentile(shuffled_values, 99),
        'Shuffled_99.9th': np.percentile(shuffled_values, 99.9),
    }
    
    stats_file = os.path.join(output_dir, f"{S}{suffix}_correlation_stats.csv")
    with open(stats_file, 'w') as f:
        f.write("Metric,Value\n")
        for metric, value in stats.items():
            f.write(f"{metric},{value}\n")
    
    # Save smoothed expression data
    smoothed_file = os.path.join(output_dir, f"{S}{suffix}_smoothed_expression.csv")
    with open(smoothed_file, 'w') as f:
        f.write("spot," + ",".join(gene_names) + "\n")
        for i in range(smoothed_expr.shape[0]):
            f.write(f"spot_{i}," + ",".join([str(x) for x in smoothed_expr[i, :]]) + "\n")
    
    # Create summary statistics
    n_communities = len(set([node[1] for node in node_labels]))
    
    summary_stats = {
        'Total_spots': gene_expr.shape[0],
        'Total_genes': gene_expr.shape[1],
        'Number_of_edges': len(edges),
        'Number_of_communities': n_communities,
        'Gaussian_SD': 46.37,
        'Min_spots_under_gaussian': 25
    }
    
    summary_file = os.path.join(output_dir, f"{S}{suffix}_summary_stats.csv")
    with open(summary_file, 'w') as f:
        f.write("Metric,Value\n")
        for metric, value in summary_stats.items():
            f.write(f"{metric},{value}\n")
    
    print(f"Smoothie analysis completed! Results saved to: {output_dir}")
    return summary_stats

def main():
    # Parse command line arguments
    if len(sys.argv) != 3:
        print("Usage: python smoothie.py <S> <suffix>")
        print("Example: python smoothie.py S1 _1")
        sys.exit(1)
    
    S = sys.argv[1]
    suffix = sys.argv[2]
    
    # Create output directory
    output_dir = "./coexpression/smoothie_Agent"
    
    # Construct file paths
    gene_expression_file = f"./Data/{S}_gene/simulated/{S}_gene_expression{suffix}.csv"
    coord_file = f"./Data/{S}_gene/simulated/{S}_coord{suffix}.csv"
    layer_file = f"./Data/{S}_gene/simulated/{S}_layer{suffix}.csv"
    
    # Check if input files exist
    for file_path in [gene_expression_file, coord_file, layer_file]:
        if not os.path.exists(file_path):
            print(f"Error: Input file not found: {file_path}")
            sys.exit(1)
    
    print(f"Starting Smoothie analysis for {S}{suffix}")
    print(f"Output directory: {output_dir}")
    
    try:
        # Load data
        gene_expr, coords, layers, gene_names, spot_names = load_data(
            gene_expression_file, coord_file, layer_file)
        
        # Quality control - filter low expression genes
        gene_totals = np.sum(gene_expr, axis=0)
        valid_genes = gene_totals >= 100  # Keep genes with at least 100 total counts
        
        gene_expr_filtered = gene_expr[:, valid_genes]
        gene_names_filtered = [gene_names[i] for i in range(len(gene_names)) if valid_genes[i]]
        
        print(f"Filtered genes: {len(gene_names)} -> {len(gene_names_filtered)}")
        
        # Normalization (CPM + log1p)
        print("Normalizing data...")
        gene_expr_norm = gene_expr_filtered.copy().astype(float)
        
        # CPM normalization
        spot_totals = np.sum(gene_expr_norm, axis=1)
        gene_expr_norm = (gene_expr_norm.T / spot_totals * 1000).T
        
        # Log1p transformation
        gene_expr_norm = np.log1p(gene_expr_norm)
        
        # Gaussian smoothing
        smoothed_expr = gaussian_smoothing(gene_expr_norm, coords)
        
        # Create shuffled control
        shuffled_coords = create_shuffled_data(gene_expr_norm, coords)
        smoothed_shuffled = gaussian_smoothing(gene_expr_norm, shuffled_coords)
        
        # Compute correlation matrices
        correlation_matrix = compute_correlation_matrix(smoothed_expr)
        shuffled_correlation_matrix = compute_correlation_matrix(smoothed_shuffled)
        
        # Calculate correlation statistics
        lower_tri_indices = np.tril_indices(correlation_matrix.shape[0], -1)
        true_values = correlation_matrix[lower_tri_indices]
        shuffled_values = shuffled_correlation_matrix[lower_tri_indices]
        
        print("Correlation distribution statistics:")
        print(f"True 99th percentile: {np.percentile(true_values, 99):.3f}")
        print(f"Shuffled 99th percentile: {np.percentile(shuffled_values, 99):.3f}")
        
        # Set PCC cutoff based on shuffled data
        pcc_cutoff = max(0.4, np.percentile(shuffled_values, 99) + 0.05)
        print(f"Selected PCC cutoff: {pcc_cutoff:.3f}")
        
        # Network analysis
        edges, node_labels = simple_clustering(correlation_matrix, gene_names_filtered, pcc_cutoff)
        
        # Save results
        summary_stats = save_results(
            output_dir, S, suffix, correlation_matrix, shuffled_correlation_matrix,
            edges, node_labels, gene_expr_filtered, smoothed_expr, coords, layers, gene_names_filtered)
        
        print(f"Smoothie analysis completed successfully for {S}{suffix}!")
        print(f"Found {summary_stats['Number_of_communities']} communities with {summary_stats['Number_of_edges']} edges")
        
    except Exception as e:
        print(f"Error during Smoothie analysis for {S}{suffix}: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
