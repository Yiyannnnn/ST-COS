
import argparse
import pandas as pd
import NaiveDE
import SpatialDE
import numpy as np
import scipy as sp

for attr in ("arange", "array", "argsort", "cumsum", "maximum", "minimum", "zeros_like"):
  if not hasattr(sp, attr):
    setattr(sp, attr, getattr(np, attr))

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run SpatialDE analysis with dynamic file paths.")
parser.add_argument('S', type=str, help="The dataset identifier (e.g., S1).")
parser.add_argument('suffix', type=str, help="The suffix for file naming (e.g., _1).")
args = parser.parse_args()

# Define file paths
data_f = f'./Data/{args.S}_gene/simulated/'
gene_f = data_f + f'{args.S}_gene_expression{args.suffix}.csv'
coord_f = data_f + f'{args.S}_coord{args.suffix}.csv'
spatialde_output_f = f'./Data/{args.S}_gene/spatialDE2{args.suffix}.csv'

# Read in data
expr = pd.read_csv(gene_f)
coords = pd.read_csv(coord_f)
coords.columns = ['x', 'y']

# Ensure indices align
sample_ids = expr.index.astype(str)
expr.index = sample_ids
coords.index = sample_ids

# Transform to genes x samples for NaiveDE utilities
expr_T = expr.T
norm_expr = NaiveDE.stabilize(expr_T)
norm_expr_df = pd.DataFrame(norm_expr, index=expr_T.index, columns=expr_T.columns)
resid_expr = NaiveDE.regress_out(coords, norm_expr_df, '1 + x + y')
resid_expr_df = pd.DataFrame(resid_expr, index=norm_expr_df.index, columns=norm_expr_df.columns)

# Transpose back to samples x genes for SpatialDE
resid_expr_cells = resid_expr_df.T

# Run SpatialDE
results = SpatialDE.run(coords[['x', 'y']].values, resid_expr_cells)

# Save results
results.to_csv(spatialde_output_f, index=False)



# """Run SpatialDE2"""

# module load python/3.7.14
# import pandas as pd
# import numpy as np
# import scanpy as sc
# import pandas as pd
# import numpy as np
# #import tracemalloc
# import SpatialDE as sd
# from tqdm import tqdm
# from functools import partialmethod
# from anndata import AnnData



# df = sd.fit(adata_st, normalized=True, control=None)
# df.set_index("gene", inplace=True)
# df = df.loc[adata_st.var_names]
        
# df.to_csv("./Data/S1_gene/spatialDE2_1.csv")
