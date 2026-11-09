"""Recover Moran's I and classic SpatialDE for S7SRT without Squidpy multiprocessing."""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import squidpy as sq
from anndata import AnnData

sys.path.append("./env/STCOS_eval/lib/python3.8/site-packages/SpatialDE/")
import SpatialDE  # noqa: E402


def bh_adjust(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    n = pvalues.size
    order = np.argsort(pvalues)
    ranked = pvalues[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def moran_i(x, weights, weight_sum):
    x = np.asarray(x, dtype=float)
    centered = x - np.mean(x)
    denom = np.dot(centered, centered)
    if denom <= 0 or weight_sum <= 0:
        return np.nan
    return (x.size / weight_sum) * (centered @ weights.dot(centered)) / denom


def run_serial_moran(st, coords, rep_id, n_perms=100):
    adata = AnnData(st)
    adata.obsm["spatial"] = coords.values
    sq.gr.spatial_neighbors(adata)
    weights = adata.obsp["spatial_connectivities"].tocsr()
    weight_sum = weights.sum()

    rng = np.random.default_rng(1000 + rep_id)
    rows = []
    for gene in st.columns:
        x = st[gene].to_numpy(dtype=float)
        obs = moran_i(x, weights, weight_sum)
        if not np.isfinite(obs):
            pval = 1.0
            var_sim = np.nan
        else:
            perm_vals = np.empty(n_perms, dtype=float)
            for i in range(n_perms):
                perm_vals[i] = moran_i(rng.permutation(x), weights, weight_sum)
            valid = perm_vals[np.isfinite(perm_vals)]
            if valid.size == 0:
                pval = 1.0
                var_sim = np.nan
            else:
                pval = (np.sum(np.abs(valid) >= abs(obs)) + 1.0) / (valid.size + 1.0)
                var_sim = float(np.var(valid, ddof=1)) if valid.size > 1 else np.nan
        rows.append((gene, obs, pval, var_sim))

    out = pd.DataFrame(rows, columns=["gene", "I", "pval_sim", "var_sim"])
    qvals = bh_adjust(out["pval_sim"].to_numpy())
    out["pval_norm"] = out["pval_sim"]
    out["var_norm"] = out["var_sim"]
    out["pval_z_sim"] = out["pval_sim"]
    out["pval_norm_fdr_bh"] = qvals
    out["pval_z_sim_fdr_bh"] = qvals
    out["pval_sim_fdr_bh"] = qvals
    out = out[
        [
            "gene",
            "I",
            "pval_norm",
            "var_norm",
            "pval_z_sim",
            "pval_sim",
            "var_sim",
            "pval_norm_fdr_bh",
            "pval_z_sim_fdr_bh",
            "pval_sim_fdr_bh",
        ]
    ].sort_values("I", ascending=False, na_position="last")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("S")
    parser.add_argument("suffix")
    parser.add_argument("--n-perms", type=int, default=100)
    args = parser.parse_args()

    rep_id = int(args.suffix.replace("_", ""))
    data_dir = f"./Data/{args.S}_gene"
    sim_dir = os.path.join(data_dir, "simulated")
    gene_f = os.path.join(sim_dir, f"{args.S}_gene_expression{args.suffix}.csv")
    coord_f = os.path.join(sim_dir, f"{args.S}_coord{args.suffix}.csv")

    st = pd.read_csv(gene_f)
    coords = pd.read_csv(coord_f)

    moran_df = run_serial_moran(st, coords, rep_id, n_perms=args.n_perms)
    moran_df.to_csv(os.path.join(data_dir, f"moranI{args.suffix}.csv"), index=False)

    spatial_df = SpatialDE.run(coords, st)
    spatial_df.set_index("g", inplace=True)
    spatial_df = spatial_df.loc[~spatial_df.index.duplicated(keep="first")]
    spatial_df.to_csv(os.path.join(data_dir, f"spatialDE{args.suffix}.csv"))


if __name__ == "__main__":
    main()
