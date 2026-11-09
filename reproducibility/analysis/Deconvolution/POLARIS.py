#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Run POLARIS deconvolution for one scenario replicate.")
    parser.add_argument("S", type=str, help="Scenario id, e.g. S1")
    parser.add_argument("suffix", type=str, help="Replicate suffix, e.g. _1")
    parser.add_argument("--epochs", type=int, default=200, help="Epochs for both ST and sc stage.")
    parser.add_argument("--batch_size", type=int, default=512, help="Batch size for both ST and sc.")
    args = parser.parse_args()

    root = "."
    sim_dir = f"{root}/Data/{args.S}_gene/simulated"
    out_csv = f"{root}/Data/{args.S}_gene/POLARIS{args.suffix}.csv"

    sc_ref_file = f"{sim_dir}/{args.S}_scref{args.suffix}.csv"
    st_cnt_file = f"{sim_dir}/{args.S}_gene_expression{args.suffix}.csv"
    layer_file = f"{sim_dir}/{args.S}_layer{args.suffix}.csv"

    if not (os.path.exists(sc_ref_file) and os.path.exists(st_cnt_file) and os.path.exists(layer_file)):
        raise FileNotFoundError("One or more required simulation files are missing.")

    run_root = f"{root}/Data/{args.S}_gene/POLARIS_input{args.suffix}"
    os.makedirs(run_root, exist_ok=True)

    # Build POLARIS-formatted files (TSV with index in col 1).
    sc_ref = pd.read_csv(sc_ref_file)
    cell_ids = [f"cell_{i}" for i in range(1, sc_ref.shape[0] + 1)]
    sc_cnt = sc_ref.drop(columns=["Cell_type"]).copy()
    sc_cnt.index = cell_ids
    sc_cnt_tsv = f"{run_root}/sc_cnt.tsv"
    sc_cnt.to_csv(sc_cnt_tsv, sep="\t", index=True)

    sc_labels = pd.DataFrame({"cell_type": sc_ref["Cell_type"].values}, index=cell_ids)
    sc_labels_tsv = f"{run_root}/sc_labels.tsv"
    sc_labels.to_csv(sc_labels_tsv, sep="\t", index=True)

    st_cnt = pd.read_csv(st_cnt_file)
    spot_ids = [f"spot_{i}" for i in range(1, st_cnt.shape[0] + 1)]
    st_cnt.index = spot_ids
    st_cnt_tsv = f"{run_root}/st_cnt.tsv"
    st_cnt.to_csv(st_cnt_tsv, sep="\t", index=True)

    layer = pd.read_csv(layer_file)
    st_label = pd.DataFrame({"layer": layer["layer"].values}, index=spot_ids)
    st_label_tsv = f"{run_root}/st_label.tsv"
    st_label.to_csv(st_label_tsv, sep="\t", index=True)

    polaris_repo = f"{root}/Deconvolution/POLARIS_repo_20260223/source_code"
    out_dir = f"{run_root}/polaris_out"
    os.makedirs(out_dir, exist_ok=True)

    polaris_python = os.environ.get(
        "POLARIS_PYTHON", f"{root}/env/AD_drug/bin/python"
    )
    if not os.path.exists(polaris_python):
        polaris_python = sys.executable

    cmd = [
        polaris_python,
        "train.py",
        "--sc_cnt", sc_cnt_tsv,
        "--sc_labels", sc_labels_tsv,
        "--st_cnt", st_cnt_tsv,
        "--st_source", "simulation",
        "--st_label", st_label_tsv,
        "--st_batch_size", str(args.batch_size),
        "--sc_batch_size", str(args.batch_size),
        "--st_epochs", str(args.epochs),
        "--sc_epochs1", str(args.epochs),
        "--out_dir", out_dir,
        "--prefix", "polaris"
    ]

    subprocess.run(cmd, cwd=polaris_repo, check=True)

    section_name = os.path.splitext(os.path.basename(st_cnt_tsv))[0]
    pred_tsv = f"{out_dir}/{section_name}/polaris.tsv"
    if not os.path.exists(pred_tsv):
        raise FileNotFoundError(f"POLARIS output not found: {pred_tsv}")

    pred = pd.read_csv(pred_tsv, sep="\t", index_col=0)
    # Convert POLARIS intensities to proportions for fair RMSE comparison.
    pred = pred.div(pred.sum(axis=1), axis=0)
    pred["spot"] = pred.index
    pred.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    main()
