#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("check_root")
    parser.add_argument("scenario")
    parser.add_argument("suffix")
    args = parser.parse_args()

    root = Path(args.check_root)
    scenario_dir = root / "Data" / f"{args.scenario}_gene"
    sim_dir = scenario_dir / "simulated"
    run_dir = scenario_dir / f"STdGCN_input{args.suffix}"
    st_dir = run_dir / "ST_path"
    sc_dir = run_dir / "sc_path"
    out_dir = run_dir / "output"
    for directory in (st_dir, sc_dir, out_dir):
        directory.mkdir(parents=True, exist_ok=True)

    expr = pd.read_csv(sim_dir / f"{args.scenario}_gene_expression{args.suffix}.csv")
    coord = pd.read_csv(sim_dir / f"{args.scenario}_coord{args.suffix}.csv")
    scref = pd.read_csv(sim_dir / f"{args.scenario}_scref{args.suffix}.csv")

    # STdGCN requires the single-cell and spatial matrices to share features.
    gene_columns = [column for column in expr.columns if column in scref.columns]
    if not gene_columns:
        raise ValueError("No shared genes between spatial counts and single-cell reference")

    spot_ids = [f"spot_{i}" for i in range(1, len(expr) + 1)]
    cell_ids = [str(i) for i in range(1, len(scref) + 1)]

    expr.index = spot_ids
    expr.to_csv(st_dir / "ST_data.tsv", sep="\t", index=True)
    pd.DataFrame(
        {"x": coord["x_coord"].to_numpy(), "y": coord["y_coord"].to_numpy()},
        index=spot_ids,
    ).to_csv(st_dir / "coordinates.csv", index=True)

    scref[gene_columns].set_axis(cell_ids).to_csv(
        sc_dir / "sc_data.tsv", sep="\t", index=False
    )
    pd.DataFrame({"Cell_type": scref["Cell_type"].to_numpy()}, index=cell_ids).to_csv(
        sc_dir / "sc_label.tsv", sep="\t", index=True, index_label="cell_id"
    )


if __name__ == "__main__":
    main()
