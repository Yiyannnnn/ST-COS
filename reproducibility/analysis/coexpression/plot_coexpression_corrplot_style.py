#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
LOCAL_PY_LIB = REPO_ROOT / ".pylibs_venn"
if LOCAL_PY_LIB.exists():
    sys.path.insert(0, str(LOCAL_PY_LIB))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd
import rdata
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform


@dataclass
class MatrixBundle:
    method: str
    matrix: np.ndarray
    labels: list[str]
    description: str
    domain: str


BASE_METHOD_ORDER = ["spMOCA", "SpaceX", "CSSN", "hdWGCNA", "Smoothie", "binSpect", "stIHC"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot co-expression matrices in a corrplot-style panel for the requested "
            "STSimulator scenario/replicate combinations."
        )
    )
    parser.add_argument(
        "--scenario",
        nargs="+",
        default=["S7"],
        help="Scenario IDs to plot. Default: S7.",
    )
    parser.add_argument(
        "--rep",
        nargs="+",
        type=int,
        help="Replicate IDs to plot. Default: all replicates with all six methods present.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parent / "corrplot_style"),
        help="Directory for output panels.",
    )
    parser.add_argument(
        "--cluster-max-size",
        type=int,
        default=400,
        help="Maximum matrix size for hierarchical reordering. Default: 400.",
    )
    return parser.parse_args()


def normalize_method_name(method: str) -> str:
    return method.lower().replace("_", "").replace("-", "")


def natural_key(text: str) -> list[object]:
    import re

    return [int(tok) if tok.isdigit() else tok.lower() for tok in re.split(r"(\d+)", text)]


def infer_labels(prefix: str, size: int) -> list[str]:
    return [f"{prefix}{idx}" for idx in range(1, size + 1)]


def truth_file_path(base_dir: Path, scenario: str) -> Path:
    return base_dir / "Data" / f"{scenario}_gene" / "simulated" / f"{scenario}_coexpression_truth.csv"


def truth_available(base_dir: Path, scenario: str, rep: int) -> bool:
    path = truth_file_path(base_dir, scenario)
    if not path.exists():
        return False
    frame = pd.read_csv(path)
    if "scenario" not in frame.columns:
        return not frame.empty
    rep_values = pd.to_numeric(frame["scenario"], errors="coerce")
    return bool((rep_values == rep).any())


def validate_square_matrix(matrix: np.ndarray, method: str) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{method} matrix is not square: shape={matrix.shape}")
    matrix = np.nan_to_num(matrix, nan=0.0, posinf=1.0, neginf=-1.0)
    return matrix


def load_spacex(base_dir: Path, scenario: str, rep: int) -> MatrixBundle:
    path = base_dir / "Data" / f"{scenario}_gene" / f"SpaceX_{rep}.csv"
    frame = pd.read_csv(path)
    matrix = validate_square_matrix(frame.to_numpy(), "SpaceX")
    labels = infer_labels("gene_", matrix.shape[0])
    return MatrixBundle("SpaceX", matrix, labels, "shared network", "gene")


def load_spmoca(base_dir: Path, scenario: str, rep: int) -> MatrixBundle:
    mat_path = base_dir / "coexpression" / "spMOCA_Agent" / f"{scenario}_{rep}_correlation_matrix.csv"
    label_path = base_dir / "coexpression" / "spMOCA_Agent" / f"{scenario}_{rep}_node_labels.csv"
    frame = pd.read_csv(mat_path)
    matrix = validate_square_matrix(frame.to_numpy(), "spMOCA")
    if label_path.exists():
        label_frame = pd.read_csv(label_path)
        if "gene_name" in label_frame.columns and len(label_frame) == matrix.shape[0]:
            labels = label_frame["gene_name"].astype(str).tolist()
        else:
            labels = infer_labels("gene_", matrix.shape[0])
    else:
        labels = infer_labels("gene_", matrix.shape[0])
    return MatrixBundle("spMOCA", matrix, labels, "estimated gene correlation", "gene")


def load_cssn(base_dir: Path, scenario: str, rep: int) -> MatrixBundle:
    mean_path = base_dir / "coexpression" / "CSSN" / f"{scenario}_{rep}_mean_coexpression_network.csv"
    path = base_dir / "coexpression" / "CSSN" / f"{scenario}_{rep}_CSSN_subset_results.rds"
    expr_path = base_dir / "Data" / f"{scenario}_gene" / "simulated" / f"{scenario}_gene_expression_{rep}.csv"
    coord_path = base_dir / "Data" / f"{scenario}_gene" / "simulated" / f"{scenario}_coord_{rep}.csv"
    expr_frame = pd.read_csv(expr_path)
    coord_frame = pd.read_csv(coord_path)
    n_spots = coord_frame.shape[0]
    n_genes = expr_frame.shape[1]

    if mean_path.exists():
        frame = pd.read_csv(mean_path)
        matrix = validate_square_matrix(frame.to_numpy(), "CSSN")
        labels = infer_labels("gene_", matrix.shape[0])
        description = "mean co-expression network"
        domain = "gene"
        return MatrixBundle("CSSN", matrix, labels, description, domain)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message='Missing constructor for R class "dsCMatrix".*',
            category=UserWarning,
        )
        obj = rdata.read_rds(path)
    corr_list = obj["Correlation Matrix"]
    matrices = [validate_square_matrix(np.asarray(mat), "CSSN") for mat in corr_list if mat is not None]
    if not matrices:
        raise ValueError(f"No CSSN correlation matrices found in {path}")
    if len(matrices) == 1:
        matrix = matrices[0]
        description = "subset correlation"
    else:
        matrix = np.mean(np.stack(matrices, axis=0), axis=0)
        description = f"mean of {len(matrices)} subset correlations"

    if matrix.shape[0] == n_genes:
        labels = infer_labels("gene_", matrix.shape[0])
        domain = "gene"
    elif matrix.shape[0] == n_spots:
        labels = infer_labels("spot_", matrix.shape[0])
        domain = "spot"
        description = f"{description}; saved CSSN output is spot-level"
    else:
        labels = infer_labels("item_", matrix.shape[0])
        domain = "unknown"
        description = f"{description}; unexpected size vs genes/spots"

    return MatrixBundle("CSSN", matrix, labels, description, domain)


def load_hdwgcna(base_dir: Path, scenario: str, rep: int) -> MatrixBundle:
    path = base_dir / "coexpression" / "hdWGCNA_Agent" / f"{scenario}_{rep}_modules.csv"
    frame = pd.read_csv(path)
    if not {"gene_name", "module"}.issubset(frame.columns):
        raise ValueError(f"Unexpected hdWGCNA columns in {path}: {list(frame.columns)}")

    labels = frame["gene_name"].astype(str).tolist()
    modules = frame["module"].astype(str).tolist()
    size = len(labels)
    matrix = np.zeros((size, size), dtype=float)
    np.fill_diagonal(matrix, 1.0)
    for i in range(size):
        for j in range(i + 1, size):
            same_module = modules[i] == modules[j] and modules[i].lower() != "grey"
            if same_module:
                matrix[i, j] = 1.0
                matrix[j, i] = 1.0
    return MatrixBundle("hdWGCNA", matrix, labels, "co-module adjacency", "gene")


def load_smoothie(base_dir: Path, scenario: str, rep: int) -> MatrixBundle:
    mat_path = base_dir / "coexpression" / "smoothie_Agent" / f"{scenario}_{rep}_correlation_matrix.npy"
    label_path = base_dir / "coexpression" / "smoothie_Agent" / f"{scenario}_{rep}_node_labels.csv"
    matrix = validate_square_matrix(np.load(mat_path), "Smoothie")
    label_frame = pd.read_csv(label_path)
    if "gene_name" in label_frame.columns and len(label_frame) == matrix.shape[0]:
        labels = label_frame["gene_name"].astype(str).tolist()
    else:
        labels = infer_labels("gene_", matrix.shape[0])
    return MatrixBundle("Smoothie", matrix, labels, "smoothed gene correlation", "gene")


def load_binspect(base_dir: Path, scenario: str, rep: int) -> MatrixBundle:
    path = base_dir / "coexpression" / "binSpect_Agent" / f"{scenario}_{rep}_correlation_matrix.csv"
    frame = pd.read_csv(path)
    matrix = validate_square_matrix(frame.to_numpy(), "binSpect")
    labels = infer_labels("gene_", matrix.shape[0])
    return MatrixBundle("binSpect", matrix, labels, "derived spatial correlation", "gene")


def load_stihc(base_dir: Path, scenario: str, rep: int) -> MatrixBundle:
    path = base_dir / "coexpression" / "stIHC_Agent" / f"{scenario}_gene_stIHC_{rep}.RDS"
    obj = rdata.read_rds(path)
    labels_arr = np.asarray(obj["label"]).astype(int)
    size = labels_arr.shape[0]
    matrix = np.equal.outer(labels_arr, labels_arr).astype(float)
    labels = infer_labels("gene_", size)
    return MatrixBundle("stIHC", matrix, labels, "co-cluster adjacency", "gene")


def load_ground_truth(base_dir: Path, scenario: str, rep: int) -> MatrixBundle:
    path = truth_file_path(base_dir, scenario)
    expr_path = base_dir / "Data" / f"{scenario}_gene" / "simulated" / f"{scenario}_gene_expression_{rep}.csv"
    frame = pd.read_csv(path)
    if "scenario" in frame.columns:
        rep_series = pd.to_numeric(frame["scenario"], errors="coerce")
        frame = frame.loc[rep_series == rep].copy()
    if frame.empty or not {"gene_id", "module_label"}.issubset(frame.columns):
        raise ValueError(f"No coexpression truth rows found for {scenario} rep {rep} in {path}")
    if not expr_path.exists():
        raise ValueError(f"Missing simulated expression file for {scenario} rep {rep}: {expr_path}")

    frame["gene_id"] = frame["gene_id"].astype(str)
    frame["module_label"] = frame["module_label"].astype(str)
    expr_frame = pd.read_csv(expr_path)
    full_labels = expr_frame.columns.astype(str).tolist()
    if "gene_index" in frame.columns:
        frame["_gene_rank"] = pd.to_numeric(frame["gene_index"], errors="coerce")
    else:
        frame["_gene_rank"] = pd.to_numeric(frame["gene_id"].str.extract(r"(\d+)$")[0], errors="coerce")
    order_cols = ["module_label", "_gene_rank", "gene_id"]
    if "module_id" in frame.columns:
        frame["_module_rank"] = pd.to_numeric(frame["module_id"], errors="coerce")
        order_cols = ["_module_rank", "module_label", "_gene_rank", "gene_id"]
    ordered_cols: list[str] = []
    for col in ["gene_id", "module_label"] + [col for col in order_cols if col in frame.columns]:
        if col not in ordered_cols:
            ordered_cols.append(col)
    ordered = frame[ordered_cols].drop_duplicates(subset=["gene_id"])
    ordered = ordered.sort_values(order_cols, kind="stable")
    truth_labels = ordered["gene_id"].astype(str).tolist()
    missing_labels = [label for label in full_labels if label not in set(truth_labels)]
    labels = truth_labels + missing_labels
    module_map = (
        frame[["gene_id", "module_label"]]
        .drop_duplicates(subset=["gene_id"])
        .set_index("gene_id")["module_label"]
        .to_dict()
    )
    size = len(labels)
    matrix = np.zeros((size, size), dtype=float)
    np.fill_diagonal(matrix, 1.0)
    for i in range(size):
        for j in range(i + 1, size):
            same_module = (
                labels[i] in module_map
                and labels[j] in module_map
                and module_map[labels[i]] == module_map[labels[j]]
            )
            if same_module:
                matrix[i, j] = 1.0
                matrix[j, i] = 1.0
    return MatrixBundle("Ground Truth", matrix, labels, "truth co-module adjacency with unassigned genes padded", "gene")


LOADERS = {
    "Ground Truth": load_ground_truth,
    "spMOCA": load_spmoca,
    "SpaceX": load_spacex,
    "CSSN": load_cssn,
    "hdWGCNA": load_hdwgcna,
    "Smoothie": load_smoothie,
    "binSpect": load_binspect,
    "stIHC": load_stihc,
}


def method_order_for_scenario(base_dir: Path, scenario: str, rep: int) -> list[str]:
    methods = list(BASE_METHOD_ORDER)
    if truth_available(base_dir, scenario, rep):
        methods.insert(0, "Ground Truth")
    return methods


def available_reps(base_dir: Path, scenario: str) -> list[int]:
    reps: list[int] = []
    for rep in range(1, 6):
        present = all(
            path_exists_for_method(base_dir, scenario, rep, method)
            for method in BASE_METHOD_ORDER
        )
        if present:
            reps.append(rep)
    return reps


def path_exists_for_method(base_dir: Path, scenario: str, rep: int, method: str) -> bool:
    if method == "SpaceX":
        return (base_dir / "Data" / f"{scenario}_gene" / f"SpaceX_{rep}.csv").exists()
    if method == "spMOCA":
        return (base_dir / "coexpression" / "spMOCA_Agent" / f"{scenario}_{rep}_correlation_matrix.csv").exists()
    if method == "CSSN":
        return (
            (base_dir / "coexpression" / "CSSN" / f"{scenario}_{rep}_mean_coexpression_network.csv").exists()
            or (base_dir / "coexpression" / "CSSN" / f"{scenario}_{rep}_CSSN_subset_results.rds").exists()
        )
    if method == "hdWGCNA":
        return (base_dir / "coexpression" / "hdWGCNA_Agent" / f"{scenario}_{rep}_modules.csv").exists()
    if method == "Smoothie":
        return (base_dir / "coexpression" / "smoothie_Agent" / f"{scenario}_{rep}_correlation_matrix.npy").exists()
    if method == "binSpect":
        return (base_dir / "coexpression" / "binSpect_Agent" / f"{scenario}_{rep}_correlation_matrix.csv").exists()
    if method == "stIHC":
        return (base_dir / "coexpression" / "stIHC_Agent" / f"{scenario}_gene_stIHC_{rep}.RDS").exists()
    return False


def reorder_to_reference(matrix: np.ndarray, labels: list[str], reference_labels: list[str]) -> tuple[np.ndarray, list[str]]:
    if len(labels) != len(reference_labels) or set(labels) != set(reference_labels):
        return matrix, labels
    index_map = {label: idx for idx, label in enumerate(labels)}
    order = [index_map[label] for label in reference_labels]
    return matrix[np.ix_(order, order)], reference_labels


def reorder_matrix(matrix: np.ndarray, labels: list[str], cluster_max_size: int) -> tuple[np.ndarray, list[str]]:
    size = matrix.shape[0]
    if size <= 2:
        return matrix, labels
    if size > cluster_max_size:
        return matrix, labels

    sym = (matrix + matrix.T) / 2.0
    np.fill_diagonal(sym, 1.0)
    clipped = np.clip(sym, -1.0, 1.0)
    dist = 1.0 - clipped
    np.fill_diagonal(dist, 0.0)
    condensed = squareform(dist, checks=False)
    linkage = hierarchy.linkage(condensed, method="average")
    order = hierarchy.leaves_list(linkage)
    reordered = matrix[np.ix_(order, order)]
    reordered_labels = [labels[idx] for idx in order]
    return reordered, reordered_labels


def style_axis(
    ax: plt.Axes,
    bundle: MatrixBundle,
    cluster_max_size: int,
    reference_labels: list[str] | None = None,
) -> None:
    if bundle.domain == "gene" and reference_labels is not None:
        ordered_matrix, _ = reorder_to_reference(bundle.matrix, bundle.labels, reference_labels)
    elif bundle.method == "Ground Truth":
        ordered_matrix = bundle.matrix
    else:
        ordered_matrix, _ = reorder_matrix(bundle.matrix, bundle.labels, cluster_max_size)
    size = ordered_matrix.shape[0]
    cmap = LinearSegmentedColormap.from_list("spdep_corr", ["#ffb6c1", "#ffffff", "#4B9CD3"])
    norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    im = ax.imshow(
        ordered_matrix,
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        aspect="equal",
        rasterized=size > 250,
    )
    ax.set_title(bundle.method, fontsize=12, fontweight="bold", pad=8)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    return im


def load_all_bundles(base_dir: Path, scenario: str, rep: int) -> list[MatrixBundle]:
    bundles: list[MatrixBundle] = []
    for method in method_order_for_scenario(base_dir, scenario, rep):
        bundles.append(LOADERS[method](base_dir, scenario, rep))
    return bundles


def plot_panel(base_dir: Path, scenario: str, rep: int, output_dir: Path, cluster_max_size: int) -> tuple[Path, Path, Path]:
    bundles = load_all_bundles(base_dir, scenario, rep)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_panels = len(bundles)
    ncols = 4 if n_panels >= 4 else n_panels
    nrows = math.ceil(n_panels / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.1 * ncols, 4.2 * nrows))
    axes = np.atleast_1d(axes).flatten()
    truth_bundle = next((bundle for bundle in bundles if bundle.method == "Ground Truth" and bundle.domain == "gene"), None)
    reference_labels = truth_bundle.labels if truth_bundle is not None else None
    images = []
    for ax, bundle in zip(axes, bundles):
        images.append(style_axis(ax, bundle, cluster_max_size, reference_labels))

    for ax in axes[len(bundles) :]:
        ax.axis("off")

    cbar_ax = fig.add_axes([0.935, 0.14, 0.018, 0.72])
    cbar = fig.colorbar(images[0], cax=cbar_ax)
    cbar.set_label("Correlation / adjacency", fontsize=10)
    cbar.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    fig.suptitle(f"{scenario} rep {rep} co-expression heatmaps", fontsize=16, fontweight="bold")
    fig.subplots_adjust(left=0.03, right=0.91, top=0.90, bottom=0.05, wspace=0.14, hspace=0.18)

    prefix = output_dir / f"{scenario}_rep{rep}_coexpression_corrplot_panel"
    png_path = prefix.with_suffix(".png")
    pdf_path = prefix.with_suffix(".pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    summary_rows = [
        {
            "method": bundle.method,
            "description": bundle.description,
            "domain": bundle.domain,
            "n_rows": bundle.matrix.shape[0],
            "n_cols": bundle.matrix.shape[1],
            "min_value": float(np.min(bundle.matrix)),
            "max_value": float(np.max(bundle.matrix)),
        }
        for bundle in bundles
    ]
    summary_path = output_dir / f"{scenario}_rep{rep}_coexpression_corrplot_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    return png_path, pdf_path, summary_path


def main() -> int:
    args = parse_args()
    base_dir = Path(__file__).resolve().parents[1]
    output_dir = Path(args.output_dir)

    all_outputs: list[tuple[str, int, Path, Path, Path]] = []
    for scenario in args.scenario:
        reps = sorted(set(args.rep)) if args.rep else available_reps(base_dir, scenario)
        if not reps:
            print(f"Skipping {scenario}: no replicates with all methods present.")
            continue

        for rep in reps:
            png_path, pdf_path, summary_path = plot_panel(
                base_dir=base_dir,
                scenario=scenario,
                rep=rep,
                output_dir=output_dir,
                cluster_max_size=args.cluster_max_size,
            )
            all_outputs.append((scenario, rep, png_path, pdf_path, summary_path))

    if not all_outputs:
        raise SystemExit("No corrplot-style co-expression panels were generated.")

    for scenario, rep, png_path, pdf_path, summary_path in all_outputs:
        print(f"{scenario} rep {rep}")
        print(f"  png: {png_path}")
        print(f"  pdf: {pdf_path}")
        print(f"  summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
