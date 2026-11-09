#!/usr/bin/env python3
"""
Pure image-driven spot clustering.

Given a histology image and a CSV of spot coordinates (e.g. Visium pixel space),
cluster each spot into K domains using only local image appearance. Outputs a
CSV with spot labels plus an overlay PNG matching the simulator-style
visualisation (colored disks on top of the slide).
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from matplotlib.tri import Triangulation
from skimage import color, exposure, filters, io, transform, util


def read_image(path: str, max_side: int | None = None) -> tuple[np.ndarray, float, float]:
    img = io.imread(path)
    if img.ndim == 2:
        img = np.stack([img] * 3, axis=-1)
    img = util.img_as_float32(img)
    h0, w0, _ = img.shape
    if max_side is not None and max(h0, w0) > max_side:
        scale = max_side / max(h0, w0)
        img_resized = transform.resize(img, (int(round(h0 * scale)), int(round(w0 * scale))), anti_aliasing=True)
        return img_resized, scale, scale
    return img, 1.0, 1.0


def extract_features(img: np.ndarray,
                     coords: np.ndarray,
                     radius: int,
                     white_threshold: float = 0.93) -> np.ndarray:
    h, w, _ = img.shape
    lab = color.rgb2lab(img)
    hed = color.rgb2hed(img)
    brightness = img.mean(axis=2)
    gradient = filters.sobel(filters.gaussian(brightness, sigma=1.0))

    feats = []
    rr = max(1, radius)
    white_mask = brightness < white_threshold
    for (x, y) in coords:
        cx = np.clip(int(round(x)), 0, w - 1)
        cy = np.clip(int(round(y)), 0, h - 1)
        xs = slice(max(0, cx - rr), min(w, cx + rr + 1))
        ys = slice(max(0, cy - rr), min(h, cy + rr + 1))
        mask = white_mask[ys, xs]
        if not np.any(mask):
            mask = np.ones_like(mask, dtype=bool)
        spot_L = lab[ys, xs, 0][mask]
        spot_a = lab[ys, xs, 1][mask]
        spot_b = lab[ys, xs, 2][mask]
        spot_H = (-hed[ys, xs, 0])[mask]
        spot_E = hed[ys, xs, 1][mask]
        spot_B = brightness[ys, xs][mask]
        spot_G = gradient[ys, xs][mask]

        def stats(arr: np.ndarray) -> tuple[float, float]:
            if arr.size == 0:
                return 0.0, 0.0
            return float(arr.mean()), float(arr.std())

        def frac(arr: np.ndarray, thr: float) -> float:
            if arr.size == 0:
                return 0.0
            return float(np.mean(arr > thr))

        feats.append([
            *stats(spot_L),
            *stats(spot_a),
            *stats(spot_b),
            *stats(spot_H),
            *stats(spot_E),
            *stats(spot_B),
            *stats(spot_G),
            frac(spot_B, 0.80),
            frac(spot_H, 0.55),
            frac(spot_E, 0.55),
        ])
    feats = np.asarray(feats, dtype=np.float32)
    # standardise
    mean = np.nanmean(feats, axis=0, keepdims=True)
    std = np.nanstd(feats, axis=0, keepdims=True)
    std[std == 0] = 1.0
    feats = (feats - mean) / std
    feats[np.isnan(feats)] = 0
    return feats


def smooth_labels(coords: np.ndarray,
                  labels: np.ndarray,
                  n_neighbors: int = 6,
                  iterations: int = 1) -> np.ndarray:
    if iterations <= 0 or n_neighbors <= 0 or len(labels) == 0:
        return labels
    n = len(labels)
    k_eff = min(n_neighbors + 1, n)
    nbrs = NearestNeighbors(n_neighbors=k_eff, metric="euclidean")
    nbrs.fit(coords)
    neighbor_idx = nbrs.kneighbors(return_distance=False)
    lab = labels.copy()
    for _ in range(iterations):
        new_lab = lab.copy()
        for i in range(n):
            neigh = neighbor_idx[i][1:]  # drop self
            if neigh.size == 0:
                continue
            counts = np.bincount(lab[neigh], minlength=lab.max() + 1)
            counts[0] = 0
            if counts.sum() == 0:
                continue
            new_lab[i] = counts.argmax()
        lab = new_lab
    return lab


def _palette(unique_labels: np.ndarray) -> dict[int, np.ndarray]:
    base_colors = np.array([
        [234, 99, 123],
        [90, 194, 120],
        [80, 130, 210],
        [240, 185, 70],
        [170, 120, 200],
        [120, 200, 200],
        [200, 160, 100],
    ], dtype=float) / 255.0
    palette = {lab: base_colors[i % len(base_colors)] for i, lab in enumerate(unique_labels)}
    return palette


def render_overlay(img: np.ndarray,
                   coords: np.ndarray,
                   labels: np.ndarray,
                   spot_radius: int,
                   out_path: str,
                   mode: str = "lines") -> None:
    h, w, _ = img.shape
    uniq = np.unique(labels)
    palette = _palette(uniq)

    fig = plt.figure(figsize=(w / 100, h / 100), dpi=100)
    ax = plt.axes([0, 0, 1, 1])
    ax.imshow(img, interpolation="bilinear")

    if mode == "spots":
        for (x, y), lab in zip(coords, labels):
            circ = plt.Circle((x, y), radius=spot_radius, facecolor=palette[lab], edgecolor="white", linewidth=0.6)
            ax.add_patch(circ)
    else:
        tri = Triangulation(coords[:, 0], coords[:, 1])
        label_list = uniq.tolist()
        label_indices = np.array([label_list.index(lab) for lab in labels])
        cmap = ListedColormap([palette[lab] for lab in label_list])
        levels = np.arange(-0.5, len(label_list) + 0.5, 1.0)
        ax.tricontourf(tri, label_indices, levels=levels, cmap=cmap, alpha=0.35, extend='neither')
        if len(label_list) > 1:
            ax.tricontour(tri, label_indices, levels=levels[1:-1], colors='white', linewidths=1.5)
        for lab in uniq:
            mask = labels == lab
            if not np.any(mask):
                continue
            ax.scatter(coords[mask, 0], coords[mask, 1], s=5, color=palette[lab], alpha=0.6)

    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")
    handles = [plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=palette[lab],
                          markersize=8, label=f"cluster {int(lab)}") for lab in uniq]
    fig.legend(handles=handles, loc="upper right", frameon=True)
    fig.savefig(out_path, dpi=100)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Spot-level image k-means clustering.")
    parser.add_argument("--image", required=True, help="Path to histology image.")
    parser.add_argument("--coords", required=True, help="CSV with x_coord,y_coord columns for spots.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix for CSV/PNG.")
    parser.add_argument("--k", type=int, default=3, help="Number of clusters.")
    parser.add_argument("--spot-radius", type=int, default=8, help="Spot radius (pixels) on processed image.")
    parser.add_argument("--patch-radius", type=int, default=10, help="Radius used for feature extraction.")
    parser.add_argument("--max-side", type=int, default=1024, help="Resize image so max side equals this value.")
    parser.add_argument("--white-threshold", type=float, default=0.93, help="Brightness threshold to mask background.")
    parser.add_argument("--seed", type=int, default=123, help="Random seed for k-means.")
    parser.add_argument("--smooth-iterations", type=int, default=1,
                        help="Number of neighbor-smoothing passes applied to labels (0 disables).")
    parser.add_argument("--smooth-k", type=int, default=6,
                        help="Number of nearest spots used for smoothing.")
    parser.add_argument("--overlay-mode", choices=("spots", "lines"), default="lines",
                        help="Visualisation style for the overlay PNG.")
    args = parser.parse_args()

    img, sx, sy = read_image(args.image, max_side=args.max_side)
    df = pd.read_csv(args.coords)
    if not {"x_coord", "y_coord"}.issubset(df.columns):
        raise ValueError("Coordinate CSV must contain x_coord and y_coord columns.")
    coords = df[["x_coord", "y_coord"]].to_numpy(dtype=float)
    h, w, _ = img.shape
    if coords[:, 0].max() > w + 1 or coords[:, 1].max() > h + 1:
        coords[:, 0] *= sx
        coords[:, 1] *= sy

    feats = extract_features(img, coords, radius=args.patch_radius, white_threshold=args.white_threshold)
    km = KMeans(n_clusters=args.k, random_state=args.seed, n_init=20)
    labels = km.fit_predict(feats) + 1  # 1-based
    labels = smooth_labels(coords, labels, n_neighbors=args.smooth_k, iterations=args.smooth_iterations)

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    df_out = df.copy()
    df_out["cluster"] = labels
    csv_path = f"{args.out_prefix}_clusters.csv"
    df_out.to_csv(csv_path, index=False)

    overlay_path = f"{args.out_prefix}_overlay.png"
    render_overlay(img, coords, labels, args.spot_radius, overlay_path, mode=args.overlay_mode)

    print(f"[output] spot clusters: {csv_path}")
    print(f"[output] overlay: {overlay_path}")


if __name__ == "__main__":
    main()
