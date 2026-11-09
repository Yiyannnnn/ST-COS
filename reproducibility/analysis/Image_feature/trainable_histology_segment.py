#!/usr/bin/env python3
"""
Learned histology segmentation driven by spot annotations.

Given:
  * an H&E slide image (JPEG/PNG/TIFF),
  * a CSV of Visium/spot coordinates with human-provided labels (layer_mapping),
this script learns a simple RandomForest classifier on spot-level image features
and then predicts a dense superpixel segmentation that mimics the provided
labels as closely as possible.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from skimage import color, exposure, filters, io, morphology, segmentation, transform, util


# Canonical colour palette keyed by lower-case class name
DEFAULT_CLASS_COLORS: Dict[str, Tuple[int, int, int]] = {
    "invasive cancer": (220, 20, 60),
    "breast glands": (0, 128, 0),
    "cancer in situ": (255, 140, 0),
    "connective tissue": (72, 61, 139),
    "immune infiltrate": (30, 144, 255),
    "adipose tissue": (255, 215, 0),
    "undetermined": (169, 169, 169),
}


@dataclass
class ImageBundle:
    rgb: np.ndarray          # float32, range [0,1]
    lab: np.ndarray          # same shape, Lab space
    hed: np.ndarray          # HED colour deconvolution
    brightness: np.ndarray   # grayscale brightness
    gradient: np.ndarray     # Sobel gradient
    tissue_mask: np.ndarray  # boolean mask
    scale_x: float           # scale factor original -> working width
    scale_y: float           # scale factor original -> working height


def read_and_prepare_image(path: str,
                           white_threshold: float = 0.92,
                           max_side: int = 768,
                           gaussian_sigma: float = 1.0) -> ImageBundle:
    raw = io.imread(path)
    if raw.ndim == 2:
        raw = np.stack([raw] * 3, axis=-1)
    raw = util.img_as_float32(raw)
    h0, w0, _ = raw.shape

    scale = 1.0
    if max(h0, w0) > max_side:
        scale = max_side / float(max(h0, w0))
        rgb = transform.resize(raw, (int(round(h0 * scale)), int(round(w0 * scale))), anti_aliasing=True)
    else:
        rgb = raw.copy()
    rgb = np.clip(rgb, 0, 1)
    h, w, _ = rgb.shape
    scale_x = w / float(w0)
    scale_y = h / float(h0)

    lab = color.rgb2lab(rgb)
    hed = color.rgb2hed(rgb)
    brightness = rgb.mean(axis=2)
    tissue_mask = brightness < white_threshold
    tissue_mask = morphology.remove_small_objects(tissue_mask, min_size=128)
    tissue_mask = morphology.remove_small_holes(tissue_mask, area_threshold=128)
    blur = filters.gaussian(brightness, sigma=gaussian_sigma)
    gradient = filters.sobel(blur)
    gradient = exposure.rescale_intensity(gradient, out_range=(0, 1))

    return ImageBundle(
        rgb=rgb,
        lab=lab,
        hed=hed,
        brightness=brightness,
        gradient=gradient,
        tissue_mask=tissue_mask,
        scale_x=scale_x,
        scale_y=scale_y,
    )


def safe_region(values: np.ndarray) -> Tuple[float, float]:
    if values.size == 0:
        return 0.0, 0.0
    return float(values.mean()), float(values.std())


def fraction_above(values: np.ndarray, threshold: float) -> float:
    if values.size == 0:
        return 0.0
    return float(np.mean(values > threshold))


def extract_circle_mask(cx: float, cy: float, radius: int, shape: Tuple[int, int]) -> np.ndarray:
    h, w = shape
    x_min = max(0, int(np.floor(cx - radius)))
    x_max = min(w - 1, int(np.ceil(cx + radius)))
    y_min = max(0, int(np.floor(cy - radius)))
    y_max = min(h - 1, int(np.ceil(cy + radius)))
    yy = np.arange(y_min, y_max + 1)[:, None]
    xx = np.arange(x_min, x_max + 1)[None, :]
    mask_local = ((xx - cx) ** 2 + (yy - cy) ** 2) <= radius ** 2
    mask = np.zeros(shape, dtype=bool)
    mask[y_min:y_max + 1, x_min:x_max + 1] = mask_local
    return mask


def compute_feature_vector(bundle: ImageBundle, mask: np.ndarray) -> np.ndarray:
    mask = mask & bundle.tissue_mask
    if not mask.any():
        mask = mask | mask  # keep empty, results will be zeros
    lab_L = bundle.lab[:, :, 0][mask]
    lab_a = bundle.lab[:, :, 1][mask]
    lab_b = bundle.lab[:, :, 2][mask]
    hed_H = (-bundle.hed[:, :, 0])[mask]  # nuclei
    hed_E = bundle.hed[:, :, 1][mask]
    brightness = bundle.brightness[mask]
    gradient = bundle.gradient[mask]

    feats = [
        *safe_region(lab_L),
        *safe_region(lab_a),
        *safe_region(lab_b),
        *safe_region(hed_H),
        *safe_region(hed_E),
        *safe_region(brightness),
        *safe_region(gradient),
        fraction_above(brightness, 0.82),
        fraction_above(hed_H, 0.55),
        fraction_above(hed_E, 0.58),
    ]
    return np.array(feats, dtype=np.float32)


def build_spot_features(bundle: ImageBundle,
                        coords: np.ndarray,
                        radius: int) -> np.ndarray:
    h, w = bundle.rgb.shape[:2]
    feats = []
    for (x, y) in coords:
        mask = extract_circle_mask(x, y, radius, (h, w))
        feats.append(compute_feature_vector(bundle, mask))
    return np.vstack(feats)


def build_spot_mask(coords: np.ndarray,
                    shape: Tuple[int, int],
                    radius: int,
                    dilation: int = 4) -> np.ndarray:
    h, w = shape
    mask = np.zeros((h, w), dtype=bool)
    for x, y in coords:
        mask |= extract_circle_mask(x, y, radius, (h, w))
    if dilation > 0:
        mask = morphology.binary_dilation(mask, morphology.disk(max(1, dilation)))
    mask = morphology.remove_small_holes(mask, area_threshold=int(radius ** 2))
    return mask

def compute_superpixel_features(bundle: ImageBundle,
                                segments: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    labels = np.unique(segments[bundle.tissue_mask])
    feats = []
    for lab in labels:
        mask = segments == lab
        feats.append(compute_feature_vector(bundle, mask))
    return labels, np.vstack(feats)


def ensure_directory(path: str) -> None:
    if path and not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def pick_color(name: str) -> Tuple[int, int, int]:
    key = name.lower()
    if key in DEFAULT_CLASS_COLORS:
        return DEFAULT_CLASS_COLORS[key]
    return (200, 200, 200)


def render_overlay(bundle: ImageBundle,
                   segments: np.ndarray,
                   segment_labels: np.ndarray,
                   class_names: Iterable[str],
                   alpha: float = 0.45) -> np.ndarray:
    solid = np.zeros_like(bundle.rgb, dtype=np.float32)
    for lab, name in zip(segment_labels, class_names):
        color_rgb = np.array(pick_color(name), dtype=np.float32) / 255.0
        solid[segments == lab] = color_rgb
    overlay = (alpha * solid + (1 - alpha) * bundle.rgb).clip(0, 1)
    # Highlight boundaries
    boundaries = segmentation.find_boundaries(segments, connectivity=1, mode="inner")
    overlay[boundaries] = 1.0
    return util.img_as_ubyte(overlay)


def main():
    parser = argparse.ArgumentParser(description="Train a spot-informed histology segmentation.")
    parser.add_argument("--image", required=True, help="Path to H&E image.")
    parser.add_argument("--labels", required=True, help="CSV with spot coordinates and labels.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix for artifacts.")
    parser.add_argument("--x-col", default="x_coord", help="Column with x coordinates (pixels).")
    parser.add_argument("--y-col", default="y_coord", help="Column with y coordinates (pixels).")
    parser.add_argument("--label-col", default="layer_mapping", help="Column with class labels.")
    parser.add_argument("--spot-radius", type=int, default=8, help="Spot radius in pixels on the downscaled image.")
    parser.add_argument("--white-threshold", type=float, default=0.94, help="Brightness cutoff for tissue mask.")
    parser.add_argument("--max-side", type=int, default=768, help="Largest side for downscaled processing image.")
    parser.add_argument("--gaussian-sigma", type=float, default=1.0, help="Gaussian smoothing before gradients.")
    parser.add_argument("--slic-segments", type=int, default=1600, help="Target number of SLIC superpixels.")
    parser.add_argument("--slic-compactness", type=float, default=18.0, help="SLIC compactness parameter.")
    parser.add_argument("--seed", type=int, default=123, help="Random seed for RandomForest and SLIC.")
    parser.add_argument("--restrict-to-spots", action="store_true",
                        help="Limit clustering to areas covered by spot coordinates.")
    parser.add_argument("--spot-mask-radius", type=int, default=15,
                        help="Radius (pixels) used to build spot coverage mask when restricting.")
    args = parser.parse_args()

    bundle = read_and_prepare_image(
        path=args.image,
        white_threshold=args.white_threshold,
        max_side=args.max_side,
        gaussian_sigma=args.gaussian_sigma,
    )
    h, w = bundle.rgb.shape[:2]

    df = pd.read_csv(args.labels)
    if args.x_col not in df.columns or args.y_col not in df.columns or args.label_col not in df.columns:
        raise ValueError("Labels file must contain columns "
                         f"{args.x_col}, {args.y_col}, {args.label_col}. Columns present: {list(df.columns)}")
    coords = df[[args.x_col, args.y_col]].to_numpy(dtype=float)
    # Auto-scale coordinates if they appear to be in original resolution
    if coords[:, 0].max() > w + 1 or coords[:, 1].max() > h + 1:
        coords[:, 0] *= bundle.scale_x
        coords[:, 1] *= bundle.scale_y

    if args.restrict_to_spots:
        spot_mask = build_spot_mask(coords, (h, w), radius=max(args.spot_mask_radius, args.spot_radius))
        bundle.tissue_mask &= spot_mask

    X_spots = build_spot_features(bundle, coords, radius=args.spot_radius)
    label_encoder = LabelEncoder()
    y_spots = label_encoder.fit_transform(df[args.label_col].astype(str))

    clf = RandomForestClassifier(
        n_estimators=300,
        max_depth=None,
        min_samples_leaf=3,
        class_weight="balanced",
        random_state=args.seed,
        n_jobs=-1,
    )
    clf.fit(X_spots, y_spots)

    # Evaluate on training spots for sanity
    spot_preds = clf.predict(X_spots)
    spot_acc = (spot_preds == y_spots).mean()
    print(f"[train] Spot classification accuracy: {spot_acc:.3f}")

    # Build SLIC superpixels and features
    segments = segmentation.slic(
        util.img_as_float(bundle.rgb),
        n_segments=args.slic_segments,
        compactness=args.slic_compactness,
        enforce_connectivity=True,
        start_label=0,
        mask=bundle.tissue_mask,
        convert2lab=False,
    )
    segment_ids, X_segments = compute_superpixel_features(bundle, segments)
    segment_preds = clf.predict(X_segments)
    segment_probs = clf.predict_proba(X_segments)
    class_names = label_encoder.inverse_transform(segment_preds)

    # Override segment labels where training spots provide direct evidence
    for idx, (x, y, label_name) in enumerate(zip(coords[:, 0], coords[:, 1], df[args.label_col].astype(str))):
        seg_label = segments[int(round(y)), int(round(x))]
        match_idx = np.where(segment_ids == seg_label)[0]
        if match_idx.size:
            segment_preds[match_idx[0]] = label_encoder.transform([label_name])[0]
            class_names[match_idx[0]] = label_name

    overlay = render_overlay(bundle, segments, segment_ids, class_names)

    out_dir = os.path.dirname(args.out_prefix)
    ensure_directory(out_dir)

    overlay_path = f"{args.out_prefix}_overlay.png"
    label_map_path = f"{args.out_prefix}_labels.npy"
    segment_csv_path = f"{args.out_prefix}_segments.csv"
    spot_csv_path = f"{args.out_prefix}_spots.csv"

    io.imsave(overlay_path, overlay)
    np.save(label_map_path, {
        "segment_ids": segment_ids,
        "segment_preds": segment_preds,
        "class_names": class_names,
        "segments_array": segments,
    }, allow_pickle=True)

    segment_df = pd.DataFrame({
        "segment_id": segment_ids,
        "predicted_class": class_names,
        "predicted_index": segment_preds,
        "probability": segment_probs.max(axis=1),
    })
    segment_df.to_csv(segment_csv_path, index=False)

    spot_df = pd.DataFrame({
        "spot_id": df.index,
        "x": coords[:, 0],
        "y": coords[:, 1],
        "label": df[args.label_col],
        "predicted": label_encoder.inverse_transform(spot_preds),
    })
    spot_df.to_csv(spot_csv_path, index=False)

    print(f"[output] overlay: {overlay_path}")
    print(f"[output] segment predictions: {segment_csv_path}")
    print(f"[output] per-spot predictions: {spot_csv_path}")


if __name__ == "__main__":
    main()
