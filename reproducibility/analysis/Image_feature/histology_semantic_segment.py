#!/usr/bin/env python3
from __future__ import annotations

"""
Rule-based semantic segmentation for H&E breast tumor slides.

The pipeline:
1. Read an RGB slide image (supports JPEG/PNG/TIFF).
2. Compute a tissue mask, color-deconvolve to Hematoxylin/Eosin space,
   derive gradient/texture statistics, and run SLIC superpixels.
3. For every superpixel collect heuristically chosen features (nuclear density,
   eosin affinity, bright lumen fraction, etc.).
4. Assign one of seven histological categories via simple rules tuned for
   HER2+ breast cancer morphology:
      1: Invasive cancer (nuclear-dense, solid)
      2: Breast glands (moderate lumen, eosin-rich)
      3: Cancer in situ (large lumen surrounded by nuclei)
      4: Connective tissue (stromal eosin signal)
      5: Immune infiltrate (hematoxylin rich, low eosin)
      6: Adipose tissue (large bright vacuoles)
      7: Undetermined / fallback
   Background (outside tissue mask) is labelled 0.
5. Save a segmentation PNG overlay and a CSV summarising superpixel statistics.

This is not a replacement for a trained deep model, but produces far more
structured maps than pure colour clustering and is tuned for slides similar to
the user-provided example.
"""

import argparse
import csv
import os
from dataclasses import dataclass

import numpy as np
from skimage import color, exposure, filters, io, morphology, segmentation, util


@dataclass
class RegionStats:
    label: int
    area: int
    mean_brightness: float
    lumen_fraction: float
    hematoxylin_mean: float
    hematoxylin_std: float
    hematoxylin_frac: float
    eosin_mean: float
    eosin_std: float
    gradient_mean: float
    gradient_std: float
    class_id: int


CLASS_INFO = {
    0: ("Background", (0, 0, 0)),
    1: ("Invasive cancer", (220, 20, 60)),       # crimson
    2: ("Breast glands", (0, 128, 0)),           # green
    3: ("Cancer in situ", (255, 140, 0)),        # orange
    4: ("Connective tissue", (72, 61, 139)),     # slate blue
    5: ("Immune infiltrate", (30, 144, 255)),    # dodger blue
    6: ("Adipose tissue", (255, 215, 0)),        # gold
    7: ("Undetermined", (169, 169, 169)),        # dark grey
}


def read_image(path: str) -> np.ndarray:
    img = io.imread(path)
    if img.ndim == 2:  # grayscale
        img = np.stack([img] * 3, axis=-1)
    if img.dtype != np.float32 and img.dtype != np.float64:
        img = util.img_as_float32(img)
    return img


def compute_tissue_mask(img: np.ndarray, white_threshold: float = 0.90) -> np.ndarray:
    brightness = img.mean(axis=2)
    tissue = brightness < white_threshold
    tissue = morphology.remove_small_objects(tissue, min_size=256)
    tissue = morphology.remove_small_holes(tissue, area_threshold=256)
    return tissue


def compute_features(img: np.ndarray, tissue_mask: np.ndarray,
                     n_segments: int, compactness: float, sigma: float):
    hed = color.rgb2hed(img)
    hematoxylin_raw = -hed[:, :, 0]  # higher = more nuclei
    eosin_raw = hed[:, :, 1]

    hematoxylin = exposure.rescale_intensity(hematoxylin_raw, in_range='image', out_range=(0, 1))
    eosin = exposure.rescale_intensity(eosin_raw, in_range='image', out_range=(0, 1))
    brightness = img.mean(axis=2)

    gradient = filters.gaussian(brightness, sigma=sigma)
    gradient = filters.sobel(gradient)
    gradient = exposure.rescale_intensity(gradient, out_range=(0, 1))

    segments = segmentation.slic(
        img,
        n_segments=n_segments,
        compactness=compactness,
        sigma=sigma,
        start_label=0,
        mask=tissue_mask,
        enforce_connectivity=True,
    )

    unique_labels = np.unique(segments[tissue_mask])
    stats = []
    luminance_vals = []
    lumen_vals = []
    h_mean_vals = []
    h_frac_vals = []
    e_mean_vals = []
    grad_mean_vals = []
    for lab in unique_labels:
        mask = segments == lab
        area = int(mask.sum())
        mean_brightness = float(brightness[mask].mean())
        lumen_fraction = float(np.mean(brightness[mask] > 0.82))
        h_vals = hematoxylin[mask]
        e_vals = eosin[mask]
        grad_vals = gradient[mask]
        h_mean = float(h_vals.mean())
        h_std = float(h_vals.std())
        h_frac = float(np.mean(h_vals > 0.55))
        e_mean = float(e_vals.mean())
        e_std = float(e_vals.std())
        grad_mean = float(grad_vals.mean())
        grad_std = float(grad_vals.std())

        region = RegionStats(
            label=int(lab),
            area=area,
            mean_brightness=mean_brightness,
            lumen_fraction=lumen_fraction,
            hematoxylin_mean=h_mean,
            hematoxylin_std=h_std,
            hematoxylin_frac=h_frac,
            eosin_mean=e_mean,
            eosin_std=e_std,
            gradient_mean=grad_mean,
            gradient_std=grad_std,
            class_id=0,  # placeholder
        )
        stats.append(region)
        luminance_vals.append(mean_brightness)
        lumen_vals.append(lumen_fraction)
        h_mean_vals.append(h_mean)
        h_frac_vals.append(h_frac)
        e_mean_vals.append(e_mean)
        grad_mean_vals.append(grad_mean)

    return segments, stats, {
        "hematoxylin": hematoxylin,
        "eosin": eosin,
        "brightness": brightness,
        "gradient": gradient,
        "brightness_list": np.array(luminance_vals),
        "lumen_list": np.array(lumen_vals),
        "hmean_list": np.array(h_mean_vals),
        "hfrac_list": np.array(h_frac_vals),
        "emean_list": np.array(e_mean_vals),
        "gmean_list": np.array(grad_mean_vals),
    }


def derive_thresholds(features: dict[str, np.ndarray]) -> dict[str, float]:
    def quant(values: np.ndarray, q: float, default: float) -> float:
        vals = values[~np.isnan(values)]
        if vals.size == 0:
            return default
        return float(np.quantile(vals, q))

    thresholds = {
        "bright_high": quant(features["brightness_list"], 0.90, 0.92),
        "bright_low": quant(features["brightness_list"], 0.30, 0.75),
        "lumen_high": quant(features["lumen_list"], 0.85, 0.55),
        "lumen_mid": quant(features["lumen_list"], 0.70, 0.35),
        "lumen_low": quant(features["lumen_list"], 0.35, 0.18),
        "nuc_mean_high": quant(features["hmean_list"], 0.80, 0.55),
        "nuc_mean_low": quant(features["hmean_list"], 0.20, 0.25),
        "nuc_frac_high": quant(features["hfrac_list"], 0.80, 0.45),
        "nuc_frac_mid": quant(features["hfrac_list"], 0.60, 0.35),
        "nuc_frac_low": quant(features["hfrac_list"], 0.25, 0.15),
        "eosin_high": quant(features["emean_list"], 0.80, 0.55),
        "eosin_low": quant(features["emean_list"], 0.20, 0.35),
        "grad_high": quant(features["gmean_list"], 0.80, 0.08),
    }
    return thresholds


def assign_category(stat: RegionStats, thresholds: dict[str, float]) -> int:
    # Background
    if stat.mean_brightness > thresholds["bright_high"] and stat.hematoxylin_mean < thresholds["nuc_mean_low"]:
        return 0
    # Adipose: bright with large vacuoles
    if stat.lumen_fraction > thresholds["lumen_high"] and stat.hematoxylin_frac < thresholds["nuc_frac_low"]:
        return 6
    # Cancer in situ: big lumen surrounded by nuclei
    if (stat.lumen_fraction > thresholds["lumen_mid"] and
            stat.hematoxylin_frac > thresholds["nuc_frac_mid"] and
            stat.gradient_mean > thresholds["grad_high"]):
        return 3
    # Breast glands: moderate lumen, eosin strong
    if (thresholds["lumen_low"] < stat.lumen_fraction <= thresholds["lumen_mid"] and
            stat.eosin_mean > thresholds["eosin_high"]):
        return 2
    # Invasive cancer: dense nuclei, solid
    if (stat.hematoxylin_mean > thresholds["nuc_mean_high"] and
            stat.hematoxylin_frac > thresholds["nuc_frac_mid"] and
            stat.lumen_fraction < thresholds["lumen_low"]):
        return 1
    # Immune infiltrate: high nuclei, low eosin
    if stat.hematoxylin_frac > thresholds["nuc_frac_mid"] and stat.eosin_mean < thresholds["eosin_low"]:
        return 5
    # Connective tissue: eosin high, nuclei sparse
    if stat.eosin_mean > thresholds["eosin_high"] and stat.hematoxylin_frac < thresholds["nuc_frac_mid"]:
        return 4
    # Default fallback
    return 7


def render_segmentation(img: np.ndarray, segments: np.ndarray, stats: list[RegionStats], tissue_mask: np.ndarray) -> np.ndarray:
    seg_map = np.zeros(segments.shape, dtype=np.uint8)
    label_to_class = {stat.label: stat.class_id for stat in stats}
    for lab, cls in label_to_class.items():
        seg_map[segments == lab] = cls
    seg_map[~tissue_mask] = 0

    overlay = (img * 255).astype(np.uint8).copy()
    boundaries = segmentation.find_boundaries(seg_map, connectivity=1, mode="inner")
    overlay[boundaries] = [255, 255, 255]

    solid = np.zeros_like(overlay)
    for cls_id, (_, color_rgb) in CLASS_INFO.items():
        mask = seg_map == cls_id
        solid[mask] = color_rgb

    alpha = 0.45
    blended = ((1 - alpha) * overlay + alpha * solid).astype(np.uint8)
    return blended, seg_map


def write_csv(path: str, stats: list[RegionStats]) -> None:
    fieldnames = [
        "superpixel_id",
        "class_id",
        "class_name",
        "area_px",
        "mean_brightness",
        "lumen_fraction",
        "hematoxylin_mean",
        "hematoxylin_std",
        "hematoxylin_fraction",
        "eosin_mean",
        "eosin_std",
        "gradient_mean",
        "gradient_std",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for stat in stats:
            class_name = CLASS_INFO.get(stat.class_id, ("Unknown", None))[0]
            writer.writerow({
                "superpixel_id": stat.label,
                "class_id": stat.class_id,
                "class_name": class_name,
                "area_px": stat.area,
                "mean_brightness": stat.mean_brightness,
                "lumen_fraction": stat.lumen_fraction,
                "hematoxylin_mean": stat.hematoxylin_mean,
                "hematoxylin_std": stat.hematoxylin_std,
                "hematoxylin_fraction": stat.hematoxylin_frac,
                "eosin_mean": stat.eosin_mean,
                "eosin_std": stat.eosin_std,
                "gradient_mean": stat.gradient_mean,
                "gradient_std": stat.gradient_std,
            })


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Heuristic semantic segmentation of H&E breast cancer slides.")
    parser.add_argument("--image", required=True, help="Path to the input slide image (PNG/JPG/TIFF).")
    parser.add_argument("--out-prefix", required=True, help="Output file prefix.")
    parser.add_argument("--slic-segments", type=int, default=1600, help="Approximate number of SLIC superpixels.")
    parser.add_argument("--slic-compactness", type=float, default=18.0, help="SLIC compactness parameter.")
    parser.add_argument("--gaussian-sigma", type=float, default=1.0, help="Gaussian smoothing prior to feature extraction.")
    parser.add_argument("--white-threshold", type=float, default=0.92, help="Brightness threshold for tissue masking.")
    return parser.parse_args()


def main():
    args = parse_args()
    img = read_image(args.image)
    tissue_mask = compute_tissue_mask(img, white_threshold=args.white_threshold)
    segments, stats, feat = compute_features(
        img,
        tissue_mask,
        n_segments=args.slic_segments,
        compactness=args.slic_compactness,
        sigma=args.gaussian_sigma,
    )
    thresholds = derive_thresholds(feat)
    enriched_stats = []
    for stat in stats:
        cls = assign_category(stat, thresholds)
        enriched_stats.append(stat.__class__(**{**stat.__dict__, "class_id": cls}))

    overlay, seg_map = render_segmentation(img, segments, enriched_stats, tissue_mask)

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    overlay_path = f"{args.out_prefix}_semantic_overlay.png"
    raw_path = f"{args.out_prefix}_semantic_labels.npy"
    csv_path = f"{args.out_prefix}_superpixels.csv"

    io.imsave(overlay_path, overlay)
    np.save(raw_path, seg_map.astype(np.uint8))
    write_csv(csv_path, enriched_stats)

    print(f"[segmentation] overlay saved to {overlay_path}")
    print(f"[segmentation] label map saved to {raw_path}")
    print(f"[segmentation] superpixel stats saved to {csv_path}")


if __name__ == "__main__":
    main()
