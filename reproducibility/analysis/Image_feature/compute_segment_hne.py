#!/usr/bin/env python3
"""Reproduce the Squidpy ``compute_segment_hne`` example on the simulator image.

The Squidpy tutorial smooths an H&E crop (sigma=4) and thresholds channel 0 at 90,
then runs watershed with ``geq=False`` so that darker nuclei are foreground. This
script mirrors those steps without the full Squidpy dependency.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage import filters, morphology, segmentation

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None


DEFAULT_IMAGE = Path(
    "./data_external/HER2/9769_C1_HE_small.jpg"
)


def load_rgb(path: Path) -> np.ndarray:
    """Load an RGB image and return an (H,W,3) uint8 array."""
    img = Image.open(path)
    img = img.convert("RGB")
    return np.array(img)


def crop_image(img: np.ndarray, crop_size: int | None) -> np.ndarray:
    if crop_size is None:
        return img
    h, w = img.shape[:2]
    size = min(crop_size, h, w)
    return img[:size, :size]


def smooth_image(img: np.ndarray, sigma: float) -> np.ndarray:
    """Apply Gaussian smoothing to preserve original scale (0-255)."""
    return filters.gaussian(img, sigma=sigma, preserve_range=True, channel_axis=-1)


def surf_segmentation(
    img: np.ndarray,
    sigma: float = 4.0,
    thresh: float = 90.0,
    crop_size: int | None = None,
    geq: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Replicate the core Squidpy pipeline: smooth, threshold channel 0, and watershed.
    """
    crop = crop_image(img, crop_size)
    smooth = smooth_image(crop, sigma)
    channel0 = smooth[..., 0]
    mask = channel0 >= thresh if geq else channel0 < thresh
    mask = morphology.remove_small_objects(mask, min_size=16)
    distance = ndi.distance_transform_edt(mask)
    local_max = morphology.local_maxima(distance)
    markers, _ = ndi.label(local_max)
    if markers.max() == 0:
        markers, _ = ndi.label(mask)
    gradient = filters.sobel(channel0)
    labels = segmentation.watershed(gradient, markers, mask=mask)
    return crop.astype(np.uint8), smooth.astype(np.float32), labels


def overlay_labels(image: np.ndarray, labels: np.ndarray) -> np.ndarray:
    boundaries = segmentation.find_boundaries(labels, connectivity=1, mode="inner")
    overlay = image.copy()
    overlay[boundaries] = [255, 0, 0]
    return overlay


def write_outputs(
    output_dir: Path,
    crop: np.ndarray,
    smooth: np.ndarray,
    labels: np.ndarray,
    overlay: np.ndarray,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    Image.fromarray(crop).save(output_dir / "crop.png")
    Image.fromarray(np.clip(smooth[..., 0], 0, 255).astype(np.uint8)).save(
        output_dir / "smoothed_channel0.png"
    )
    Image.fromarray(overlay).save(output_dir / "segmentation_overlay.png")
    np.savez(
        output_dir / "labels.npz",
        labels=labels,
        num_segments=int(np.max(labels)),
        unique=int(len(np.unique(labels)) - (1 if np.any(labels == 0) else 0)),
    )


def plot_summary(output_dir: Path, crop: np.ndarray, labels: np.ndarray) -> bool:
    if plt is None:
        return False
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    axes[0].imshow(crop)
    axes[0].set_title("Cropped H&E")
    axes[0].axis("off")
    axes[1].imshow(labels, cmap="tab20")
    axes[1].set_title("Watershed segmentation")
    axes[1].axis("off")
    plt.tight_layout()
    fig.savefig(output_dir / "segmentation_summary.png", dpi=150)
    plt.close(fig)
    return True


def export_segment_metadata(output_dir: Path, labels: np.ndarray) -> None:
    rows: list[tuple[int, float, float, int]] = []
    for seg in np.unique(labels):
        if seg == 0:
            continue
        coords = np.nonzero(labels == seg)
        if coords[0].size == 0:
            continue
        y_mean = float(coords[0].mean())
        x_mean = float(coords[1].mean())
        rows.append((int(seg), x_mean, y_mean, int(coords[0].size)))
    if not rows:
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "segments.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["segment_id", "x", "y", "pixel_count"])
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the Squidpy compute_segment_hne example.")
    parser.add_argument(
        "--image-path",
        type=Path,
        default=DEFAULT_IMAGE,
        help="H&E image to segment (defaults to the simulator example).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("STSimulator/Image_feature/segment_hne_output"),
        help="Directory for summary images/labels.",
    )
    parser.add_argument(
        "--crop-size",
        type=int,
        default=None,
        help="Crop size applied to the image; omit this flag to process the full image.",
    )
    parser.add_argument("--sigma", type=float, default=4.0, help="Gaussian smoothing sigma.")
    parser.add_argument("--threshold", type=float, default=90.0, help="Channel-0 threshold.")
    parser.add_argument(
        "--geq",
        action="store_true",
        help="If set, treat values >= threshold as foreground (not the default).",
    )
    args = parser.parse_args()

    image = load_rgb(args.image_path)
    crop, smooth, labels = surf_segmentation(
        image, sigma=args.sigma, thresh=args.threshold, crop_size=args.crop_size, geq=args.geq
    )
    overlay = overlay_labels(crop, labels)
    write_outputs(args.output_dir, crop, smooth, labels, overlay)
    summary_written = plot_summary(args.output_dir, crop, labels)
    export_segment_metadata(args.output_dir, labels)

    num_segments = len(np.unique(labels)) - (1 if np.any(labels == 0) else 0)
    print(f"Computed watershed segmentation with {num_segments} unique segments.")
    print(f"Outputs written to: {args.output_dir}")
    if not summary_written:
        print("Optional matplotlib dependency not found; skipped segmentation_summary.png.")


if __name__ == "__main__":
    main()
