#!/usr/bin/env python3
"""Overlay segment centroids (from `segments.csv`) onto the H&E crop."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from PIL import Image, ImageDraw


def read_segments(path: Path) -> list[tuple[int, float, float, int]]:
    rows: list[tuple[int, float, float, int]] = []
    with path.open() as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split(",")
            if len(parts) < 4:
                continue
            try:
                seg_id = int(parts[0])
                x = float(parts[1])
                y = float(parts[2])
                count = int(parts[3])
            except ValueError:
                continue
            rows.append((seg_id, x, y, count))
    return rows


def draw_segments(
    image_path: Path,
    segments: list[tuple[int, float, float, int]],
    output_path: Path,
    radius: float = 6.0,
    size_scale: float = 0.0,
    color: tuple[int, int, int, int] = (255, 0, 0, 200),
    annotate: bool = False,
    max_segments: int | None = None,
) -> None:
    image = Image.open(image_path).convert("RGBA")
    overlay = Image.new("RGBA", image.size)
    draw = ImageDraw.Draw(overlay)
    to_draw = segments if max_segments is None else segments[:max_segments]
    for seg_id, x, y, count in to_draw:
        scaled_radius = radius + size_scale * math.sqrt(count)
        bbox = (x - scaled_radius, y - scaled_radius, x + scaled_radius, y + scaled_radius)
        draw.ellipse(bbox, outline=color, width=2)
        if annotate:
            text = str(seg_id)
            draw.text((x + scaled_radius, y), text, fill=color)
    combined = Image.alpha_composite(image, overlay)
    combined.save(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Draw segment centroids on the H&E crop.")
    parser.add_argument(
        "--segments",
        type=Path,
        default=Path("STSimulator/Image_feature/segment_hne_output/segments.csv"),
        help="CSV produced by compute_segment_hne.py",
    )
    parser.add_argument(
        "--image",
        type=Path,
        default=Path("STSimulator/Image_feature/segment_hne_output/crop.png"),
        help="Base H&E image (crop or full slide).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("STSimulator/Image_feature/segment_hne_output/crop_segments.png"),
        help="Overlay image to write.",
    )
    parser.add_argument("--radius", type=float, default=6.0, help="Base radius (pixels) of the spots.")
    parser.add_argument(
        "--size-scale",
        type=float,
        default=0.0,
        help="Additional radius added per sqrt(pixel_count).",
    )
    parser.add_argument(
        "--alpha",
        type=int,
        default=200,
        help="Alpha value for spot color (0-255).",
    )
    parser.add_argument(
        "--annotate",
        action="store_true",
        help="Draw segment IDs next to each spot.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit to the first N segments (useful for very dense overlays).",
    )

    args = parser.parse_args()
    color = (255, 0, 0, max(0, min(255, args.alpha)))
    segments = read_segments(args.segments)
    draw_segments(
        args.image,
        segments,
        args.output,
        radius=args.radius,
        size_scale=args.size_scale,
        color=color,
        annotate=args.annotate,
        max_segments=args.limit,
    )
    print("Wrote overlay to", args.output)


if __name__ == "__main__":
    main()
