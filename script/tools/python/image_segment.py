#!/usr/bin/env python3
"""Create a tissue mask and full-resolution DBiT spot positions."""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path

import cv2
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a tissue mask and tissue_positions.tsv.gz from a full-resolution "
            "image and the adjacent mask.png frame mask."
        )
    )
    parser.add_argument("--image_path", required=True, type=Path)
    parser.add_argument("--result_path", required=True, type=Path)
    parser.add_argument("--barcodeA_whitelist", required=True, type=Path)
    parser.add_argument("--barcodeB_whitelist", required=True, type=Path)
    parser.add_argument("--x_spots_number", type=int, required=True)
    parser.add_argument("--y_spots_number", type=int, required=True)
    parser.add_argument("--length_spot", type=int, required=True)
    parser.add_argument("--interval", type=int, required=True)
    parser.add_argument("--pixel_length", type=float, required=True)
    return parser.parse_args()


BACKGROUND_THRESHOLD = 15
LOCAL_DENSITY_KERNEL = 201
MIN_LOCAL_SIGNAL_FRACTION = 0.02
MORPHOLOGY_KERNEL = 11
MIN_COMPONENT_FRACTION = 0.3
MIN_COMPONENT_RELATIVE_TO_LARGEST = 0.1
MIN_TISSUE_FRACTION = 0.03


@dataclass(frozen=True)
class GridConfig:
    row_count: int
    col_count: int
    length_spot: int
    interval: int
    pixel_length: float


@dataclass(frozen=True)
class FramePlacement:
    x: int
    y: int
    visible_x: int
    visible_y: int
    visible_width: int
    visible_height: int
    visible_offset_x: int
    visible_offset_y: int


def read_barcode_components(whitelist_path: Path) -> list[str]:
    path = str(whitelist_path)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        components = [line.strip() for line in handle if line.strip()]
    if len(components) != len(set(components)):
        raise ValueError(f"Whitelist contains duplicate barcodes: {path}")
    return components


def frame_region_from_mask(mask: np.ndarray) -> np.ndarray:
    if mask.ndim == 2:
        return mask > 0
    if mask.ndim != 3 or mask.shape[2] not in (3, 4):
        raise ValueError(f"Unsupported frame-mask shape: {mask.shape}")
    if mask.shape[2] == 4:
        alpha = mask[:, :, 3]
        if np.any(alpha == 0) and np.any(alpha > 0):
            return alpha > 0
    grayscale = cv2.cvtColor(mask[:, :, :3], cv2.COLOR_BGR2GRAY)
    return grayscale > 0


def infer_axis_placement(
    start: int,
    extent: int,
    canvas_extent: int,
    expected_extent: int,
    axis_name: str,
) -> tuple[int, int]:
    if extent >= expected_extent:
        return start, 0

    touches_low = start == 0
    touches_high = start + extent == canvas_extent
    if touches_low == touches_high:
        raise ValueError(
            f"Visible frame {axis_name} extent {extent} is smaller than the "
            f"configured extent {expected_extent}, but it does not touch exactly "
            f"one image boundary; the clipped-frame position is ambiguous"
        )
    if touches_low:
        return start + extent - expected_extent, expected_extent - extent
    return start, 0


def locate_frame(
    mask_path: Path,
    image_shape,
    expected_size: tuple[int, int],
) -> tuple[np.ndarray, FramePlacement]:
    mask = cv2.imread(str(mask_path), cv2.IMREAD_UNCHANGED)
    if mask is None:
        raise ValueError(f"Unable to read frame mask: {mask_path}")
    if mask.shape[:2] != image_shape[:2]:
        raise ValueError(
            f"Frame mask size {mask.shape[:2]} does not match image size "
            f"{image_shape[:2]}"
        )
    frame_region = frame_region_from_mask(mask)
    nonzero = cv2.findNonZero(frame_region.astype(np.uint8))
    if nonzero is None:
        raise ValueError(f"Frame mask contains no nonzero region: {mask_path}")
    visible_x, visible_y, visible_width, visible_height = cv2.boundingRect(nonzero)
    expected_width, expected_height = expected_size
    frame_x, visible_offset_x = infer_axis_placement(
        visible_x,
        visible_width,
        image_shape[1],
        expected_width,
        "horizontal",
    )
    frame_y, visible_offset_y = infer_axis_placement(
        visible_y,
        visible_height,
        image_shape[0],
        expected_height,
        "vertical",
    )
    return frame_region, FramePlacement(
        x=frame_x,
        y=frame_y,
        visible_x=visible_x,
        visible_y=visible_y,
        visible_width=visible_width,
        visible_height=visible_height,
        visible_offset_x=visible_offset_x,
        visible_offset_y=visible_offset_y,
    )


def grayscale_uint8(image: np.ndarray) -> np.ndarray:
    if image.ndim == 2:
        grayscale = image
    elif image.ndim == 3 and image.shape[2] == 1:
        grayscale = image[:, :, 0]
    elif image.ndim == 3 and image.shape[2] >= 3:
        grayscale = cv2.cvtColor(image[:, :, :3], cv2.COLOR_BGR2GRAY)
    else:
        raise ValueError(f"Unsupported image shape: {image.shape}")
    if grayscale.dtype == np.uint8:
        return grayscale
    return cv2.normalize(grayscale, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)


def save_fullres_grayscale(image: np.ndarray, output_path: Path) -> Path:
    grayscale = grayscale_uint8(image)
    if grayscale.shape[:2] != image.shape[:2]:
        raise RuntimeError("Grayscale conversion changed the image resolution")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not cv2.imwrite(str(output_path), grayscale):
        raise RuntimeError(f"Failed to write full-resolution grayscale image: {output_path}")
    return output_path


def remove_small_regions(mask: np.ndarray) -> np.ndarray:
    contours, _ = cv2.findContours(
        mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
    )
    if not contours:
        return np.zeros_like(mask)
    areas = [cv2.contourArea(contour) for contour in contours]
    min_area = max(
        64,
        mask.size * MIN_COMPONENT_FRACTION,
        max(areas) * MIN_COMPONENT_RELATIVE_TO_LARGEST,
    )
    retained = [
        contour for contour, area in zip(contours, areas) if area >= min_area
    ]
    cleaned = np.zeros_like(mask)
    cv2.drawContours(cleaned, retained, -1, 255, thickness=cv2.FILLED)
    return cleaned


def generate_tissue_mask(
    frame_image: np.ndarray,
    frame_region: np.ndarray,
    output_path: Path,
) -> np.ndarray:
    intensity = grayscale_uint8(frame_image)
    signal = np.where(intensity > BACKGROUND_THRESHOLD, 255, 0).astype(np.uint8)
    local_density = cv2.boxFilter(
        signal,
        ddepth=-1,
        ksize=(LOCAL_DENSITY_KERNEL, LOCAL_DENSITY_KERNEL),
        normalize=True,
    )
    otsu_threshold, _ = cv2.threshold(
        local_density, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU
    )
    density_threshold = max(
        round(255 * MIN_LOCAL_SIGNAL_FRACTION),
        round(otsu_threshold / 2),
    )
    tissue_mask = np.where(
        local_density >= density_threshold, 255, 0
    ).astype(np.uint8)
    kernel = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE, (MORPHOLOGY_KERNEL, MORPHOLOGY_KERNEL)
    )
    tissue_mask = cv2.morphologyEx(tissue_mask, cv2.MORPH_CLOSE, kernel)
    tissue_mask = remove_small_regions(tissue_mask)
    tissue_mask[~frame_region] = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not cv2.imwrite(str(output_path), tissue_mask):
        raise RuntimeError(f"Failed to write tissue mask: {output_path}")
    return tissue_mask


def spot_bounds(row: int, col: int, grid: GridConfig) -> tuple[int, int, int, int]:
    pitch = grid.length_spot + grid.interval
    display_col = grid.col_count - 1 - col
    x_start = int(display_col * pitch / grid.pixel_length)
    y_start = int(row * pitch / grid.pixel_length)
    spot_pixels = int(grid.length_spot / grid.pixel_length)
    return x_start, y_start, x_start + spot_pixels, y_start + spot_pixels


def build_tissue_positions(
    tissue_mask: np.ndarray,
    frame_x: int,
    frame_y: int,
    visible_offset_x: int,
    visible_offset_y: int,
    image_width: int,
    image_height: int,
    grid: GridConfig,
    barcode_as: list[str],
    barcode_bs: list[str],
) -> pd.DataFrame:
    if len(barcode_as) != grid.row_count:
        raise ValueError(
            f"Barcode A whitelist has {len(barcode_as)} entries; "
            f"expected {grid.row_count}"
        )
    if len(barcode_bs) != grid.col_count:
        raise ValueError(
            f"Barcode B whitelist has {len(barcode_bs)} entries; "
            f"expected {grid.col_count}"
        )

    records = []
    for row in range(grid.row_count):
        for col in range(grid.col_count):
            x_start, y_start, x_end, y_end = spot_bounds(row, col, grid)
            visible_x_start = max(0, x_start - visible_offset_x)
            visible_y_start = max(0, y_start - visible_offset_y)
            visible_x_end = min(tissue_mask.shape[1], x_end - visible_offset_x)
            visible_y_end = min(tissue_mask.shape[0], y_end - visible_offset_y)
            tissue_pixels = 0
            if (
                visible_x_start < visible_x_end
                and visible_y_start < visible_y_end
            ):
                spot_mask = tissue_mask[
                    visible_y_start:visible_y_end,
                    visible_x_start:visible_x_end,
                ]
                tissue_pixels = int(np.count_nonzero(spot_mask))
            spot_area = (x_end - x_start) * (y_end - y_start)
            tissue_fraction = tissue_pixels / spot_area
            center_x = frame_x + (x_start + x_end) // 2
            center_y = frame_y + (y_start + y_end) // 2
            center_in_image = (
                0 <= center_x < image_width and 0 <= center_y < image_height
            )
            records.append(
                {
                    "barcode": barcode_bs[col] + barcode_as[row],
                    "in_tissue": int(
                        center_in_image
                        and tissue_fraction >= MIN_TISSUE_FRACTION
                    ),
                    "array_row": f"{row:02d}",
                    "array_col": f"{col:02d}",
                    "pxl_row_in_fullres": center_y,
                    "pxl_col_in_fullres": center_x,
                }
            )
    return pd.DataFrame.from_records(
        records,
        columns=[
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ],
    )


def process_image(
    image_path: Path,
    result_path: Path,
    grid: GridConfig,
    barcode_a_whitelist: Path,
    barcode_b_whitelist: Path,
) -> pd.DataFrame:
    image = cv2.imread(str(image_path), cv2.IMREAD_UNCHANGED)
    if image is None:
        raise ValueError(f"Unable to read image: {image_path}")
    frame_mask_path = image_path.with_name("mask.png")
    if not frame_mask_path.is_file():
        raise FileNotFoundError(f"Adjacent frame mask not found: {frame_mask_path}")

    expected_width = int(
        (
            grid.col_count * grid.length_spot
            + (grid.col_count - 1) * grid.interval
        )
        / grid.pixel_length
    )
    expected_height = int(
        (
            grid.row_count * grid.length_spot
            + (grid.row_count - 1) * grid.interval
        )
        / grid.pixel_length
    )
    frame_region_full, placement = locate_frame(
        frame_mask_path,
        image.shape,
        (expected_width, expected_height),
    )
    frame_image = image[
        placement.visible_y : placement.visible_y + placement.visible_height,
        placement.visible_x : placement.visible_x + placement.visible_width,
    ]
    frame_region = frame_region_full[
        placement.visible_y : placement.visible_y + placement.visible_height,
        placement.visible_x : placement.visible_x + placement.visible_width,
    ]

    result_path.mkdir(parents=True, exist_ok=True)
    grayscale_path = save_fullres_grayscale(
        image,
        result_path / "fullres_grayscale.png",
    )
    tissue_mask = generate_tissue_mask(
        frame_image,
        frame_region,
        result_path / "tissue_mask.png",
    )
    positions = build_tissue_positions(
        tissue_mask,
        placement.x,
        placement.y,
        placement.visible_offset_x,
        placement.visible_offset_y,
        image.shape[1],
        image.shape[0],
        grid,
        read_barcode_components(barcode_a_whitelist),
        read_barcode_components(barcode_b_whitelist),
    )
    output_path = result_path / "tissue_positions.tsv.gz"
    positions.to_csv(output_path, sep="\t", index=False, compression="gzip")
    print(
        f"Frame bbox: x={placement.x}, y={placement.y}, "
        f"width={expected_width}, height={expected_height}"
    )
    print(
        f"Visible frame bbox: x={placement.visible_x}, y={placement.visible_y}, "
        f"width={placement.visible_width}, height={placement.visible_height}"
    )
    tissue_spots = int((positions["in_tissue"] == 1).sum())
    print(f"Tissue spots: {tissue_spots}/{len(positions)}")
    print(f"Wrote full-resolution grayscale image: {grayscale_path.resolve()}")
    print(f"Wrote tissue mask: {(result_path / 'tissue_mask.png').resolve()}")
    print(f"Wrote tissue positions: {output_path.resolve()}")
    return positions


def main() -> None:
    args = parse_args()
    process_image(
        image_path=args.image_path,
        result_path=args.result_path,
        grid=GridConfig(
            args.x_spots_number,
            args.y_spots_number,
            args.length_spot,
            args.interval,
            args.pixel_length,
        ),
        barcode_a_whitelist=args.barcodeA_whitelist,
        barcode_b_whitelist=args.barcodeB_whitelist,
    )


if __name__ == "__main__":
    main()
