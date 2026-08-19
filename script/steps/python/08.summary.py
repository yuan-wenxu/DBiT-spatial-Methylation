#!/usr/bin/env python3
"""Summarize DBiT methylation coverage, read metrics, and spatial heatmaps."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Optional

if "MPLCONFIGDIR" not in os.environ:
    _MPL_CACHE = tempfile.TemporaryDirectory(prefix="dbitm-matplotlib-")
    os.environ["MPLCONFIGDIR"] = _MPL_CACHE.name

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam


VALID_FLAGS = {83, 99, 147, 163}
CH_CONTEXTS = ("ca", "cc", "ct")
CONTEXT_SPECS = {
    "cg": {
        "suffix": "CG",
        "mean_field": "mean_methylation",
        "site_field": "cpg_site_count",
        "site_heatmap": "cpg_site_count_heatmap.png",
        "site_title": "CpG Sites per Spot",
        "site_label": "CpG sites",
        "mean_heatmap": "mean_methylation_heatmap.png",
        "mean_violin": "mean_methylation_violin.png",
        "mean_title": "Mean CpG Methylation per Spot",
        "violin_title": "Spot Mean CpG Methylation Distribution",
        "host_mean_field": "host_spot_mean_methylation",
        "host_median_field": "host_spot_median_cpg_sites",
        "mito_mean_field": "host_mito_mean_methylation",
        "spike_mean_suffix": "mean_methylation",
    },
    **{
        context: {
            "suffix": context.upper(),
            "mean_field": f"mean_{context}_methylation",
            "site_field": f"{context}_site_count",
            "site_heatmap": f"{context}_site_count_heatmap.png",
            "site_title": f"{context.upper()} Sites per Spot",
            "site_label": f"{context.upper()} sites",
            "mean_heatmap": f"mean_{context}_methylation_heatmap.png",
            "mean_violin": f"mean_{context}_methylation_violin.png",
            "mean_title": f"Mean {context.upper()} Methylation per Spot",
            "violin_title": f"Spot Mean {context.upper()} Methylation Distribution",
            "host_mean_field": f"host_spot_mean_{context}_methylation",
            "host_median_field": f"host_spot_median_{context}_sites",
            "mito_mean_field": f"host_mito_mean_{context}_methylation",
            "spike_mean_suffix": f"mean_{context}_methylation",
        }
        for context in CH_CONTEXTS
    },
}


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate per-spot/sample summary TSVs and spatial heatmaps."
    )
    parser.add_argument("--work-dir", required=True, help="Path to the dbitm directory.")
    parser.add_argument("--output-dir")
    parser.add_argument("--cb-tag", default="CB")
    parser.add_argument("--min-mapping-quality", type=int, default=10)
    parser.add_argument(
        "--context-mode",
        choices=("cg", "ch", "both"),
        default="cg",
        help="Coverage contexts to summarize; ch expands to CA, CC, and CT (default: cg).",
    )
    parser.add_argument("--spike-in-name", action="append", default=[])
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def to_int(value: object) -> Optional[int]:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)) and math.isfinite(value):
        return int(value)
    if isinstance(value, str):
        try:
            return int(float(value.strip()))
        except ValueError:
            return None
    return None


def format_float(value: Optional[float]) -> str:
    return "NA" if value is None else f"{value:.6f}"


def format_int(value: Optional[int]) -> str:
    return "NA" if value is None else str(value)


def format_percentage(numerator: Optional[int], denominator: Optional[int]) -> str:
    if numerator is None or denominator is None or denominator <= 0:
        return "NA"
    return f"{100 * numerator / denominator:.2f}%"


def parse_cov_stats(path: Path) -> Optional[tuple[float, int]]:
    if not path.is_file():
        return None
    methylation_sum = 0.0
    site_count = 0
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            try:
                methylation = float(fields[3])
                methylated = int(fields[4])
                unmethylated = int(fields[5])
            except ValueError:
                continue
            if methylated + unmethylated <= 0:
                continue
            methylation_sum += methylation
            site_count += 1
    if site_count == 0:
        return None
    return methylation_sum / site_count, site_count


def read_saturation_rate(path: Path) -> Optional[str]:
    if not path.is_file():
        return None
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            value = (row.get("saturation_rate") or "").strip()
            return value or None
    return None


def load_spot_manifest(
    path: Path,
) -> tuple[dict[str, tuple[str, str]], dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"spot manifest not found: {path}")
    coordinates: dict[str, tuple[str, str]] = {}
    cb_to_spot: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"spot_id", "barcode1_index", "barcode2_index", "raw_cb"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"spot manifest lacks required columns: {path}")
        for row in reader:
            spot = (row.get("spot_id") or "").strip()
            raw_cb = (row.get("raw_cb") or "").strip().upper()
            if not spot or not raw_cb:
                continue
            x_value = to_int(row.get("barcode1_index"))
            y_value = to_int(row.get("barcode2_index"))
            if x_value is None or y_value is None:
                continue
            coordinates[spot] = (str(x_value), str(y_value))
            cb_to_spot[raw_cb] = spot
            cb_to_spot[raw_cb.replace("+", "")] = spot
    if not coordinates:
        raise ValueError(f"spot manifest contains no usable entries: {path}")
    return coordinates, cb_to_spot


def is_unique_mapped(record: pysam.AlignedSegment, minimum_mapq: int) -> bool:
    if record.is_unmapped or record.mapping_quality < minimum_mapq:
        return False
    try:
        nh = record.get_tag("NH")
    except KeyError:
        return True
    return not isinstance(nh, int) or nh <= 1


def count_host_metrics(
    bam_path: Path,
    cb_tag: str,
    cb_to_spot: dict[str, str],
    minimum_mapq: int,
) -> tuple[dict[str, int], Optional[int], Optional[int]]:
    if not bam_path.is_file():
        return {}, None, None
    spot_counts: dict[str, int] = defaultdict(int)
    mapped_reads = 0
    valid_reads = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for record in bam.fetch(until_eof=True):
            try:
                cb_value = str(record.get_tag(cb_tag)).upper()
            except KeyError:
                cb_value = ""
            spot = cb_to_spot.get(cb_value)
            if spot is not None:
                spot_counts[spot] += 1
            if not is_unique_mapped(record, minimum_mapq):
                continue
            mapped_reads += 1
            if record.flag in VALID_FLAGS:
                valid_reads += 1
    return dict(spot_counts), mapped_reads, valid_reads


def count_mapped_reads(path: Path, minimum_mapq: int) -> Optional[int]:
    if not path.is_file():
        return None
    count = 0
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for record in bam.fetch(until_eof=True):
            if is_unique_mapped(record, minimum_mapq):
                count += 1
    return count


def parse_fastp_reads(path: Path) -> Optional[int]:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    summary = payload.get("summary")
    if isinstance(summary, dict):
        before = summary.get("before_filtering")
        if isinstance(before, dict):
            value = to_int(before.get("total_reads"))
            if value is not None:
                return value
        value = to_int(summary.get("total_reads"))
        if value is not None:
            return value
    return to_int(payload.get("total_reads"))


def parse_barcoded_reads(path: Path) -> Optional[int]:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    kept_pairs = to_int(payload.get("kept_reads"))
    return None if kept_pairs is None else kept_pairs * 2


def summarize_spots(
    context_cov_paths: dict[str, list[Path]],
    spot_counts: dict[str, int],
    coordinates: dict[str, tuple[str, str]],
) -> list[dict[str, str]]:
    context_stats: dict[str, dict[str, Optional[tuple[float, int]]]] = {}
    observed_spots = set(coordinates) | set(spot_counts)
    for context, cov_paths in context_cov_paths.items():
        suffix = str(CONTEXT_SPECS[context]["suffix"])
        cov_stats: dict[str, Optional[tuple[float, int]]] = {}
        for path in cov_paths:
            spot = path.name.removesuffix(f".{suffix}.cov")
            cov_stats[spot] = parse_cov_stats(path)
        context_stats[context] = cov_stats
        observed_spots.update(cov_stats)

    rows: list[dict[str, str]] = []
    for spot in sorted(observed_spots):
        x_index, y_index = coordinates.get(spot, ("NA", "NA"))
        row = {"X_index": x_index, "Y_index": y_index, "spot": spot}
        for context in context_cov_paths:
            spec = CONTEXT_SPECS[context]
            stats = context_stats[context].get(spot)
            row[str(spec["mean_field"])] = format_float(stats[0] if stats else None)
            row[str(spec["site_field"])] = format_int(stats[1] if stats else None)
        row["reads"] = str(spot_counts.get(spot, 0))
        rows.append(row)
    return rows


def mean_from_rows(
    rows: list[dict[str, str]], mean_field: str, site_field: str
) -> tuple[Optional[float], Optional[float]]:
    weighted_sum = 0.0
    total_sites = 0
    counts: list[int] = []
    for row in rows:
        if row[mean_field] == "NA" or row[site_field] == "NA":
            continue
        site_count = int(row[site_field])
        if site_count <= 0:
            continue
        weighted_sum += float(row[mean_field]) * site_count
        total_sites += site_count
        counts.append(site_count)
    mean = weighted_sum / total_sites if total_sites else None
    median = statistics.median(counts) if counts else None
    return mean, float(median) if median is not None else None


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def heatmap_ticks(maximum: int) -> list[int]:
    if maximum <= 10:
        return list(range(maximum + 1))
    step = max(1, math.ceil((maximum + 1) / 10))
    ticks = list(range(0, maximum + 1, step))
    if ticks[-1] != maximum:
        ticks.append(maximum)
    return ticks


def write_heatmap(
    rows: list[dict[str, str]],
    field: str,
    output: Path,
    title: str,
    colorbar_label: str,
    cmap_name: str,
    vmax: Optional[float],
) -> None:
    points: list[tuple[int, int, float]] = []
    for row in rows:
        if row[field] == "NA" or row["X_index"] == "NA" or row["Y_index"] == "NA":
            continue
        points.append((int(row["X_index"]), int(row["Y_index"]), float(row[field])))
    figure, axis = plt.subplots(figsize=(6, 5))
    if not points:
        axis.axis("off")
        axis.text(0.5, 0.5, "No valid data", ha="center", va="center")
        axis.set_title(title)
    else:
        max_x = max(point[0] for point in points)
        max_y = max(point[1] for point in points)
        matrix = [[math.nan for _ in range(max_x + 1)] for _ in range(max_y + 1)]
        for x_index, y_index, value in points:
            matrix[y_index][x_index] = value
        cmap = plt.get_cmap(cmap_name).copy()
        cmap.set_bad("#f2f2f2")
        image = axis.imshow(
            matrix,
            origin="lower",
            interpolation="nearest",
            aspect="equal",
            cmap=cmap,
            vmin=0.0,
            vmax=vmax,
            extent=(-0.5, max_x + 0.5, -0.5, max_y + 0.5),
        )
        axis.set_xticks(heatmap_ticks(max_x))
        axis.set_yticks(heatmap_ticks(max_y))
        axis.set_xlabel("X index", fontsize=8)
        axis.set_ylabel("Y index", fontsize=8)
        axis.set_title(title, fontsize=12)
        colorbar = figure.colorbar(image, ax=axis, shrink=0.9)
        colorbar.set_label(colorbar_label)
    figure.tight_layout()
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    try:
        figure.savefig(temporary, format="png", dpi=200)
        os.replace(temporary, output)
    finally:
        plt.close(figure)
        temporary.unlink(missing_ok=True)


def write_methylation_violin(
    rows: list[dict[str, str]], output: Path, field: str, title: str
) -> None:
    values: list[float] = []
    for row in rows:
        text = row[field]
        if text == "NA":
            continue
        value = float(text)
        if math.isfinite(value):
            values.append(value)

    figure, axis = plt.subplots(figsize=(4, 5))
    if not values:
        axis.axis("off")
        axis.text(0.5, 0.5, "No valid data", ha="center", va="center")
        axis.set_title(title)
    else:
        if len(set(values)) > 1:
            violin = axis.violinplot(
                values,
                positions=[1],
                widths=0.7,
                showmeans=True,
                showmedians=True,
                showextrema=True,
            )
            for body in violin["bodies"]:
                body.set_facecolor("tab:blue")
                body.set_edgecolor("black")
                body.set_alpha(0.7)
        else:
            axis.scatter([1], [values[0]], color="tab:blue", s=35, zorder=3)
            axis.hlines(values[0], 0.8, 1.2, color="black", linewidth=1.5)
        axis.set_xticks([1], [f"All spots\n(n={len(values)})"])
        axis.set_xlim(0.5, 1.5)
        axis.set_ylabel("Mean methylation (%)", fontsize=10)
        axis.set_title(title, fontsize=12)
        axis.grid(axis="y", alpha=0.25)

    figure.tight_layout()
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    try:
        figure.savefig(temporary, format="png", dpi=200)
        os.replace(temporary, output)
    finally:
        plt.close(figure)
        temporary.unlink(missing_ok=True)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    if args.min_mapping_quality < 0:
        print("[summary] error: --min-mapping-quality must be >= 0", file=sys.stderr)
        return 1
    work_dir = Path(args.work_dir)
    output_dir = Path(args.output_dir) if args.output_dir else work_dir / "summary"
    coverage_dir = work_dir / "coverage"
    manifest_path = coverage_dir / "spot_manifest.tsv"
    host_bam = work_dir / "pooled" / "pooled.cb.bam"
    saturation_summary = work_dir / "saturation" / "saturation_summary.tsv"
    if args.context_mode == "cg":
        contexts = ["cg"]
    elif args.context_mode == "ch":
        contexts = list(CH_CONTEXTS)
    else:
        contexts = ["cg", *CH_CONTEXTS]
    host_cov_paths = {
        context: sorted(
            (coverage_dir / "host").rglob(
                f"*.{CONTEXT_SPECS[context]['suffix']}.cov"
            )
        )
        for context in contexts
    }
    spike_names = list(
        dict.fromkeys(name.strip() for name in args.spike_in_name if name.strip())
    )
    if not spike_names:
        discovered_spikes: set[str] = set()
        for context in contexts:
            suffix = str(CONTEXT_SPECS[context]["suffix"])
            discovered_spikes.update(
                path.name.removesuffix(f".{suffix}.cov")
                for path in coverage_dir.glob(f"*.{suffix}.cov")
                if path.name != f"host_mito.{suffix}.cov"
            )
        spike_names = sorted(discovered_spikes)

    print(f"[summary] work-dir={work_dir}")
    print(f"[summary] context-mode={args.context_mode}")
    for context in contexts:
        print(f"[summary] host-{context}-cov-count={len(host_cov_paths[context])}")
    print(f"[summary] spike-names={','.join(spike_names) if spike_names else 'none'}")
    print(f"[summary] saturation-summary={saturation_summary}")
    print(f"[summary] output-dir={output_dir}")
    if args.dry_run:
        print("[summary] dry-run=1")
        return 0

    try:
        coordinates, cb_to_spot = load_spot_manifest(manifest_path)
        spot_counts, host_mapped, host_valid = count_host_metrics(
            host_bam, args.cb_tag, cb_to_spot, args.min_mapping_quality
        )
        per_spot_rows = summarize_spots(host_cov_paths, spot_counts, coordinates)
        if not per_spot_rows:
            raise ValueError("no per-spot coverage or CB-tagged reads were found")
        output_dir.mkdir(parents=True, exist_ok=True)
        per_spot_path = output_dir / "per_spot_summary.tsv"
        per_spot_fields = ["X_index", "Y_index", "spot"]
        for context in contexts:
            spec = CONTEXT_SPECS[context]
            per_spot_fields.extend([str(spec["mean_field"]), str(spec["site_field"])])
        per_spot_fields.append("reads")
        write_tsv(
            per_spot_path,
            per_spot_fields,
            per_spot_rows,
        )

        sample_row: dict[str, str] = {
            "saturation_rate": read_saturation_rate(saturation_summary) or "NA"
        }
        for context in contexts:
            spec = CONTEXT_SPECS[context]
            host_mean, host_median_sites = mean_from_rows(
                per_spot_rows,
                str(spec["mean_field"]),
                str(spec["site_field"]),
            )
            suffix = str(spec["suffix"])
            mito_stats = parse_cov_stats(coverage_dir / f"host_mito.{suffix}.cov")
            sample_row[str(spec["host_mean_field"])] = format_float(host_mean)
            sample_row[str(spec["host_median_field"])] = format_float(host_median_sites)
            sample_row[str(spec["mito_mean_field"])] = format_float(
                mito_stats[0] if mito_stats else None
            )
            for spike_name in spike_names:
                stats = parse_cov_stats(coverage_dir / f"{spike_name}.{suffix}.cov")
                sample_row[f"{spike_name}_{spec['spike_mean_suffix']}"] = format_float(
                    stats[0] if stats else None
                )
        raw_reads = parse_fastp_reads(work_dir / "fastp" / "fastp.json")
        barcoded_reads = parse_barcoded_reads(work_dir / "barcode" / "stats.json")
        sample_row.update(
            {
                "raw_reads": format_int(raw_reads),
                "barcoded_reads": format_int(barcoded_reads),
                "barcoded_reads_rate": format_percentage(barcoded_reads, raw_reads),
                "host_mapped_reads": format_int(host_mapped),
            }
        )
        for spike_name in spike_names:
            mapped = count_mapped_reads(
                work_dir / "pooled" / f"pooled.{spike_name}.bam",
                args.min_mapping_quality,
            )
            sample_row[f"{spike_name}_mapped_reads"] = format_int(mapped)
        sample_row["host_valid_reads"] = format_int(host_valid)
        sample_row["valid_reads_rate"] = format_percentage(host_valid, raw_reads)
        write_tsv(
            output_dir / "sample_summary.tsv",
            list(sample_row),
            [sample_row],
        )
        write_heatmap(
            per_spot_rows,
            "reads",
            output_dir / "reads_heatmap.png",
            "Reads per Spot",
            "Reads",
            "Reds",
            None,
        )
        for context in contexts:
            spec = CONTEXT_SPECS[context]
            write_heatmap(
                per_spot_rows,
                str(spec["site_field"]),
                output_dir / str(spec["site_heatmap"]),
                str(spec["site_title"]),
                str(spec["site_label"]),
                "magma",
                None,
            )
            write_heatmap(
                per_spot_rows,
                str(spec["mean_field"]),
                output_dir / str(spec["mean_heatmap"]),
                str(spec["mean_title"]),
                f"Mean {context.upper()} methylation (%)",
                "coolwarm",
                100.0,
            )
            write_methylation_violin(
                per_spot_rows,
                output_dir / str(spec["mean_violin"]),
                str(spec["mean_field"]),
                str(spec["violin_title"]),
            )
    except (OSError, ValueError) as error:
        print(f"[summary] error: {error}", file=sys.stderr)
        return 1

    print(f"[summary] per-spot-rows={len(per_spot_rows)}")
    print(f"[summary] output={output_dir}")
    print("[summary] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
