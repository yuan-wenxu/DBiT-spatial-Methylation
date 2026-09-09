#!/usr/bin/env python3
"""Estimate host CpG saturation from per-spot coverage and BAM read counts."""

from __future__ import annotations

import argparse
import csv
import gzip
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


DEFAULT_FRACTIONS = (
    0.01,
    0.02,
    0.05,
    0.10,
    0.20,
    0.30,
    0.40,
    0.50,
    0.60,
    0.70,
    0.80,
    0.90,
    1.00,
)

SMC_EXCLUDED_SPOTS = frozenset({"00_01"})
MIN_THRESHOLD_INFERENCE_SPOTS = 100
MIN_THRESHOLD_CLASS_SPOTS = 10
MIN_THRESHOLD_RETAINED_SPOTS = 100
MIN_THRESHOLD_CLASS_FRACTION = 0.01
MIN_THRESHOLD_SEPARATION = 0.60


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Estimate host CpG saturation across spatial spots."
    )
    parser.add_argument("--work-dir", required=True, help="Path to the dbitm directory.")
    parser.add_argument(
        "--assay",
        choices=("taps", "taps-v2", "emseq", "cabernet", "smc"),
        default="taps",
        help="Assay type used to interpret spatial barcodes. Default: taps.",
    )
    parser.add_argument("--output-dir")
    parser.add_argument("--cb-tag", default="CB")
    parser.add_argument(
        "--reads-threshold",
        type=float,
        default=1_000_000.0,
        help=(
            "Fallback reads threshold used when automatic inference fails. "
            "Default: 1000000."
        ),
    )
    parser.add_argument("--pred-fraction", type=float, default=2.0)
    parser.add_argument("--linear-r2-threshold", type=float, default=0.99)
    parser.add_argument("--fastp-json")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if args.reads_threshold <= 0:
        raise ValueError("--reads-threshold must be > 0")
    if args.pred_fraction <= 0:
        raise ValueError("--pred-fraction must be > 0")
    if not 0 < args.linear_r2_threshold <= 1:
        raise ValueError("--linear-r2-threshold must be in (0, 1]")


def format_optional_int(value: Optional[float]) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return str(int(round(value)))


def format_optional_float(value: Optional[float]) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.6f}"


def load_cb_to_spot(path: Path, c_to_t: bool = False) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"spot manifest not found: {path}")
    cb_to_spot: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"spot_id", "raw_cb"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"spot manifest lacks required columns: {path}")
        for row in reader:
            spot = (row.get("spot_id") or "").strip()
            raw_cb = (row.get("raw_cb") or "").strip().upper()
            if not spot or not raw_cb:
                continue
            lookup_cb = raw_cb.replace("C", "T") if c_to_t else raw_cb
            cb_to_spot[lookup_cb] = spot
            cb_to_spot[lookup_cb.replace("+", "")] = spot
    if not cb_to_spot:
        raise ValueError(f"spot manifest contains no usable entries: {path}")
    return cb_to_spot


def count_spot_reads(
    bam_paths: list[Path], cb_tag: str, cb_to_spot: dict[str, str]
) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for bam_path in bam_paths:
        if not bam_path.is_file():
            raise FileNotFoundError(f"host BAM not found: {bam_path}")
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for record in bam.fetch(until_eof=True):
                try:
                    cb_value = str(record.get_tag(cb_tag)).upper()
                except KeyError:
                    continue
                spot = cb_to_spot.get(cb_value)
                if spot is not None:
                    counts[spot] += 1
    return dict(counts)


def infer_reads_threshold(
    spot_reads: dict[str, int],
    excluded_spots: frozenset[str] = frozenset(),
) -> tuple[Optional[int], Optional[float]]:
    log_reads = sorted(
        math.log10(reads + 1.0)
        for spot, reads in spot_reads.items()
        if reads > 0 and spot not in excluded_spots
    )
    value_count = len(log_reads)
    if value_count < MIN_THRESHOLD_INFERENCE_SPOTS:
        return None, None

    total = sum(log_reads)
    mean = total / value_count
    total_sum_squares = sum((value - mean) ** 2 for value in log_reads)
    if total_sum_squares <= 0:
        return None, None

    minimum_class_size = max(
        MIN_THRESHOLD_CLASS_SPOTS,
        math.ceil(value_count * MIN_THRESHOLD_CLASS_FRACTION),
    )
    left_sum = 0.0
    best_between_sum_squares = -1.0
    best_split: Optional[int] = None
    for index in range(value_count - 1):
        left_sum += log_reads[index]
        left_count = index + 1
        right_count = value_count - left_count
        if (
            left_count < minimum_class_size
            or right_count < MIN_THRESHOLD_RETAINED_SPOTS
        ):
            continue
        if log_reads[index] == log_reads[index + 1]:
            continue
        left_mean = left_sum / left_count
        right_mean = (total - left_sum) / right_count
        between_sum_squares = (
            left_count
            * right_count
            / value_count
            * (left_mean - right_mean) ** 2
        )
        if between_sum_squares > best_between_sum_squares:
            best_between_sum_squares = between_sum_squares
            best_split = index

    if best_split is None:
        return None, None
    separation = best_between_sum_squares / total_sum_squares
    if separation < MIN_THRESHOLD_SEPARATION:
        return None, separation

    log_threshold = (
        log_reads[best_split] + log_reads[best_split + 1]
    ) / 2.0
    threshold = max(1, int(round(10**log_threshold - 1.0)))
    return threshold, separation


def parse_barcoded_cov_histograms(
    path: Path, selected_spots: set[str]
) -> dict[str, dict[int, int]]:
    histograms: dict[str, dict[int, int]] = {}
    if not path.is_file():
        return histograms
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            spot = fields[6]
            if spot not in selected_spots:
                continue
            try:
                depth = int(fields[4]) + int(fields[5])
            except ValueError:
                continue
            if depth > 0:
                histogram = histograms.setdefault(spot, {})
                histogram[depth] = histogram.get(depth, 0) + 1
    return histograms


def expected_unique(histogram: dict[int, int], fraction: float) -> float:
    return sum(
        count * (1.0 - (1.0 - fraction) ** depth)
        for depth, count in histogram.items()
    )


def median_and_iqr(values: list[float]) -> tuple[float, float, float]:
    median = float(statistics.median(values))
    if len(values) == 1:
        return median, 0.0, 0.0
    q1, _, q3 = statistics.quantiles(values, n=4, method="inclusive")
    return median, max(0.0, median - q1), max(0.0, q3 - median)


def saturation_function(fraction: float, maximum: float, rate: float) -> float:
    return maximum * (1.0 - math.exp(-rate * fraction))


def fit_linear_through_origin(fractions: list[float], values: list[float]) -> float:
    denominator = sum(fraction * fraction for fraction in fractions)
    if denominator <= 0:
        return 0.0
    return sum(
        fraction * value for fraction, value in zip(fractions, values)
    ) / denominator


def r_squared(values: list[float], predictions: list[float]) -> float:
    mean_value = sum(values) / len(values)
    total = sum((value - mean_value) ** 2 for value in values)
    residual = sum(
        (value - prediction) ** 2
        for value, prediction in zip(values, predictions)
    )
    if total <= 0:
        return 1.0 if residual <= 0 else 0.0
    return 1.0 - residual / total


def fit_saturation_curve(
    fractions: list[float], values: list[float]
) -> tuple[float, float]:
    best_maximum = 0.0
    best_rate = 1.0
    best_error = math.inf

    def search(log10_low: float, log10_high: float, count: int) -> None:
        nonlocal best_maximum, best_rate, best_error
        step = (log10_high - log10_low) / float(count - 1)
        for index in range(count):
            rate = 10 ** (log10_low + step * index)
            transformed = [1.0 - math.exp(-rate * fraction) for fraction in fractions]
            denominator = sum(value * value for value in transformed)
            if denominator <= 0:
                continue
            maximum = max(
                0.0,
                sum(
                    value * transformed_value
                    for value, transformed_value in zip(values, transformed)
                )
                / denominator,
            )
            error = sum(
                (value - maximum * transformed_value) ** 2
                for value, transformed_value in zip(values, transformed)
            )
            if error < best_error:
                best_maximum = maximum
                best_rate = rate
                best_error = error

    search(-4.0, 2.0, 800)
    for _ in range(3):
        center = math.log10(best_rate if best_rate > 0 else 1.0)
        search(center - 0.7, center + 0.7, 240)
    return best_maximum, best_rate


def load_sequencing_gbp(path: Path) -> Optional[float]:
    if not path.is_file():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    summary = payload.get("summary", {})
    for section in ("after_filtering", "before_filtering"):
        value = summary.get(section, {}).get("total_bases")
        if value is not None and int(value) > 0:
            return int(value) / 1e9
    return None


SUMMARY_FIELDS = [
    "observed_median_unique_cpgs",
    "theoretical_max_median_unique_cpgs",
    "prediction_fraction",
    "predicted_median_unique_cpgs",
    "saturation_rate",
    "extrapolation_model",
    "hq_spot_count",
]

HISTOGRAM_FIELDS = ["spot", "reads", "depth", "site_count"]
THRESHOLD_FIELDS = ["reads_threshold"]


def write_summary(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerow(row)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_reads_threshold(path: Path, threshold: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=THRESHOLD_FIELDS,
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerow({"reads_threshold": f"{threshold:.12g}"})
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def read_reads_threshold(path: Path) -> float:
    if not path.is_file():
        raise FileNotFoundError(f"reads threshold file not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "reads_threshold" not in reader.fieldnames:
            raise ValueError(
                f"reads threshold file lacks reads_threshold column: {path}"
            )
        row = next(reader, None)
    if row is None:
        raise ValueError(f"reads threshold file contains no data row: {path}")
    raw_threshold = (row.get("reads_threshold") or "").strip()
    try:
        threshold = float(raw_threshold)
    except ValueError as error:
        raise ValueError(
            f"invalid reads threshold in {path}: {raw_threshold or 'empty'}"
        ) from error
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError(f"reads threshold must be > 0 in {path}: {raw_threshold}")
    return threshold


def write_spot_histograms(
    path: Path,
    spot_reads: dict[str, int],
    histograms: dict[str, dict[int, int]],
    excluded_spots: frozenset[str] = frozenset(),
) -> tuple[int, int]:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    spot_count = 0
    row_count = 0
    try:
        with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=HISTOGRAM_FIELDS,
                delimiter="\t",
            )
            writer.writeheader()
            for spot in sorted(histograms):
                if spot in excluded_spots:
                    continue
                histogram = histograms[spot]
                if not histogram:
                    continue
                spot_count += 1
                for depth in sorted(histogram):
                    writer.writerow(
                        {
                            "spot": spot,
                            "reads": spot_reads.get(spot, 0),
                            "depth": depth,
                            "site_count": histogram[depth],
                        }
                    )
                    row_count += 1
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return spot_count, row_count


def read_spot_histograms(
    path: Path,
) -> tuple[dict[str, int], dict[str, dict[int, int]]]:
    if not path.is_file():
        raise FileNotFoundError(f"spot depth histogram not found: {path}")
    spot_reads: dict[str, int] = {}
    histograms: dict[str, dict[int, int]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = set(HISTOGRAM_FIELDS)
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"spot depth histogram lacks required columns: {path}")
        for line_number, row in enumerate(reader, start=2):
            spot = (row.get("spot") or "").strip()
            try:
                reads = int(row.get("reads") or "")
                depth = int(row.get("depth") or "")
                site_count = int(row.get("site_count") or "")
            except ValueError as error:
                raise ValueError(
                    f"invalid histogram value at {path}:{line_number}"
                ) from error
            if not spot or reads < 0 or depth <= 0 or site_count <= 0:
                raise ValueError(
                    f"invalid histogram row at {path}:{line_number}"
                )
            previous_reads = spot_reads.setdefault(spot, reads)
            if previous_reads != reads:
                raise ValueError(
                    f"inconsistent reads for spot {spot} at {path}:{line_number}"
                )
            histogram = histograms.setdefault(spot, {})
            histogram[depth] = histogram.get(depth, 0) + site_count
    return spot_reads, histograms


def save_figure(figure: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        figure.savefig(temporary, format="png", dpi=300, bbox_inches="tight")
        os.replace(temporary, path)
    finally:
        plt.close(figure)
        temporary.unlink(missing_ok=True)


def write_reads_rank_plot(
    path: Path,
    spot_reads: dict[str, int],
    threshold: float,
    threshold_source: str,
    excluded_spots: frozenset[str] = frozenset(),
) -> None:
    ranked_reads = sorted(
        (
            reads
            for spot, reads in spot_reads.items()
            if reads > 0 and spot not in excluded_spots
        ),
        reverse=True,
    )
    figure, axis = plt.subplots(figsize=(5, 4))
    if not ranked_reads:
        axis.axis("off")
        axis.text(
            0.5,
            0.5,
            "No CB-tagged reads were assigned to spots",
            ha="center",
            va="center",
            fontsize=12,
        )
    else:
        ranks = list(range(1, len(ranked_reads) + 1))
        retained_count = sum(reads > threshold for reads in ranked_reads)
        axis.plot(ranks, ranked_reads, color="steelblue", linewidth=1.8)
        axis.axhline(
            threshold,
            color="firebrick",
            linestyle="--",
            linewidth=1.6,
            label=f"Threshold: {threshold:g}",
        )
        axis.set_yscale("log")
        axis.set_xlabel("Spot rank (highest reads first)", fontsize=10)
        axis.set_ylabel("Reads per spot", fontsize=10)
        axis.grid(True, alpha=0.3)
        axis.legend(loc="best")
        axis.text(
            0.98,
            0.04,
            f"Spots above threshold: {retained_count}",
            transform=axis.transAxes,
            ha="right",
            va="bottom",
            fontsize=10,
        )
    axis.set_title(f"Reads per spot rank ({threshold_source})", fontsize=12)
    figure.tight_layout()
    save_figure(figure, path)


def write_empty_plot(path: Path, message: str) -> None:
    figure, axis = plt.subplots(figsize=(5, 4))
    axis.axis("off")
    axis.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
    axis.set_title(f"Saturation analysis\nSaturation rate: NA")
    figure.tight_layout()
    save_figure(figure, path)


def write_plot(
    path: Path,
    fractions: list[float],
    medians: list[float],
    errors_low: list[float],
    errors_high: list[float],
    model: str,
    maximum: float,
    rate: float,
    slope: float,
    pred_fraction: float,
    predicted: Optional[float],
    saturation_rate: Optional[float],
    sequencing_gbp: Optional[float],
) -> None:
    scale = 1e4
    x_scale = sequencing_gbp if sequencing_gbp else 1.0
    x_label = "Sequencing depth (Gbp)" if sequencing_gbp else "Coverage fraction"
    x_values = [fraction * x_scale for fraction in fractions]
    max_fraction = max(max(fractions), pred_fraction, 1.0)
    fit_fractions = [max_fraction * index / 250.0 for index in range(251)]
    fit_x = [fraction * x_scale for fraction in fit_fractions]
    if model == "linear":
        fit_y = [slope * fraction / scale for fraction in fit_fractions]
        fit_label = "Fitted line (linear)"
    else:
        fit_y = [
            saturation_function(fraction, maximum, rate) / scale
            for fraction in fit_fractions
        ]
        fit_label = "Fitted curve (saturation)"

    figure, axis = plt.subplots(figsize=(5, 4))
    axis.errorbar(
        x_values,
        [value / scale for value in medians],
        yerr=(
            [value / scale for value in errors_low],
            [value / scale for value in errors_high],
        ),
        fmt="o-",
        color="blue",
        linewidth=2,
        markersize=5,
        capsize=4,
        label="Observed median (IQR)",
    )
    axis.plot(fit_x, fit_y, "r--", linewidth=2, label=fit_label)
    axis.scatter(
        [pred_fraction * x_scale],
        [(predicted or 0.0) / scale],
        color="green",
        s=55,
        zorder=4,
        label=f"Prediction at {pred_fraction:g}x",
    )
    if model == "saturation" and maximum > 0:
        axis.axhline(
            maximum / scale,
            color="purple",
            linestyle="--",
            linewidth=1.8,
            label="Max CpG sites",
        )
    if saturation_rate is not None:
        rate_text = f"{saturation_rate:.2f}%"
    elif model == "linear":
        rate_text = "NA (linear, unsaturated)"
    else:
        rate_text = "NA"
    axis.set_title(f"Saturation analysis\nSaturation rate: {rate_text}")
    axis.set_xlabel(x_label)
    axis.set_ylabel("Median unique CpGs per spot (x10^4)")
    axis.grid(True, alpha=0.35)
    axis.legend(loc="lower right")
    figure.tight_layout()
    save_figure(figure, path)


def empty_row(hq_spot_count: int, prediction_fraction: float) -> dict[str, str]:
    return {
        "observed_median_unique_cpgs": "NA",
        "theoretical_max_median_unique_cpgs": "NA",
        "prediction_fraction": format_optional_float(prediction_fraction),
        "predicted_median_unique_cpgs": "NA",
        "saturation_rate": "NA",
        "extrapolation_model": "NA",
        "hq_spot_count": str(hq_spot_count),
    }


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        work_dir = Path(args.work_dir)
        manifest_path = work_dir / "coverage" / "spot_manifest.tsv"
        if args.assay == "smc":
            bam_paths = [
                work_dir / "pooled" / "pooled.watson.cb.bam",
                work_dir / "pooled" / "pooled.crick.cb.bam",
            ]
        else:
            bam_paths = [work_dir / "pooled" / "pooled.cb.bam"]
        coverage_path = work_dir / "coverage" / "host" / "host.CG.cov"
        fastp_path = Path(args.fastp_json) if args.fastp_json else work_dir / "fastp" / "fastp.json"
        output_dir = Path(args.output_dir) if args.output_dir else work_dir / "saturation"
        summary_path = output_dir / "saturation_summary.tsv"
        plot_path = output_dir / "saturation_curve.png"
        histogram_path = output_dir / "per_spot_cpg_depth_histogram.tsv.gz"
        threshold_path = output_dir / "reads_threshold.tsv"
        reads_rank_plot_path = output_dir / "reads_per_spot_rank.png"

        barcode_c_to_t = args.assay in {"emseq", "cabernet"}
        excluded_spots = SMC_EXCLUDED_SPOTS if args.assay == "smc" else frozenset()
        print(f"[saturation] assay={args.assay}")
        print(f"[saturation] barcode-c-to-t={str(barcode_c_to_t).lower()}")
        print(f"[saturation] work-dir={work_dir}")
        print(
            "[saturation] host-bams="
            + ",".join(str(path) for path in bam_paths)
        )
        print(f"[saturation] spot-manifest={manifest_path}")
        print(f"[saturation] host-cg-cov={coverage_path}")
        print(f"[saturation] configured-reads-threshold={args.reads_threshold}")
        print(f"[saturation] reads-threshold-file={threshold_path}")
        print(
            "[saturation] excluded-spots="
            + (",".join(sorted(excluded_spots)) if excluded_spots else "none")
        )
        print(f"[saturation] pred-fraction={args.pred_fraction}")
        print(f"[saturation] linear-r2-threshold={args.linear_r2_threshold}")
        print(f"[saturation] fastp-json={fastp_path}")
        print(f"[saturation] spot-depth-histogram={histogram_path}")
        print(f"[saturation] reads-rank-plot={reads_rank_plot_path}")
        print(f"[saturation] output-dir={output_dir}")
        if args.dry_run:
            if threshold_path.is_file():
                print(
                    "[saturation] reads-threshold="
                    f"{read_reads_threshold(threshold_path):g}"
                )
                print("[saturation] reads-threshold-source=file")
            else:
                print("[saturation] reads-threshold-source=pending-auto-inference")
            print(
                "[saturation] spot-depth-histogram-source="
                + ("cache" if histogram_path.is_file() else "pending-generation")
            )
            print("[saturation] dry-run=1")
            return 0

        if histogram_path.is_file():
            spot_reads, spot_histograms = read_spot_histograms(histogram_path)
            histogram_source = "cache"
            histogram_spot_count = len(spot_histograms)
            histogram_row_count = sum(
                len(histogram) for histogram in spot_histograms.values()
            )
        else:
            cb_to_spot = load_cb_to_spot(manifest_path, c_to_t=barcode_c_to_t)
            all_spot_reads = count_spot_reads(bam_paths, args.cb_tag, cb_to_spot)
            if not coverage_path.is_file():
                raise FileNotFoundError(f"host CG coverage not found: {coverage_path}")
            included_spots = set(cb_to_spot.values()) - excluded_spots
            spot_histograms = parse_barcoded_cov_histograms(
                coverage_path, included_spots
            )
            spot_reads = {
                spot: all_spot_reads.get(spot, 0) for spot in spot_histograms
            }
            histogram_spot_count, histogram_row_count = write_spot_histograms(
                histogram_path,
                spot_reads,
                spot_histograms,
            )
            histogram_source = "generated"
        print(f"[saturation] spot-depth-histogram-source={histogram_source}")
        print(f"[saturation] histogram-spot-count={histogram_spot_count}")
        print(f"[saturation] histogram-row-count={histogram_row_count}")
        threshold_source = "file"
        if not threshold_path.is_file():
            inferred_threshold, separation = infer_reads_threshold(
                spot_reads,
                excluded_spots,
            )
            if inferred_threshold is None:
                selected_threshold = args.reads_threshold
                threshold_source = "config fallback"
                print("[saturation] warning=reads_threshold_inference_failed")
                if separation is not None:
                    print(f"[saturation] threshold-separation={separation:.6f}")
            else:
                selected_threshold = float(inferred_threshold)
                threshold_source = "automatic inference"
                print(f"[saturation] threshold-separation={separation:.6f}")
            write_reads_threshold(threshold_path, selected_threshold)
        reads_threshold = read_reads_threshold(threshold_path)
        print(f"[saturation] reads-threshold={reads_threshold:g}")
        print(f"[saturation] reads-threshold-source={threshold_source}")
        write_reads_rank_plot(
            reads_rank_plot_path,
            spot_reads,
            reads_threshold,
            threshold_source,
            excluded_spots,
        )
        hq_candidates = {
            spot
            for spot, reads in spot_reads.items()
            if reads > reads_threshold and spot not in excluded_spots
        }
        hq_spots = sorted(hq_candidates & spot_histograms.keys())
        sequencing_gbp = load_sequencing_gbp(fastp_path)
        if sequencing_gbp is None:
            print("[saturation] warning=fastp_depth_unavailable")

        if not hq_spots:
            print("[saturation] warning=no_hq_spots_after_filter")
            write_summary(summary_path, empty_row(0, args.pred_fraction))
            write_empty_plot(plot_path, "No HQ spots after reads filter")
            print("[saturation] done")
            return 0

        fraction_values = {fraction: [] for fraction in DEFAULT_FRACTIONS}
        usable_spots: list[str] = []
        for spot in hq_spots:
            histogram = spot_histograms[spot]
            if not histogram:
                continue
            usable_spots.append(spot)
            for fraction in DEFAULT_FRACTIONS:
                fraction_values[fraction].append(expected_unique(histogram, fraction))

        fractions: list[float] = []
        medians: list[float] = []
        errors_low: list[float] = []
        errors_high: list[float] = []
        for fraction in DEFAULT_FRACTIONS:
            values = fraction_values[fraction]
            if values:
                median, low, high = median_and_iqr(values)
                fractions.append(fraction)
                medians.append(median)
                errors_low.append(low)
                errors_high.append(high)

        if not fractions:
            print("[saturation] warning=no_coverage_values_for_hq_spots")
            write_summary(summary_path, empty_row(0, args.pred_fraction))
            write_empty_plot(plot_path, "No valid CG coverage rows")
            print("[saturation] done")
            return 0

        observed = medians[-1]
        maximum, rate = fit_saturation_curve(fractions, medians)
        slope = fit_linear_through_origin(fractions, medians)
        linear_r2 = r_squared(medians, [slope * fraction for fraction in fractions])
        exponential_r2 = r_squared(
            medians,
            [saturation_function(fraction, maximum, rate) for fraction in fractions],
        )
        if linear_r2 >= args.linear_r2_threshold:
            model = "linear"
            theoretical = None
            predicted = slope * args.pred_fraction
            saturation_rate = None
        else:
            model = "saturation"
            theoretical = maximum if maximum > 0 else None
            predicted = (
                saturation_function(args.pred_fraction, maximum, rate)
                if theoretical is not None
                else None
            )
            saturation_rate = (
                observed / theoretical * 100.0
                if theoretical is not None and theoretical > 0
                else None
            )

        print(f"[saturation] linear-r2={linear_r2:.6f}")
        print(f"[saturation] exponential-r2={exponential_r2:.6f}")
        print(f"[saturation] extrapolation-model={model}")
        print(f"[saturation] saturation-rate={format_optional_float(saturation_rate)}")
        write_plot(
            plot_path,
            fractions,
            medians,
            errors_low,
            errors_high,
            model,
            maximum,
            rate,
            slope,
            args.pred_fraction,
            predicted,
            saturation_rate,
            sequencing_gbp,
        )
        write_summary(
            summary_path,
            {
                "observed_median_unique_cpgs": format_optional_int(observed),
                "theoretical_max_median_unique_cpgs": format_optional_int(theoretical),
                "prediction_fraction": format_optional_float(args.pred_fraction),
                "predicted_median_unique_cpgs": format_optional_int(predicted),
                "saturation_rate": format_optional_float(saturation_rate),
                "extrapolation_model": model,
                "hq_spot_count": str(len(usable_spots)),
            },
        )
        print(f"[saturation] hq-spot-count={len(usable_spots)}")
        print("[saturation] done")
        return 0
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"[saturation] error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
