#!/usr/bin/env python3
"""Infer read-end trimming from an M-bias TSV."""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class CyclePoint:
    rate: float
    coverage: int


@dataclass(frozen=True)
class Cutoff:
    read: str
    read_length: int
    axis_start: int
    axis_end: int
    stable_start: int
    stable_end: int
    left_trim: int
    right_trim: int
    baseline_rate: float
    coverage_threshold: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Infer R1/R2 trimming directly from cycle-level methylation rates."
        )
    )
    parser.add_argument("--mbias-tsv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--r1-original-length", type=int, default=150)
    parser.add_argument(
        "--rate-tolerance",
        type=float,
        default=0.05,
        help="Maximum absolute methylation-rate deviation from baseline. Default: 0.05.",
    )
    parser.add_argument(
        "--coverage-fraction",
        type=float,
        default=0.10,
        help="Minimum coverage as a fraction of central-cycle median coverage. Default: 0.10.",
    )
    parser.add_argument("--min-coverage", type=int, default=100)
    parser.add_argument("--smoothing-window", type=int, default=5)
    parser.add_argument("--stable-cycles", type=int, default=5)
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.r1_original_length < 0:
        raise ValueError("r1-original-length must be >= 0")
    if not 0 < args.rate_tolerance < 1:
        raise ValueError("rate-tolerance must be in (0, 1)")
    if not 0 < args.coverage_fraction <= 1:
        raise ValueError("coverage-fraction must be in (0, 1]")
    if args.min_coverage <= 0:
        raise ValueError("min-coverage must be > 0")
    if args.smoothing_window <= 0 or args.smoothing_window % 2 == 0:
        raise ValueError("smoothing-window must be a positive odd integer")
    if args.stable_cycles <= 0:
        raise ValueError("stable-cycles must be > 0")
    if not Path(args.mbias_tsv).is_file():
        raise FileNotFoundError(args.mbias_tsv)


def load_mbias(path: str) -> dict[str, dict[int, CyclePoint]]:
    points: dict[str, dict[int, CyclePoint]] = {"R1": {}, "R2": {}}
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"read", "cycle", "coverage", "methylation_rate"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"M-bias TSV lacks required columns: {path}")
        for row in reader:
            read = row["read"]
            if read not in points or row["methylation_rate"] == "NA":
                continue
            cycle = int(row["cycle"])
            points[read][cycle] = CyclePoint(
                rate=float(row["methylation_rate"]),
                coverage=int(row["coverage"]),
            )
    return points


def find_run(stable: dict[int, bool], start: int, end: int, length: int) -> int:
    run_start: Optional[int] = None
    run_length = 0
    for cycle in range(start, end + 1):
        if stable.get(cycle, False):
            if run_start is None:
                run_start = cycle
            run_length += 1
            if run_length >= length:
                return run_start
        else:
            run_start = None
            run_length = 0
    raise ValueError(f"no run of {length} stable cycles between {start} and {end}")


def infer_cutoff(
    read: str,
    points: dict[int, CyclePoint],
    r1_original_length: int,
    rate_tolerance: float,
    coverage_fraction: float,
    min_coverage: int,
    smoothing_window: int,
    stable_cycles: int,
) -> Cutoff:
    covered_points = [
        point for point in points.values() if point.coverage >= min_coverage
    ]
    if not covered_points:
        raise ValueError(f"no well-covered cycles available for {read}")
    preliminary_coverage = max(
        min_coverage,
        round(
            statistics.median(point.coverage for point in covered_points)
            * coverage_fraction
        ),
    )
    supported_cycles = [
        cycle
        for cycle, point in points.items()
        if point.coverage >= preliminary_coverage
    ]
    if not supported_cycles:
        raise ValueError(f"no coverage-supported cycle range available for {read}")

    if read == "R1" and r1_original_length:
        axis_start = min(supported_cycles)
        axis_end = r1_original_length
        if axis_start > axis_end:
            raise ValueError(
                f"R1 cycle range starts after original length: {axis_start} > {axis_end}"
            )
    else:
        axis_start = 1
        axis_end = max(points)
    read_length = axis_end - axis_start + 1

    central_start = axis_start + read_length // 4
    central_end = axis_end - read_length // 4
    central_points = [
        point
        for cycle, point in points.items()
        if central_start <= cycle <= central_end and point.coverage >= min_coverage
    ]
    if len(central_points) < stable_cycles:
        raise ValueError(
            f"insufficient well-covered central cycles for {read}: {len(central_points)}"
        )
    baseline_rate = statistics.median(point.rate for point in central_points)
    median_coverage = statistics.median(point.coverage for point in central_points)
    coverage_threshold = max(min_coverage, round(median_coverage * coverage_fraction))

    half_window = smoothing_window // 2
    stable: dict[int, bool] = {}
    for cycle in range(axis_start, axis_end + 1):
        point = points.get(cycle)
        if point is None or point.coverage < coverage_threshold:
            stable[cycle] = False
            continue
        local_rates = [
            local_point.rate
            for local_cycle in range(cycle - half_window, cycle + half_window + 1)
            if (local_point := points.get(local_cycle)) is not None
            and local_point.coverage >= coverage_threshold
        ]
        stable[cycle] = (
            len(local_rates) >= half_window + 1
            and abs(statistics.median(local_rates) - baseline_rate) <= rate_tolerance
        )

    stable_start = find_run(stable, axis_start, axis_end, stable_cycles)
    reversed_stable = {
        axis_start + axis_end - cycle: value for cycle, value in stable.items()
    }
    reversed_start = find_run(
        reversed_stable, axis_start, axis_end, stable_cycles
    )
    stable_end = axis_start + axis_end - reversed_start
    if stable_start > stable_end:
        raise ValueError(f"stable {read} boundaries cross: {stable_start} > {stable_end}")

    return Cutoff(
        read=read,
        read_length=read_length,
        axis_start=axis_start,
        axis_end=axis_end,
        stable_start=stable_start,
        stable_end=stable_end,
        left_trim=stable_start - axis_start,
        right_trim=axis_end - stable_end,
        baseline_rate=baseline_rate,
        coverage_threshold=coverage_threshold,
    )


def write_cutoffs(path: Path, cutoffs: list[Cutoff]) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "read",
                    "read_length",
                    "axis_start",
                    "axis_end",
                    "stable_start",
                    "stable_end",
                    "left_trim",
                    "right_trim",
                    "baseline_methylation_rate",
                    "coverage_threshold",
                ]
            )
            for cutoff in cutoffs:
                writer.writerow(
                    [
                        cutoff.read,
                        cutoff.read_length,
                        cutoff.axis_start,
                        cutoff.axis_end,
                        cutoff.stable_start,
                        cutoff.stable_end,
                        cutoff.left_trim,
                        cutoff.right_trim,
                        f"{cutoff.baseline_rate:.6f}",
                        cutoff.coverage_threshold,
                    ]
                )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_empty(path: Path) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        temporary.touch()
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> int:
    args = parse_args()
    validate_args(args)
    points = load_mbias(args.mbias_tsv)
    reads_without_data = [read for read in ("R1", "R2") if not points[read]]
    if reads_without_data:
        output = Path(args.output)
        write_empty(output)
        print(f"[mbias-cutoff] no usable cycles for: {', '.join(reads_without_data)}")
        print(f"[mbias-cutoff] output={output}")
        return 0
    cutoffs: list[Cutoff] = []
    failed_reads: list[str] = []
    for read in ("R1", "R2"):
        try:
            cutoffs.append(
                infer_cutoff(
                    read,
                    points[read],
                    args.r1_original_length,
                    args.rate_tolerance,
                    args.coverage_fraction,
                    args.min_coverage,
                    args.smoothing_window,
                    args.stable_cycles,
                )
            )
        except ValueError as error:
            failed_reads.append(read)
            print(f"[mbias-cutoff] warning: {error}")
    if failed_reads:
        output = Path(args.output)
        write_empty(output)
        print(f"[mbias-cutoff] no suitable cutoff for: {', '.join(failed_reads)}")
        print(f"[mbias-cutoff] output={output} (empty; no trimming will be applied)")
        return 0
    output = Path(args.output)
    write_cutoffs(output, cutoffs)
    for cutoff in cutoffs:
        print(
            f"[mbias-cutoff] read={cutoff.read} length={cutoff.read_length} "
            f"stable-cycles={cutoff.stable_start}-{cutoff.stable_end} "
            f"left-trim={cutoff.left_trim} right-trim={cutoff.right_trim} "
            f"baseline={cutoff.baseline_rate:.4f}"
        )
    print(f"[mbias-cutoff] output={output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, OSError, ValueError) as error:
        print(f"[mbias-cutoff] error: {error}", file=sys.stderr)
        raise SystemExit(1) from error
