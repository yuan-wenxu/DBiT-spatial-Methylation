#!/usr/bin/env python3
"""Count-preserving quantification of paired TAPS/TAPS-beta data."""

from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CountKey = tuple[str, int]
Aggregate = dict[CountKey, list[int]]
COV_COLUMN_INDEXES = {
    "chrom": 0,
    "position": 1,
    "barcode": 6,
    "methylated": 4,
    "unmethylated": 5,
}


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Coverage-aware paired TAPS/TAPS-beta methylation matrices."
    )
    p.add_argument("--taps-cov", required=True, type=Path)
    p.add_argument("--taps-beta-cov", required=True, type=Path)
    p.add_argument("--spot-map", required=True, type=Path)
    p.add_argument("--bin-size", type=int, default=2000)
    p.add_argument("--min-coverage-per-site", type=int, default=1)
    p.add_argument("--max-coverage-per-site", type=int)
    p.add_argument("--min-coverage-per-bin", type=int, default=1)
    p.add_argument("--min-site-per-bin", type=int, default=1)
    p.add_argument(
        "--min-valid-spots-per-bin",
        type=int,
        default=6,
        help="Minimum number of valid spots required before variability filtering.",
    )
    p.add_argument(
        "--variable-fraction",
        type=float,
        default=0.02,
        help=(
            "Fraction of coverage-eligible bins retained by residual variance; "
            "set to 1 to retain every eligible bin."
        ),
    )
    p.add_argument("--shrinkage-lambda", type=float, default=1.0)
    p.add_argument("--chunksize", type=int, default=1_000_000)
    p.add_argument(
        "--threads",
        type=int,
        default=1,
        help=(
            "Worker processes used after indexing each uncompressed cov file by "
            "chromosome; one thread streams without building an index."
        ),
    )
    return p


def check_args(args: argparse.Namespace) -> None:
    for name in (
        "min_coverage_per_site",
        "min_coverage_per_bin",
        "min_site_per_bin",
        "min_valid_spots_per_bin",
        "chunksize",
        "bin_size",
        "threads",
    ):
        if getattr(args, name) <= 0:
            raise ValueError(f"--{name.replace('_', '-')} must be greater than zero")
    if (
        args.max_coverage_per_site is not None
        and args.max_coverage_per_site < args.min_coverage_per_site
    ):
        raise ValueError(
            "--max-coverage-per-site must be >= --min-coverage-per-site"
        )
    if args.shrinkage_lambda < 0:
        raise ValueError("--shrinkage-lambda must be non-negative")
    if not 0 < args.variable_fraction <= 1:
        raise ValueError("--variable-fraction must be greater than zero and <= 1")
    for path in (args.taps_cov, args.taps_beta_cov, args.spot_map):
        if not path.is_file():
            raise FileNotFoundError(path)


@dataclass
class AssayData:
    chromosome_files: dict[str, Path]
    observed_barcodes: set[str]
    coverage_counts: dict[int, int]


@dataclass
class AssayPartition:
    chromosome: str
    aggregate_file: Path
    observed_barcodes: set[str]
    coverage_counter: CoverageCounter


@dataclass(frozen=True)
class ChromosomeRange:
    chromosome: str
    start: int
    end: int


@dataclass
class ChromosomeAnalysis:
    chromosome: str
    result_file: Path
    taps_coverage_counts: dict[int, int]
    taps_site_counts: dict[int, int]
    beta_coverage_counts: dict[int, int]
    beta_site_counts: dict[int, int]


class CoverageCounter:
    """Count observations for each integer site coverage."""

    def __init__(self) -> None:
        self.counts: dict[int, int] = defaultdict(int)

    def update(self, values: np.ndarray) -> None:
        unique, counts = np.unique(values, return_counts=True)
        for value, count in zip(unique.tolist(), counts.tolist()):
            self.counts[int(value)] += int(count)

    def merge(self, other: CoverageCounter) -> None:
        for coverage, count in other.counts.items():
            self.counts[coverage] += count


def read_spot_map(path: Path) -> pd.DataFrame:
    spot_map = pd.read_csv(path, sep="\t")
    spot_map = spot_map.loc[spot_map["overlap"]].copy()
    taps_spots = (
        spot_map["moving_array_row"].astype(str).str.zfill(2)
        + "_"
        + spot_map["moving_array_col"].astype(str).str.zfill(2)
    )
    beta_spots = (
        spot_map["nearest_taps_beta_array_row"].astype(str).str.zfill(2)
        + "_"
        + spot_map["nearest_taps_beta_array_col"].astype(str).str.zfill(2)
    )
    matched = pd.DataFrame(
        {
            "taps_barcode": taps_spots,
            "taps_beta_barcode": beta_spots,
        }
    ).reset_index(drop=True)
    return matched


class FileSlice(io.RawIOBase):
    """Expose one byte range of a file as a readable binary stream."""

    def __init__(self, path: Path, start: int, end: int) -> None:
        self._handle = path.open("rb")
        self._handle.seek(start)
        self._end = end

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: bytearray) -> int:
        remaining = self._end - self._handle.tell()
        if remaining <= 0:
            return 0
        return self._handle.readinto(memoryview(buffer)[: min(len(buffer), remaining)])

    def close(self) -> None:
        self._handle.close()
        super().close()


def chromosome_filename(chromosome: str) -> str:
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", chromosome).strip("._")
    return f"{safe_name or 'chromosome'}.npz"


def write_chromosome_index(
    path: Path, index_path: Path
) -> list[ChromosomeRange]:
    """Index contiguous chromosome byte ranges in an uncompressed cov file."""
    ranges: list[ChromosomeRange] = []
    completed: set[str] = set()
    current_chromosome: str | None = None
    current_start = 0
    with path.open("rb") as handle:
        while True:
            line_start = handle.tell()
            line = handle.readline()
            if not line:
                file_end = handle.tell()
                break
            if not line.strip() or line.startswith(b"#"):
                continue
            chromosome = line.split(b"\t", 1)[0].decode("utf-8")
            if current_chromosome is None:
                current_chromosome = chromosome
                current_start = line_start
            elif chromosome != current_chromosome:
                ranges.append(
                    ChromosomeRange(current_chromosome, current_start, line_start)
                )
                completed.add(current_chromosome)
                if chromosome in completed:
                    raise ValueError(
                        f"{path}: chromosome {chromosome!r} is not contiguous"
                    )
                current_chromosome = chromosome
                current_start = line_start
    if current_chromosome is not None:
        ranges.append(ChromosomeRange(current_chromosome, current_start, file_end))

    with index_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("chromosome\tstart_byte\tend_byte\n")
        for record in ranges:
            handle.write(
                f"{record.chromosome}\t{record.start}\t{record.end}\n"
            )
    return ranges


def cov_chunks(
    path: Path,
    args: argparse.Namespace,
    byte_range: tuple[int, int] | None = None,
) -> Iterator[pd.DataFrame]:
    usecols = sorted(COV_COLUMN_INDEXES.values())
    reverse_names = {value: key for key, value in COV_COLUMN_INDEXES.items()}
    source: Path | io.BufferedReader
    source_context: io.BufferedReader | None = None
    if byte_range is None:
        source = path
    else:
        source_context = io.BufferedReader(FileSlice(path, *byte_range))
        source = source_context
    try:
        for chunk in pd.read_csv(
            source,
            sep="\t",
            header=None,
            comment="#",
            usecols=usecols,
            dtype=str,
            chunksize=args.chunksize,
        ):
            chunk = chunk.rename(columns=reverse_names)[list(COV_COLUMN_INDEXES)]
            chunk["chrom"] = chunk["chrom"].astype(str)
            chunk["barcode"] = chunk["barcode"].astype(str)
            for column in ("position", "methylated", "unmethylated"):
                chunk[column] = pd.to_numeric(
                    chunk[column], errors="raise"
                ).astype(np.int64)
            chunk["coverage"] = chunk["methylated"] + chunk["unmethylated"]
            yield chunk
    finally:
        if source_context is not None:
            source_context.close()


def update_aggregates(
    aggregates: Aggregate,
    frame: pd.DataFrame,
    args: argparse.Namespace,
    barcode_lookup: dict[str, str] | None,
    valid_barcodes: set[str],
) -> None:
    passing = frame["coverage"].to_numpy() >= args.min_coverage_per_site
    if args.max_coverage_per_site is not None:
        passing &= frame["coverage"].to_numpy() <= args.max_coverage_per_site
    filtered = frame.loc[passing].copy()
    if barcode_lookup is None:
        filtered = filtered[filtered["barcode"].isin(valid_barcodes)]
        filtered["spatial_id"] = filtered["barcode"]
    else:
        filtered["spatial_id"] = filtered["barcode"].map(barcode_lookup)
        filtered = filtered[filtered["spatial_id"].notna()]
    if filtered.empty:
        return
    clear_site = (filtered["methylated"] == 0) | (
        filtered["unmethylated"] == 0
    )
    filtered = filtered.loc[clear_site].copy()
    if filtered.empty:
        return
    filtered["bin_start"] = (
        filtered["position"].to_numpy() // args.bin_size
    ) * args.bin_size
    filtered["methylated_sites"] = (filtered["methylated"] > 0).astype(
        np.int64
    )
    bins_in_chunk = (
        filtered.groupby(
            ["spatial_id", "bin_start"],
            sort=False,
            observed=True,
            as_index=False,
        )
        .agg(
            methylated_sites=("methylated_sites", "sum"),
            coverage=("coverage", "sum"),
            n_sites=("position", "size"),
        )
    )
    for row in bins_in_chunk.itertuples(index=False):
        key = (str(row.spatial_id), int(row.bin_start))
        counts = [
            int(row.methylated_sites),
            int(row.coverage),
            int(row.n_sites),
        ]
        existing = aggregates.get(key)
        if existing is None:
            aggregates[key] = counts
        else:
            for index, value in enumerate(counts):
                existing[index] += value


def write_chromosome_aggregates(path: Path, aggregates: Aggregate) -> None:
    """Write one chromosome's aggregated spot-bin counts."""
    items = sorted(aggregates.items())
    spatial_ids = np.asarray([key[0] for key, _ in items], dtype=str)
    barcodes, barcode_index = np.unique(spatial_ids, return_inverse=True)
    with path.open("wb") as handle:
        np.savez_compressed(
            handle,
            barcodes=barcodes,
            barcode_index=barcode_index.astype(np.int32, copy=False),
            bin_start=np.fromiter(
                (key[1] for key, _ in items),
                dtype=np.int64,
                count=len(items),
            ),
            methylated_sites=np.fromiter(
                (counts[0] for _, counts in items),
                dtype=np.int64,
                count=len(items),
            ),
            coverage=np.fromiter(
                (counts[1] for _, counts in items),
                dtype=np.int64,
                count=len(items),
            ),
            n_sites=np.fromiter(
                (counts[2] for _, counts in items),
                dtype=np.int64,
                count=len(items),
            ),
        )


def read_assay_partition(
    path: Path,
    args: argparse.Namespace,
    barcode_lookup: dict[str, str] | None,
    valid_barcodes: set[str],
    chromosome: str,
    byte_range: tuple[int, int],
    aggregate_file: Path,
) -> AssayPartition:
    aggregates: Aggregate = {}
    observed_barcodes: set[str] = set()
    coverage_counter = CoverageCounter()
    for chunk in cov_chunks(path, args, byte_range):
        chromosomes = chunk["chrom"].unique().tolist()
        if chromosomes != [chromosome]:
            raise ValueError(
                f"{path}: indexed range for {chromosome!r} contains "
                f"chromosomes {chromosomes}"
            )
        observed_barcodes.update(chunk["barcode"].unique().tolist())
        coverage_counter.update(chunk["coverage"].to_numpy())
        update_aggregates(
            aggregates, chunk, args, barcode_lookup, valid_barcodes
        )
    write_chromosome_aggregates(aggregate_file, aggregates)
    return AssayPartition(
        chromosome=chromosome,
        aggregate_file=aggregate_file,
        observed_barcodes=observed_barcodes,
        coverage_counter=coverage_counter,
    )


def read_assay(
    path: Path,
    args: argparse.Namespace,
    barcode_lookup: dict[str, str] | None,
    valid_barcodes: set[str],
    intermediate_dir: Path,
) -> AssayData:
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    if args.threads == 1:
        chromosome_files: dict[str, Path] = {}
        observed_barcodes: set[str] = set()
        coverage_counter = CoverageCounter()
        completed: set[str] = set()
        current_chromosome: str | None = None
        aggregates: Aggregate = {}

        def flush_current() -> None:
            if current_chromosome is None:
                return
            aggregate_file = intermediate_dir / chromosome_filename(
                current_chromosome
            )
            write_chromosome_aggregates(aggregate_file, aggregates)
            chromosome_files[current_chromosome] = aggregate_file

        for chunk in cov_chunks(path, args):
            observed_barcodes.update(chunk["barcode"].unique().tolist())
            coverage_counter.update(chunk["coverage"].to_numpy())
            run_ids = chunk["chrom"].ne(chunk["chrom"].shift()).cumsum()
            for _, chromosome_chunk in chunk.groupby(run_ids, sort=False):
                chromosome = str(chromosome_chunk["chrom"].iloc[0])
                if current_chromosome is None:
                    current_chromosome = chromosome
                elif chromosome != current_chromosome:
                    flush_current()
                    completed.add(current_chromosome)
                    if chromosome in completed:
                        raise ValueError(
                            f"{path}: chromosome {chromosome!r} is not contiguous"
                        )
                    current_chromosome = chromosome
                    aggregates = {}
                update_aggregates(
                    aggregates,
                    chromosome_chunk,
                    args,
                    barcode_lookup,
                    valid_barcodes,
                )
        flush_current()
        return AssayData(
            chromosome_files=chromosome_files,
            observed_barcodes=observed_barcodes,
            coverage_counts=dict(coverage_counter.counts),
        )

    if path.suffix.lower() in {
        ".gz",
        ".gzip",
        ".bgz",
        ".bz2",
        ".xz",
        ".zip",
        ".zst",
    }:
        raise ValueError(
            f"{path}: --threads > 1 requires an uncompressed cov file"
        )
    chromosome_ranges = write_chromosome_index(
        path, intermediate_dir / "chromosome_index.tsv"
    )
    print(
        f"[dbitm] paired-taps: indexed {len(chromosome_ranges)} chromosomes "
        f"in {path}"
    )
    if not chromosome_ranges:
        return AssayData({}, set(), {})
    with ProcessPoolExecutor(
        max_workers=min(args.threads, len(chromosome_ranges))
    ) as executor:
        futures = [
            executor.submit(
                read_assay_partition,
                path,
                args,
                barcode_lookup,
                valid_barcodes,
                record.chromosome,
                (record.start, record.end),
                intermediate_dir / chromosome_filename(record.chromosome),
            )
            for record in chromosome_ranges
        ]
        partitions = [future.result() for future in futures]

    chromosome_files: dict[str, Path] = {}
    observed_barcodes: set[str] = set()
    coverage_counter = CoverageCounter()
    for partition in partitions:
        chromosome_files[partition.chromosome] = partition.aggregate_file
        observed_barcodes.update(partition.observed_barcodes)
        coverage_counter.merge(partition.coverage_counter)
    return AssayData(
        chromosome_files=chromosome_files,
        observed_barcodes=observed_barcodes,
        coverage_counts=dict(coverage_counter.counts),
    )


def natural_key(value: str) -> tuple[Any, ...]:
    return tuple(
        int(part) if part.isdigit() else part for part in re.split(r"(\d+)", value)
    )


def load_chromosome_matrices(
    aggregate_file: Path | None,
    feature_starts: np.ndarray,
    output_row_by_barcode: dict[str, int],
    n_spots: int,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix, sparse.csr_matrix]:
    """Load one chromosome's aggregate file into spot-by-bin matrices."""
    shape = (n_spots, len(feature_starts))
    if aggregate_file is None:
        empty = sparse.csr_matrix(shape, dtype=np.int64)
        return empty, empty.copy(), empty.copy()
    with np.load(aggregate_file) as archive:
        barcodes = archive["barcodes"].astype(str)
        barcode_index = archive["barcode_index"]
        starts = archive["bin_start"]
        methylated_sites = archive["methylated_sites"]
        coverage = archive["coverage"]
        n_sites = archive["n_sites"]
    barcode_rows = np.asarray(
        [output_row_by_barcode.get(barcode, -1) for barcode in barcodes],
        dtype=np.int64,
    )
    rows = barcode_rows[barcode_index]
    columns = np.searchsorted(feature_starts, starts)
    in_layout = columns < len(feature_starts)
    if in_layout.any():
        in_layout[in_layout] &= (
            feature_starts[columns[in_layout]] == starts[in_layout]
        )
    keep = (rows >= 0) & in_layout
    rows = rows[keep]
    columns = columns[keep]

    def make_matrix(values: np.ndarray) -> sparse.csr_matrix:
        result = sparse.csr_matrix(
            (values[keep], (rows, columns)), shape=shape, dtype=np.int64
        )
        result.eliminate_zeros()
        return result

    return (
        make_matrix(methylated_sites),
        make_matrix(coverage),
        make_matrix(n_sites),
    )


def valid_mask(
    read_coverage: sparse.csr_matrix,
    n_sites: sparse.csr_matrix,
    min_coverage: int,
    min_sites: int,
) -> sparse.csr_matrix:
    result = (read_coverage >= min_coverage).multiply(n_sites >= min_sites)
    return result.astype(np.bool_).tocsr()


def ratio_matrix(
    modified: sparse.csr_matrix,
    total: sparse.csr_matrix,
    valid: sparse.csr_matrix,
) -> sparse.csr_matrix:
    rows, columns = valid.nonzero()
    if not rows.size:
        return sparse.csr_matrix(valid.shape, dtype=np.float64)
    modified_values = np.asarray(modified[rows, columns]).ravel().astype(float)
    total_values = np.asarray(total[rows, columns]).ravel().astype(float)
    fractions = modified_values / total_values
    # Preserve explicit zeros at valid spot-bin coordinates. An absent coordinate
    # represents a missing or QC-filtered observation, not a measured zero.
    return sparse.csr_matrix(
        (fractions, (rows, columns)), shape=valid.shape, dtype=np.float64
    )


def paired_values(
    matrix: sparse.csr_matrix, mask: sparse.csr_matrix
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows, columns = mask.nonzero()
    if not rows.size:
        return rows, columns, np.asarray([], dtype=float)
    values = np.asarray(matrix[rows, columns]).ravel().astype(float)
    return rows, columns, values


def baseline(
    modified: sparse.csr_matrix, total: sparse.csr_matrix
) -> np.ndarray:
    numerator = np.asarray(modified.sum(axis=0)).ravel().astype(float)
    denominator = np.asarray(total.sum(axis=0)).ravel().astype(float)
    result = np.full(total.shape[1], np.nan, dtype=float)
    np.divide(numerator, denominator, out=result, where=denominator > 0)
    return result


def residual_matrix(
    fraction: sparse.csr_matrix,
    valid: sparse.csr_matrix,
    baseline_values: np.ndarray,
    weight: sparse.csr_matrix,
) -> sparse.csr_matrix:
    rows, columns, values = paired_values(fraction, valid)
    if not rows.size:
        return sparse.csr_matrix(fraction.shape, dtype=np.float64)
    raw = values - baseline_values[columns]
    weight_values = np.asarray(weight[rows, columns]).ravel().astype(float)
    residuals = raw * weight_values
    # Keep explicit zeros at valid spot-bin coordinates. They represent observed
    # values equal to the baseline, whereas absent coordinates represent missing
    # or QC-filtered observations.
    return sparse.csr_matrix(
        (residuals, (rows, columns)), shape=fraction.shape, dtype=np.float64
    )


def feature_variance(
    residual: sparse.csr_matrix,
    valid: sparse.csr_matrix,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate residual variance and valid-spot counts by feature."""
    valid_spots = np.asarray(valid.sum(axis=0)).ravel().astype(np.int64)
    residual_sum = np.asarray(residual.sum(axis=0)).ravel().astype(float)
    residual_square_sum = np.asarray(residual.power(2).sum(axis=0)).ravel().astype(
        float
    )
    variance = np.full(valid.shape[1], np.nan, dtype=float)
    observed = valid_spots > 0
    variance[observed] = (
        residual_square_sum[observed] / valid_spots[observed]
        - (residual_sum[observed] / valid_spots[observed]) ** 2
    )
    # Avoid tiny negative values caused by floating-point cancellation.
    variance[observed] = np.maximum(variance[observed], 0.0)

    return variance, valid_spots


def select_variable_features(
    variance: np.ndarray,
    eligible: np.ndarray,
    variable_fraction: float,
) -> np.ndarray:
    """Select a global fraction of eligible features by residual variance."""
    eligible_indices = np.flatnonzero(eligible)
    if not eligible_indices.size:
        return np.asarray([], dtype=np.int64)

    selected_count = max(
        1, int(np.ceil(eligible_indices.size * variable_fraction))
    )
    ranked = eligible_indices[
        np.argsort(variance[eligible_indices], kind="stable")
    ]
    selected = np.sort(ranked[-selected_count:]).astype(np.int64, copy=False)
    return selected


def shrinkage_weight(
    n_observations: sparse.csr_matrix, shrinkage: float
) -> sparse.csr_matrix:
    result = n_observations.astype(np.float64).copy()
    result.data = result.data / (result.data + shrinkage)
    return result


def effective_observations(
    taps_sites: sparse.csr_matrix,
    beta_sites: sparse.csr_matrix,
    valid: sparse.csr_matrix,
) -> sparse.csr_matrix:
    rows, columns = valid.nonzero()
    if not rows.size:
        return sparse.csr_matrix(valid.shape, dtype=np.float64)
    taps_values = np.asarray(taps_sites[rows, columns]).ravel().astype(float)
    beta_values = np.asarray(beta_sites[rows, columns]).ravel().astype(float)
    values = taps_values * beta_values / (taps_values + beta_values)
    return sparse.csr_matrix(
        (values, (rows, columns)), shape=valid.shape, dtype=np.float64
    )


def sparse_payload(name: str, matrix: sparse.csr_matrix) -> dict[str, np.ndarray]:
    matrix = matrix.tocsr()
    return {
        f"{name}_data": matrix.data,
        f"{name}_indices": matrix.indices,
        f"{name}_indptr": matrix.indptr,
        f"{name}_shape": np.asarray(matrix.shape, dtype=np.int64),
    }


def sparse_from_archive(archive: Any, name: str) -> sparse.csr_matrix:
    shape = tuple(int(value) for value in archive[f"{name}_shape"])
    return sparse.csr_matrix(
        (
            archive[f"{name}_data"],
            archive[f"{name}_indices"],
            archive[f"{name}_indptr"],
        ),
        shape=shape,
    )


def count_values(values: np.ndarray) -> dict[int, int]:
    if not values.size:
        return {}
    unique, counts = np.unique(values, return_counts=True)
    return {
        int(value): int(count)
        for value, count in zip(unique.tolist(), counts.tolist())
    }


def analyse_chromosome(
    args: argparse.Namespace,
    chromosome: str,
    taps_file: Path | None,
    beta_file: Path,
    spots: list[str],
    result_file: Path,
) -> ChromosomeAnalysis:
    """Calculate and persist all matrices and statistics for one chromosome."""
    with np.load(beta_file) as archive:
        feature_starts = np.unique(archive["bin_start"])
    output_row_by_barcode = {
        barcode: output_row for output_row, barcode in enumerate(spots)
    }
    taps_methylated, taps_coverage, taps_sites = load_chromosome_matrices(
        taps_file, feature_starts, output_row_by_barcode, len(spots)
    )
    beta_methylated, beta_coverage, beta_sites = load_chromosome_matrices(
        beta_file, feature_starts, output_row_by_barcode, len(spots)
    )
    taps_valid = valid_mask(
        taps_coverage,
        taps_sites,
        args.min_coverage_per_bin,
        args.min_site_per_bin,
    )
    beta_valid = valid_mask(
        beta_coverage,
        beta_sites,
        args.min_coverage_per_bin,
        args.min_site_per_bin,
    )
    hmc_valid = taps_valid.multiply(beta_valid).astype(np.bool_).tocsr()
    taps_fraction = ratio_matrix(taps_methylated, taps_sites, taps_valid)
    beta_fraction = ratio_matrix(beta_methylated, beta_sites, beta_valid)
    hmc_rows, hmc_columns, taps_values = paired_values(
        taps_fraction, hmc_valid
    )
    beta_values = np.asarray(
        beta_fraction[hmc_rows, hmc_columns]
    ).ravel().astype(float)
    hmc_fraction = sparse.csr_matrix(
        (taps_values - beta_values, (hmc_rows, hmc_columns)),
        shape=hmc_valid.shape,
        dtype=np.float64,
    )

    taps_mu = baseline(taps_methylated, taps_sites)
    beta_mu = baseline(beta_methylated, beta_sites)
    hmc_mu = taps_mu - beta_mu
    mc_weight = shrinkage_weight(beta_sites, args.shrinkage_lambda).multiply(
        beta_valid
    ).tocsr()
    n_eff = effective_observations(taps_sites, beta_sites, hmc_valid)
    hmc_weight = shrinkage_weight(n_eff, args.shrinkage_lambda)
    mc_residual = residual_matrix(beta_fraction, beta_valid, beta_mu, mc_weight)
    hmc_residual = residual_matrix(hmc_fraction, hmc_valid, hmc_mu, hmc_weight)
    mc_variance, mc_valid_spots = feature_variance(mc_residual, beta_valid)
    hmc_variance, hmc_valid_spots = feature_variance(hmc_residual, hmc_valid)

    payload: dict[str, np.ndarray] = {
        "feature_starts": feature_starts,
        "mc_variance": mc_variance,
        "hmc_variance": hmc_variance,
        "mc_valid_spots": mc_valid_spots,
        "hmc_valid_spots": hmc_valid_spots,
    }
    for name, matrix in (
        ("mc_fraction", beta_fraction),
        ("hmc_fraction", hmc_fraction),
        ("mc_residual", mc_residual),
        ("hmc_residual", hmc_residual),
    ):
        payload.update(sparse_payload(name, matrix))

    with result_file.open("wb") as handle:
        np.savez_compressed(handle, **payload)

    return ChromosomeAnalysis(
        chromosome=chromosome,
        result_file=result_file,
        taps_coverage_counts=count_values(taps_coverage.data),
        taps_site_counts=count_values(taps_sites.data),
        beta_coverage_counts=count_values(beta_coverage.data),
        beta_site_counts=count_values(beta_sites.data),
    )


def write_10x(
    output: Path,
    name: str,
    matrix: sparse.csr_matrix,
    spots: list[str],
    features: list[str],
    feature_type: str,
) -> None:
    """Write a feature-by-spot matrix in 10x-style format."""
    bundle = output / name
    bundle.mkdir(parents=True, exist_ok=True)
    with (bundle / "matrix.mtx.gz").open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, mtime=0
        ) as matrix_handle:
            mmwrite(
                matrix_handle,
                matrix.transpose().tocoo(),
                comment="generated by DBiT-spatial-Methylation",
                field="real",
                symmetry="general",
            )
    with (bundle / "barcodes.tsv.gz").open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, mtime=0
        ) as gzip_handle:
            with io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="") as handle:
                for spot in spots:
                    handle.write(f"{spot}\n")
    with (bundle / "features.tsv.gz").open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, mtime=0
        ) as gzip_handle:
            with io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="") as handle:
                for feature in features:
                    handle.write(f"{feature}\t{feature}\t{feature_type}\n")


def plot_site_coverage(
    output: Path,
    taps_counts: dict[int, int],
    beta_counts: dict[int, int],
) -> None:
    rows = [
        {"assay": assay, "coverage": coverage, "observations": observations}
        for assay, counts in (("TAPS", taps_counts), ("TAPS-beta", beta_counts))
        for coverage, observations in sorted(counts.items())
    ]
    pd.DataFrame(rows).to_csv(
        output / "site_coverage_counts.tsv", sep="\t", index=False
    )
    figure, axes = plt.subplots(1, 2, figsize=(10, 4))
    for axis, title, counts, color in zip(
        axes,
        ("TAPS", "TAPS-beta"),
        (taps_counts, beta_counts),
        ("#4DBBD5", "#E64B35"),
    ):
        coverage = np.arange(1, 21, dtype=int)
        frequency = np.asarray(
            [counts.get(value, 0) for value in coverage[:-1]]
            + [sum(count for value, count in counts.items() if value >= 20)],
            dtype=np.int64,
        )
        tick_labels = [str(value) for value in coverage[:-1]] + ["≥20"]
        axis.bar(coverage, frequency, width=0.8, color=color)
        axis.set_xticks(coverage, tick_labels, fontsize=10)
        axis.set_title(title, fontsize=12)
        axis.set_xlabel("Site coverage", fontsize=10)
        axis.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Number of spot-CpG observations", fontsize=10)
    figure.subplots_adjust(left=0.11, right=0.98, bottom=0.14, top=0.90, wspace=0.08)
    figure.savefig(output / "site_coverage_histogram.png", dpi=300)
    plt.close(figure)


def plot_bin_qc(
    output: Path,
    taps_coverage: dict[int, int],
    taps_sites: dict[int, int],
    beta_coverage: dict[int, int],
    beta_sites: dict[int, int],
    min_coverage: int,
    min_sites: int,
) -> None:
    distributions = (
        ("TAPS", "coverage", taps_coverage),
        ("TAPS", "sites", taps_sites),
        ("TAPS-beta", "coverage", beta_coverage),
        ("TAPS-beta", "sites", beta_sites),
    )
    rows = [
        {
            "assay": assay,
            "metric": metric,
            "value": value,
            "spot_bin_pairs": pairs,
        }
        for assay, metric, counts in distributions
        for value, pairs in sorted(counts.items())
    ]
    pd.DataFrame(
        rows, columns=["assay", "metric", "value", "spot_bin_pairs"]
    ).to_csv(output / "bin_qc_counts.tsv", sep="\t", index=False)

    values = np.arange(1, 21, dtype=int)
    tick_labels = [str(value) for value in values[:-1]] + ["≥20"]

    def grouped_frequencies(counts: dict[int, int]) -> np.ndarray:
        return np.asarray(
            [counts.get(value, 0) for value in values[:-1]]
            + [sum(count for value, count in counts.items() if value >= 20)],
            dtype=np.int64,
        )

    figure, axes = plt.subplots(2, 2, figsize=(10, 8))
    for column, (title, coverage, sites, color) in enumerate(
        (
            ("TAPS", taps_coverage, taps_sites, "#4DBBD5"),
            ("TAPS-beta", beta_coverage, beta_sites, "#E64B35"),
        )
    ):
        axes[0, column].bar(
            values, grouped_frequencies(coverage), width=0.8, color=color
        )
        axes[0, column].axvline(min_coverage, color="black", linestyle="--")
        axes[0, column].set_xticks(values, tick_labels, fontsize=10)
        axes[0, column].set_title(f"{title} coverage", fontsize=12)
        axes[0, column].set_xlabel("Bin coverage", fontsize=10)
        axes[0, column].grid(axis="y", alpha=0.25)

        axes[1, column].bar(
            values, grouped_frequencies(sites), width=0.8, color=color
        )
        axes[1, column].axvline(min_sites, color="black", linestyle="--")
        axes[1, column].set_xticks(values, tick_labels, fontsize=10)
        axes[1, column].set_title(f"{title} sites", fontsize=12)
        axes[1, column].set_xlabel("Sites per bin", fontsize=10)
        axes[1, column].grid(axis="y", alpha=0.25)

    axes[0, 0].set_ylabel("Number of observed spot-bin pairs", fontsize=10)
    axes[1, 0].set_ylabel("Number of observed spot-bin pairs", fontsize=10)
    figure.subplots_adjust(
        left=0.10, right=0.98, bottom=0.08, top=0.94, hspace=0.30, wspace=0.16
    )
    figure.savefig(output / "bin_qc_histogram.png", dpi=300)
    plt.close(figure)


def analyse_assays(
    args: argparse.Namespace,
    output: Path,
    taps_data: AssayData,
    beta_data: AssayData,
    registered_beta: pd.DataFrame,
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    spots = registered_beta["taps_beta_barcode"].tolist()
    chromosomes = sorted(beta_data.chromosome_files, key=natural_key)
    if not chromosomes:
        raise ValueError("no TAPS-beta CpGs remain after site-level QC")
    analysis_dir = output / "intermediate" / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    jobs = [
        (
            args,
            chromosome,
            taps_data.chromosome_files.get(chromosome),
            beta_data.chromosome_files[chromosome],
            spots,
            analysis_dir / chromosome_filename(chromosome),
        )
        for chromosome in chromosomes
    ]
    if args.threads == 1:
        results = [analyse_chromosome(*job) for job in jobs]
    else:
        with ProcessPoolExecutor(
            max_workers=min(args.threads, len(jobs))
        ) as executor:
            futures = [executor.submit(analyse_chromosome, *job) for job in jobs]
            results = [future.result() for future in futures]

    def merged_counts(attribute: str) -> dict[int, int]:
        merged: dict[int, int] = defaultdict(int)
        for result in results:
            for value, count in getattr(result, attribute).items():
                merged[value] += count
        return dict(merged)

    plot_bin_qc(
        output,
        merged_counts("taps_coverage_counts"),
        merged_counts("taps_site_counts"),
        merged_counts("beta_coverage_counts"),
        merged_counts("beta_site_counts"),
        args.min_coverage_per_bin,
        args.min_site_per_bin,
    )

    feature_starts_by_chromosome: list[np.ndarray] = []
    mc_variance_parts: list[np.ndarray] = []
    hmc_variance_parts: list[np.ndarray] = []
    mc_valid_spot_parts: list[np.ndarray] = []
    hmc_valid_spot_parts: list[np.ndarray] = []
    for result in results:
        with np.load(result.result_file) as archive:
            feature_starts_by_chromosome.append(archive["feature_starts"])
            mc_variance_parts.append(archive["mc_variance"])
            hmc_variance_parts.append(archive["hmc_variance"])
            mc_valid_spot_parts.append(archive["mc_valid_spots"])
            hmc_valid_spot_parts.append(archive["hmc_valid_spots"])
    n_features = sum(len(starts) for starts in feature_starts_by_chromosome)
    if not n_features:
        raise ValueError("no TAPS-beta CpGs remain after site-level QC")
    mc_variance = np.concatenate(mc_variance_parts)
    hmc_variance = np.concatenate(hmc_variance_parts)
    mc_valid_spots = np.concatenate(mc_valid_spot_parts)
    hmc_valid_spots = np.concatenate(hmc_valid_spot_parts)
    mc_eligible = mc_valid_spots >= args.min_valid_spots_per_bin
    hmc_eligible = hmc_valid_spots >= args.min_valid_spots_per_bin
    mc_selected = select_variable_features(
        mc_variance, mc_eligible, args.variable_fraction
    )
    hmc_selected = select_variable_features(
        hmc_variance, hmc_eligible, args.variable_fraction
    )

    if not mc_selected.size:
        raise ValueError(
            "no 5mC bins pass --min-valid-spots-per-bin; lower the threshold"
        )

    mc_output = output / "5mc"
    hmc_output = output / "5hmc"
    mc_output.mkdir(parents=True, exist_ok=True)
    hmc_output.mkdir(parents=True, exist_ok=True)
    mc_local_selected: list[np.ndarray] = []
    hmc_local_selected: list[np.ndarray] = []
    offset = 0
    for starts in feature_starts_by_chromosome:
        next_offset = offset + len(starts)
        mc_left, mc_right = np.searchsorted(
            mc_selected, [offset, next_offset]
        )
        hmc_left, hmc_right = np.searchsorted(
            hmc_selected, [offset, next_offset]
        )
        mc_local_selected.append(mc_selected[mc_left:mc_right] - offset)
        hmc_local_selected.append(hmc_selected[hmc_left:hmc_right] - offset)
        offset = next_offset

    mc_fraction_parts: list[sparse.csr_matrix] = []
    hmc_fraction_parts: list[sparse.csr_matrix] = []
    mc_residual_parts: list[sparse.csr_matrix] = []
    hmc_residual_parts: list[sparse.csr_matrix] = []
    mc_features: list[str] = []
    hmc_features: list[str] = []
    for index, result in enumerate(results):
        chromosome = result.chromosome
        starts = feature_starts_by_chromosome[index]
        mc_local = mc_local_selected[index]
        hmc_local = hmc_local_selected[index]
        features = np.asarray(
            [
                f"{chromosome}:{int(start)}-{int(start) + args.bin_size}"
                for start in starts
            ],
            dtype=str,
        )
        mc_features.extend(str(features[local_index]) for local_index in mc_local)
        hmc_features.extend(
            str(features[local_index]) for local_index in hmc_local
        )

        with np.load(result.result_file) as archive:
            if mc_local.size:
                mc_fraction_parts.append(
                    sparse_from_archive(archive, "mc_fraction")[
                        :, mc_local
                    ].tocsr()
                )
                mc_residual_parts.append(
                    sparse_from_archive(archive, "mc_residual")[
                        :, mc_local
                    ].tocsr()
                )
            if hmc_local.size:
                hmc_fraction_parts.append(
                    sparse_from_archive(archive, "hmc_fraction")[
                        :, hmc_local
                    ].tocsr()
                )
                hmc_residual_parts.append(
                    sparse_from_archive(archive, "hmc_residual")[
                        :, hmc_local
                    ].tocsr()
                )

    def combine(parts: list[sparse.csr_matrix]) -> sparse.csr_matrix:
        if not parts:
            return sparse.csr_matrix((len(spots), 0), dtype=np.float64)
        return sparse.hstack(parts, format="csr")

    write_10x(
        mc_output,
        "methylation_fractions",
        combine(mc_fraction_parts),
        spots,
        mc_features,
        "5mC Fraction",
    )
    write_10x(
        mc_output,
        "mean_shrunken_residuals",
        combine(mc_residual_parts),
        spots,
        mc_features,
        "5mC Residual",
    )
    write_10x(
        hmc_output,
        "methylation_fractions",
        combine(hmc_fraction_parts),
        spots,
        hmc_features,
        "5hmC Fraction",
    )
    write_10x(
        hmc_output,
        "mean_shrunken_residuals",
        combine(hmc_residual_parts),
        spots,
        hmc_features,
        "5hmC Residual",
    )


def run(args: argparse.Namespace) -> None:
    check_args(args)
    output_dir = args.spot_map.resolve().parent / "paired-taps"
    if output_dir.exists():
        if not output_dir.is_dir():
            raise ValueError(f"output path is not a directory: {output_dir}")
        if any(output_dir.iterdir()):
            raise ValueError(f"output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    spot_map = read_spot_map(args.spot_map)
    beta_spots = set(spot_map["taps_beta_barcode"])
    taps_to_beta = dict(
        zip(spot_map["taps_barcode"], spot_map["taps_beta_barcode"])
    )
    intermediate = output_dir / "intermediate"

    taps_data = read_assay(
        args.taps_cov,
        args,
        taps_to_beta,
        beta_spots,
        intermediate / "taps",
    )
    beta_data = read_assay(
        args.taps_beta_cov,
        args,
        None,
        beta_spots,
        intermediate / "taps_beta",
    )

    plot_site_coverage(
        output_dir,
        taps_data.coverage_counts,
        beta_data.coverage_counts,
    )

    beta_mask = spot_map["taps_beta_barcode"].isin(beta_data.observed_barcodes)
    registered_beta = spot_map.loc[beta_mask].copy()
    if registered_beta.empty:
        raise ValueError("no spot-map rows have barcodes observed in TAPS-beta")

    analyse_assays(
        args,
        output_dir,
        taps_data,
        beta_data,
        registered_beta,
    )


def main() -> int:
    args = parser().parse_args()
    try:
        run(args)
    except (FileNotFoundError, OSError, ValueError, RuntimeError, ImportError) as error:
        print(f"[dbitm] paired-taps: error: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
