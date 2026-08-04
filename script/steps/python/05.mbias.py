#!/usr/bin/env python3
"""Generate read-cycle M-bias QC from a pooled methylation BAM."""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import sys
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Optional, Tuple

import matplotlib
import pysam

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PAIRED_FLAGS = {83, 99, 147, 163}
TOP_FLAGS = {99, 147}
BOT_FLAGS = {83, 163}

# key: (read, cycle from 5' end), value: [methylated, unmethylated]
MbiasCounts = DefaultDict[Tuple[str, int], list[int]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate aggregate CG M-bias for TAPS/TAPS-v2, EM-seq, or "
            "Cabernet from "
            "a pooled coordinate-sorted BAM."
        )
    )
    parser.add_argument(
        "--assay",
        choices=("taps", "taps-v2", "emseq", "cabernet"),
        required=True,
    )
    parser.add_argument("--bam", required=True, help="Input pooled BAM.")
    parser.add_argument("--reference", required=True, help="Reference FASTA with .fai index.")
    parser.add_argument("--label", required=True, help="Output label, for example host or lambda or puc19.")
    parser.add_argument("--output-dir", required=True, help="Directory for TSV and PNG outputs.")
    parser.add_argument(
        "--subsample-fraction",
        type=float,
        default=0.1,
        help="Deterministic record sampling fraction in (0, 1]. Default: 0.1.",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=10_000_000,
        help="Maximum selected alignment records; 0 means unlimited. Default: 10,000,000.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Sampling seed. Default: 42.")
    parser.add_argument(
        "--max-cycle", type=int, default=150, help="Maximum cycle to report. Default: 150."
    )
    parser.add_argument(
        "--r1-original-length",
        type=int,
        default=150,
        help=(
            "Original R1 length before barcode removal; R1 cycles are shifted "
            "back to these original read coordinates. Default: 150."
        ),
    )
    parser.add_argument(
        "--min-base-quality",
        type=int,
        default=30,
        help="Minimum quality for bases used in a call. Default: 30.",
    )
    parser.add_argument(
        "--min-mapping-quality",
        type=int,
        default=10,
        help="Minimum read mapping quality. Default: 10.",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Validate and print resolved inputs only."
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not 0 < args.subsample_fraction <= 1:
        raise ValueError("subsample-fraction must be in (0, 1]")
    if args.max_records < 0:
        raise ValueError("max-records must be >= 0")
    if args.max_cycle <= 0:
        raise ValueError("max-cycle must be > 0")
    if args.r1_original_length < 0:
        raise ValueError("r1-original-length must be >= 0")
    if args.min_base_quality < 0:
        raise ValueError("min-base-quality must be >= 0")
    if args.min_mapping_quality < 0:
        raise ValueError("min-mapping-quality must be >= 0")
    if not args.label or any(char in args.label for char in "/\\"):
        raise ValueError("label must be a non-empty filename component")

    bam_path = Path(args.bam)
    reference_path = Path(args.reference)
    if not bam_path.is_file():
        raise FileNotFoundError(f"input BAM not found: {bam_path}")
    if not reference_path.is_file():
        raise FileNotFoundError(f"reference FASTA not found: {reference_path}")


def base_quality_passes(
    qualities: Optional[list[int]], query_positions: Tuple[int, ...], minimum: int
) -> bool:
    if minimum == 0:
        return True
    if qualities is None:
        return False
    return all(qualities[position] >= minimum for position in query_positions)


def classify_taps(
    record: pysam.AlignedSegment,
    query_pos: int,
    ref_pos: int,
    reference_sequence: str,
    min_base_quality: int,
) -> Optional[Tuple[bool, int]]:
    """Return TAPS methylation status and cytosine-derived query position."""
    sequence = record.query_sequence
    if not sequence or query_pos + 1 >= len(sequence):
        return None
    if ref_pos < 0 or ref_pos + 1 >= len(reference_sequence):
        return None
    if reference_sequence[ref_pos : ref_pos + 2] != "CG":
        return None
    if not base_quality_passes(
        record.query_qualities, (query_pos, query_pos + 1), min_base_quality
    ):
        return None

    dinucleotide = sequence[query_pos : query_pos + 2].upper()
    if dinucleotide == "TG":
        return True, query_pos
    if dinucleotide == "CA":
        return True, query_pos + 1
    if dinucleotide == "CG":
        if record.flag in TOP_FLAGS:
            return False, query_pos
        if record.flag in BOT_FLAGS:
            return False, query_pos + 1
    return None


def classify_emseq(
    record: pysam.AlignedSegment,
    query_pos: int,
    ref_pos: int,
    reference_sequence: str,
    min_base_quality: int,
) -> Optional[Tuple[bool, int]]:
    """Return EM-seq methylation status and cytosine-derived query position."""
    sequence = record.query_sequence
    if not sequence or query_pos + 1 >= len(sequence):
        return None
    if ref_pos < 0 or ref_pos + 1 >= len(reference_sequence):
        return None
    if reference_sequence[ref_pos : ref_pos + 2] != "CG":
        return None
    if not base_quality_passes(
        record.query_qualities, (query_pos, query_pos + 1), min_base_quality
    ):
        return None

    dinucleotide = sequence[query_pos : query_pos + 2].upper()
    if dinucleotide == "TG":
        return False, query_pos
    if dinucleotide == "CA":
        return False, query_pos + 1
    if dinucleotide == "CG":
        if record.flag in TOP_FLAGS:
            return True, query_pos
        if record.flag in BOT_FLAGS:
            return True, query_pos + 1
    return None


def add_observation(
    counts: MbiasCounts,
    record: pysam.AlignedSegment,
    query_pos: int,
    read_length: int,
    methylated: bool,
    max_cycle: int,
    r1_original_length: int,
) -> None:
    read_label = "R1" if record.is_read1 else "R2"
    cycle_offset = 0
    if record.is_read1 and r1_original_length:
        if read_length > r1_original_length:
            raise ValueError(
                "trimmed R1 length exceeds configured original R1 length: "
                f"{read_length} > {r1_original_length}"
            )
        cycle_offset = r1_original_length - read_length
    if record.is_reverse:
        cycle_from_5p = cycle_offset + read_length - query_pos
    else:
        cycle_from_5p = cycle_offset + query_pos + 1

    value_index = 0 if methylated else 1
    if 1 <= cycle_from_5p <= max_cycle:
        counts[(read_label, cycle_from_5p)][value_index] += 1


def count_mbias(args: argparse.Namespace) -> Tuple[MbiasCounts, int, int, int]:
    counts: MbiasCounts = defaultdict(lambda: [0, 0])
    rng = random.Random(args.seed)
    records_scanned = 0
    records_selected = 0
    observations = 0

    with (
        pysam.AlignmentFile(args.bam, "rb") as bam_file,
        pysam.FastaFile(args.reference) as fasta,
    ):
        current_contig: Optional[str] = None
        reference_sequence = ""
        for record in bam_file.fetch(until_eof=True):
            records_scanned += 1
            if record.is_unmapped or record.is_secondary or record.is_supplementary:
                continue
            if record.flag not in PAIRED_FLAGS:
                continue
            if record.mapping_quality < args.min_mapping_quality:
                continue
            if args.subsample_fraction < 1 and rng.random() >= args.subsample_fraction:
                continue
            if args.max_records and records_selected >= args.max_records:
                break

            sequence = record.query_sequence
            if not sequence or record.reference_name is None:
                continue
            records_selected += 1
            if record.reference_name != current_contig:
                current_contig = record.reference_name
                reference_sequence = fasta.fetch(current_contig).upper()

            for query_pos, ref_pos in record.get_aligned_pairs(matches_only=False):
                if query_pos is None or ref_pos is None:
                    continue
                if args.assay in {"emseq", "cabernet"}:
                    observation = classify_emseq(
                        record,
                        query_pos,
                        ref_pos,
                        reference_sequence,
                        args.min_base_quality,
                    )
                else:
                    observation = classify_taps(
                        record,
                        query_pos,
                        ref_pos,
                        reference_sequence,
                        args.min_base_quality,
                    )
                if observation is None:
                    continue
                methylated, call_query_pos = observation
                add_observation(
                    counts,
                    record,
                    call_query_pos,
                    len(sequence),
                    methylated,
                    args.max_cycle,
                    args.r1_original_length,
                )
                observations += 1

    return counts, records_scanned, records_selected, observations


def write_tsv(path: Path, counts: MbiasCounts) -> None:
    read_order = {"R1": 0, "R2": 1}
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "read",
                "cycle",
                "context",
                "methylated_count",
                "unmethylated_count",
                "coverage",
                "methylation_rate",
            ]
        )
        for (read_label, cycle), values in sorted(
            counts.items(), key=lambda item: (read_order[item[0][0]], item[0][1])
        ):
            methylated, unmethylated = values
            coverage = methylated + unmethylated
            rate = methylated / coverage if coverage else math.nan
            writer.writerow(
                [
                    read_label,
                    cycle,
                    "CG",
                    methylated,
                    unmethylated,
                    coverage,
                    f"{rate:.6f}" if coverage else "NA",
                ]
            )


def write_png(
    path: Path,
    counts: MbiasCounts,
    label: str,
    assay: str,
    max_cycle: int,
    r1_original_length: int,
) -> None:
    colors = {"R1": "#FEA830", "R2": "#36B7C9"}
    fig, axis = plt.subplots(figsize=(5, 4), dpi=300)
    for read_label in ("R1", "R2"):
        points = []
        for cycle in range(1, max_cycle + 1):
            methylated, unmethylated = counts.get((read_label, cycle), (0, 0))
            coverage = methylated + unmethylated
            if coverage:
                points.append((cycle, 100 * methylated / coverage))
        if points:
            axis.plot(
                [point[0] for point in points],
                [point[1] for point in points],
                color=colors[read_label],
                label=read_label,
                linewidth=1.5,
            )
    axis.set_xlabel(
        "Cycle from original 5' end"
        if r1_original_length
        else "Cycle from trimmed 5' end",
        fontsize=10
    )
    axis.set_xlim(1, max_cycle)
    axis.set_ylim(0, 100)
    axis.set_ylabel("Methylation rate (%)", fontsize=10)
    axis.grid(alpha=0.25)
    if axis.lines:
        axis.legend(loc="best")
    fig.suptitle(f" M-bias: {label} ({assay})", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, format="png")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    validate_args(args)
    output_dir = Path(args.output_dir)
    output_tsv = output_dir / f"{args.label}.mbias.tsv"
    output_png = output_dir / f"{args.label}.mbias.png"

    print(f"[mbias] assay={args.assay}")
    print(f"[mbias] label={args.label}")
    print(f"[mbias] bam={args.bam}")
    print(f"[mbias] reference={args.reference}")
    print(f"[mbias] subsample-fraction={args.subsample_fraction}")
    print(f"[mbias] max-records={args.max_records or 'unlimited'}")
    print(f"[mbias] max-cycle={args.max_cycle}")
    print(f"[mbias] r1-original-length={args.r1_original_length or 'trimmed-relative'}")
    print(f"[mbias] output-tsv={output_tsv}")
    print(f"[mbias] output-png={output_png}")
    if args.dry_run:
        print("[mbias] dry-run=1")
        return 0

    if output_tsv.exists() and output_png.exists():
        print(f"[mbias] skip-existing-outputs={args.label}")
        return 0
    if output_tsv.exists() or output_png.exists():
        raise FileExistsError(
            f"[mbias] partial M-bias output exists for {args.label}; remove or rename it before retrying"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    counts, scanned, selected, observations = count_mbias(args)
    tsv_temporary = output_tsv.with_name(f".{output_tsv.name}.tmp.{os.getpid()}")
    png_temporary = output_png.with_name(f".{output_png.name}.tmp.{os.getpid()}")
    try:
        write_tsv(tsv_temporary, counts)
        write_png(
            png_temporary,
            counts,
            args.label,
            args.assay,
            args.max_cycle,
            args.r1_original_length,
        )
        os.replace(tsv_temporary, output_tsv)
        os.replace(png_temporary, output_png)
    finally:
        tsv_temporary.unlink(missing_ok=True)
        png_temporary.unlink(missing_ok=True)

    print(f"[mbias] records-scanned={scanned}")
    print(f"[mbias] records-selected={selected}")
    print(f"[mbias] cpg-observations={observations}")
    print("[mbias] done")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, FileExistsError, ValueError) as error:
        print(f"[mbias] error: {error}", file=sys.stderr)
        raise SystemExit(1) from error
