#!/usr/bin/env python3
"""Standalone filter for SmC alignments with three successive methylated C calls."""

from __future__ import annotations

import argparse
import math
import multiprocessing
import os
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import pysam

from bam_intervals import interval_records, reference_offsets


TOP_FLAGS = {99, 147}
BOT_FLAGS = {83, 163}
PAIRED_FLAGS = TOP_FLAGS | BOT_FLAGS


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply the SmC-seq Bismark-style XXX non-conversion filter to "
            "a coordinate-sorted BAM and index the filtered output."
        )
    )
    parser.add_argument("--bam", required=True, help="Input BAM.")
    parser.add_argument("--reference", required=True, help="Reference FASTA.")
    parser.add_argument("--output-bam", required=True, help="Filtered BAM output.")
    parser.add_argument("--r1-left-trim", type=int, default=0)
    parser.add_argument("--r1-right-trim", type=int, default=0)
    parser.add_argument("--r2-left-trim", type=int, default=0)
    parser.add_argument("--r2-right-trim", type=int, default=0)
    parser.add_argument(
        "--chromosomes",
        default="",
        help=(
            "Comma-separated chromosomes to apply the filter to. Records on "
            "other chromosomes are copied through unchanged. Default: filter "
            "every chromosome."
        ),
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="BAM coordinate-interval filter processes. Default: 1.",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def is_trimmed(
    record: pysam.AlignedSegment,
    query_pos: int,
    read_len: int,
    r1_left: int,
    r1_right: int,
    r2_left: int,
    r2_right: int,
) -> bool:
    """Return whether a query position is outside the retained cycle range."""
    if record.is_read1:
        left_trim, right_trim = r1_left, r1_right
    elif record.is_read2:
        left_trim, right_trim = r2_left, r2_right
    else:
        return False
    if record.is_reverse:
        left_cycle = read_len - query_pos
        right_cycle = query_pos + 1
    else:
        left_cycle = query_pos + 1
        right_cycle = read_len - query_pos
    return left_cycle <= left_trim or right_cycle <= right_trim


def has_unclear_conversion_strand(record: pysam.AlignedSegment) -> bool:
    """Return whether BISCUIT marked the conversion strand as unknown."""
    try:
        return record.get_tag("YD") == "u"
    except KeyError:
        return False


def has_xxx_pattern(
    record: pysam.AlignedSegment,
    reference_sequence: str,
    r1_left: int,
    r1_right: int,
    r2_left: int,
    r2_right: int,
    reference_start: int = 0,
) -> bool:
    """Return whether a read has three successive methylated C calls."""
    sequence = record.query_sequence
    if not sequence or record.flag not in PAIRED_FLAGS:
        return False

    query_sequence = sequence.upper()
    read_length = len(query_sequence)
    top_strand = record.flag in TOP_FLAGS
    methylated_run = 0
    for query_pos, reference_pos in record.get_aligned_pairs(matches_only=True):
        if query_pos is None or reference_pos is None:
            continue
        reference_offset = reference_pos - reference_start
        if reference_offset < 0 or reference_offset >= len(reference_sequence):
            continue
        if is_trimmed(
            record,
            query_pos,
            read_length,
            r1_left,
            r1_right,
            r2_left,
            r2_right,
        ):
            continue

        reference_base = reference_sequence[reference_offset]
        observed_base = query_sequence[query_pos]
        if top_strand:
            if reference_base != "C" or observed_base not in "CT":
                continue
            methylated = observed_base == "C"
        else:
            if reference_base != "G" or observed_base not in "GA":
                continue
            methylated = observed_base == "G"

        methylated_run = methylated_run + 1 if methylated else 0
        if methylated_run >= 3:
            return True
    return False


def selected_chromosomes(raw_list: str) -> Optional[set[str]]:
    """Return the chromosome set to filter, or None when everything filters."""
    names = [name.strip() for name in raw_list.split(",") if name.strip()]
    return set(names) if names else None


def validate_args(args: argparse.Namespace) -> None:
    if args.jobs <= 0:
        raise ValueError("--jobs must be greater than zero")
    for name in (
        "r1_left_trim",
        "r1_right_trim",
        "r2_left_trim",
        "r2_right_trim",
    ):
        if getattr(args, name) < 0:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 0")

    bam_path = Path(args.bam).resolve()
    reference_path = Path(args.reference).resolve()
    output_path = Path(args.output_bam).resolve()
    if bam_path == output_path:
        raise ValueError("--output-bam must differ from --bam")
    if not bam_path.is_file():
        raise ValueError(f"input BAM not found: {bam_path}")
    bam_indexes = (
        Path(str(bam_path) + ".bai"),
        bam_path.with_suffix(".bai"),
    )
    if not any(path.is_file() for path in bam_indexes):
        raise ValueError(f"BAI index not found for: {bam_path}")
    if not reference_path.is_file():
        raise ValueError(f"reference FASTA not found: {reference_path}")
    if not Path(str(reference_path) + ".fai").is_file():
        raise ValueError(f"reference FASTA index not found: {reference_path}.fai")


PartitionResult = tuple[int, str, int, int, int]


def filter_partition(
    partition_index: int,
    chromosome: Optional[str],
    start: int,
    end: int,
    evaluate: bool,
    input_path: str,
    reference_path: str,
    part_path: str,
    r1_left: int,
    r1_right: int,
    r2_left: int,
    r2_right: int,
    first_record_offset: Optional[int] = None,
) -> PartitionResult:
    """Filter records starting in one coordinate interval into a BAM part."""
    total = 0
    evaluated = 0
    filtered = 0
    label = f"{chromosome}:{start}-{end}" if chromosome is not None else "unplaced"

    with pysam.AlignmentFile(input_path, "rb") as input_bam, pysam.AlignmentFile(
        part_path, "wb", template=input_bam
    ) as output_bam:
        records = (
            input_bam.fetch("*")
            if chromosome is None
            else interval_records(input_bam, chromosome, start, end, first_record_offset)
        )
        if evaluate and chromosome is not None:
            with pysam.FastaFile(reference_path) as fasta:
                reference_sequence = fasta.fetch(chromosome, start, end).upper()
                for record in records:
                    total += 1
                    should_evaluate = (
                        not record.is_unmapped
                        and record.flag in PAIRED_FLAGS
                        and not has_unclear_conversion_strand(record)
                    )
                    if should_evaluate:
                        evaluated += 1
                        alignment_reference = reference_sequence
                        alignment_reference_start = start
                        if (
                            record.reference_end is not None
                            and record.reference_end > end
                        ):
                            alignment_reference_start = record.reference_start
                            alignment_reference = fasta.fetch(
                                chromosome,
                                alignment_reference_start,
                                record.reference_end,
                            ).upper()
                        if has_xxx_pattern(
                            record,
                            alignment_reference,
                            r1_left,
                            r1_right,
                            r2_left,
                            r2_right,
                            reference_start=alignment_reference_start,
                        ):
                            filtered += 1
                            continue
                    output_bam.write(record)
        else:
            for record in records:
                total += 1
                output_bam.write(record)
    return partition_index, label, total, evaluated, filtered


def filter_bam(args: argparse.Namespace) -> tuple[int, int, int]:
    """Filter BAM partitions in parallel and atomically replace the output."""
    input_path = Path(args.bam).resolve()
    output_path = Path(args.output_bam).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(
        f".{output_path.name}.tmp.{os.getpid()}.bam"
    )
    temporary_index = Path(str(temporary_path) + ".bai")
    output_index = Path(str(output_path) + ".bai")
    part_dir = Path(
        tempfile.mkdtemp(
            prefix=f".{output_path.name}.parts.", dir=output_path.parent
        )
    )

    try:
        selected = selected_chromosomes(args.chromosomes)
        chromosome_offsets = reference_offsets(input_path)
        with pysam.AlignmentFile(input_path, "rb") as input_bam, pysam.FastaFile(
            args.reference
        ) as fasta:
            reference_names = set(fasta.references)
            indexed_counts = {
                stats.contig: stats.mapped + stats.unmapped
                for stats in input_bam.get_index_statistics()
            }
            chromosomes = [
                chromosome
                for chromosome in input_bam.references
                if indexed_counts.get(chromosome, 0) > 0
            ]
            chromosome_lengths = dict(zip(input_bam.references, input_bam.lengths))
            missing = [
                chromosome
                for chromosome in chromosomes
                if chromosome not in reference_names
            ]
            if missing:
                raise ValueError(
                    "BAM chromosomes absent from reference: " + ",".join(missing)
                )
            unknown = (
                sorted(selected - reference_names) if selected is not None else []
            )
            if unknown:
                raise ValueError(
                    "--chromosomes absent from reference: " + ",".join(unknown)
                )

        total_bases = sum(chromosome_lengths[name] for name in chromosomes)
        target_bases = (
            max(1, math.ceil(total_bases / (args.jobs * 4)))
            if args.jobs > 1
            else max(chromosome_lengths.values(), default=1)
        )
        intervals = [
            (
                chromosome,
                start,
                min(start + target_bases, chromosome_lengths[chromosome]),
            )
            for chromosome in chromosomes
            for start in range(0, chromosome_lengths[chromosome], target_bases)
        ]
        partition_specs: list[tuple[int, Optional[str], int, int, bool, Path]] = [
            (
                index,
                chromosome,
                start,
                end,
                selected is None or chromosome in selected,
                part_dir / f"{index:06d}.bam",
            )
            for index, (chromosome, start, end) in enumerate(intervals)
        ]
        partition_specs.append(
            (
                len(partition_specs),
                None,
                0,
                0,
                False,
                part_dir / f"{len(partition_specs):06d}.bam",
            )
        )
        worker_args = (
            str(input_path),
            str(Path(args.reference).resolve()),
            args.r1_left_trim,
            args.r1_right_trim,
            args.r2_left_trim,
            args.r2_right_trim,
        )
        results: list[PartitionResult] = []
        max_workers = min(args.jobs, len(partition_specs))
        if max_workers == 1:
            for (
                partition_index,
                chromosome,
                start,
                end,
                evaluate,
                part_path,
            ) in partition_specs:
                result = filter_partition(
                    partition_index,
                    chromosome,
                    start,
                    end,
                    evaluate,
                    worker_args[0],
                    worker_args[1],
                    str(part_path),
                    *worker_args[2:],
                    first_record_offset=(
                        chromosome_offsets[chromosome]
                        if chromosome is not None and start == 0
                        else None
                    ),
                )
                results.append(result)
                print(
                    f"[smc-xxx-filter] partition={result[1]} "
                    f"records={result[2]} evaluated={result[3]} "
                    f"filtered={result[4]}",
                    flush=True,
                )
        else:
            fork_context = multiprocessing.get_context("fork")
            with ProcessPoolExecutor(
                max_workers=max_workers, mp_context=fork_context
            ) as executor:
                future_to_label = {
                    executor.submit(
                        filter_partition,
                        partition_index,
                        chromosome,
                        start,
                        end,
                        evaluate,
                        worker_args[0],
                        worker_args[1],
                        str(part_path),
                        *worker_args[2:],
                        first_record_offset=(
                            chromosome_offsets[chromosome]
                            if chromosome is not None and start == 0
                            else None
                        ),
                    ): (
                        f"{chromosome}:{start}-{end}"
                        if chromosome is not None
                        else "unplaced"
                    )
                    for (
                        partition_index,
                        chromosome,
                        start,
                        end,
                        evaluate,
                        part_path,
                    ) in partition_specs
                }
                for future in as_completed(future_to_label):
                    label = future_to_label[future]
                    try:
                        result = future.result()
                    except Exception as error:
                        raise RuntimeError(
                            f"partition {label} failed: {error}"
                        ) from error
                    results.append(result)
                    print(
                        f"[smc-xxx-filter] partition={result[1]} "
                        f"records={result[2]} evaluated={result[3]} "
                        f"filtered={result[4]}",
                        flush=True,
                    )

        ordered_parts = [
            str(partition_specs[result[0]][5]) for result in sorted(results)
        ]
        pysam.cat(
            "--no-PG",
            "-o",
            str(temporary_path),
            *ordered_parts,
            catch_stdout=False,
        )
        pysam.index("-@", str(args.jobs), str(temporary_path))
        os.replace(temporary_path, output_path)
        os.replace(temporary_index, output_index)
    finally:
        temporary_path.unlink(missing_ok=True)
        temporary_index.unlink(missing_ok=True)
        shutil.rmtree(part_dir, ignore_errors=True)
    total = sum(result[2] for result in results)
    evaluated = sum(result[3] for result in results)
    filtered = sum(result[4] for result in results)
    return total, evaluated, filtered


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
    except ValueError as error:
        print(f"[smc-xxx-filter] error: {error}", file=sys.stderr)
        return 1

    print(f"[smc-xxx-filter] bam={Path(args.bam).resolve()}")
    print(f"[smc-xxx-filter] reference={Path(args.reference).resolve()}")
    print(f"[smc-xxx-filter] output-bam={Path(args.output_bam).resolve()}")
    print(
        f"[smc-xxx-filter] trim R1={args.r1_left_trim},"
        f"{args.r1_right_trim} R2={args.r2_left_trim},{args.r2_right_trim}"
    )
    print(f"[smc-xxx-filter] jobs={args.jobs}")
    selected = selected_chromosomes(args.chromosomes)
    if selected is None:
        print("[smc-xxx-filter] chromosomes=all")
    else:
        print("[smc-xxx-filter] chromosomes=" + ",".join(sorted(selected)))
    if args.dry_run:
        print("[smc-xxx-filter] dry-run=1")
        return 0

    try:
        total, evaluated, filtered = filter_bam(args)
    except (OSError, RuntimeError, ValueError, pysam.SamtoolsError) as error:
        print(f"[smc-xxx-filter] error: {error}", file=sys.stderr)
        return 1
    filtered_pct = 100 * filtered / evaluated if evaluated else 0.0
    print(
        f"[smc-xxx-filter] done records={total} evaluated={evaluated} "
        f"filtered={filtered} filtered-percent={filtered_pct:.2f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
