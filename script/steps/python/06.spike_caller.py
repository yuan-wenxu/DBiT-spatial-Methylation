#!/usr/bin/env python3
"""Aggregate CpG/CH methylation calling for mitochondrial and spike-in BAMs."""

from __future__ import annotations

import argparse
import importlib.util
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from types import ModuleType
from typing import Optional

import pysam


def load_shared_caller() -> ModuleType:
    caller_path = Path(__file__).with_name("06.methy_caller.py")
    spec = importlib.util.spec_from_file_location("dbitm_methy_caller", caller_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load shared caller: {caller_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SHARED = load_shared_caller()


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate methylation caller for DBiT mito/spike-in BAMs."
    )
    parser.add_argument(
        "--assay",
        choices=("taps", "taps-v2", "emseq", "cabernet"),
        required=True,
    )
    parser.add_argument("--bam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-name", required=True)
    parser.add_argument("--chromosomes", required=True)
    parser.add_argument(
        "--context-mode", choices=("cg", "ch", "both"), default="cg"
    )
    parser.add_argument(
        "--cg-strand-mode", choices=("separate", "merged"), default="separate"
    )
    parser.add_argument("--min-base-quality", type=int, default=30)
    parser.add_argument("--min-mapping-quality", type=int, default=10)
    parser.add_argument("--max-depth", type=int, default=50_000)
    parser.add_argument("--r1-left-trim", type=int, default=0)
    parser.add_argument("--r1-right-trim", type=int, default=0)
    parser.add_argument("--r2-left-trim", type=int, default=0)
    parser.add_argument("--r2-right-trim", type=int, default=0)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> list[str]:
    if args.assay == "cabernet":
        args.cg_strand_mode = "separate"
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", args.output_name) is None:
        raise ValueError("--output-name must be filename-safe")
    chromosomes = SHARED.parse_chromosomes(args.chromosomes)
    if not chromosomes:
        raise ValueError("--chromosomes must not be empty")
    for name in (
        "min_base_quality",
        "min_mapping_quality",
        "r1_left_trim",
        "r1_right_trim",
        "r2_left_trim",
        "r2_right_trim",
    ):
        if getattr(args, name) < 0:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 0")
    if args.max_depth <= 0:
        raise ValueError("--max-depth must be > 0")
    return chromosomes


def bam_index_exists(path: Path) -> bool:
    candidates = [
        Path(str(path) + ".bai"),
        path.with_suffix(".bai"),
        Path(str(path) + ".csi"),
    ]
    return any(candidate.is_file() for candidate in candidates)


def count_cg_column(pileup_col, args: argparse.Namespace, counts) -> None:
    pos = pileup_col.pos
    for pileup_read in pileup_col.pileups:
        record = pileup_read.alignment
        if record.flag not in SHARED.PAIRED_FLAGS:
            continue
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        query_pos = pileup_read.query_position
        if query_pos is None:
            continue
        sequence = record.query_sequence
        if not sequence or query_pos + 1 >= len(sequence):
            continue
        next_pos = query_pos + 1
        if SHARED.is_trimmed(
            record,
            query_pos,
            len(sequence),
            args.r1_left_trim,
            args.r1_right_trim,
            args.r2_left_trim,
            args.r2_right_trim,
        ) or SHARED.is_trimmed(
            record,
            next_pos,
            len(sequence),
            args.r1_left_trim,
            args.r1_right_trim,
            args.r2_left_trim,
            args.r2_right_trim,
        ):
            continue
        observation = SHARED.classify_cg_observation(
            sequence[query_pos : next_pos + 1].upper(), record.flag, args.assay
        )
        if observation is None:
            continue
        methylated, strand = observation
        output_pos = pos if strand == "+" or args.cg_strand_mode == "merged" else pos + 1
        output_strand = strand if args.cg_strand_mode == "separate" else "."
        values = counts[(output_pos, output_strand)]
        values[0 if methylated else 1] += 1


def count_ch_column(
    pileup_col,
    context: str,
    strand: str,
    args: argparse.Namespace,
    counts,
) -> None:
    pos = pileup_col.pos
    for pileup_read in pileup_col.pileups:
        record = pileup_read.alignment
        if record.flag not in SHARED.PAIRED_FLAGS:
            continue
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        query_pos = pileup_read.query_position
        if query_pos is None:
            continue
        forward_sequence = record.get_forward_sequence()
        if not forward_sequence:
            continue
        read_length = len(forward_sequence)
        forward_index = SHARED.get_forward_index(record, query_pos, read_length)
        if forward_index < 0 or forward_index >= read_length:
            continue
        if SHARED.is_trimmed(
            record,
            query_pos,
            read_length,
            args.r1_left_trim,
            args.r1_right_trim,
            args.r2_left_trim,
            args.r2_right_trim,
        ):
            continue
        base = forward_sequence[forward_index].upper()
        if strand == "+":
            if forward_index + 1 >= read_length:
                continue
            if forward_sequence[forward_index + 1].upper() != context[1]:
                continue
        observation = SHARED.classify_ch_observation(base, strand, args.assay)
        if observation is None:
            continue
        methylated, unmethylated = observation
        values = counts[(pos, context, strand)]
        values[0] += methylated
        values[1] += unmethylated


def call_chromosome(
    bam: pysam.AlignmentFile,
    fasta: pysam.FastaFile,
    chromosome: str,
    args: argparse.Namespace,
):
    sequence = fasta.fetch(chromosome).upper()
    cg_counts = defaultdict(lambda: [0, 0])
    ch_counts = defaultdict(lambda: [0, 0])
    for pileup_col in bam.pileup(
        chromosome,
        0,
        len(sequence),
        truncate=True,
        ignore_overlaps=True,
        min_base_quality=args.min_base_quality,
        stepper="samtools",
        redo_baq=False,
        max_depth=args.max_depth,
        ignore_orphans=True,
        compute_baq=True,
        min_mapping_quality=args.min_mapping_quality,
    ):
        is_cpg, ch_target = SHARED.classify_reference_position(
            sequence, pileup_col.pos
        )
        if is_cpg and args.context_mode in {"cg", "both"}:
            count_cg_column(pileup_col, args, cg_counts)
        if ch_target is not None and args.context_mode in {"ch", "both"}:
            context, strand = ch_target
            count_ch_column(pileup_col, context, strand, args, ch_counts)
    return cg_counts, ch_counts


def write_outputs(args: argparse.Namespace, chromosomes: list[str]) -> tuple[int, int]:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cg_path = output_dir / f"{args.output_name}.CG.cov"
    ch_path = output_dir / f"{args.output_name}.CH.cov"
    cg_temp = cg_path.with_name(f".{cg_path.name}.tmp.{os.getpid()}")
    ch_temp = ch_path.with_name(f".{ch_path.name}.tmp.{os.getpid()}")
    cg_lines = 0
    ch_lines = 0
    try:
        cg_handle = cg_temp.open("w", encoding="utf-8") if args.context_mode in {"cg", "both"} else None
        ch_handle = ch_temp.open("w", encoding="utf-8") if args.context_mode in {"ch", "both"} else None
        try:
            with pysam.AlignmentFile(args.bam, "rb") as bam, pysam.FastaFile(
                args.reference
            ) as fasta:
                for chromosome in chromosomes:
                    cg_counts, ch_counts = call_chromosome(
                        bam, fasta, chromosome, args
                    )
                    if cg_handle is not None:
                        for (pos, strand), (methylated, unmethylated) in sorted(cg_counts.items()):
                            output_strand = strand if args.cg_strand_mode == "separate" else None
                            cg_handle.write(
                                SHARED.format_cg_line(
                                    chromosome, pos, methylated, unmethylated, output_strand
                                )
                            )
                            cg_lines += 1
                    if ch_handle is not None:
                        for (pos, context, strand), (methylated, unmethylated) in sorted(ch_counts.items()):
                            ch_handle.write(
                                SHARED.format_ch_line(
                                    chromosome,
                                    pos,
                                    methylated,
                                    unmethylated,
                                    context,
                                    strand,
                                )
                            )
                            ch_lines += 1
        finally:
            if cg_handle is not None:
                cg_handle.close()
            if ch_handle is not None:
                ch_handle.close()
        if args.context_mode in {"cg", "both"}:
            os.replace(cg_temp, cg_path)
        if args.context_mode in {"ch", "both"}:
            os.replace(ch_temp, ch_path)
    finally:
        cg_temp.unlink(missing_ok=True)
        ch_temp.unlink(missing_ok=True)
    return cg_lines, ch_lines


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    try:
        chromosomes = validate_args(args)
    except ValueError as error:
        print(f"[spike-caller] error: {error}", file=sys.stderr)
        return 1

    print(f"[spike-caller] assay={args.assay}")
    print(f"[spike-caller] bam={args.bam}")
    print(f"[spike-caller] reference={args.reference}")
    print(f"[spike-caller] output-name={args.output_name}")
    print(f"[spike-caller] chromosomes={','.join(chromosomes)}")
    print(f"[spike-caller] context-mode={args.context_mode}")
    print(f"[spike-caller] cg-strand-mode={args.cg_strand_mode}")
    if args.dry_run:
        print("[spike-caller] dry-run=1")
        return 0

    bam_path = Path(args.bam)
    reference_path = Path(args.reference)
    if not bam_path.is_file():
        print(f"[spike-caller] error: BAM not found: {bam_path}", file=sys.stderr)
        return 1
    if not bam_index_exists(bam_path):
        print(f"[spike-caller] error: BAM index not found: {bam_path}", file=sys.stderr)
        return 1
    if not reference_path.is_file() or not Path(str(reference_path) + ".fai").is_file():
        print(
            f"[spike-caller] error: reference FASTA/.fai not found: {reference_path}",
            file=sys.stderr,
        )
        return 1
    try:
        SHARED.validate_chromosomes(args.reference, chromosomes)
        cg_lines, ch_lines = write_outputs(args, chromosomes)
    except (OSError, ValueError) as error:
        print(f"[spike-caller] error: {error}", file=sys.stderr)
        return 1
    print(f"[spike-caller] done cg-sites={cg_lines} ch-sites={ch_lines}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
