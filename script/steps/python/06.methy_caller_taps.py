#!/usr/bin/env python3
"""Methylation caller for DBiT TAPS/TAPS-v2 data.

Reads a coordinate-sorted pooled BAM, performs pileup over reference CpG and
CH positions, groups counts by spot (CB tag), and writes per-spot .CG.cov /
.CH.cov files. Chromosomes are processed serially; genomic intervals within
one chromosome are processed in parallel and merged by the parent process.
"""

from __future__ import annotations

import argparse
import multiprocessing
import os
import re
import shutil
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, TextIO, Tuple

import pysam

# properly-paired primary alignment flags for paired-end data
TOP_FLAGS: Set[int] = {99, 147}
BOT_FLAGS: Set[int] = {83, 163}
PAIRED_FLAGS: Set[int] = TOP_FLAGS | BOT_FLAGS

# forward-strand CH contexts and their reverse-complement mapping
CH_FORWARD_CONTEXTS: Set[str] = {"CA", "CC", "CT"}
CH_REVERSE_MAP: Dict[str, str] = {"T": "CA", "G": "CC", "A": "CT"}

# Set in the parent immediately before forking interval workers. The immutable
# chromosome string is then shared copy-on-write instead of being serialized
# once per interval.
WORKER_CHROMOSOME_SEQUENCE: Optional[str] = None


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Methylation caller for DBiT TAPS/TAPS-v2 data."
    )
    parser.add_argument(
        "-b", "--bam",
        required=True,
        help="Path to a coordinate-sorted and indexed pooled BAM.",
    )
    parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Reference FASTA file.",
    )
    parser.add_argument(
        "-o", "--out-dir",
        required=True,
        help="Output directory for per-spot .CG.cov / .CH.cov files.",
    )
    parser.add_argument(
        "--barcode-whitelist",
        required=True,
        help=(
            "Ordered barcode whitelist. Use '<index> <sequence>' or one "
            "sequence per line."
        ),
    )
    parser.add_argument(
        "-c", "--chromosomes",
        required=True,
        help="Comma-separated chromosome names to call.",
    )
    parser.add_argument(
        "--cb-tag",
        default="CB",
        help="BAM tag for spot barcode. Default: CB.",
    )
    parser.add_argument(
        "--context-mode",
        choices=("cg", "ch", "both"),
        default="cg",
        help="Methylation context to call. Default: cg.",
    )
    parser.add_argument(
        "--min-base-quality",
        type=int,
        default=30,
        help="Minimum base quality. Default: 30.",
    )
    parser.add_argument(
        "--min-mapping-quality",
        type=int,
        default=10,
        help="Minimum mapping quality. Default: 10.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=50_000,
        help="Maximum pileup depth. Default: 50,000.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10_000_000,
        help="Genomic batch size in bp. Default: 10,000,000.",
    )
    parser.add_argument(
        "--r1-left-trim", type=int, default=0,
        help="Trim N bp from R1 5' end. Default: 0.",
    )
    parser.add_argument(
        "--r1-right-trim", type=int, default=0,
        help="Trim N bp from R1 3' end. Default: 0.",
    )
    parser.add_argument(
        "--r2-left-trim", type=int, default=0,
        help="Trim N bp from R2 5' end. Default: 0.",
    )
    parser.add_argument(
        "--r2-right-trim", type=int, default=0,
        help="Trim N bp from R2 3' end. Default: 0.",
    )
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=8,
        help="Worker processes for intervals within one chromosome. Default: 8.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print parameters and exit.",
    )
    return parser.parse_args(argv)


def parse_chromosomes(csv: str) -> List[str]:
    return [c.strip() for c in csv.split(",") if c.strip()]


def validate_chromosomes(reference: str, chromosomes: List[str]) -> None:
    with pysam.FastaFile(reference) as fa:
        available = set(fa.references)
    missing = [c for c in chromosomes if c not in available]
    if missing:
        raise ValueError(
            f"chromosomes not found in reference '{reference}': {','.join(missing)}"
        )


BarcodeMap = Dict[str, str]


def load_barcode_whitelist(path: str) -> Tuple[BarcodeMap, List[Tuple[str, str]], int]:
    """Return sequence-to-index mapping, ordered entries, and barcode length."""
    whitelist_path = Path(path)
    if not whitelist_path.is_file():
        raise ValueError(f"barcode whitelist not found: {path}")

    raw_entries: List[Tuple[Optional[str], str]] = []
    with whitelist_path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) == 1:
                barcode_id: Optional[str] = None
                sequence = fields[0].upper()
            else:
                barcode_id = fields[0]
                sequence = fields[-1].upper()
            if not sequence or re.fullmatch(r"[ACGTN]+", sequence) is None:
                raise ValueError(
                    f"invalid barcode sequence at {path}:{line_number}: {sequence}"
                )
            raw_entries.append((barcode_id, sequence))

    if not raw_entries:
        raise ValueError(f"barcode whitelist is empty: {path}")

    generated_width = max(2, len(str(len(raw_entries) - 1)))
    entries: List[Tuple[str, str]] = []
    for position, (barcode_id, sequence) in enumerate(raw_entries):
        resolved_id = (
            barcode_id
            if barcode_id is not None
            else f"{position:0{generated_width}d}"
        )
        if re.fullmatch(r"[A-Za-z0-9.-]+", resolved_id) is None:
            raise ValueError(
                f"barcode index is not filename-safe: {resolved_id}"
            )
        entries.append((resolved_id, sequence))

    barcode_lengths = {len(sequence) for _, sequence in entries}
    if len(barcode_lengths) != 1:
        raise ValueError("all whitelist barcode sequences must have equal length")
    ids = [barcode_id for barcode_id, _ in entries]
    sequences = [sequence for _, sequence in entries]
    if len(ids) != len(set(ids)):
        raise ValueError("barcode whitelist contains duplicate indices")
    if len(sequences) != len(set(sequences)):
        raise ValueError("barcode whitelist contains duplicate sequences")

    return (
        {sequence: barcode_id for barcode_id, sequence in entries},
        entries,
        barcode_lengths.pop(),
    )


def spot_id_from_cb(
    cb: str, barcode_map: BarcodeMap, barcode_length: int
) -> str:
    if "+" in cb:
        barcodes = cb.split("+")
        if len(barcodes) != 2:
            raise ValueError(
                f"CB tag '{cb}' must contain exactly two barcodes"
            )
        barcode1, barcode2 = barcodes
        if len(barcode1) != barcode_length or len(barcode2) != barcode_length:
            raise ValueError(
                f"CB tag '{cb}' contains barcodes with unexpected lengths; "
                f"expected {barcode_length}+{barcode_length} bases"
            )
    else:
        expected_length = barcode_length * 2
        if len(cb) != expected_length:
            raise ValueError(
                f"CB tag '{cb}' has length {len(cb)}; "
                f"expected {expected_length} without a separator"
            )
        barcode1 = cb[:barcode_length]
        barcode2 = cb[barcode_length:]

    barcode1 = barcode1.upper()
    barcode2 = barcode2.upper()
    try:
        index1 = barcode_map[barcode1]
        index2 = barcode_map[barcode2]
    except KeyError as exc:
        raise ValueError(
            f"CB tag '{cb}' contains a barcode absent from the whitelist"
        ) from exc
    return f"{index1}_{index2}"


def write_spot_manifest(
    path: Path, entries: List[Tuple[str, str]]
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "spot_id\tbarcode1_index\tbarcode2_index\t"
            "barcode1\tbarcode2\traw_cb\n"
        )
        for index1, barcode1 in entries:
            for index2, barcode2 in entries:
                handle.write(
                    f"{index1}_{index2}\t{index1}\t{index2}\t"
                    f"{barcode1}\t{barcode2}\t{barcode1}+{barcode2}\n"
                )


def classify_reference_position(
    sequence: str, pos: int
) -> Tuple[bool, Optional[Tuple[str, str]]]:
    """Return whether pos is CpG-C and an optional CH (context, strand)."""
    base = sequence[pos].upper()
    next_base = (
        sequence[pos + 1].upper() if pos + 1 < len(sequence) else ""
    )
    if base == "C":
        dinucleotide = base + next_base
        if dinucleotide == "CG":
            return True, None
        if dinucleotide in CH_FORWARD_CONTEXTS:
            return False, (dinucleotide, "+")
        return False, None
    if base == "G" and pos > 0:
        context = CH_REVERSE_MAP.get(sequence[pos - 1].upper())
        if context is not None:
            return False, (context, "-")
    return False, None


def is_trimmed(
    record: pysam.AlignedSegment,
    query_pos: int,
    read_len: int,
    r1_left: int, r1_right: int,
    r2_left: int, r2_right: int,
) -> bool:
    """Return True if query_pos falls within trimmed region."""
    if record.is_read1:
        left_trim, right_trim = r1_left, r1_right
    elif record.is_read2:
        left_trim, right_trim = r2_left, r2_right
    else:
        return False
    if left_trim == 0 and right_trim == 0:
        return False
    if record.is_reverse:
        left_cycle = read_len - query_pos
        right_cycle = query_pos + 1
    else:
        left_cycle = query_pos + 1
        right_cycle = read_len - query_pos
    return left_cycle <= left_trim or right_cycle <= right_trim


# CG counters: per-spot → position → [methylated, unmethylated]
CgCounts = Dict[str, Dict[int, Tuple[int, int]]]
# CH counters: per-spot → position → [methylated, unmethylated, context, strand]
ChCounts = Dict[str, Dict[int, Tuple[int, int, str, str]]]


def count_cg_column(
    pileup_col,
    cb_tag: str,
    cg_counts: CgCounts,
    r1_left: int, r1_right: int,
    r2_left: int, r2_right: int,
) -> None:
    pos = pileup_col.pos
    for pileup_read in pileup_col.pileups:
        record = pileup_read.alignment
        if record.flag not in PAIRED_FLAGS:
            continue
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        qpos = pileup_read.query_position
        if qpos is None:
            continue
        seq = record.query_sequence
        if not seq:
            continue
        read_len = len(seq)
        if qpos + 1 >= read_len:
            continue
        if is_trimmed(record, qpos, read_len,
                       r1_left, r1_right, r2_left, r2_right):
            continue
        next_qpos = qpos + 1
        if is_trimmed(record, next_qpos, read_len,
                       r1_left, r1_right, r2_left, r2_right):
            continue
        dinuc = (seq[qpos] + seq[next_qpos]).upper()
        # pos is the forward-strand C; pos + 1 is the reference G that
        # represents the complementary-strand C.
        if dinuc == "TG":
            target_pos = pos
            value_index = 0
        elif dinuc == "CA":
            target_pos = pos + 1
            value_index = 0
        elif dinuc == "CG":
            target_pos = pos if record.flag in TOP_FLAGS else pos + 1
            value_index = 1
        else:
            continue
        try:
            cb = record.get_tag(cb_tag)
        except KeyError:
            continue

        spot = cg_counts.setdefault(cb, {}).setdefault(target_pos, [0, 0])
        spot[value_index] += 1


def get_forward_index(record: pysam.AlignedSegment, query_pos: int,
                      read_len: int) -> int:
    if record.is_reverse:
        return read_len - query_pos - 1
    return query_pos


def count_ch_column(
    pileup_col,
    ctx: str,
    strand: str,
    cb_tag: str,
    ch_counts: ChCounts,
    r1_left: int, r1_right: int,
    r2_left: int, r2_right: int,
) -> None:
    pos = pileup_col.pos
    for pileup_read in pileup_col.pileups:
        record = pileup_read.alignment
        if record.flag not in PAIRED_FLAGS:
            continue
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        qpos = pileup_read.query_position
        if qpos is None:
            continue
        fwd_seq = record.get_forward_sequence()
        if not fwd_seq:
            continue
        read_len = len(fwd_seq)
        fwd_idx = get_forward_index(record, qpos, read_len)
        if fwd_idx < 0 or fwd_idx >= read_len:
            continue
        if is_trimmed(record, qpos, read_len,
                       r1_left, r1_right, r2_left, r2_right):
            continue

        base = fwd_seq[fwd_idx].upper()
        if strand == "+":
            if fwd_idx + 1 >= read_len:
                continue
            if fwd_seq[fwd_idx + 1].upper() != ctx[1]:
                continue
            if base == "T":
                meth, unmeth = 1, 0
            elif base == "C":
                meth, unmeth = 0, 1
            else:
                continue
        else:  # strand == "-"
            if base == "A":
                meth, unmeth = 1, 0
            elif base == "G":
                meth, unmeth = 0, 1
            else:
                continue

        try:
            cb = record.get_tag(cb_tag)
        except KeyError:
            continue
        spot = ch_counts.setdefault(cb, {}).setdefault(
            pos, [0, 0, ctx, strand]
        )
        spot[0] += meth
        spot[1] += unmeth


def format_cg_line(chrom: str, pos: int, meth: int, unmeth: int) -> str:
    total = meth + unmeth
    pct = round((meth / total) * 100, 2) if total > 0 else 0.0
    return f"{chrom}\t{pos + 1}\t{pos + 1}\t{pct:.2f}\t{meth}\t{unmeth}\n"


def format_ch_line(chrom: str, pos: int, meth: int, unmeth: int,
                   ctx: str, strand: str) -> str:
    total = meth + unmeth
    pct = round((meth / total) * 100, 2) if total > 0 else 0.0
    return (
        f"{chrom}\t{pos + 1}\t{pos + 1}\t{pct:.2f}\t{meth}\t{unmeth}\t"
        f"{ctx}\t{strand}\n"
    )


def write_cg_part(
    cg_counts: CgCounts,
    chrom: str,
    part_path: Path,
    barcode_map: BarcodeMap,
    barcode_length: int,
) -> Tuple[int, int]:
    processed = 0
    covered = 0
    spot_maps = [
        (spot_id_from_cb(cb, barcode_map, barcode_length), pos_map)
        for cb, pos_map in cg_counts.items()
        if pos_map
    ]
    with part_path.open("w", encoding="utf-8") as handle:
        for spot_id, pos_map in sorted(spot_maps):
            for pos in sorted(pos_map):
                meth, unmeth = pos_map[pos]
                handle.write(
                    f"{spot_id}\t{format_cg_line(chrom, pos, meth, unmeth)}"
                )
                processed += 1
                covered += int(meth + unmeth > 0)
    return processed, covered


def write_ch_part(
    ch_counts: ChCounts,
    chrom: str,
    part_path: Path,
    barcode_map: BarcodeMap,
    barcode_length: int,
) -> Tuple[int, int]:
    processed = 0
    covered = 0
    spot_maps = [
        (spot_id_from_cb(cb, barcode_map, barcode_length), pos_map)
        for cb, pos_map in ch_counts.items()
        if pos_map
    ]
    with part_path.open("w", encoding="utf-8") as handle:
        for spot_id, pos_map in sorted(spot_maps):
            for pos in sorted(pos_map):
                meth, unmeth, ctx, strand = pos_map[pos]
                handle.write(
                    f"{spot_id}\t"
                    f"{format_ch_line(chrom, pos, meth, unmeth, ctx, strand)}"
                )
                processed += 1
                covered += int(meth + unmeth > 0)
    return processed, covered


BatchResult = Tuple[int, str, str, int, int, int, int]


def process_interval(
    bam_path: str,
    chrom: str,
    batch_index: int,
    start: int,
    stop: int,
    part_dir: str,
    cb_tag: str,
    context_mode: str,
    min_base_quality: int,
    min_mapping_quality: int,
    max_depth: int,
    r1_left: int, r1_right: int,
    r2_left: int, r2_right: int,
    barcode_map: BarcodeMap,
    barcode_length: int,
) -> BatchResult:
    sequence = WORKER_CHROMOSOME_SEQUENCE
    if sequence is None:
        raise RuntimeError("worker chromosome sequence was not initialized")

    cg_counts: CgCounts = defaultdict(dict)
    ch_counts: ChCounts = defaultdict(dict)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        pileup_iter = bam.pileup(
            chrom, start, stop,
            truncate=True,
            ignore_overlaps=True,
            min_base_quality=min_base_quality,
            stepper="samtools",
            redo_baq=False,
            max_depth=max_depth,
            ignore_orphans=True,
            compute_baq=True,
            min_mapping_quality=min_mapping_quality,
        )
        for pileup_col in pileup_iter:
            pos = pileup_col.pos
            is_cpg, ch_target = classify_reference_position(sequence, pos)
            if is_cpg and context_mode in ("cg", "both"):
                count_cg_column(
                    pileup_col, cb_tag, cg_counts,
                    r1_left, r1_right, r2_left, r2_right,
                )
            if ch_target is not None and context_mode in ("ch", "both"):
                ctx, strand = ch_target
                count_ch_column(
                    pileup_col, ctx, strand, cb_tag, ch_counts,
                    r1_left, r1_right, r2_left, r2_right,
                )

    part_root = Path(part_dir)
    cg_part = part_root / f"batch_{batch_index:06d}.CG.part"
    ch_part = part_root / f"batch_{batch_index:06d}.CH.part"
    cpg_processed, cpg_covered = write_cg_part(
        cg_counts, chrom, cg_part, barcode_map, barcode_length
    )
    ch_processed, ch_covered = write_ch_part(
        ch_counts, chrom, ch_part, barcode_map, barcode_length
    )
    return (
        batch_index,
        str(cg_part),
        str(ch_part),
        cpg_processed,
        cpg_covered,
        ch_processed,
        ch_covered,
    )


def append_part_to_outputs(
    part_path: Path,
    context: str,
    out_dir: Path,
    output_paths: Dict[Path, Path],
    temp_suffix: str,
) -> None:
    active_spot: Optional[str] = None
    output_handle: Optional[TextIO] = None
    try:
        with part_path.open(encoding="utf-8") as part_handle:
            for raw_line in part_handle:
                spot_id, payload = raw_line.split("\t", 1)
                if spot_id != active_spot:
                    if output_handle is not None:
                        output_handle.close()
                    prefix = spot_id.split("_", 1)[0]
                    final_path = (
                        out_dir / "host" / prefix / f"{spot_id}.{context}.cov"
                    )
                    final_path.parent.mkdir(parents=True, exist_ok=True)
                    if final_path in output_paths:
                        temp_path = output_paths[final_path]
                        mode = "a"
                    else:
                        temp_path = final_path.with_name(
                            final_path.name + temp_suffix
                        )
                        output_paths[final_path] = temp_path
                        mode = "w"
                    output_handle = temp_path.open(mode, encoding="utf-8")
                    active_spot = spot_id
                output_handle.write(payload)
    finally:
        if output_handle is not None:
            output_handle.close()


def process_chromosome(
    args: argparse.Namespace,
    chrom: str,
    chromosome_index: int,
    part_root: Path,
    barcode_map: BarcodeMap,
    barcode_length: int,
    output_paths: Dict[Path, Path],
    temp_suffix: str,
) -> Tuple[int, int, int, int]:
    global WORKER_CHROMOSOME_SEQUENCE

    with pysam.FastaFile(args.reference) as fasta:
        sequence = fasta.fetch(chrom)
    if not sequence.isupper():
        sequence = sequence.upper()
    WORKER_CHROMOSOME_SEQUENCE = sequence

    chromosome_length = len(sequence)
    intervals = [
        (batch_index, start, min(start + args.batch_size, chromosome_length))
        for batch_index, start in enumerate(
            range(0, chromosome_length, args.batch_size)
        )
    ]
    chromosome_part_dir = part_root / f"{chromosome_index:03d}"
    chromosome_part_dir.mkdir(parents=True)
    results: List[BatchResult] = []

    print(
        f"[methy-caller] chrom={chrom} length={chromosome_length} "
        f"batches={len(intervals)} workers={min(args.jobs, len(intervals))}"
    )

    worker_args = (
        args.bam,
        chrom,
        chromosome_part_dir,
        args.cb_tag,
        args.context_mode,
        args.min_base_quality,
        args.min_mapping_quality,
        args.max_depth,
        args.r1_left_trim,
        args.r1_right_trim,
        args.r2_left_trim,
        args.r2_right_trim,
        barcode_map,
        barcode_length,
    )
    if args.jobs == 1 or len(intervals) == 1:
        for batch_index, start, stop in intervals:
            results.append(
                process_interval(
                    worker_args[0], worker_args[1],
                    batch_index, start, stop,
                    str(worker_args[2]), *worker_args[3:],
                )
            )
    else:
        fork_context = multiprocessing.get_context("fork")
        with ProcessPoolExecutor(
            max_workers=min(args.jobs, len(intervals)),
            mp_context=fork_context,
        ) as executor:
            future_to_batch = {
                executor.submit(
                    process_interval,
                    worker_args[0], worker_args[1],
                    batch_index, start, stop,
                    str(worker_args[2]), *worker_args[3:],
                ): batch_index
                for batch_index, start, stop in intervals
            }
            for future in as_completed(future_to_batch):
                batch_index = future_to_batch[future]
                try:
                    results.append(future.result())
                except Exception as exc:
                    raise RuntimeError(
                        f"chrom={chrom} batch={batch_index + 1} failed: {exc}"
                    ) from exc

    totals = [0, 0, 0, 0]
    for result in sorted(results):
        _, cg_part, ch_part, cp, cc, hp, hc = result
        if args.context_mode in ("cg", "both"):
            append_part_to_outputs(
                Path(cg_part), "CG", Path(args.out_dir),
                output_paths, temp_suffix,
            )
        if args.context_mode in ("ch", "both"):
            append_part_to_outputs(
                Path(ch_part), "CH", Path(args.out_dir),
                output_paths, temp_suffix,
            )
        totals[0] += cp
        totals[1] += cc
        totals[2] += hp
        totals[3] += hc

    shutil.rmtree(chromosome_part_dir)
    WORKER_CHROMOSOME_SEQUENCE = None
    return totals[0], totals[1], totals[2], totals[3]


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    if args.dry_run:
        print(f"[methy-caller] bam={args.bam}")
        print(f"[methy-caller] reference={args.reference}")
        print(f"[methy-caller] out-dir={args.out_dir}")
        print(f"[methy-caller] barcode-whitelist={args.barcode_whitelist}")
        print(f"[methy-caller] chromosomes={args.chromosomes}")
        print(f"[methy-caller] context-mode={args.context_mode}")
        print(f"[methy-caller] dry-run=1")
        return 0

    bam_path = Path(args.bam)
    if not bam_path.exists():
        print(f"[methy-caller] error: BAM not found: {args.bam}", file=sys.stderr)
        return 1
    bam_index = Path(str(args.bam) + ".bai")
    if not bam_index.exists():
        bam_index = Path(str(args.bam)[:-4] + ".bai")
    if not bam_index.exists() and not Path(str(args.bam) + ".csi").exists():
        print(f"[methy-caller] error: BAM index not found for {args.bam}",
              file=sys.stderr)
        return 1

    if args.jobs <= 0:
        print("[methy-caller] error: --jobs must be greater than zero",
              file=sys.stderr)
        return 1
    if args.batch_size <= 0:
        print("[methy-caller] error: --batch-size must be greater than zero",
              file=sys.stderr)
        return 1
    if args.max_depth <= 0:
        print("[methy-caller] error: --max-depth must be greater than zero",
              file=sys.stderr)
        return 1

    try:
        chromosomes = parse_chromosomes(args.chromosomes)
        validate_chromosomes(args.reference, chromosomes)
        barcode_map, barcode_entries, barcode_length = load_barcode_whitelist(
            args.barcode_whitelist
        )
    except (OSError, ValueError) as exc:
        print(f"[methy-caller] error: {exc}", file=sys.stderr)
        return 1

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_suffix = f".tmp.{os.getpid()}"
    manifest_path = out_dir / "spot_manifest.tsv"
    manifest_temp = manifest_path.with_name(manifest_path.name + temp_suffix)
    write_spot_manifest(manifest_temp, barcode_entries)

    print(f"[methy-caller] bam={args.bam}")
    print(f"[methy-caller] reference={args.reference}")
    print(f"[methy-caller] barcode-whitelist={args.barcode_whitelist}")
    print(f"[methy-caller] spots={len(barcode_entries) ** 2}")
    print(f"[methy-caller] chromosomes={len(chromosomes)}")
    print(f"[methy-caller] context-mode={args.context_mode}")
    print(f"[methy-caller] interval-workers={args.jobs}")
    print(f"[methy-caller] out-dir={args.out_dir}")

    total_cpg_proc = 0
    total_cpg_cov = 0
    total_ch_proc = 0
    total_ch_cov = 0

    output_paths: Dict[Path, Path] = {}
    try:
        with tempfile.TemporaryDirectory(
            prefix=".taps-call-parts.", dir=out_dir
        ) as part_directory:
            for chromosome_index, chrom in enumerate(chromosomes):
                cp, cc, hp, hc = process_chromosome(
                    args,
                    chrom,
                    chromosome_index,
                    Path(part_directory),
                    barcode_map,
                    barcode_length,
                    output_paths,
                    temp_suffix,
                )
                total_cpg_proc += cp
                total_cpg_cov += cc
                total_ch_proc += hp
                total_ch_cov += hc
                print(
                    f"[methy-caller] chrom={chrom} merged "
                    f"cpg-processed={cp} cpg-covered={cc} "
                    f"ch-processed={hp} ch-covered={hc}"
                )
        for final_path, temp_path in output_paths.items():
            os.replace(temp_path, final_path)
        os.replace(manifest_temp, manifest_path)
    except Exception as exc:
        for temp_path in output_paths.values():
            try:
                temp_path.unlink()
            except FileNotFoundError:
                pass
        try:
            manifest_temp.unlink()
        except FileNotFoundError:
            pass
        print(f"[methy-caller] error: {exc}", file=sys.stderr)
        return 1

    print(
        f"[methy-caller] done "
        f"cpg-processed={total_cpg_proc} cpg-covered={total_cpg_cov} "
        f"ch-processed={total_ch_proc} ch-covered={total_ch_cov}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
