#!/usr/bin/env python3
"""Methylation caller for DBiT TAPS/TAPS-v2, EM-seq, and Cabernet data.
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
from typing import BinaryIO, Dict, List, Optional, Set, Tuple

import pysam

# properly-paired primary alignment flags grouped by read alignment strand
TOP_FLAGS: Set[int] = {99, 163}
BOT_FLAGS: Set[int] = {83, 147}
PAIRED_FLAGS: Set[int] = TOP_FLAGS | BOT_FLAGS

# forward-strand CH contexts and their reverse-complement mapping
CH_FORWARD_CONTEXTS: Set[str] = {"CA", "CC", "CT"}
CH_REVERSE_MAP: Dict[str, str] = {"T": "CA", "G": "CC", "A": "CT"}
CH_REVERSE_NEIGHBOR: Dict[str, str] = {"CA": "T", "CC": "G", "CT": "A"}

# Set in the parent immediately before forking interval workers. The immutable
# chromosome string is then shared copy-on-write instead of being serialized
# once per interval.
WORKER_CHROMOSOME_SEQUENCE: Optional[str] = None


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Methylation caller for DBiT TAPS/TAPS-v2, EM-seq, and Cabernet data."
        )
    )
    parser.add_argument(
        "--assay",
        choices=("taps", "taps-v2", "emseq", "cabernet"),
        required=True,
        help="Methylation assay used to interpret converted bases.",
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
        help=(
            "Output directory for combined host.CG.cov / host.CA.cov / "
            "host.CC.cov / host.CT.cov files."
        ),
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
        help=(
            "Maximum genomic batch size in bp. It is reduced per chromosome "
            "when needed to keep interval workers occupied. Default: "
            "10,000,000."
        ),
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


def effective_batch_size(
    chromosome_length: int,
    configured_batch_size: int,
    jobs: int,
) -> int:
    """Return a batch size that creates at least one interval per worker."""
    return min(
        configured_batch_size,
        max(1, chromosome_length // jobs),
    )


BarcodeMap = Dict[str, str]


def load_barcode_whitelist(
    path: str, c_to_t: bool = False
) -> Tuple[BarcodeMap, List[Tuple[str, str]], int]:
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
    lookup_sequences = [
        sequence.replace("C", "T") if c_to_t else sequence
        for sequence in sequences
    ]
    if len(lookup_sequences) != len(set(lookup_sequences)):
        raise ValueError(
            "barcode whitelist contains duplicate sequences after C-to-T conversion"
        )

    return (
        {
            lookup_sequence: barcode_id
            for (barcode_id, _), lookup_sequence in zip(entries, lookup_sequences)
        },
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


# CG counters: per-spot → forward-C position → [methylated, unmethylated]
CgCounts = Dict[str, Dict[int, List[int]]]
# CH counters: per-spot → position → [methylated, unmethylated, context, strand]
ChCounts = Dict[str, Dict[int, Tuple[int, int, str, str]]]


def classify_cg_observation(
    dinucleotide: str,
    flag: int,
    assay: str,
) -> Optional[Tuple[bool, str]]:
    """Return methylation status and cytosine strand for a CpG observation."""
    if dinucleotide == "TG":
        strand = "+"
    elif dinucleotide == "CA":
        strand = "-"
    elif dinucleotide == "CG":
        if flag in TOP_FLAGS:
            strand = "+"
        elif flag in BOT_FLAGS:
            strand = "-"
        else:
            return None
    else:
        return None

    converted = dinucleotide in {"TG", "CA"}
    methylated = not converted if assay in {"emseq", "cabernet"} else converted
    return methylated, strand


def count_cg_column(
    pileup_col,
    cb_tag: str,
    cg_counts: CgCounts,
    assay: str,
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
        # pos is the forward-strand C used for evidence from both strands.
        observation = classify_cg_observation(dinuc, record.flag, assay)
        if observation is None:
            continue
        methylated, _ = observation
        value_index = 0 if methylated else 1
        try:
            cb = record.get_tag(cb_tag)
        except KeyError:
            continue

        spot = cg_counts.setdefault(cb, {}).setdefault(pos, [0, 0])
        spot[value_index] = int(spot[value_index]) + 1


def classify_ch_observation(
    base: str,
    strand: str,
    assay: str,
) -> Optional[Tuple[int, int]]:
    """Return methylated and unmethylated increments for a CH observation."""
    if strand == "+":
        if base not in {"C", "T"}:
            return None
        retained = base == "C"
    else:
        if base not in {"G", "A"}:
            return None
        retained = base == "G"

    methylated = retained if assay in {"emseq", "cabernet"} else not retained
    return (1, 0) if methylated else (0, 1)


def ch_neighbor_matches(context: str, strand: str, base: str) -> bool:
    """Return whether an observed neighbor is compatible with a CH context."""
    observed = base.upper()
    if context == "CC":
        return observed in (("C", "T") if strand == "+" else ("G", "A"))
    expected = (
        context[1] if strand == "+" else CH_REVERSE_NEIGHBOR[context]
    )
    return observed == expected


def count_ch_column(
    pileup_col,
    ctx: str,
    strand: str,
    cb_tag: str,
    ch_counts: ChCounts,
    assay: str,
    r1_left: int, r1_right: int,
    r2_left: int, r2_right: int,
) -> None:
    pos = pileup_col.pos
    strand_flags = TOP_FLAGS if strand == "+" else BOT_FLAGS
    for pileup_read in pileup_col.pileups:
        record = pileup_read.alignment
        if record.flag not in strand_flags:
            continue
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        qpos = pileup_read.query_position
        if qpos is None:
            continue
        seq = record.query_sequence
        if not seq or qpos >= len(seq):
            continue
        if is_trimmed(record, qpos, len(seq),
                      r1_left, r1_right, r2_left, r2_right):
            continue

        base = seq[qpos].upper()
        if strand == "+":
            if qpos + 1 >= len(seq):
                continue
            if not ch_neighbor_matches(ctx, strand, seq[qpos + 1]):
                continue
        else:
            if qpos == 0:
                continue
            if not ch_neighbor_matches(ctx, strand, seq[qpos - 1]):
                continue
        observation = classify_ch_observation(base, strand, assay)
        if observation is None:
            continue
        meth, unmeth = observation

        try:
            cb = record.get_tag(cb_tag)
        except KeyError:
            continue
        spot = ch_counts.setdefault(cb, {}).setdefault(
            pos, [0, 0, ctx, strand]
        )
        spot[0] += meth
        spot[1] += unmeth


def format_cg_line(
    chrom: str,
    pos: int,
    meth: int,
    unmeth: int,
) -> str:
    total = meth + unmeth
    pct = round((meth / total) * 100, 2) if total > 0 else 0.0
    fields = [chrom, str(pos), str(pos), f"{pct:.2f}", str(meth), str(unmeth)]
    return "\t".join(fields) + "\n"


def format_ch_line(
    chrom: str,
    pos: int,
    meth: int,
    unmeth: int,
    spot_id: str,
    strand: str,
) -> str:
    total = meth + unmeth
    pct = round((meth / total) * 100, 2) if total > 0 else 0.0
    return (
        f"{chrom}\t{pos}\t{pos}\t{pct:.2f}\t{meth}\t{unmeth}\t"
        f"{spot_id}\t{strand}\n"
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
                meth_count = int(meth)
                unmeth_count = int(unmeth)
                cg_line = format_cg_line(
                    chrom, pos, meth_count, unmeth_count
                ).rstrip("\n")
                handle.write(f"{cg_line}\t{spot_id}\n")
                processed += 1
                covered += int(meth_count + unmeth_count > 0)
    return processed, covered


def write_ch_parts(
    ch_counts: ChCounts,
    chrom: str,
    part_paths: Dict[str, Path],
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
    handles = {
        context: path.open("w", encoding="utf-8")
        for context, path in part_paths.items()
    }
    try:
        for spot_id, pos_map in sorted(spot_maps):
            for pos, values in sorted(
                pos_map.items(), key=lambda item: (item[1][2], item[0])
            ):
                meth, unmeth, ctx, strand = values
                handles[ctx].write(
                    format_ch_line(
                        chrom, pos, meth, unmeth, spot_id, strand
                    )
                )
                processed += 1
                covered += int(meth + unmeth > 0)
    finally:
        for handle in handles.values():
            handle.close()
    return processed, covered


BatchResult = Tuple[int, Dict[str, str], int, int, int, int]


def process_interval(
    bam_path: str,
    chrom: str,
    batch_index: int,
    start: int,
    stop: int,
    part_dir: str,
    cb_tag: str,
    assay: str,
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
                    assay,
                    r1_left, r1_right, r2_left, r2_right,
                )
            if ch_target is not None and context_mode in ("ch", "both"):
                ctx, strand = ch_target
                count_ch_column(
                    pileup_col, ctx, strand, cb_tag, ch_counts,
                    assay,
                    r1_left, r1_right, r2_left, r2_right,
                )

    part_root = Path(part_dir)
    part_paths: Dict[str, Path] = {}
    cpg_processed = 0
    cpg_covered = 0
    ch_processed = 0
    ch_covered = 0
    if context_mode in ("cg", "both"):
        part_paths["CG"] = part_root / f"batch_{batch_index:06d}.CG.part"
        cpg_processed, cpg_covered = write_cg_part(
            cg_counts, chrom, part_paths["CG"], barcode_map, barcode_length
        )
    if context_mode in ("ch", "both"):
        ch_part_paths = {
            context: part_root / f"batch_{batch_index:06d}.{context}.part"
            for context in sorted(CH_FORWARD_CONTEXTS)
        }
        part_paths.update(ch_part_paths)
        ch_processed, ch_covered = write_ch_parts(
            ch_counts, chrom, ch_part_paths, barcode_map, barcode_length
        )
    return (
        batch_index,
        {context: str(path) for context, path in part_paths.items()},
        cpg_processed,
        cpg_covered,
        ch_processed,
        ch_covered,
    )


def append_part_to_output(
    part_path: Path,
    context: str,
    out_dir: Path,
    output_paths: Dict[Path, Path],
    output_handles: Dict[Path, BinaryIO],
    temp_suffix: str,
) -> None:
    if part_path.stat().st_size == 0:
        return
    final_path = out_dir / "host" / f"host.{context}.cov"
    if final_path not in output_handles:
        final_path.parent.mkdir(parents=True, exist_ok=True)
        temp_path = final_path.with_name(final_path.name + temp_suffix)
        output_paths[final_path] = temp_path
        output_handles[final_path] = temp_path.open("wb")
    with part_path.open("rb") as part_handle:
        shutil.copyfileobj(
            part_handle,
            output_handles[final_path],
            length=16 * 1024 * 1024,
        )


def process_chromosome(
    args: argparse.Namespace,
    chrom: str,
    chromosome_index: int,
    part_root: Path,
    barcode_map: BarcodeMap,
    barcode_length: int,
    output_paths: Dict[Path, Path],
    output_handles: Dict[Path, BinaryIO],
    temp_suffix: str,
) -> Tuple[int, int, int, int]:
    global WORKER_CHROMOSOME_SEQUENCE

    with pysam.FastaFile(args.reference) as fasta:
        sequence = fasta.fetch(chrom)
    if not sequence.isupper():
        sequence = sequence.upper()
    WORKER_CHROMOSOME_SEQUENCE = sequence

    chromosome_length = len(sequence)
    chrom_batch_size = effective_batch_size(
        chromosome_length,
        args.batch_size,
        args.jobs,
    )
    intervals = [
        (batch_index, start, min(start + chrom_batch_size, chromosome_length))
        for batch_index, start in enumerate(
            range(0, chromosome_length, chrom_batch_size)
        )
    ]
    chromosome_part_dir = part_root / f"{chromosome_index:03d}"
    chromosome_part_dir.mkdir(parents=True)
    results: List[BatchResult] = []

    print(
        f"[methy-caller] chrom={chrom} length={chromosome_length} "
        f"batch-size={chrom_batch_size} batches={len(intervals)} "
        f"workers={min(args.jobs, len(intervals))}"
    )

    worker_args = (
        args.bam,
        chrom,
        chromosome_part_dir,
        args.cb_tag,
        args.assay,
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
        _, part_paths, cp, cc, hp, hc = result
        for context in ("CG", *sorted(CH_FORWARD_CONTEXTS)):
            part_path = part_paths.get(context)
            if part_path is None:
                continue
            append_part_to_output(
                Path(part_path), context, Path(args.out_dir),
                output_paths, output_handles, temp_suffix,
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
        print(f"[methy-caller] assay={args.assay}")
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
        barcode_c_to_t = args.assay in ("emseq", "cabernet")
        barcode_map, barcode_entries, barcode_length = load_barcode_whitelist(
            args.barcode_whitelist, c_to_t=barcode_c_to_t
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

    print(f"[methy-caller] assay={args.assay}")
    print(f"[methy-caller] bam={args.bam}")
    print(f"[methy-caller] reference={args.reference}")
    print(f"[methy-caller] barcode-whitelist={args.barcode_whitelist}")
    print(f"[methy-caller] barcode-c-to-t={str(barcode_c_to_t).lower()}")
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
    output_handles: Dict[Path, BinaryIO] = {}
    try:
        with tempfile.TemporaryDirectory(
            prefix=".methylation-call-parts.", dir=out_dir
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
                    output_handles,
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
        for output_handle in output_handles.values():
            output_handle.close()
        output_handles.clear()
        for final_path, temp_path in output_paths.items():
            os.replace(temp_path, final_path)
        os.replace(manifest_temp, manifest_path)
    except Exception as exc:
        for output_handle in output_handles.values():
            output_handle.close()
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
