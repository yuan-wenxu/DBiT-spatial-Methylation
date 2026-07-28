#!/usr/bin/env python3
"""Methylation caller for DBiT TAPS/TAPS-v2 data.

Reads a coordinate-sorted pooled BAM, performs pileup over reference CpG and
CH positions, groups counts by spot (CB tag), and writes per-spot .CG.cov /
.CH.cov files in one BAM pass.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pysam

# properly-paired primary alignment flags for paired-end data
PAIRED_FLAGS: Set[int] = {99, 147, 83, 163}

# forward-strand CH contexts and their reverse-complement mapping
CH_FORWARD_CONTEXTS: Set[str] = {"CA", "CC", "CT"}
CH_REVERSE_MAP: Dict[str, str] = {"T": "CA", "G": "CC", "A": "CT"}


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
        default=250,
        help="Maximum pileup depth. Default: 250.",
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
        help="Worker processes for parallel chromosome processing. Default: 8.",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Verbose logging.",
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


def find_cpg_positions(reference: str, chromosome: str) -> List[int]:
    """Return 0-based positions of C in CpG contexts."""
    with pysam.FastaFile(reference) as fa:
        seq = fa.fetch(chromosome).upper()
    positions: List[int] = []
    for i in range(len(seq) - 1):
        if seq[i] == "C" and seq[i + 1] == "G":
            positions.append(i)
    return positions


def find_ch_positions(
    reference: str, chromosome: str
) -> List[Tuple[int, str, str]]:
    """Return (0-based C position, context, strand) for CH sites."""
    with pysam.FastaFile(reference) as fa:
        seq = fa.fetch(chromosome).upper()
    results: List[Tuple[int, str, str]] = []
    # forward strand: C followed by A/C/T
    for i in range(len(seq) - 1):
        ctx = seq[i : i + 2]
        if ctx in CH_FORWARD_CONTEXTS:
            results.append((i, ctx, "+"))
    # reverse strand: G preceded by T/G/A → complement → forward context
    for i in range(1, len(seq) - 1):
        if seq[i] != "G":
            continue
        prev = seq[i - 1]
        ctx = CH_REVERSE_MAP.get(prev)
        if ctx is not None:
            results.append((i, ctx, "-"))
    results.sort(key=lambda x: x[0])
    return results


def build_batches(
    cpg_positions: List[int],
    ch_positions: List[Tuple[int, str, str]],
    batch_size: int,
) -> List[Tuple[int, int, Set[int], Dict[int, Tuple[str, str]]]]:
    """Yield (start, end, cpg_set, ch_map) batches covering positions."""
    cpg_set = set(cpg_positions)
    ch_map = {
        pos: (context, strand)
        for pos, context, strand in ch_positions
    }
    all_positions: Set[int] = cpg_set | set(ch_map)
    if not all_positions:
        return []
    sorted_positions = sorted(all_positions)
    batches: List[
        Tuple[int, int, Set[int], Dict[int, Tuple[str, str]]]
    ] = []
    batch_start = sorted_positions[0]
    current_cpg: Set[int] = set()
    current_ch: Dict[int, Tuple[str, str]] = {}
    for pos in sorted_positions:
        if pos - batch_start > batch_size:
            batches.append((batch_start, pos - 1, current_cpg, current_ch))
            batch_start = pos
            current_cpg = set()
            current_ch = {}
        if pos in cpg_set:
            current_cpg.add(pos)
        ch_context = ch_map.get(pos)
        if ch_context is not None:
            current_ch[pos] = ch_context
    if current_cpg or current_ch:
        batches.append(
            (batch_start, sorted_positions[-1], current_cpg, current_ch)
        )
    return batches


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
    cpg_set: Set[int],
    cb_tag: str,
    cg_counts: CgCounts,
    r1_left: int, r1_right: int,
    r2_left: int, r2_right: int,
) -> None:
    pos = pileup_col.pos
    if pos not in cpg_set:
        return
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
        try:
            cb = record.get_tag(cb_tag)
        except KeyError:
            continue

        spot = cg_counts.setdefault(cb, {}).setdefault(pos, [0, 0])
        if dinuc in ("TG", "CA"):
            spot[0] += 1  # methylated
        elif dinuc == "CG":
            spot[1] += 1  # unmethylated


def get_forward_index(record: pysam.AlignedSegment, query_pos: int,
                      read_len: int) -> int:
    if record.is_reverse:
        return read_len - query_pos - 1
    return query_pos


def count_ch_column(
    pileup_col,
    ch_map: Dict[int, Tuple[str, str]],
    cb_tag: str,
    ch_counts: ChCounts,
    r1_left: int, r1_right: int,
    r2_left: int, r2_right: int,
) -> None:
    pos = pileup_col.pos
    target = ch_map.get(pos)
    if target is None:
        return
    ctx, strand = target
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


def flush_batch_cg(
    cg_counts: CgCounts, chrom: str, out_dir: Path,
) -> None:
    for cb, pos_map in cg_counts.items():
        if not pos_map:
            continue
        subdir = out_dir / "host" / cb[:2]
        subdir.mkdir(parents=True, exist_ok=True)
        out_path = subdir / f"{cb}.CG.cov"
        with out_path.open("a", encoding="utf-8") as fh:
            for pos in sorted(pos_map):
                meth, unmeth = pos_map[pos]
                fh.write(format_cg_line(chrom, pos, meth, unmeth))


def flush_batch_ch(
    ch_counts: ChCounts, chrom: str, out_dir: Path,
) -> None:
    for cb, pos_map in ch_counts.items():
        if not pos_map:
            continue
        subdir = out_dir / "host" / cb[:2]
        subdir.mkdir(parents=True, exist_ok=True)
        out_path = subdir / f"{cb}.CH.cov"
        with out_path.open("a", encoding="utf-8") as fh:
            for pos in sorted(pos_map):
                meth, unmeth, ctx, strand = pos_map[pos]
                fh.write(format_ch_line(chrom, pos, meth, unmeth, ctx, strand))


def process_chromosome(
    bam_path: str,
    reference: str,
    chrom: str,
    out_dir: Path,
    cb_tag: str,
    context_mode: str,
    min_base_quality: int,
    min_mapping_quality: int,
    max_depth: int,
    batch_size: int,
    r1_left: int, r1_right: int,
    r2_left: int, r2_right: int,
    verbose: bool,
) -> Tuple[int, int, int, int]:
    cpg_processed = 0
    cpg_covered = 0
    ch_processed = 0
    ch_covered = 0

    cpg_positions: List[int] = []
    ch_positions: List[Tuple[int, str, str]] = []
    if context_mode in ("cg", "both"):
        cpg_positions = find_cpg_positions(reference, chrom)
    if context_mode in ("ch", "both"):
        ch_positions = find_ch_positions(reference, chrom)

    if not cpg_positions and not ch_positions:
        if verbose:
            print(f"[methy-caller] chrom={chrom} no positions, skipping")
        return 0, 0, 0, 0

    batches = build_batches(cpg_positions, ch_positions, batch_size)
    if verbose:
        print(
            f"[methy-caller] chrom={chrom} cpg-sites={len(cpg_positions)} "
            f"ch-sites={len(ch_positions)} batches={len(batches)}"
        )

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for batch_idx, (start, end, cpg_set, ch_map) in enumerate(batches):
            cg_counts: CgCounts = defaultdict(dict)
            ch_counts: ChCounts = defaultdict(dict)

            pileup_iter = bam.pileup(
                chrom, start, end,
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
                if cpg_set:
                    count_cg_column(pileup_col, cpg_set, cb_tag,
                                    cg_counts, r1_left, r1_right,
                                    r2_left, r2_right)
                if ch_map:
                    count_ch_column(pileup_col, ch_map, cb_tag,
                                    ch_counts, r1_left, r1_right,
                                    r2_left, r2_right)

            if cg_counts:
                flush_batch_cg(cg_counts, chrom, out_dir)
                for _, pos_map in cg_counts.items():
                    cpg_processed += len(pos_map)
                    cpg_covered += sum(
                        1 for v in pos_map.values()
                        if v[0] + v[1] > 0
                    )
            if ch_counts:
                flush_batch_ch(ch_counts, chrom, out_dir)
                for _, pos_map in ch_counts.items():
                    ch_processed += len(pos_map)
                    ch_covered += sum(
                        1 for v in pos_map.values()
                        if v[0] + v[1] > 0
                    )

            if verbose and batch_idx % 10 == 0:
                print(
                    f"[methy-caller] chrom={chrom} batch={batch_idx + 1}/{len(batches)} "
                    f"region={start}-{end}"
                )

    return cpg_processed, cpg_covered, ch_processed, ch_covered


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    if args.dry_run:
        print(f"[methy-caller] bam={args.bam}")
        print(f"[methy-caller] reference={args.reference}")
        print(f"[methy-caller] out-dir={args.out_dir}")
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

    chromosomes = parse_chromosomes(args.chromosomes)
    validate_chromosomes(args.reference, chromosomes)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[methy-caller] bam={args.bam}")
    print(f"[methy-caller] reference={args.reference}")
    print(f"[methy-caller] chromosomes={len(chromosomes)}")
    print(f"[methy-caller] context-mode={args.context_mode}")
    print(f"[methy-caller] out-dir={args.out_dir}")

    total_cpg_proc = 0
    total_cpg_cov = 0
    total_ch_proc = 0
    total_ch_cov = 0

    from concurrent.futures import ProcessPoolExecutor, as_completed

    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = {}
        for chrom in chromosomes:
            future = executor.submit(
                process_chromosome,
                args.bam, args.reference, chrom, out_dir, args.cb_tag,
                args.context_mode,
                args.min_base_quality, args.min_mapping_quality,
                args.max_depth, args.batch_size,
                args.r1_left_trim, args.r1_right_trim,
                args.r2_left_trim, args.r2_right_trim,
                args.verbose,
            )
            futures[future] = chrom

        for future in as_completed(futures):
            chrom = futures[future]
            try:
                cp, cc, hp, hc = future.result()
            except Exception as exc:
                print(f"[methy-caller] chrom={chrom} failed: {exc}",
                      file=sys.stderr)
                raise
            total_cpg_proc += cp
            total_cpg_cov += cc
            total_ch_proc += hp
            total_ch_cov += hc
            if args.verbose:
                print(
                    f"[methy-caller] chrom={chrom} done "
                    f"cpg-processed={cp} cpg-covered={cc} "
                    f"ch-processed={hp} ch-covered={hc}"
                )

    print(
        f"[methy-caller] done "
        f"cpg-processed={total_cpg_proc} cpg-covered={total_cpg_cov} "
        f"ch-processed={total_ch_proc} ch-covered={total_ch_cov}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
