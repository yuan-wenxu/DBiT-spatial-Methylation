#!/usr/bin/env python3
"""Extract SmC barcodes and split pairs by conversion lineage."""

from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing as mp
import queue
import time
import traceback
from collections import Counter
from contextlib import ExitStack
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
from typing import Any, Iterator, TextIO


FastqRecord = tuple[str, str, str, str]


@dataclass(frozen=True)
class ConversionScore:
    c_matches: int
    t_matches: int
    other_bases: int

    @property
    def watson_mismatches(self) -> int:
        return self.t_matches + self.other_bases

    @property
    def crick_mismatches(self) -> int:
        return self.c_matches + self.other_bases


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Locate SmC-seq barcodes with linker-BC, trim the barcode mate "
            "after insert-left, and split pairs by conversion lineage."
        )
    )
    parser.add_argument(
        "--barcode-fastq",
        required=True,
        help="Barcode-bearing FASTQ(.gz).",
    )
    parser.add_argument(
        "--genomic-fastq",
        required=True,
        help="Long genomic mate FASTQ(.gz).",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--chunks",
        type=int,
        required=True,
        help="Number of parallel output chunks.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        required=True,
        help="Read pairs dispatched to one worker at a time.",
    )
    parser.add_argument(
        "--compression",
        choices=("gzip", "none"),
        required=True,
        help="Write gzip-compressed or uncompressed FASTQ output.",
    )
    parser.add_argument(
        "--linker-bc",
        required=True,
        help="Sequence between barcode2 and barcode1.",
    )
    parser.add_argument(
        "--insert-left",
        required=True,
        help="Anchor immediately before genome.",
    )
    parser.add_argument(
        "--anchor-edit-distance",
        type=int,
        required=True,
        help="Allowed mismatches outside expected conversion states.",
    )
    parser.add_argument(
        "--max-conversion-mismatches",
        type=int,
        required=True,
        help="Maximum C/T-position mismatches for one lineage.",
    )
    parser.add_argument(
        "--minimum-score-margin",
        type=int,
        required=True,
        help="Required C-versus-T evidence margin in each anchor.",
    )
    parser.add_argument(
        "--barcode-whitelist",
        required=True,
        help=(
            "Barcode whitelist, one sequence per line. Each barcode must "
            "match exactly or correct to the unique nearest whitelist "
            "entry; otherwise the pair is written to discarded output."
        ),
    )
    parser.add_argument(
        "--barcode-hamming-distance",
        type=int,
        required=True,
        help="Maximum Hamming distance for barcode correction.",
    )
    parser.add_argument(
        "--gzip-level", type=int, required=True, help="Output gzip level."
    )
    parser.add_argument(
        "--progress-reads",
        type=int,
        required=True,
        help="Progress interval per chunk; 0 disables it.",
    )
    return parser.parse_args(argv)


def open_text(path: str | Path, mode: str, gzip_level: int = 1) -> TextIO:
    path_text = str(path)
    if path_text.endswith(".gz"):
        if "w" in mode:
            return gzip.open(
                path_text,
                mode + "t",
                encoding="utf-8",
                compresslevel=gzip_level,
            )
        return gzip.open(path_text, mode + "t", encoding="utf-8")
    return open(path_text, mode, encoding="utf-8")


def fastq_iter(handle: TextIO) -> Iterator[FastqRecord]:
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        if not (sequence and plus and quality):
            raise ValueError("incomplete FASTQ record")
        record = tuple(
            value.rstrip("\r\n")
            for value in (header, sequence, plus, quality)
        )
        if not record[0].startswith("@") or not record[2].startswith("+"):
            raise ValueError(f"invalid FASTQ record: {record[0]!r}")
        if len(record[1]) != len(record[3]):
            raise ValueError(f"sequence/quality length mismatch: {record[0]!r}")
        yield record  # type: ignore[misc]


def read_name(header: str) -> str:
    return (
        header.split(maxsplit=1)[0]
        .removeprefix("@")
        .removesuffix("/1")
        .removesuffix("/2")
    )


def read_insert_left(path: str | Path) -> str:
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.split()
            if len(fields) >= 2 and fields[0].lower() == "insert-left":
                return fields[1].upper()
    raise ValueError(f"insert-left row not found in linker table: {path}")


def read_whitelist(path: str) -> tuple[set[str], int]:
    barcodes: set[str] = set()
    with open_text(path, "r") as handle:
        for raw_line in handle:
            barcode = raw_line.strip().upper()
            if barcode:
                barcodes.add(barcode)
    if not barcodes:
        raise ValueError(f"empty barcode whitelist: {path}")
    barcode_lengths = sorted({len(barcode) for barcode in barcodes})
    if len(barcode_lengths) != 1:
        raise ValueError(
            "whitelist contains inconsistent barcode lengths: "
            f"{barcode_lengths}"
        )
    return barcodes, barcode_lengths[0]


def hamming_distance(a: str, b: str, stop_at: int | None = None) -> int:
    if len(a) != len(b):
        raise ValueError("hamming distance requires equal-length strings")
    distance = 0
    for left, right in zip(a, b):
        if left != right:
            distance += 1
            if stop_at is not None and distance > stop_at:
                return distance
    return distance


def build_correction_table(whitelist: set[str]) -> dict[str, str | None]:
    table: dict[str, str | None] = {seq: seq for seq in whitelist}
    for seq in whitelist:
        for pos, original in enumerate(seq):
            prefix, suffix = seq[:pos], seq[pos + 1 :]
            for base in "ACGT":
                if base == original:
                    continue
                variant = prefix + base + suffix
                if variant in whitelist:
                    continue
                if variant in table:
                    if table[variant] is not None:
                        table[variant] = None
                else:
                    table[variant] = seq
    return table


def match_whitelist(
    observed: str,
    whitelist: set[str],
    correction_table: dict[str, str | None] | None,
    max_hamming: int,
) -> str | None:
    if observed in whitelist:
        return observed
    if max_hamming <= 0:
        return None
    if max_hamming == 1:
        if correction_table is None:
            raise RuntimeError(
                "distance-1 matching requires a correction table"
            )
        return correction_table.get(observed)
    best: str | None = None
    best_distance = max_hamming + 1
    tied = False
    for candidate in whitelist:
        if len(candidate) != len(observed):
            continue
        distance = hamming_distance(observed, candidate, stop_at=best_distance)
        if distance < best_distance:
            best = candidate
            best_distance = distance
            tied = False
        elif distance == best_distance:
            tied = True
    if best is None or best_distance > max_hamming or tied:
        return None
    return best


def allowed_bases(expected_base: str) -> str:
    if expected_base == "C":
        return "CT"
    if expected_base == "G":
        return "GA"
    return expected_base


def anchor_mismatches(observed: str, expected: str) -> int:
    if len(observed) != len(expected):
        raise ValueError("anchor lengths differ")
    return sum(
        observed_base not in allowed_bases(expected_base)
        for observed_base, expected_base in zip(observed, expected)
    )


def find_conversion_anchor(
    sequence: str,
    expected: str,
    window_start: int,
    window_end: int,
    max_mismatches: int,
) -> tuple[int, int] | None:
    """Find a unique best C/T- and G/A-compatible fixed-length anchor."""
    window_start = max(0, window_start)
    window_end = min(len(sequence), window_end)
    length = len(expected)
    if not expected or window_start + length > window_end:
        return None

    best_distance = max_mismatches + 1
    best_positions: list[int] = []
    for position in range(window_start, window_end - length + 1):
        distance = anchor_mismatches(
            sequence[position : position + length], expected
        )
        if distance < best_distance:
            best_distance = distance
            best_positions = [position]
        elif distance == best_distance:
            best_positions.append(position)
    if best_distance > max_mismatches or len(best_positions) != 1:
        return None
    start = best_positions[0]
    return start, start + length


def conversion_score(observed: str, expected: str) -> ConversionScore:
    c_matches = 0
    t_matches = 0
    other_bases = 0
    for observed_base, expected_base in zip(observed, expected):
        if expected_base != "C":
            continue
        if observed_base == "C":
            c_matches += 1
        elif observed_base == "T":
            t_matches += 1
        else:
            other_bases += 1
    return ConversionScore(c_matches, t_matches, other_bases)


def classify_score(
    score: ConversionScore,
    max_mismatches: int,
    minimum_margin: int,
) -> str:
    if (
        score.watson_mismatches <= max_mismatches
        and score.c_matches - score.t_matches >= minimum_margin
    ):
        return "watson"
    if (
        score.crick_mismatches <= max_mismatches
        and score.t_matches - score.c_matches >= minimum_margin
    ):
        return "crick"
    return "ambiguous"


def consensus_class(linker_class: str, insert_class: str) -> str:
    informative = {linker_class, insert_class} - {"ambiguous"}
    if len(informative) == 1:
        return informative.pop()
    return "ambiguous"


def annotate_header(
    header: str, barcode1: str, barcode2: str, strand: str
) -> str:
    parts = header.split(" ", 1)
    first = parts[0].removeprefix("@")
    annotated = f"@{barcode1}+{barcode2}:{strand}:{first}"
    return annotated if len(parts) == 1 else f"{annotated} {parts[1]}"


def write_record(handle: TextIO, record: FastqRecord) -> None:
    handle.write("\n".join(record) + "\n")


GROUPS = ("watson", "crick", "ambiguous", "discarded")
REJECT_NAMES = (
    "short_r1",
    "structure_mismatch",
    "linker2_not_found",
    "tn5_not_found",
    "barcode1_not_in_whitelist",
    "barcode2_not_in_whitelist",
)
STRAND_AMBIGUITY_NAMES = (
    "ambiguous_conversion",
    "anchor_class_conflict",
)
COUNT_KEYS = (
    "input_pairs",
    "watson_pairs",
    "crick_pairs",
    "ambiguous_pairs",
    "discarded_pairs",
    *REJECT_NAMES,
    *STRAND_AMBIGUITY_NAMES,
)


@dataclass(frozen=True)
class SplitConfig:
    linker_bc: str
    insert_left: str
    barcode_length: int
    anchor_edit_distance: int
    max_conversion_mismatches: int
    minimum_score_margin: int
    minimum_length: int
    whitelist: frozenset[str]
    correction_table: dict[str, str | None] | None
    barcode_hamming_distance: int


@dataclass(frozen=True)
class SplitResult:
    group: str
    first_record: FastqRecord
    genomic_record: FastqRecord
    reject_reasons: tuple[str, ...]


def chunk_output_paths(
    output_dir: Path, chunk_index: int, compression: str
) -> dict[str, Path]:
    chunk = f"{chunk_index:04d}"
    suffix = ".fastq.gz" if compression == "gzip" else ".fastq"
    return {
        "watson_first": output_dir / f"{chunk}.watson.short-genomic{suffix}",
        "watson_genomic": output_dir / f"{chunk}.watson.genomic{suffix}",
        "crick_first": output_dir / f"{chunk}.crick.short-genomic{suffix}",
        "crick_genomic": output_dir / f"{chunk}.crick.genomic{suffix}",
        "ambiguous_first": output_dir / f"{chunk}.ambiguous.barcode-raw{suffix}",
        "ambiguous_genomic": output_dir / f"{chunk}.ambiguous.genomic{suffix}",
        "discarded_first": output_dir / f"{chunk}.discarded.barcode-raw{suffix}",
        "discarded_genomic": output_dir / f"{chunk}.discarded.genomic{suffix}",
        "stats": output_dir / f"{chunk}.stats.json",
    }


def open_chunk_outputs(
    paths: dict[str, Path], gzip_level: int
) -> tuple[ExitStack, dict[str, TextIO]]:
    stack = ExitStack()
    try:
        handles = {
            key: stack.enter_context(open_text(path, "w", gzip_level))
            for key, path in paths.items()
            if key != "stats"
        }
    except BaseException:
        stack.close()
        raise
    return stack, handles


def rejected_result(
    reasons: tuple[str, ...],
    barcode_record: FastqRecord,
    genomic_record: FastqRecord,
    barcode1: str | None = None,
    barcode2: str | None = None,
) -> SplitResult:
    if barcode1 is None or barcode2 is None:
        return SplitResult(
            "discarded",
            barcode_record,
            genomic_record,
            reasons,
        )
    first: FastqRecord = (
        annotate_header(barcode_record[0], barcode1, barcode2, "D"),
        barcode_record[1],
        barcode_record[2],
        barcode_record[3],
    )
    genomic: FastqRecord = (
        annotate_header(genomic_record[0], barcode1, barcode2, "D"),
        genomic_record[1],
        genomic_record[2],
        genomic_record[3],
    )
    return SplitResult("discarded", first, genomic, reasons)


def split_pair(
    records: tuple[FastqRecord, FastqRecord], config: SplitConfig
) -> SplitResult:
    barcode_record, genomic_record = records
    sequence = barcode_record[1].upper()
    group = "ambiguous"
    reject_reasons: tuple[str, ...] = ()
    insert_end: int | None = None

    if len(sequence) < config.minimum_length:
        return rejected_result(
            ("short_r1",), barcode_record, genomic_record
        )
    linker_position = find_conversion_anchor(
        sequence,
        config.linker_bc,
        0,
        min(
            len(sequence),
            config.barcode_length + len(config.linker_bc) + 50,
        ),
        config.anchor_edit_distance,
    )
    if linker_position is None:
        return rejected_result(
            ("structure_mismatch", "linker2_not_found"),
            barcode_record,
            genomic_record,
        )
    linker_start, linker_end = linker_position
    barcode2_start = linker_start - config.barcode_length
    barcode1_end = linker_end + config.barcode_length
    if barcode2_start < 0 or barcode1_end > len(sequence):
        return rejected_result(
            ("structure_mismatch",), barcode_record, genomic_record
        )
    barcode2 = sequence[barcode2_start:linker_start]
    barcode1 = sequence[linker_end:barcode1_end]
    barcode1_fixed = match_whitelist(
        barcode1,
        config.whitelist,
        config.correction_table,
        config.barcode_hamming_distance,
    )
    if barcode1_fixed is None:
        return rejected_result(
            ("barcode1_not_in_whitelist",),
            barcode_record,
            genomic_record,
            barcode1,
            barcode2,
        )
    barcode2_fixed = match_whitelist(
        barcode2,
        config.whitelist,
        config.correction_table,
        config.barcode_hamming_distance,
    )
    if barcode2_fixed is None:
        return rejected_result(
            ("barcode2_not_in_whitelist",),
            barcode_record,
            genomic_record,
            barcode1_fixed,
            barcode2,
        )
    barcode1 = barcode1_fixed
    barcode2 = barcode2_fixed
    linker_class = classify_score(
        conversion_score(
            sequence[linker_start:linker_end], config.linker_bc
        ),
        config.max_conversion_mismatches,
        config.minimum_score_margin,
    )
    insert_position = find_conversion_anchor(
        sequence,
        config.insert_left,
        barcode1_end,
        barcode1_end + len(config.insert_left) + 40,
        config.anchor_edit_distance,
    )
    if insert_position is None:
        reject_reasons = (
            "structure_mismatch",
            "tn5_not_found",
        )
    else:
        insert_start, insert_end = insert_position
        if insert_end >= len(sequence):
            reject_reasons = (
                "structure_mismatch",
                "tn5_not_found",
            )
        else:
            insert_class = classify_score(
                conversion_score(
                    sequence[insert_start:insert_end],
                    config.insert_left,
                ),
                config.max_conversion_mismatches,
                config.minimum_score_margin,
            )
            group = consensus_class(linker_class, insert_class)
            if (
                linker_class != "ambiguous"
                and insert_class != "ambiguous"
                and linker_class != insert_class
            ):
                reject_reasons = ("anchor_class_conflict",)
            elif group == "ambiguous":
                reject_reasons = ("ambiguous_conversion",)

    class_code = {
        "watson": "W",
        "crick": "C",
        "ambiguous": "U",
        "discarded": "D",
    }[group]
    output_first: FastqRecord = (
        annotate_header(barcode_record[0], barcode1, barcode2, class_code),
        barcode_record[1],
        barcode_record[2],
        barcode_record[3],
    )
    output_genomic: FastqRecord = (
        annotate_header(genomic_record[0], barcode1, barcode2, class_code),
        genomic_record[1],
        genomic_record[2],
        genomic_record[3],
    )
    if group in {"watson", "crick"}:
        if insert_end is None:
            raise RuntimeError("classified read lacks insert-left end")
        output_first = (
            output_first[0],
            barcode_record[1][insert_end:],
            output_first[2],
            barcode_record[3][insert_end:],
        )
    return SplitResult(
        group,
        output_first,
        output_genomic,
        reject_reasons,
    )


def run_chunk_worker(
    chunk_index: int,
    batch_queue: Any,
    result_queue: Any,
    output_dir_text: str,
    compression: str,
    gzip_level: int,
    progress_reads: int,
    config: SplitConfig,
) -> None:
    output_dir = Path(output_dir_text)
    paths = chunk_output_paths(output_dir, chunk_index, compression)
    counts = Counter({key: 0 for key in COUNT_KEYS})
    group_reason_counts = {group: Counter() for group in GROUPS}
    stack: ExitStack | None = None
    started_at = time.monotonic()
    last_report_at = started_at
    reads_at_last_report = 0
    next_progress_report = progress_reads

    try:
        stack, handles = open_chunk_outputs(paths, gzip_level)
        while True:
            batch = batch_queue.get()
            if batch is None:
                break
            for records in batch:
                result = split_pair(records, config)
                counts["input_pairs"] += 1
                for reason in result.reject_reasons:
                    counts[reason] += 1
                    group_reason_counts[result.group][reason] += 1
                counts[f"{result.group}_pairs"] += 1
                write_record(handles[f"{result.group}_first"], result.first_record)
                write_record(
                    handles[f"{result.group}_genomic"], result.genomic_record
                )

                if progress_reads > 0 and counts["input_pairs"] >= next_progress_report:
                    now = time.monotonic()
                    elapsed = max(now - last_report_at, 1e-9)
                    speed = (
                        counts["input_pairs"] - reads_at_last_report
                    ) / elapsed
                    print(
                        f"[extract-bc] mode=smc chunk={chunk_index:04d} "
                        f"reads={counts['input_pairs']} "
                        f"watson={counts['watson_pairs']} "
                        f"crick={counts['crick_pairs']} "
                        f"ambiguous={counts['ambiguous_pairs']} "
                        f"discarded={counts['discarded_pairs']} "
                        f"speed={speed:.0f} reads/s",
                        flush=True,
                    )
                    last_report_at = now
                    reads_at_last_report = counts["input_pairs"]
                    next_progress_report += progress_reads
    except BaseException:
        result_queue.put(
            {
                "ok": False,
                "chunk": chunk_index,
                "error": traceback.format_exc(),
            }
        )
        raise
    finally:
        if stack is not None:
            stack.close()

    total = counts["input_pairs"]
    informative = counts["watson_pairs"] + counts["crick_pairs"]
    kept = total - counts["discarded_pairs"]
    row = {
        "chunk": f"{chunk_index:04d}",
        "total_reads": total,
        "kept_reads": kept,
        "discarded_reads": counts["discarded_pairs"],
        "kept_fraction": kept / total if total else 0.0,
        "reject_counts": {
            name: group_reason_counts["discarded"][name]
            for name in REJECT_NAMES
        },
        "watson_reads": counts["watson_pairs"],
        "crick_reads": counts["crick_pairs"],
        "informative_reads": informative,
        "ambiguous_reads": counts["ambiguous_pairs"],
        "informative_fraction": informative / total if total else 0.0,
        "strand_ambiguity_counts": {
            name: group_reason_counts["ambiguous"][name]
            for name in (
                "structure_mismatch",
                "tn5_not_found",
                *STRAND_AMBIGUITY_NAMES,
            )
        },
        "output_r1_watson": paths["watson_first"].name,
        "output_r2_watson": paths["watson_genomic"].name,
        "output_r1_crick": paths["crick_first"].name,
        "output_r2_crick": paths["crick_genomic"].name,
        "output_r1_ambiguous": paths["ambiguous_first"].name,
        "output_r2_ambiguous": paths["ambiguous_genomic"].name,
        "output_r1_discarded": paths["discarded_first"].name,
        "output_r2_discarded": paths["discarded_genomic"].name,
        "stats": paths["stats"].name,
    }
    paths["stats"].write_text(json.dumps(row, indent=2), encoding="utf-8")
    result_queue.put({"ok": True, "row": row, "counts": dict(counts)})


def put_batch(
    batch_queues: list[Any],
    chunk_processes: list[mp.Process],
    chunk_index: int,
    batch: object,
) -> None:
    while True:
        try:
            batch_queues[chunk_index].put(batch, timeout=1)
            return
        except queue.Full:
            failed = [
                process
                for process in chunk_processes
                if process.exitcode not in (None, 0)
            ]
            if failed:
                raise RuntimeError(
                    f"chunk worker failed: {failed[0].name} "
                    f"(exit code {failed[0].exitcode})"
                )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    for name in (
        "anchor_edit_distance",
        "max_conversion_mismatches",
        "minimum_score_margin",
        "barcode_hamming_distance",
        "progress_reads",
    ):
        if getattr(args, name) < 0:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 0")
    if args.chunks < 1:
        raise ValueError("--chunks must be >= 1")
    if args.batch_size < 1:
        raise ValueError("--batch-size must be >= 1")
    if not 0 <= args.gzip_level <= 9:
        raise ValueError("--gzip-level must be between 0 and 9")

    linker_bc = args.linker_bc.upper()
    insert_left = (args.insert_left.upper())

    if not linker_bc or not insert_left:
        raise ValueError("linker-BC and insert-left must not be empty")

    whitelist, barcode_length = read_whitelist(args.barcode_whitelist)
    correction_table: dict[str, str | None] | None = None
    if args.barcode_hamming_distance == 1:
        correction_table = build_correction_table(whitelist)

    config = SplitConfig(
        linker_bc=linker_bc,
        insert_left=insert_left,
        barcode_length=barcode_length,
        anchor_edit_distance=args.anchor_edit_distance,
        max_conversion_mismatches=args.max_conversion_mismatches,
        minimum_score_margin=args.minimum_score_margin,
        minimum_length=2 * barcode_length + len(linker_bc),
        whitelist=frozenset(whitelist),
        correction_table=correction_table,
        barcode_hamming_distance=args.barcode_hamming_distance,
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    context = mp.get_context("fork")
    result_queue = context.Queue()
    batch_queues = [context.Queue(maxsize=2) for _ in range(args.chunks)]
    chunk_processes = [
        context.Process(
            target=run_chunk_worker,
            args=(
                index + 1,
                batch_queues[index],
                result_queue,
                str(output_dir),
                args.compression,
                args.gzip_level,
                args.progress_reads,
                config,
            ),
            name=f"smc-split-chunk-{index + 1:04d}",
        )
        for index in range(args.chunks)
    ]
    started_at = time.monotonic()
    input_pairs = 0
    batch_count = 0
    batch: list[tuple[FastqRecord, FastqRecord]] = []

    for process in chunk_processes:
        process.start()
    try:
        with open_text(args.barcode_fastq, "r") as barcode_input, open_text(
            args.genomic_fastq, "r"
        ) as genomic_input:
            for barcode_record, genomic_record in zip_longest(
                fastq_iter(barcode_input), fastq_iter(genomic_input)
            ):
                if barcode_record is None or genomic_record is None:
                    raise ValueError("paired FASTQs have different record counts")
                if read_name(barcode_record[0]) != read_name(genomic_record[0]):
                    raise ValueError(
                        "paired FASTQ names differ: "
                        f"{barcode_record[0]!r} vs {genomic_record[0]!r}"
                    )
                batch.append((barcode_record, genomic_record))
                input_pairs += 1
                if len(batch) == args.batch_size:
                    put_batch(
                        batch_queues,
                        chunk_processes,
                        batch_count % args.chunks,
                        batch,
                    )
                    batch = []
                    batch_count += 1
            if batch:
                put_batch(
                    batch_queues,
                    chunk_processes,
                    batch_count % args.chunks,
                    batch,
                )
                batch_count += 1
        for index in range(args.chunks):
            put_batch(batch_queues, chunk_processes, index, None)
        for process in chunk_processes:
            process.join()
    except BaseException:
        for process in chunk_processes:
            if process.is_alive():
                process.terminate()
        for process in chunk_processes:
            process.join()
        raise
    finally:
        for batch_queue in batch_queues:
            batch_queue.close()

    messages: list[dict[str, object]] = []
    while len(messages) < args.chunks:
        try:
            messages.append(result_queue.get(timeout=1))
        except queue.Empty:
            break
    result_queue.close()
    failures = [message for message in messages if not message.get("ok")]
    failed_processes = [
        process for process in chunk_processes if process.exitcode != 0
    ]
    if failures or failed_processes or len(messages) != args.chunks:
        details = "\n".join(
            str(message.get("error", "")) for message in failures
        )
        raise RuntimeError(
            f"strand chunk processing failed; received {len(messages)}/"
            f"{args.chunks} results\n{details}"
        )

    messages.sort(key=lambda message: str(dict(message["row"])["chunk"]))
    chunk_rows = [dict(message["row"]) for message in messages]
    total_counts = {
        key: sum(int(dict(message["counts"])[key]) for message in messages)
        for key in COUNT_KEYS
    }
    if total_counts["input_pairs"] != input_pairs:
        raise RuntimeError(
            "processed read count mismatch: "
            f"input={input_pairs}, processed={total_counts['input_pairs']}"
        )
    classified_pairs = sum(
        total_counts[f"{group}_pairs"] for group in GROUPS
    )
    if classified_pairs != input_pairs:
        raise RuntimeError(
            "output group count mismatch: "
            f"input={input_pairs}, grouped={classified_pairs}"
        )
    informative = total_counts["watson_pairs"] + total_counts["crick_pairs"]
    kept = input_pairs - total_counts["discarded_pairs"]
    stats = {
        "mode": "smc",
        "input_r1": Path(args.barcode_fastq).name,
        "input_r2": Path(args.genomic_fastq).name,
        "output_dir": str(output_dir),
        "compression": args.compression,
        "chunk_count": args.chunks,
        "batch_size": args.batch_size,
        "batch_count": batch_count,
        "chunks": chunk_rows,
        "total_reads": input_pairs,
        "kept_reads": kept,
        "discarded_reads": total_counts["discarded_pairs"],
        "kept_fraction": kept / input_pairs if input_pairs else 0.0,
        "reject_counts": {
            name: sum(
                int(dict(row["reject_counts"])[name]) for row in chunk_rows
            )
            for name in REJECT_NAMES
        },
        "watson_reads": total_counts["watson_pairs"],
        "crick_reads": total_counts["crick_pairs"],
        "informative_reads": informative,
        "ambiguous_reads": total_counts["ambiguous_pairs"],
        "informative_fraction": informative / input_pairs if input_pairs else 0.0,
        "strand_ambiguity_counts": {
            name: sum(
                int(dict(row["strand_ambiguity_counts"])[name])
                for row in chunk_rows
            )
            for name in (
                "structure_mismatch",
                "tn5_not_found",
                *STRAND_AMBIGUITY_NAMES,
            )
        },
    }
    stats_path = output_dir / "stats.json"
    stats_path.write_text(json.dumps(stats, indent=2), encoding="utf-8")

    elapsed = max(time.monotonic() - started_at, 1e-9)
    print("[extract-bc] mode=smc")
    print(f"[extract-bc] input-r1={args.barcode_fastq}")
    print(f"[extract-bc] input-r2={args.genomic_fastq}")
    print(f"[extract-bc] output-dir={output_dir}")
    print(f"[extract-bc] chunks={args.chunks}")
    print(f"[extract-bc] batch-size={args.batch_size}")
    print(f"[extract-bc] batch-count={batch_count}")
    print(f"[extract-bc] compression={args.compression}")
    print(f"[extract-bc] stats={stats_path}")
    print(
        f"[extract-bc] informative={informative}/{input_pairs} "
        f"informative-rate={informative / input_pairs if input_pairs else 0.0:.4f} "
        f"watson={total_counts['watson_pairs']} "
        f"crick={total_counts['crick_pairs']} "
        f"ambiguous={total_counts['ambiguous_pairs']} "
        f"discarded={total_counts['discarded_pairs']} "
        f"avg-speed={input_pairs / elapsed:.1f} reads/s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
