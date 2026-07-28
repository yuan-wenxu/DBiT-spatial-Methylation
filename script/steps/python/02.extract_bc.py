#!/usr/bin/env python3
"""Extract DBiT barcodes in TAPS, TAPS-v2, or EMSeq mode."""

from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing as mp
import queue
import re
import time
import traceback
from contextlib import ExitStack
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
from typing import Any, Iterator, TextIO

from fuzzysearch import find_near_matches


FastqRecord = tuple[str, str, str, str]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract DBiT barcodes with standard TAPS matching, TAPS-v2 "
            "C/T-aware matching, or EMSeq C/T-aware matching."
        )
    )
    parser.add_argument(
        "--mode",
        "--assay",
        dest="mode",
        required=True,
        choices=("taps", "taps-v2", "emseq"),
        help="Barcode extraction mode.",
    )
    parser.add_argument("r1", help="Input R1 FASTQ(.gz).")
    parser.add_argument("r2", help="Input R2 FASTQ(.gz).")
    parser.add_argument(
        "-b1",
        "--barcode1-whitelist",
        required=True,
        help="Whitelist for barcode1, one barcode per line.",
    )
    parser.add_argument(
        "-b2",
        "--barcode2-whitelist",
        required=True,
        help="Whitelist for barcode2, one barcode per line.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help=(
            "Output directory for numbered demux/spike-in FASTQ chunks and stats.json."
        ),
    )
    parser.add_argument(
        "--chunks",
        type=int,
        default=1,
        help="Number of output chunks. Default: 1.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50_000,
        help="Read pairs dispatched to one chunk at a time. Default: 50000.",
    )
    parser.add_argument(
        "--compression",
        choices=("gzip", "none"),
        default="gzip",
        help="Compress FASTQ while writing. Default: gzip.",
    )
    parser.add_argument(
        "--linker-bc",
        required=True,
        help="Expected linker sequence between barcode2 and barcode1.",
    )
    parser.add_argument(
        "--insert-left",
        required=True,
        help="Expected insert-left anchor (methylated linker in taps-v2 mode).",
    )
    parser.add_argument(
        "--linker-edit-distance",
        type=int,
        default=1,
        help="Maximum edit distance for linker matching. Default: 1.",
    )
    parser.add_argument(
        "--barcode-hamming-distance",
        type=int,
        default=1,
        help="Maximum Hamming distance for barcode correction. Default: 1.",
    )
    parser.add_argument(
        "--insert-left-edit-distance",
        type=int,
        default=1,
        help=(
            "Maximum mismatches at non-C positions when locating the C/T-aware "
            "insert-left in taps-v2 or emseq mode. Original C positions must "
            "be C or T. Default: 1."
        ),
    )
    parser.add_argument(
        "--methylated-c-positions",
        default="3,6,10",
        help=(
            "Comma-separated 0-based methylated C positions in insert-left "
            "for taps-v2 conversion QC. Default: 3,6,10."
        ),
    )
    parser.add_argument(
        "--progress-reads",
        type=int,
        default=1_000_000,
        help="Report progress every N processed read pairs per chunk; 0 disables it. Default: 1000000.",
    )
    parser.add_argument(
        "--gzip-level",
        type=int,
        default=6,
        help="Python gzip compression level, 0-9. Default: 6.",
    )
    return parser.parse_args(argv)


def open_text(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def open_output_text(path: Path, compression: str, gzip_level: int) -> TextIO:
    if compression == "gzip":
        return gzip.open(path, "wt", encoding="utf-8", compresslevel=gzip_level)
    return path.open("w", encoding="utf-8")


def read_whitelist(path: str) -> set[str]:
    barcodes: set[str] = set()
    with open_text(path, "r") as handle:
        for raw_line in handle:
            barcode = raw_line.strip().upper()
            if barcode:
                barcodes.add(barcode)
    if not barcodes:
        raise ValueError(f"empty barcode whitelist: {path}")
    return barcodes


def fastq_iter(handle: TextIO) -> Iterator[FastqRecord]:
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        if not (sequence and plus and quality):
            raise ValueError("incomplete FASTQ record detected")
        yield (
            header.rstrip("\n"),
            sequence.rstrip("\n"),
            plus.rstrip("\n"),
            quality.rstrip("\n"),
        )


def annotate_header(header: str, barcode1: str, barcode2: str) -> str:
    if not header.startswith("@"):
        raise ValueError(f"invalid FASTQ header: {header}")
    parts = header.split(" ", 1)
    parts[0] = f"@{barcode1}+{barcode2}:{parts[0][1:]}"
    return parts[0] if len(parts) == 1 else f"{parts[0]} {parts[1]}"


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


def find_exact_positions(
    text: str, pattern: str, start_pos: int, end_pos: int
) -> list[int]:
    positions: list[int] = []
    position = text.find(pattern, start_pos, end_pos)
    while position >= 0:
        positions.append(position)
        position = text.find(pattern, position + 1, end_pos)
    return positions


def find_unique_pattern_in_window(
    text: str,
    pattern: str,
    start_pos: int,
    end_pos: int,
    max_edit_distance: int,
) -> int | None:
    start_pos = max(0, start_pos)
    end_pos = min(len(text), end_pos)
    if start_pos >= end_pos:
        return None
    exact_positions = find_exact_positions(text, pattern, start_pos, end_pos)
    if len(exact_positions) == 1:
        return exact_positions[0]
    if len(exact_positions) > 1:
        return None
    if max_edit_distance <= 0:
        return None
    near_matches = find_near_matches(
        pattern,
        text[start_pos:end_pos],
        max_l_dist=max_edit_distance,
    )
    if not near_matches:
        return None
    best_distance = min(match.dist for match in near_matches)
    best_matches = [
        match for match in near_matches if match.dist == best_distance
    ]
    return start_pos + best_matches[0].start if len(best_matches) == 1 else None


def match_whitelist(
    observed: str, whitelist: set[str], max_hamming: int
) -> str | None:
    if observed in whitelist:
        return observed
    if max_hamming <= 0:
        return None
    best: str | None = None
    best_distance = max_hamming + 1
    tied = False
    for candidate in whitelist:
        if len(candidate) != len(observed):
            continue
        distance = hamming_distance(observed, candidate, stop_at=best_distance - 1)
        if distance < best_distance:
            best = candidate
            best_distance = distance
            tied = False
        elif distance == best_distance:
            tied = True
    if best is None or best_distance > max_hamming or tied:
        return None
    return best


def find_standard_insert(
    sequence: str,
    insert_left: str,
    after_pos: int,
    max_edit_distance: int,
    window_span: int = 40,
) -> tuple[int, int] | None:
    seed = insert_left[-15:] if len(insert_left) > 15 else insert_left
    offset = len(insert_left) - len(seed)
    window_start = after_pos + offset
    window_end = min(len(sequence), window_start + window_span)
    if window_start >= window_end:
        return None
    seed_position = find_unique_pattern_in_window(
        sequence,
        seed,
        window_start,
        window_end,
        max_edit_distance,
    )
    if seed_position is None:
        return None
    insert_start = seed_position - offset
    insert_end = insert_start + len(insert_left)
    if insert_start < 0 or insert_end > len(sequence):
        return None
    return insert_start, insert_end


def c_t_pattern(sequence: str) -> str:
    return "".join("[CT]" if base == "C" else base for base in sequence.upper())


def convert_c_to_t(sequence: str) -> str:
    return sequence.upper().replace("C", "T")


def c_t_mismatch_count(observed: str, expected: str) -> int | None:
    if len(observed) != len(expected):
        raise ValueError("observed and expected sequence length mismatch")
    mismatches = 0
    for observed_base, expected_base in zip(observed, expected):
        if expected_base == "C":
            if observed_base not in "CT":
                return None
        elif observed_base != expected_base:
            mismatches += 1
    return mismatches


def find_c_t_pattern_in_window(
    sequence: str,
    expected: str,
    window_start: int,
    window_end: int,
    max_mismatch: int,
) -> tuple[int, int] | None:
    pattern_length = len(expected)
    window_start = max(0, window_start)
    window_end = min(len(sequence), window_end)
    if pattern_length == 0 or window_start + pattern_length > window_end:
        return None

    regex = re.compile(
        "".join("[CT]" if base == "C" else re.escape(base) for base in expected)
    )
    window = sequence[window_start:window_end]
    exact_hits = [window_start + match.start() for match in regex.finditer(window)]
    exact_hits = [
        position
        for position in exact_hits
        if position + pattern_length <= window_end
    ]
    if len(exact_hits) == 1:
        position = exact_hits[0]
        return position, position + pattern_length
    if len(exact_hits) > 1 or max_mismatch <= 0:
        return None

    best_positions: list[int] = []
    best_distance = max_mismatch + 1
    for position in range(window_start, window_end - pattern_length + 1):
        distance = c_t_mismatch_count(
            sequence[position : position + pattern_length], expected
        )
        if distance is None:
            continue
        if distance < best_distance:
            best_distance = distance
            best_positions = [position]
        elif distance == best_distance:
            best_positions.append(position)
    if best_distance > max_mismatch or len(best_positions) != 1:
        return None
    position = best_positions[0]
    return position, position + pattern_length


def c_positions(insert_left: str) -> list[int]:
    return [index for index, base in enumerate(insert_left) if base == "C"]


def parse_methylated_c_positions(raw: str, insert_left: str) -> tuple[int, ...]:
    positions: list[int] = []
    seen: set[int] = set()
    for item in raw.split(","):
        token = item.strip()
        if not token:
            continue
        try:
            position = int(token)
        except ValueError as error:
            raise ValueError(
                f"invalid methylated C position: {token!r}"
            ) from error
        if position < 0 or position >= len(insert_left):
            raise ValueError(
                f"methylated C position {position} is outside insert-left "
                f"length {len(insert_left)}"
            )
        if insert_left[position] != "C":
            raise ValueError(
                f"methylated C position {position} is "
                f"{insert_left[position]!r} in insert-left, expected 'C'"
            )
        if position not in seen:
            positions.append(position)
            seen.add(position)
    if not positions:
        raise ValueError("--methylated-c-positions must contain at least one position")
    return tuple(positions)


def score_conversion(
    observed_insert: str, positions: list[int]
) -> tuple[int, int, list[str]]:
    converted = 0
    retained = 0
    site_bases: list[str] = []
    for position in positions:
        base = observed_insert[position] if position < len(observed_insert) else "."
        site_bases.append(base)
        if base == "T":
            converted += 1
        elif base == "C":
            retained += 1
    return converted, retained, site_bases


def write_record(handle: TextIO, record: FastqRecord) -> None:
    handle.write("\n".join(record) + "\n")


def chunk_output_paths(
    output_dir: Path, chunk_index: int, compression: str
) -> dict[str, Path]:
    chunk = f"{chunk_index:04d}"
    suffix = ".fastq.gz" if compression == "gzip" else ".fastq"
    return {
        "r1_demux": output_dir / f"{chunk}.R1.demux{suffix}",
        "r1_spike": output_dir / f"{chunk}.R1.spike-in{suffix}",
        "r2_demux": output_dir / f"{chunk}.R2.demux{suffix}",
        "r2_spike": output_dir / f"{chunk}.R2.spike-in{suffix}",
        "stats": output_dir / f"{chunk}.stats.json",
    }


def open_chunk_outputs(
    paths: dict[str, Path],
    compression: str,
    gzip_level: int,
) -> tuple[ExitStack, dict[str, TextIO]]:
    stack = ExitStack()
    try:
        handles = {
            key: stack.enter_context(
                open_output_text(paths[key], compression, gzip_level)
            )
            for key in ("r1_demux", "r1_spike", "r2_demux", "r2_spike")
        }
    except BaseException:
        stack.close()
        raise
    return stack, handles


@dataclass(frozen=True)
class ExtractionConfig:
    mode: str
    linker_bc: str
    insert_left: str
    linker_edit_distance: int
    barcode_hamming_distance: int
    insert_left_edit_distance: int
    barcode1_whitelist: frozenset[str]
    barcode2_whitelist: frozenset[str]
    barcode1_length: int
    barcode2_length: int
    minimum_length: int
    all_c_positions: tuple[int, ...]
    methylated_c_positions: tuple[int, ...]


@dataclass(frozen=True)
class ExtractionResult:
    kept: bool
    reject_reasons: tuple[str, ...]
    record1: FastqRecord
    record2: FastqRecord
    observed_insert: str | None = None


_EXTRACTION_CONFIG: ExtractionConfig | None = None


def initialize_extraction_worker(config: ExtractionConfig) -> None:
    global _EXTRACTION_CONFIG
    _EXTRACTION_CONFIG = config


def extract_record(
    records: tuple[FastqRecord, FastqRecord],
) -> ExtractionResult:
    config = _EXTRACTION_CONFIG
    if config is None:
        raise RuntimeError("extraction worker was not initialized")

    record1, record2 = records
    header1, sequence1, plus1, quality1 = record1
    header2, sequence2, plus2, quality2 = record2
    sequence1_upper = sequence1.upper()

    def rejected(*reasons: str, observed: str | None = None) -> ExtractionResult:
        return ExtractionResult(False, reasons, record1, record2, observed)

    if (
        len(sequence1_upper) < config.minimum_length
        or len(quality1) < config.minimum_length
    ):
        return rejected("short_r1")

    linker_window_end = min(
        len(sequence1_upper),
        config.barcode2_length + len(config.linker_bc) + 12 + 50,
    )
    if config.mode == "emseq":
        linker_position = find_c_t_pattern_in_window(
            sequence1_upper,
            config.linker_bc,
            0,
            linker_window_end,
            config.linker_edit_distance,
        )
    else:
        linker_start = find_unique_pattern_in_window(
            sequence1_upper,
            config.linker_bc,
            0,
            linker_window_end,
            config.linker_edit_distance,
        )
        linker_position = (
            None
            if linker_start is None
            else (linker_start, linker_start + len(config.linker_bc))
        )
    if linker_position is None:
        return rejected("structure_mismatch", "linker2_not_found")

    linker_start, linker_end = linker_position
    barcode2_start = linker_start - config.barcode2_length
    barcode1_end = linker_end + config.barcode1_length
    if barcode2_start < 0 or barcode1_end > len(sequence1_upper):
        return rejected("structure_mismatch")

    barcode2_observed = sequence1_upper[barcode2_start:linker_start]
    barcode1_observed = sequence1_upper[linker_end:barcode1_end]
    if config.mode == "taps":
        insert_position = find_standard_insert(
            sequence1_upper,
            config.insert_left,
            barcode1_end,
            config.linker_edit_distance,
        )
    else:
        insert_position = find_c_t_pattern_in_window(
            sequence1_upper,
            config.insert_left,
            barcode1_end,
            barcode1_end + len(config.insert_left) + 40,
            config.insert_left_edit_distance,
        )
    if insert_position is None:
        return rejected("structure_mismatch", "tn5_not_found")

    insert_start, insert_end = insert_position
    observed_insert = sequence1_upper[insert_start:insert_end]
    barcode1 = match_whitelist(
        barcode1_observed,
        config.barcode1_whitelist,
        config.barcode_hamming_distance,
    )
    barcode2 = match_whitelist(
        barcode2_observed,
        config.barcode2_whitelist,
        config.barcode_hamming_distance,
    )
    if barcode1 is None:
        return rejected("barcode1_not_in_whitelist", observed=observed_insert)
    if barcode2 is None:
        return rejected("barcode2_not_in_whitelist", observed=observed_insert)

    return ExtractionResult(
        True,
        (),
        (
            annotate_header(header1, barcode1, barcode2),
            sequence1[insert_end:],
            plus1,
            quality1[insert_end:],
        ),
        (
            annotate_header(header2, barcode1, barcode2),
            sequence2,
            plus2,
            quality2,
        ),
        observed_insert,
    )


def run_chunk_worker(
    chunk_index: int,
    batch_queue: Any,
    result_queue: Any,
    output_dir_text: str,
    compression: str,
    gzip_level: int,
    progress_reads: int,
    config: ExtractionConfig,
) -> None:
    reject_names = (
        "short_r1",
        "structure_mismatch",
        "linker2_not_found",
        "tn5_not_found",
        "barcode1_not_in_whitelist",
        "barcode2_not_in_whitelist",
    )
    reject_counts = {name: 0 for name in reject_names}
    total = 0
    kept = 0
    scored_fully_converted = 0
    scored_inserts = 0
    all_t_methylated_reads = 0
    methylated_rate_sum = 0.0
    methylated_rate_count = 0
    site_converted = [0] * len(config.all_c_positions)
    site_retained = [0] * len(config.all_c_positions)
    site_other = [0] * len(config.all_c_positions)
    output_dir = Path(output_dir_text)
    paths = chunk_output_paths(output_dir, chunk_index, compression)
    stack: ExitStack | None = None
    started_at = time.monotonic()
    last_report_at = started_at
    reads_at_last_report = 0
    next_progress_report = progress_reads
    all_c_position_indexes = {
        position: index for index, position in enumerate(config.all_c_positions)
    }
    methylated_c_indexes = tuple(
        all_c_position_indexes[position]
        for position in config.methylated_c_positions
    )

    try:
        stack, handles = open_chunk_outputs(paths, compression, gzip_level)
        initialize_extraction_worker(config)
        while True:
            batch = batch_queue.get()
            if batch is None:
                break
            for records in batch:
                result = extract_record(records)
                total += 1
                for reason in result.reject_reasons:
                    reject_counts[reason] += 1

                if result.kept:
                    kept += 1
                    if (
                        config.mode == "taps-v2"
                        and result.observed_insert is not None
                    ):
                        scored_inserts += 1
                        _, _, bases = score_conversion(
                            result.observed_insert, list(config.all_c_positions)
                        )
                        for index, base in enumerate(bases):
                            if base == "T":
                                site_converted[index] += 1
                            elif base == "C":
                                site_retained[index] += 1
                            else:
                                site_other[index] += 1
                        if methylated_c_indexes:
                            methylated_bases = [
                                bases[index] for index in methylated_c_indexes
                            ]
                            converted = methylated_bases.count("T")
                            retained = methylated_bases.count("C")
                            informative = converted + retained
                            if informative:
                                methylated_rate_sum += converted / informative
                                methylated_rate_count += 1
                            if all(base == "T" for base in methylated_bases):
                                all_t_methylated_reads += 1
                                scored_fully_converted += 1
                    write_record(handles["r1_demux"], result.record1)
                    write_record(handles["r2_demux"], result.record2)
                else:
                    write_record(handles["r1_spike"], result.record1)
                    write_record(handles["r2_spike"], result.record2)

                if progress_reads > 0 and total >= next_progress_report:
                    now = time.monotonic()
                    elapsed = max(now - last_report_at, 1e-9)
                    speed = (total - reads_at_last_report) / elapsed
                    print(
                        f"[extract-bc] chunk={chunk_index:04d} reads={total} "
                        f"kept={kept} speed={speed:.0f} reads/s",
                        flush=True,
                    )
                    last_report_at = now
                    reads_at_last_report = total
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

    row: dict[str, object] = {
        "chunk": f"{chunk_index:04d}",
        "total_reads": total,
        "kept_reads": kept,
        "spike_in_reads": total - kept,
        "kept_fraction": kept / total if total else 0.0,
        "output_r1_demux": paths["r1_demux"].name,
        "output_r1_spike_in": paths["r1_spike"].name,
        "output_r2_demux": paths["r2_demux"].name,
        "output_r2_spike_in": paths["r2_spike"].name,
        "stats": paths["stats"].name,
        "reject_counts": reject_counts,
    }
    paths["stats"].write_text(json.dumps(row, indent=2), encoding="utf-8")
    result_queue.put(
        {
            "ok": True,
            "row": row,
            "scored_fully_converted": scored_fully_converted,
            "scored_inserts": scored_inserts,
            "all_t_methylated_reads": all_t_methylated_reads,
            "methylated_rate_sum": methylated_rate_sum,
            "methylated_rate_count": methylated_rate_count,
            "site_converted": site_converted,
            "site_retained": site_retained,
            "site_other": site_other,
        }
    )


def build_extraction_config(args: argparse.Namespace) -> ExtractionConfig:
    linker_bc_input = args.linker_bc.upper()
    insert_left_input = args.insert_left.upper()
    linker_bc = linker_bc_input
    insert_left = insert_left_input
    if not linker_bc:
        raise ValueError("--linker-bc must not be empty")
    if not insert_left:
        raise ValueError("--insert-left must not be empty")

    barcode1_whitelist = read_whitelist(args.barcode1_whitelist)
    barcode2_whitelist = read_whitelist(args.barcode2_whitelist)
    if args.mode == "emseq":
        barcode1_whitelist = {
            convert_c_to_t(barcode) for barcode in barcode1_whitelist
        }
        barcode2_whitelist = {
            convert_c_to_t(barcode) for barcode in barcode2_whitelist
        }
    barcode1_length = len(next(iter(barcode1_whitelist)))
    barcode2_length = len(next(iter(barcode2_whitelist)))
    if any(len(barcode) != barcode1_length for barcode in barcode1_whitelist):
        raise ValueError("barcode1 whitelist contains inconsistent barcode lengths")
    if any(len(barcode) != barcode2_length for barcode in barcode2_whitelist):
        raise ValueError("barcode2 whitelist contains inconsistent barcode lengths")

    all_c_positions = (
        tuple(c_positions(insert_left)) if args.mode == "taps-v2" else ()
    )
    methylated_c_positions = (
        parse_methylated_c_positions(args.methylated_c_positions, insert_left)
        if args.mode == "taps-v2"
        else ()
    )
    return ExtractionConfig(
        mode=args.mode,
        linker_bc=linker_bc,
        insert_left=insert_left,
        linker_edit_distance=args.linker_edit_distance,
        barcode_hamming_distance=args.barcode_hamming_distance,
        insert_left_edit_distance=args.insert_left_edit_distance,
        barcode1_whitelist=frozenset(barcode1_whitelist),
        barcode2_whitelist=frozenset(barcode2_whitelist),
        barcode1_length=barcode1_length,
        barcode2_length=barcode2_length,
        minimum_length=(
            barcode2_length + len(linker_bc) + barcode1_length + len(insert_left)
        ),
        all_c_positions=all_c_positions,
        methylated_c_positions=methylated_c_positions,
    )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not 0 <= args.gzip_level <= 9:
        raise ValueError("--gzip-level must be between 0 and 9")
    if args.linker_edit_distance < 0:
        raise ValueError("--linker-edit-distance must be >= 0")
    if args.barcode_hamming_distance < 0:
        raise ValueError("--barcode-hamming-distance must be >= 0")
    if args.insert_left_edit_distance < 0:
        raise ValueError("--insert-left-edit-distance must be >= 0")
    if args.chunks < 1:
        raise ValueError("--chunks must be >= 1")
    if args.batch_size < 1:
        raise ValueError("--batch-size must be >= 1")
    if args.progress_reads < 0:
        raise ValueError("--progress-reads must be >= 0")

    config = build_extraction_config(args)
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
            name=f"barcode-chunk-{index + 1:04d}",
        )
        for index in range(args.chunks)
    ]
    started_at = time.monotonic()

    def put_batch(chunk_index: int, batch: object) -> None:
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

    for process in chunk_processes:
        process.start()
    batch: list[tuple[FastqRecord, FastqRecord]] = []
    batch_index = 0
    input_pairs = 0
    try:
        with open_text(args.r1, "r") as r1_input, open_text(args.r2, "r") as r2_input:
            for record1, record2 in zip_longest(
                fastq_iter(r1_input), fastq_iter(r2_input)
            ):
                if record1 is None or record2 is None:
                    raise ValueError("R1 and R2 FASTQ record counts are inconsistent")
                batch.append((record1, record2))
                input_pairs += 1
                if len(batch) == args.batch_size:
                    put_batch(batch_index % args.chunks, batch)
                    batch = []
                    batch_index += 1
            if batch:
                put_batch(batch_index % args.chunks, batch)
                batch_index += 1
        for index in range(args.chunks):
            put_batch(index, None)
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
    failed_processes = [
        process for process in chunk_processes if process.exitcode != 0
    ]
    failures = [message for message in messages if not message.get("ok")]
    if failed_processes or failures or len(messages) != args.chunks:
        details = "\n".join(str(message.get("error", "")) for message in failures)
        raise RuntimeError(
            f"barcode chunk processing failed; received {len(messages)}/"
            f"{args.chunks} results\n{details}"
        )

    messages.sort(key=lambda message: str(dict(message["row"])["chunk"]))
    chunk_rows = [dict(message["row"]) for message in messages]
    total = sum(int(row["total_reads"]) for row in chunk_rows)
    kept = sum(int(row["kept_reads"]) for row in chunk_rows)
    if total != input_pairs:
        raise RuntimeError(
            f"processed read count mismatch: input={input_pairs}, processed={total}"
        )
    reject_names = (
        "short_r1",
        "structure_mismatch",
        "linker2_not_found",
        "tn5_not_found",
        "barcode1_not_in_whitelist",
        "barcode2_not_in_whitelist",
    )
    reject_counts = {
        name: sum(int(dict(row["reject_counts"])[name]) for row in chunk_rows)
        for name in reject_names
    }
    stats: dict[str, object] = {
        "mode": args.mode,
        "input_r1": Path(args.r1).name,
        "input_r2": Path(args.r2).name,
        "output_dir": str(output_dir),
        "compression": args.compression,
        "chunk_count": args.chunks,
        "batch_size": args.batch_size,
        "batch_count": batch_index,
        "chunks": chunk_rows,
        "total_reads": total,
        "kept_reads": kept,
        "spike_in_reads": total - kept,
        "kept_fraction": kept / total if total else 0.0,
        "reject_counts": reject_counts,
    }
    if args.mode == "taps-v2":
        scored_inserts = sum(int(message["scored_inserts"]) for message in messages)
        scored_fully_converted = sum(
            int(message["scored_fully_converted"]) for message in messages
        )
        all_t_methylated_reads = sum(
            int(message["all_t_methylated_reads"]) for message in messages
        )
        methylated_rate_sum = sum(
            float(message["methylated_rate_sum"]) for message in messages
        )
        methylated_rate_count = sum(
            int(message["methylated_rate_count"]) for message in messages
        )
        site_converted = [
            sum(int(list(message["site_converted"])[index]) for message in messages)
            for index in range(len(config.all_c_positions))
        ]
        site_retained = [
            sum(int(list(message["site_retained"])[index]) for message in messages)
            for index in range(len(config.all_c_positions))
        ]
        site_other = [
            sum(int(list(message["site_other"])[index]) for message in messages)
            for index in range(len(config.all_c_positions))
        ]
        site_counts: dict[str, dict[str, int]] = {}
        rates_per_position: dict[str, float | None] = {}
        for index, position in enumerate(config.all_c_positions):
            converted = site_converted[index]
            retained = site_retained[index]
            other = site_other[index]
            site_counts[str(position)] = {
                "T": converted,
                "C": retained,
                "other": other,
            }
            informative = converted + retained
            rates_per_position[str(position)] = (
                converted / informative if informative else None
            )
        stats.update(
            {
                "insert_left_c_t_pattern": c_t_pattern(config.insert_left),
                "insert_left_edit_distance": args.insert_left_edit_distance,
                "scored_fully_converted": scored_fully_converted,
                "insert_left_scored_reads": scored_inserts,
                "insert_left_methylated_c_positions": list(
                    config.methylated_c_positions
                ),
                "insert_left_site_counts": site_counts,
                "insert_left_rate_per_c_pos": rates_per_position,
                "insert_left_mean_methylation_rate": (
                    methylated_rate_sum / methylated_rate_count
                    if methylated_rate_count
                    else None
                ),
                "insert_left_reads_fraction_all_T": (
                    all_t_methylated_reads / scored_inserts
                    if scored_inserts
                    else None
                ),
            }
        )
    elif args.mode == "emseq":
        stats.update(
            {
                "barcode_whitelist_c_to_t_conversion": True,
                "effective_linker_bc": config.linker_bc,
                "linker_c_t_pattern": c_t_pattern(config.linker_bc),
                "effective_insert_left": config.insert_left,
                "insert_left_c_t_pattern": c_t_pattern(config.insert_left),
            }
        )
    stats_path = output_dir / "stats.json"
    stats_path.write_text(json.dumps(stats, indent=2), encoding="utf-8")

    elapsed = max(time.monotonic() - started_at, 1e-9)
    print(f"[extract-bc] mode={args.mode}")
    print(f"[extract-bc] input-r1={args.r1}")
    print(f"[extract-bc] input-r2={args.r2}")
    print(f"[extract-bc] output-dir={output_dir}")
    print(f"[extract-bc] chunks={args.chunks}")
    print(f"[extract-bc] batch-size={args.batch_size}")
    print(f"[extract-bc] batch-count={batch_index}")
    print(f"[extract-bc] stats={stats_path}")
    print(
        f"[extract-bc] kept={kept}/{total} "
        f"keep-rate={kept / total if total else 0.0:.4f} "
        f"spike-in={total - kept} avg-speed={total / elapsed:.1f} reads/s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
