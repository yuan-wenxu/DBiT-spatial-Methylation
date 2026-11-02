#!/usr/bin/env python3
"""Find an R1 motif and measure C/T fractions at selected middle C sites."""

from __future__ import annotations

import argparse
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


Anchor = tuple[int, re.Pattern[str]]


@dataclass(frozen=True)
class MotifConfig:
    left: str
    middle: str
    right: str
    c_positions_1based: tuple[int, ...]
    c_positions_0based: tuple[int, ...]
    max_mismatches: int
    left_anchors: tuple[Anchor, ...]
    right_anchors: tuple[Anchor, ...]
    assay: str | None = None

    @property
    def length(self) -> int:
        return len(self.left) + len(self.middle) + len(self.right)


def parse_inputs(
    argv: list[str] | None = None,
) -> tuple[Path, MotifConfig, argparse.ArgumentParser]:
    """Collect the R1 path and motif settings from command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("r1", type=Path, help="R1 FASTQ or FASTQ.gz file")
    parser.add_argument(
        "--assay",
        choices=("taps", "taps-v2", "cabernet"),
        help=(
            "Report 5hmC as C for Cabernet or T for TAPS/TAPS-v2; "
            "Cabernet also converts unselected middle C to T. "
            "Omit to report only raw C/T statistics."
        ),
    )
    parser.add_argument(
        "--left-flank", required=True,
        help="Fixed left sequence; reference C accepts C/T for Cabernet only",
    )
    parser.add_argument("--middle", required=True, help="Middle sequence; N matches A/C/G/T")
    parser.add_argument(
        "--right-flank", required=True,
        help="Fixed right sequence; reference C accepts C/T for Cabernet only",
    )
    parser.add_argument(
        "--c-positions",
        "--c-indices",
        dest="c_positions",
        required=True,
        help="1-based positions of C bases within --middle, e.g. 2,6",
    )
    args = parser.parse_args(argv)
    return args.r1, build_config(args, parser), parser


def make_anchors(
    flank: str, motif_offset: int, max_mismatches: int,
    c_to_t: bool = False,
) -> tuple[Anchor, ...]:
    """Partition a flank into anchors, allowing T only at reference C sites."""
    anchor_count = max_mismatches + 1
    block_length = len(flank) // anchor_count
    if block_length == 0:
        return ()
    anchors = []
    for block in range(anchor_count):
        start = block * block_length
        expected = flank[start : start + block_length]
        pattern = "".join(
            "[CT]" if c_to_t and base == "C" else re.escape(base)
            for base in expected
        )
        # Lookahead keeps overlapping occurrences, as in the original search.
        anchors.append((motif_offset + start, re.compile(f"(?={pattern})")))
    return tuple(anchors)


def build_config(args: argparse.Namespace, parser: argparse.ArgumentParser) -> MotifConfig:
    left = args.left_flank.upper()
    middle = args.middle.upper()
    right = args.right_flank.upper()
    if not left or set(left) - set("ACGT"):
        parser.error("--left-flank must contain only A/C/G/T")
    if not right or set(right) - set("ACGT"):
        parser.error("--right-flank must contain only A/C/G/T")
    if not middle or set(middle) - set("ACGTN"):
        parser.error("--middle must contain only A/C/G/T/N")

    try:
        c_positions = tuple(int(value) for value in args.c_positions.split(","))
    except ValueError:
        parser.error("--c-positions must be comma-separated positions, such as 2,6")
    if (
        not c_positions
        or len(c_positions) != len(set(c_positions))
        or any(position < 1 or position > len(middle) for position in c_positions)
    ):
        parser.error(
            f"--c-positions must contain distinct positions from 1 to {len(middle)}"
        )
    if any(middle[position - 1] != "C" for position in c_positions):
        parser.error("each --c-positions value must point to C in --middle")

    if args.assay == "cabernet":
        selected_positions = {position - 1 for position in c_positions}
        middle = "".join(
            "T" if base == "C" and index not in selected_positions else base
            for index, base in enumerate(middle)
        )

    max_mismatches = (len(left) + len(right)) // 10
    return MotifConfig(
        left=left,
        middle=middle,
        right=right,
        c_positions_1based=c_positions,
        c_positions_0based=tuple(position - 1 for position in c_positions),
        max_mismatches=max_mismatches,
        left_anchors=make_anchors(
            left, 0, max_mismatches, c_to_t=args.assay == "cabernet"
        ),
        right_anchors=make_anchors(
            right, len(left) + len(middle), max_mismatches,
            c_to_t=args.assay == "cabernet",
        ),
        assay=args.assay,
    )


def candidate_starts(
    sequence: str, anchors: tuple[Anchor, ...], last_start: int
) -> set[int]:
    candidate_starts: set[int] = set()
    for offset, anchor in anchors:
        for match in anchor.finditer(sequence):
            start = match.start() - offset
            if 0 <= start <= last_start:
                candidate_starts.add(start)
    return candidate_starts


def flank_mismatch_count(
    observed: str, expected: str, c_to_t: bool
) -> int | None:
    """Use barcode extraction's reference-C rule for Cabernet flanks."""
    if len(observed) != len(expected):
        raise ValueError("observed and expected sequence length mismatch")
    mismatches = 0
    for observed_base, expected_base in zip(observed, expected):
        if c_to_t and expected_base == "C":
            if observed_base not in "CT":
                return None
        elif observed_base != expected_base:
            mismatches += 1
    return mismatches


def middle_matches(observed: str, config: MotifConfig) -> bool:
    selected_positions = set(config.c_positions_0based)
    for position, (base, expected) in enumerate(zip(observed, config.middle)):
        if position in selected_positions:
            if base not in "CT":
                return False
        elif expected == "N":
            if base not in "ACGT":
                return False
        elif base != expected:
            return False
    return True


def find_matches(
    sequence: str, config: MotifConfig
) -> Iterator[tuple[int, str, int, int]]:
    sequence = sequence.upper()
    last_start = len(sequence) - config.length
    if last_start < 0:
        return

    all_starts = set(range(last_start + 1))
    left_starts = (
        candidate_starts(sequence, config.left_anchors, last_start)
        if config.left_anchors
        else all_starts
    )
    right_starts = (
        candidate_starts(sequence, config.right_anchors, last_start)
        if config.right_anchors
        else all_starts
    )

    for start in sorted(left_starts & right_starts):
        left_read = sequence[start : start + len(config.left)]
        right_start = start + len(config.left) + len(config.middle)
        right_read = sequence[
            right_start : right_start + len(config.right)
        ]
        left_mismatches = flank_mismatch_count(
            left_read, config.left, c_to_t=config.assay == "cabernet"
        )
        if left_mismatches is None or left_mismatches > config.max_mismatches:
            continue
        right_mismatches = flank_mismatch_count(
            right_read, config.right, c_to_t=config.assay == "cabernet"
        )
        if (
            right_mismatches is None
            or left_mismatches + right_mismatches > config.max_mismatches
        ):
            continue

        middle_read = sequence[
            start + len(config.left) : start + len(config.left) + len(config.middle)
        ]
        if not middle_matches(middle_read, config):
            continue
        matched = sequence[start : start + config.length]
        yield start, matched, left_mismatches, right_mismatches


def fastq_reads(handle: TextIO) -> Iterator[tuple[str, str]]:
    while header := handle.readline():
        sequence = handle.readline()
        separator = handle.readline()
        quality = handle.readline()
        if not sequence or not separator or not quality:
            raise ValueError("Incomplete FASTQ record")
        if not header.startswith("@") or not separator.startswith("+"):
            raise ValueError(f"Invalid FASTQ record at {header.rstrip()!r}")
        sequence = sequence.rstrip("\r\n")
        if len(sequence) != len(quality.rstrip("\r\n")):
            raise ValueError(f"Sequence/quality lengths differ at {header.rstrip()!r}")
        yield header[1:].split()[0], sequence


def main() -> None:
    r1, config, parser = parse_inputs()

    opener = gzip.open if r1.suffix == ".gz" else open
    total_reads = 0
    matched_reads = 0
    matches = 0
    site_counts = {
        position: {"C": 0, "T": 0} for position in config.c_positions_1based
    }

    try:
        with opener(r1, "rt") as handle:
            for read_id, sequence in fastq_reads(handle):
                total_reads += 1
                found_in_read = False
                for start, matched, left_mismatches, right_mismatches in find_matches(
                    sequence, config
                ):
                    middle_read = matched[
                        len(config.left) : len(config.left) + len(config.middle)
                    ]
                    site_bases = tuple(
                        middle_read[position] for position in config.c_positions_0based
                    )

                    for position, base in zip(config.c_positions_1based, site_bases):
                        site_counts[position][base] += 1
                    matches += 1
                    found_in_read = True
                matched_reads += found_in_read
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    hmc_base = "C" if config.assay == "cabernet" else "T"
    assay_summary = (
        f"assay={config.assay} hmc_base={hmc_base}"
        if config.assay is not None else ""
    )
    print(
        f"[find-r1-5hmc-motif] reads={total_reads} "
        f"matched_reads={matched_reads} matches={matches} {assay_summary}"
    )
    for position in config.c_positions_1based:
        c_count = site_counts[position]["C"]
        t_count = site_counts[position]["T"]
        total = c_count + t_count
        hmc_summary = ""
        if config.assay is not None:
            hmc_count = site_counts[position][hmc_base]
            hmc_fraction = hmc_count / total if total else 0.0
            hmc_summary = (
                f" hmc_count={hmc_count} hmc_fraction={hmc_fraction:.4f}"
            )
        print(
            f"[find-r1-5hmc-motif] middle_position_1based={position} "
            f"C={c_count} T={t_count} "
            f"C_fraction={c_count / total if total else 0.0:.4f} "
            f"T_fraction={t_count / total if total else 0.0:.4f}"
            f"{hmc_summary}"
        )


if __name__ == "__main__":
    main()
