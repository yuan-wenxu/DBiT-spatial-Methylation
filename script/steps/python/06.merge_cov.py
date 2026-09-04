#!/usr/bin/env python3
"""Merge independently called Watson/Crick coverage files."""

from __future__ import annotations

import argparse
import heapq
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


@dataclass(frozen=True)
class CoverageRow:
    key: tuple[object, ...]
    chrom: str
    start: int
    end: int
    methylated: int
    unmethylated: int
    trailing: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge coordinate-compatible coverage calls by summing methylated "
            "and unmethylated counts."
        )
    )
    parser.add_argument("--input", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--chromosomes", required=True)
    parser.add_argument(
        "--spatial",
        action="store_true",
        help="Input is host coverage ordered by spot and then position.",
    )
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--batch-size", type=int)
    parser.add_argument("--jobs", type=int)
    return parser.parse_args()


def chromosome_order(raw: str) -> dict[str, int]:
    names = [name.strip() for name in raw.split(",") if name.strip()]
    if not names or len(names) != len(set(names)):
        raise ValueError("--chromosomes must contain unique chromosome names")
    return {name: index for index, name in enumerate(names)}


def spatial_batch_sizes(
    reference: Path | None,
    configured_batch_size: int | None,
    jobs: int | None,
    order: dict[str, int],
) -> dict[str, int]:
    if reference is None or configured_batch_size is None or jobs is None:
        raise ValueError("--spatial requires --reference, --batch-size, and --jobs")
    if configured_batch_size <= 0 or jobs <= 0:
        raise ValueError("--batch-size and --jobs must be greater than zero")
    fai = Path(str(reference) + ".fai")
    if not fai.is_file():
        raise FileNotFoundError(fai)
    lengths: dict[str, int] = {}
    with fai.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[0] in order:
                lengths[fields[0]] = int(fields[1])
    missing = [name for name in order if name not in lengths]
    if missing:
        raise ValueError(
            "selected chromosomes absent from FASTA index: " + ",".join(missing)
        )
    return {
        name: min(configured_batch_size, max(1, length // jobs))
        for name, length in lengths.items()
    }


def rows(
    handle: TextIO,
    path: Path,
    order: dict[str, int],
    spatial: bool,
    batch_sizes: dict[str, int],
) -> Iterator[CoverageRow]:
    previous_key: tuple[object, ...] | None = None
    for line_number, line in enumerate(handle, start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            raise ValueError(f"{path}:{line_number}: expected at least 6 columns")
        chrom = fields[0]
        try:
            chrom_index = order[chrom]
            start = int(fields[1])
            end = int(fields[2])
            methylated = int(fields[4])
            unmethylated = int(fields[5])
        except KeyError as error:
            raise ValueError(
                f"{path}:{line_number}: chromosome {chrom!r} is not selected"
            ) from error
        except ValueError as error:
            raise ValueError(
                f"{path}:{line_number}: invalid coordinate or count"
            ) from error
        if methylated < 0 or unmethylated < 0:
            raise ValueError(f"{path}:{line_number}: counts must be non-negative")
        trailing = tuple(fields[6:])
        if spatial:
            if not trailing:
                raise ValueError(
                    f"{path}:{line_number}: spatial coverage lacks a spot column"
                )
            batch_index = start // batch_sizes[chrom]
            key = (
                chrom_index,
                batch_index,
                trailing[0],
                start,
                end,
                *trailing[1:],
            )
        else:
            key = (chrom_index, start, end, *trailing)
        if previous_key is not None and key < previous_key:
            raise ValueError(f"{path}:{line_number}: input is not in caller order")
        previous_key = key
        yield CoverageRow(
            key, chrom, start, end, methylated, unmethylated, trailing
        )


def write_row(handle: TextIO, row: CoverageRow, methylated: int, unmethylated: int) -> None:
    total = methylated + unmethylated
    percentage = round(100.0 * methylated / total, 2) if total else 0.0
    fields = [
        row.chrom,
        str(row.start),
        str(row.end),
        f"{percentage:.2f}",
        str(methylated),
        str(unmethylated),
        *row.trailing,
    ]
    handle.write("\t".join(fields) + "\n")


def merge_files(
    inputs: list[Path],
    output: Path,
    order: dict[str, int],
    spatial: bool,
    batch_sizes: dict[str, int],
) -> tuple[int, int]:
    handles: list[TextIO] = []
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    input_rows = 0
    output_rows = 0
    try:
        for path in inputs:
            if not path.is_file():
                raise FileNotFoundError(path)
            handles.append(path.open(encoding="utf-8"))
        iterators = [
            rows(handle, path, order, spatial, batch_sizes)
            for handle, path in zip(handles, inputs)
        ]
        merged = heapq.merge(*iterators, key=lambda row: row.key)
        output.parent.mkdir(parents=True, exist_ok=True)
        with temporary.open("w", encoding="utf-8", newline="") as destination:
            current: CoverageRow | None = None
            methylated = 0
            unmethylated = 0
            for row in merged:
                input_rows += 1
                if current is not None and row.key != current.key:
                    write_row(destination, current, methylated, unmethylated)
                    output_rows += 1
                    methylated = 0
                    unmethylated = 0
                current = row
                methylated += row.methylated
                unmethylated += row.unmethylated
            if current is not None:
                write_row(destination, current, methylated, unmethylated)
                output_rows += 1
        os.replace(temporary, output)
    finally:
        for handle in handles:
            handle.close()
        temporary.unlink(missing_ok=True)
    return input_rows, output_rows


def main() -> int:
    args = parse_args()
    try:
        order = chromosome_order(args.chromosomes)
        batch_sizes = (
            spatial_batch_sizes(args.reference, args.batch_size, args.jobs, order)
            if args.spatial
            else {}
        )
        input_rows, output_rows = merge_files(
            args.input, args.output, order, args.spatial, batch_sizes
        )
    except (FileNotFoundError, OSError, ValueError) as error:
        print(f"[merge-cov] error: {error}", file=sys.stderr)
        return 1
    print(
        f"[merge-cov] inputs={len(args.input)} input-rows={input_rows} "
        f"output-rows={output_rows} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
