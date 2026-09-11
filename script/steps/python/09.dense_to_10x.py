#!/usr/bin/env python3
"""Convert a MethSCAn barcode-by-VMR CSV matrix to 10x-style sparse files."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import shutil
import tempfile
from contextlib import contextmanager
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterator, TextIO


MISSING_VALUES = {"", "na", "nan", "null", "none"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert a MethSCAn dense CSV/CSV.GZ matrix (barcodes by VMRs) "
            "to a 10x-style sparse matrix (VMRs by barcodes)."
        )
    )
    parser.add_argument("input", type=Path, help="Input .csv or .csv.gz matrix")
    parser.add_argument("output", type=Path, help="Output 10x-style directory")
    parser.add_argument(
        "--feature-type",
        default="Methylation",
        help="Third column in features.tsv.gz (default: Methylation)",
    )
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("r", newline="")


@contextmanager
def open_reproducible_gzip(path: Path) -> Iterator[TextIO]:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            fileobj=raw_handle, mode="wb", filename="", mtime=0
        ) as gzip_handle:
            with io.TextIOWrapper(
                gzip_handle, encoding="utf-8", newline=""
            ) as text_handle:
                yield text_handle


def temporary_path(output_dir: Path, suffix: str) -> Path:
    descriptor, name = tempfile.mkstemp(
        prefix=".dense_to_10x.", suffix=suffix, dir=output_dir
    )
    os.close(descriptor)
    return Path(name)


def validate_unique(values: list[str], label: str) -> None:
    seen: set[str] = set()
    for value in values:
        if not value:
            raise ValueError(f"empty {label}")
        if value in seen:
            raise ValueError(f"duplicate {label}: {value}")
        seen.add(value)


def nonzero_numeric(value: str, row_number: int, column_number: int) -> bool:
    stripped = value.strip()
    if stripped.lower() in MISSING_VALUES:
        return False
    try:
        number = Decimal(stripped)
    except InvalidOperation as error:
        raise ValueError(
            f"invalid numeric value at CSV row {row_number}, column {column_number}: "
            f"{value!r}"
        ) from error
    if not number.is_finite():
        raise ValueError(
            f"non-finite numeric value at CSV row {row_number}, column {column_number}: "
            f"{value!r}"
        )
    return not number.is_zero()


def convert(
    input_path: Path, output_dir: Path, feature_type: str
) -> tuple[int, int, int]:
    if not input_path.is_file():
        raise FileNotFoundError(f"input matrix not found: {input_path}")
    if not feature_type or "\t" in feature_type or "\n" in feature_type:
        raise ValueError(
            "feature type must be non-empty and cannot contain tabs or newlines"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    entries_path = temporary_path(output_dir, ".entries")
    matrix_tmp = temporary_path(output_dir, ".matrix.mtx.gz")
    barcodes_tmp = temporary_path(output_dir, ".barcodes.tsv.gz")
    features_tmp = temporary_path(output_dir, ".features.tsv.gz")
    temporary_files = (entries_path, matrix_tmp, barcodes_tmp, features_tmp)

    try:
        with open_text(input_path) as input_handle:
            reader = csv.reader(input_handle)
            try:
                header = next(reader)
            except StopIteration as error:
                raise ValueError("input matrix is empty") from error

            if len(header) < 2:
                raise ValueError(
                    "input header must contain a row-label column and at least one VMR"
                )
            features = [feature.strip() for feature in header[1:]]
            validate_unique(features, "feature")

            barcode_count = 0
            nonzero_count = 0
            seen_barcodes: set[str] = set()
            with entries_path.open(
                "w", encoding="utf-8", newline=""
            ) as entries_handle, open_reproducible_gzip(
                barcodes_tmp
            ) as barcodes_handle:
                for csv_row_number, row in enumerate(reader, start=2):
                    if len(row) != len(header):
                        raise ValueError(
                            f"CSV row {csv_row_number} has {len(row)} columns; "
                            f"expected {len(header)}"
                        )
                    barcode = row[0].strip()
                    if not barcode:
                        raise ValueError(f"empty barcode at CSV row {csv_row_number}")
                    if barcode in seen_barcodes:
                        raise ValueError(
                            f"duplicate barcode at CSV row {csv_row_number}: {barcode}"
                        )
                    if "\t" in barcode or "\n" in barcode:
                        raise ValueError(
                            f"barcode contains a tab or newline: {barcode!r}"
                        )
                    seen_barcodes.add(barcode)
                    barcode_count += 1
                    barcodes_handle.write(f"{barcode}\n")

                    for feature_index, value in enumerate(row[1:], start=1):
                        if nonzero_numeric(
                            value, csv_row_number, feature_index + 1
                        ):
                            entries_handle.write(
                                f"{feature_index} {barcode_count} {value.strip()}\n"
                            )
                            nonzero_count += 1

        if barcode_count == 0:
            raise ValueError("input matrix contains no barcode rows")

        with open_reproducible_gzip(features_tmp) as features_handle:
            for feature in features:
                if "\t" in feature or "\n" in feature:
                    raise ValueError(
                        f"feature contains a tab or newline: {feature!r}"
                    )
                features_handle.write(
                    f"{feature}\t{feature}\t{feature_type}\n"
                )

        with open_reproducible_gzip(matrix_tmp) as matrix_handle:
            matrix_handle.write("%%MatrixMarket matrix coordinate real general\n")
            matrix_handle.write("% generated by DBiT-spatial-Methylation\n")
            matrix_handle.write(
                f"{len(features)} {barcode_count} {nonzero_count}\n"
            )
            with entries_path.open("r", encoding="utf-8") as entries_handle:
                shutil.copyfileobj(entries_handle, matrix_handle)

        os.replace(matrix_tmp, output_dir / "matrix.mtx.gz")
        os.replace(barcodes_tmp, output_dir / "barcodes.tsv.gz")
        os.replace(features_tmp, output_dir / "features.tsv.gz")
        return len(features), barcode_count, nonzero_count
    finally:
        for path in temporary_files:
            path.unlink(missing_ok=True)


def main() -> None:
    args = parse_args()
    try:
        feature_count, barcode_count, nonzero_count = convert(
            args.input, args.output, args.feature_type
        )
    except (OSError, ValueError, csv.Error) as error:
        raise SystemExit(f"dense-to-10x: error: {error}") from error

    density = nonzero_count / (feature_count * barcode_count)
    print(f"[dbitm] input: {args.input}")
    print(f"[dbitm] output: {args.output}")
    print(
        f"[dbitm] shape: {feature_count} features x {barcode_count} barcodes; "
        f"nonzero: {nonzero_count}; density: {density:.6%}"
    )


if __name__ == "__main__":
    main()
