"""Read indexed BAM intervals including reference-assigned unknown positions."""

from __future__ import annotations

import struct
from pathlib import Path
from typing import BinaryIO, Iterator, Optional

import pysam


BAI_METADATA_BIN = 37450


def read_exact(handle: BinaryIO, size: int) -> bytes:
    data = handle.read(size)
    if len(data) != size:
        raise ValueError("truncated BAI index")
    return data


def read_uint32(handle: BinaryIO) -> int:
    return struct.unpack("<I", read_exact(handle, 4))[0]


def reference_offsets(bam_path: Path) -> dict[str, int]:
    """Read chromosome start virtual offsets from a local BAI index."""
    candidates = (
        Path(str(bam_path) + ".bai"),
        bam_path.with_suffix(".bai"),
    )
    index_path = next((path for path in candidates if path.is_file()), None)
    if index_path is None:
        raise FileNotFoundError(f"BAI index not found for: {bam_path}")
    with index_path.open("rb") as handle:
        if read_exact(handle, 4) != b"BAI\x01":
            raise ValueError(f"invalid BAI index format: {index_path}")
        offsets: list[Optional[int]] = []
        for _ in range(read_uint32(handle)):
            first_offset: Optional[int] = None
            for _ in range(read_uint32(handle)):
                bin_id = read_uint32(handle)
                chunk_count = read_uint32(handle)
                for chunk_index in range(chunk_count):
                    begin, _ = struct.unpack("<QQ", read_exact(handle, 16))
                    # The second pseudo-bin pair contains counts, not offsets.
                    if bin_id != BAI_METADATA_BIN or chunk_index == 0:
                        first_offset = (
                            begin if first_offset is None else min(first_offset, begin)
                        )
            read_exact(handle, read_uint32(handle) * 8)
            offsets.append(first_offset)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if len(offsets) > len(bam.references):
            raise ValueError(f"BAM index reference count exceeds BAM header: {index_path}")
        return {
            name: offset
            for name, offset in zip(bam.references, offsets)
            if offset is not None
        }


def interval_records(
    bam: pysam.AlignmentFile,
    reference: str,
    start: int,
    end: int,
    first_record_offset: Optional[int] = None,
) -> Iterator[pysam.AlignedSegment]:
    """Assign unknown positions to the first interval of a sorted reference."""
    if start == 0 and first_record_offset is not None:
        reference_id = bam.get_tid(reference)
        bam.seek(first_record_offset)
        for record in bam.fetch(until_eof=True):
            if record.reference_id != reference_id or record.reference_start >= end:
                break
            yield record
    else:
        for record in bam.fetch(reference, start, end):
            if start <= record.reference_start < end:
                yield record
