"""
Posting-list decoder for presence.sqlite (the API must not import ``gene_corpus``).

Mirror of ``gene_corpus/sql/postings.py``; kept in step by ``tests/test_presence_api.py``.
Tag 0x01: LEB128 deltas of ascending ordinals. Tag 0x02: little-endian bitset.
"""
from __future__ import annotations

from collections.abc import Iterator

TAG_VARINT = 0x01
TAG_BITSET = 0x02


def iter_bits(mask: int) -> Iterator[int]:
    """Yield the set bit positions of a non-negative int, ascending."""
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def to_bitset(blob: bytes | None) -> int:
    """Return the blob as a Python-int bitset (bit ``o`` set = ordinal ``o``)."""
    if not blob:
        return 0
    tag = blob[0]
    if tag == TAG_BITSET:
        return int.from_bytes(blob[1:], "little")
    if tag != TAG_VARINT:
        raise ValueError(f"unknown posting tag {tag:#x}")
    mask = 0
    value = shift = prev = 0
    for byte in blob[1:]:
        value |= (byte & 0x7F) << shift
        if byte & 0x80:
            shift += 7
        else:
            prev += value
            mask |= 1 << prev
            value = shift = 0
    return mask
