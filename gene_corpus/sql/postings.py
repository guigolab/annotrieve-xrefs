"""
Posting-list codec for the served presence index.

A posting list is the strictly increasing set of annotation ordinals (``ord``,
equal to ``annotation.annotation_key``) in which one gene name occurs.

Blob format (first byte is a tag):

* ``0x01`` varint:  LEB128 deltas of the sorted ordinals. The first value is
  stored as a delta from 0. Used when it is smaller than the bitset.
* ``0x02`` bitset:  ``ceil((max_ord + 1) / 8)`` bytes, little-endian bit order
  (bit ``o`` of the integer is set when ordinal ``o`` is present), so
  ``int.from_bytes(blob[1:], "little")`` is a Python-int bitset. Used for
  dense names, where it is smaller than the varint list.

The API has its own decoder (it must not import ``gene_corpus``); keep the two
in step through ``tests/test_serve_presence.py`` fixtures and the format above.
"""
from __future__ import annotations

from collections.abc import Iterator, Sequence

TAG_VARINT = 0x01
TAG_BITSET = 0x02


def bitset_nbytes(max_ord: int) -> int:
    """Bytes needed for a bitset holding ordinals ``0..max_ord``."""
    if max_ord < 0:
        raise ValueError("max_ord must be >= 0")
    return (max_ord + 8) // 8


def _varint_encode(ords: Sequence[int]) -> bytearray:
    out = bytearray()
    prev = 0
    first = True
    for o in ords:
        delta = o - prev
        if delta < 0 or (not first and delta == 0):
            raise ValueError("ordinals must be strictly increasing")
        first = False
        prev = o
        while delta >= 0x80:
            out.append((delta & 0x7F) | 0x80)
            delta >>= 7
        out.append(delta)
    return out


def _bitset_encode(ords: Sequence[int], nbytes: int) -> bytes:
    buf = bytearray(nbytes)
    for o in ords:
        buf[o >> 3] |= 1 << (o & 7)
    return bytes([TAG_BITSET]) + bytes(buf)


def encode_postings(ords: Sequence[int], max_ord: int) -> bytes:
    """Encode a strictly increasing, non-empty ordinal list as a blob."""
    n = len(ords)
    if n == 0:
        raise ValueError("empty posting list")
    if ords[-1] > max_ord or ords[0] < 0:
        raise ValueError("ordinal outside 0..max_ord")
    nbytes = bitset_nbytes(max_ord)
    # A varint list needs at least one byte per entry, so a bitset always wins here.
    if n >= nbytes:
        _check_increasing(ords)
        return _bitset_encode(ords, nbytes)
    varint = _varint_encode(ords)
    if len(varint) < nbytes:
        return bytes([TAG_VARINT]) + bytes(varint)
    return _bitset_encode(ords, nbytes)


def _check_increasing(ords: Sequence[int]) -> None:
    prev = -1
    for o in ords:
        if o <= prev:
            raise ValueError("ordinals must be strictly increasing")
        prev = o


def decode_postings(blob: bytes) -> list[int]:
    """Return the sorted ordinal list of a blob."""
    return list(iter_postings(blob))


def iter_postings(blob: bytes) -> Iterator[int]:
    if not blob:
        return
    tag = blob[0]
    if tag == TAG_VARINT:
        value = 0
        shift = 0
        prev = 0
        for byte in blob[1:]:
            value |= (byte & 0x7F) << shift
            if byte & 0x80:
                shift += 7
            else:
                prev += value
                yield prev
                value = 0
                shift = 0
    elif tag == TAG_BITSET:
        for index, byte in enumerate(blob[1:]):
            if byte:
                base = index * 8
                for bit in range(8):
                    if byte >> bit & 1:
                        yield base + bit
    else:
        raise ValueError(f"unknown posting tag {tag:#x}")


def decode_bitset(blob: bytes) -> int:
    """Return the postings as a Python-int bitset (bit ``o`` set = ordinal ``o``)."""
    if not blob:
        return 0
    if blob[0] == TAG_BITSET:
        return int.from_bytes(blob[1:], "little")
    value = 0
    for o in iter_postings(blob):
        value |= 1 << o
    return value
