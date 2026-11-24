"""Tests for the postings codec (gene_corpus.sql.postings)."""
from __future__ import annotations

import pytest

from gene_corpus.sql.postings import (
    TAG_BITSET,
    TAG_VARINT,
    bitset_nbytes,
    decode_bitset,
    decode_postings,
    encode_postings,
)


# --------------------------------------------------------------------------
# codec
# --------------------------------------------------------------------------
def test_roundtrip_small_uses_varint():
    ords = [1, 5, 300, 20000]
    blob = encode_postings(ords, 20000)
    assert blob[0] == TAG_VARINT
    assert decode_postings(blob) == ords
    assert decode_bitset(blob) == sum(1 << o for o in ords)


def test_roundtrip_dense_uses_bitset():
    ords = list(range(1, 3000, 2))
    blob = encode_postings(ords, 3000)
    assert blob[0] == TAG_BITSET
    assert len(blob) == 1 + bitset_nbytes(3000)
    assert decode_postings(blob) == ords
    assert decode_bitset(blob) == sum(1 << o for o in ords)


def test_chooses_smaller_encoding():
    # 5 entries spread widely: varint is 3 bytes each = 15 > nbytes(1000)=126? no -> varint
    blob = encode_postings([10, 400, 700, 900, 1000], 1000)
    assert blob[0] == TAG_VARINT
    assert len(blob) <= 1 + bitset_nbytes(1000)


def test_single_and_first_ordinal_zero():
    assert decode_postings(encode_postings([0], 10)) == [0]
    assert decode_postings(encode_postings([7], 7)) == [7]


def test_rejects_bad_input():
    with pytest.raises(ValueError):
        encode_postings([], 10)
    with pytest.raises(ValueError):
        encode_postings([3, 3], 10)
    with pytest.raises(ValueError):
        encode_postings([5, 2], 10)
    with pytest.raises(ValueError):
        encode_postings([11], 10)
    with pytest.raises(ValueError):
        decode_postings(bytes([0x7F, 1]))


def test_publish_version_switches_the_current_link(tmp_path):
    import json

    from gene_corpus.sql.serve_presence import publish_version

    out = tmp_path / "serve"
    out.mkdir()
    for version in ("aaa", "bbb"):
        built = out / ".building.sqlite"
        built.write_bytes(b"x")
        publish_version(out, built, {"data_version": version})
        assert (out / "current").resolve().name == f"v{version}"
        assert json.loads((out / f"v{version}" / "build_report.json").read_text())["data_version"] == version
    assert (out / "vaaa" / "presence.sqlite").is_file()          # older versions stay in place
