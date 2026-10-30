"""Unit tests for genomic window coalescing."""
from __future__ import annotations

import pytest

from helpers.windows import GenomicWindow, coalesce_windows, slice_sequence


def test_coalesce_overlapping_same_accession():
    windows = [
        GenomicWindow("NC_1.1", 100, 200),
        GenomicWindow("NC_1.1", 180, 300),
    ]
    merged = coalesce_windows(windows, gap_bp=0, max_merged_bp=1_000_000)
    assert len(merged) == 1
    assert merged[0].start == 100
    assert merged[0].stop == 300
    assert len(merged[0].members) == 2


def test_coalesce_within_gap():
    windows = [
        GenomicWindow("NC_1.1", 1, 100),
        GenomicWindow("NC_1.1", 150, 200),
    ]
    merged = coalesce_windows(windows, gap_bp=50, max_merged_bp=1_000_000)
    assert len(merged) == 1
    assert merged[0].start == 1
    assert merged[0].stop == 200


def test_coalesce_beyond_gap_stays_separate():
    windows = [
        GenomicWindow("NC_1.1", 1, 100),
        GenomicWindow("NC_1.1", 200, 300),
    ]
    merged = coalesce_windows(windows, gap_bp=50, max_merged_bp=1_000_000)
    assert len(merged) == 2


def test_coalesce_respects_max_merged():
    windows = [
        GenomicWindow("NC_1.1", 1, 600),
        GenomicWindow("NC_1.1", 601, 1200),
    ]
    merged = coalesce_windows(windows, gap_bp=10, max_merged_bp=1000)
    assert len(merged) == 2


def test_coalesce_groups_by_accession():
    windows = [
        GenomicWindow("NC_1.1", 1, 50),
        GenomicWindow("NC_2.1", 1, 50),
        GenomicWindow("NC_1.1", 40, 80),
    ]
    merged = coalesce_windows(windows, gap_bp=0, max_merged_bp=1_000_000)
    assert len(merged) == 2
    by_acc = {m.accession: m for m in merged}
    assert by_acc["NC_1.1"].stop == 80
    assert by_acc["NC_2.1"].stop == 50


def test_slice_sequence():
    seq = "ACGTACGTAA"
    assert slice_sequence(seq, fetched_start=100, slice_start=100, slice_stop=104) == "ACGTA"
    assert slice_sequence(seq, fetched_start=100, slice_start=103, slice_stop=106) == "TACG"


def test_slice_sequence_bounds():
    with pytest.raises(ValueError):
        slice_sequence("ACGT", fetched_start=10, slice_start=9, slice_stop=12)
