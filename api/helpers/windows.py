"""Interval coalescing for NCBI genomic windows."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GenomicWindow:
    accession: str
    start: int  # 1-based inclusive
    stop: int  # 1-based inclusive

    def __post_init__(self) -> None:
        if self.start < 1:
            raise ValueError("start must be >= 1")
        if self.stop < self.start:
            raise ValueError("stop must be >= start")

    @property
    def length(self) -> int:
        return self.stop - self.start + 1


@dataclass(frozen=True)
class CoalescedWindow:
    """Merged fetch window covering one or more group windows."""

    accession: str
    start: int
    stop: int
    members: tuple[GenomicWindow, ...]

    @property
    def length(self) -> int:
        return self.stop - self.start + 1


def coalesce_windows(
    windows: list[GenomicWindow],
    *,
    gap_bp: int = 10_000,
    max_merged_bp: int = 1_000_000,
) -> list[CoalescedWindow]:
    """
    Merge windows on the same accession that overlap or are within ``gap_bp``.

    Merges stop when the combined length would exceed ``max_merged_bp``.
    Result order follows sorted (accession, start, stop).
    """
    if not windows:
        return []
    if gap_bp < 0:
        raise ValueError("gap_bp must be >= 0")
    if max_merged_bp < 1:
        raise ValueError("max_merged_bp must be >= 1")

    by_acc: dict[str, list[GenomicWindow]] = {}
    for w in windows:
        by_acc.setdefault(w.accession, []).append(w)

    out: list[CoalescedWindow] = []
    for accession in sorted(by_acc):
        ordered = sorted(by_acc[accession], key=lambda w: (w.start, w.stop))
        cur_start = ordered[0].start
        cur_stop = ordered[0].stop
        members: list[GenomicWindow] = [ordered[0]]
        for w in ordered[1:]:
            can_bridge = w.start <= cur_stop + gap_bp + 1
            merged_stop = max(cur_stop, w.stop)
            merged_len = merged_stop - cur_start + 1
            if can_bridge and merged_len <= max_merged_bp:
                cur_stop = merged_stop
                members.append(w)
            else:
                out.append(
                    CoalescedWindow(
                        accession=accession,
                        start=cur_start,
                        stop=cur_stop,
                        members=tuple(members),
                    )
                )
                cur_start = w.start
                cur_stop = w.stop
                members = [w]
        out.append(
            CoalescedWindow(
                accession=accession,
                start=cur_start,
                stop=cur_stop,
                members=tuple(members),
            )
        )
    return out


def slice_sequence(
    sequence: str,
    *,
    fetched_start: int,
    slice_start: int,
    slice_stop: int,
) -> str:
    """
    Cut ``[slice_start, slice_stop]`` (1-based) out of a sequence that
    itself starts at ``fetched_start`` (1-based).
    """
    if slice_start < fetched_start or slice_stop < slice_start:
        raise ValueError("invalid slice bounds")
    offset = slice_start - fetched_start
    length = slice_stop - slice_start + 1
    end = offset + length
    if end > len(sequence):
        raise ValueError(
            f"sequence too short for slice: need {end}, have {len(sequence)}"
        )
    return sequence[offset:end]
