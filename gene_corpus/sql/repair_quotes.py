"""Repair quote-wrapped symbol/alias accessions in gene_corpus.sqlite."""
from __future__ import annotations

import sqlite3
import sys
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from gene_corpus.sql.attach import rebuild_namespace_stats
from gene_corpus.sql.schema import checkpoint_wal
from gene_corpus.stream.symbol_gate import strip_wrapping_quotes

QUOTED_STATE_NAME = "quoted.txt"
AFFECTED_BARE_STATE_NAME = "affected_bare.txt"
DEFAULT_NAMESPACES: tuple[str, ...] = ("symbol", "alias")
_SAMPLE_LIMIT = 15
_COMMIT_EVERY = 200
_CHECKPOINT_EVERY = 1000


@dataclass(frozen=True)
class RepairResult:
    """Summary of a repair (or dry-run) pass."""

    quoted_keys: int
    affected_bare: int
    moved_rows: int
    dup_loci_skipped: int
    keys_repaired: int
    wrapped_remaining: int
    applied: bool


def is_wrapped_accession(accession: str) -> bool:
    """True when accession has one matching wrapping ASCII quote pair."""
    text = accession or ""
    return len(text) >= 3 and text[0] == text[-1] and text[0] in "'\""


def iter_wrapped_accessions(
    conn: sqlite3.Connection,
    namespaces: Sequence[str],
) -> Iterable[tuple[str, str]]:
    """
    Yield ``(namespace, accession)`` for quote-wrapped keys in ``xref_meta``.

    Streams via a cursor (no full-namespace ``fetchall``).
    """
    for ns in namespaces:
        cur = conn.execute(
            """
            SELECT namespace, accession
            FROM xref_meta
            WHERE namespace = ?
              AND length(accession) >= 3
              AND substr(accession, 1, 1) = substr(accession, -1, 1)
              AND substr(accession, 1, 1) IN ('''', '"')
            """,
            (ns,),
        )
        yield from cur


def count_wrapped_accessions(
    conn: sqlite3.Connection,
    namespaces: Sequence[str],
) -> int:
    n = 0
    for _ in iter_wrapped_accessions(conn, namespaces):
        n += 1
    return n


def _write_state_files(
    state_dir: Path,
    quoted: Sequence[tuple[str, str]],
    affected_bare: Sequence[tuple[str, str]],
) -> None:
    state_dir.mkdir(parents=True, exist_ok=True)
    quoted_path = state_dir / QUOTED_STATE_NAME
    bare_path = state_dir / AFFECTED_BARE_STATE_NAME
    with quoted_path.open("w", encoding="utf-8") as fh:
        for ns, acc in quoted:
            fh.write(f"{ns}\t{acc}\n")
    with bare_path.open("w", encoding="utf-8") as fh:
        for ns, acc in affected_bare:
            fh.write(f"{ns}\t{acc}\n")


def _read_ns_acc_file(path: Path) -> list[tuple[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"repair state file missing: {path}")
    rows: list[tuple[str, str]] = []
    with path.open(encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            text = line.rstrip("\n")
            if not text.strip():
                continue
            parts = text.split("\t", 1)
            if len(parts) != 2 or not parts[0] or not parts[1]:
                raise ValueError(
                    f"invalid repair state line {line_no} in {path}: {text!r}"
                )
            rows.append((parts[0], parts[1]))
    return rows


def _move_quoted_key(
    conn: sqlite3.Connection,
    namespace: str,
    quoted_acc: str,
    bare_acc: str,
) -> tuple[int, int]:
    """
    Copy gene_hit rows to bare accession; delete quoted gene_hit + xref_meta.

    Returns ``(moved_rows, dup_loci_skipped)``.
    """
    moved = 0
    dups = 0
    rows = conn.execute(
        """
        SELECT annotation_key, local_id, seqid, start, end, strand,
               feature_type, biotype, primary_name
        FROM gene_hit
        WHERE namespace = ? AND accession = ?
        """,
        (namespace, quoted_acc),
    )
    for ak, lid, seqid, start, end, strand, ft, bt, pname in rows:
        try:
            conn.execute(
                """
                INSERT INTO gene_hit (
                    namespace, accession, annotation_key, local_id,
                    seqid, start, end, strand,
                    feature_type, biotype, primary_name
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    namespace,
                    bare_acc,
                    ak,
                    lid,
                    seqid,
                    start,
                    end,
                    strand,
                    ft,
                    bt,
                    pname,
                ),
            )
            moved += 1
        except sqlite3.IntegrityError:
            dups += 1

    conn.execute(
        "DELETE FROM gene_hit WHERE namespace = ? AND accession = ?",
        (namespace, quoted_acc),
    )
    conn.execute(
        "DELETE FROM xref_meta WHERE namespace = ? AND accession = ?",
        (namespace, quoted_acc),
    )
    return moved, dups


def _recompute_xref_meta_for_keys(
    conn: sqlite3.Connection,
    keys: Sequence[tuple[str, str]],
    *,
    progress_every: int = 500,
) -> None:
    for i, (namespace, accession) in enumerate(keys, start=1):
        n_ann, n_loci = conn.execute(
            """
            SELECT COUNT(DISTINCT annotation_key), COUNT(*)
            FROM gene_hit
            WHERE namespace = ? AND accession = ?
            """,
            (namespace, accession),
        ).fetchone()
        if int(n_loci) == 0:
            conn.execute(
                "DELETE FROM xref_meta WHERE namespace = ? AND accession = ?",
                (namespace, accession),
            )
        else:
            conn.execute(
                """
                INSERT INTO xref_meta (
                    namespace, accession, n_annotations, n_loci
                ) VALUES (?, ?, ?, ?)
                ON CONFLICT(namespace, accession) DO UPDATE SET
                    n_annotations = excluded.n_annotations,
                    n_loci = excluded.n_loci
                """,
                (namespace, accession, int(n_ann), int(n_loci)),
            )
        if progress_every > 0 and i % progress_every == 0:
            conn.commit()
            print(
                f"repair meta: {i}/{len(keys)}",
                file=sys.stderr,
                flush=True,
            )
    conn.commit()


def recompute_xref_meta_from_state(
    conn: sqlite3.Connection,
    state_dir: Path,
) -> RepairResult:
    """
    Recompute ``xref_meta`` + ``namespace_stats`` from ``affected_bare.txt``.

    Recovery path when ``--apply`` died after gene_hit moves.
    """
    bare_path = state_dir / AFFECTED_BARE_STATE_NAME
    affected = _read_ns_acc_file(bare_path)
    print(
        f"repair recompute: {len(affected)} bare accession(s) from {bare_path}",
        file=sys.stderr,
        flush=True,
    )
    _recompute_xref_meta_for_keys(conn, affected)
    rebuild_namespace_stats(conn)
    checkpoint_wal(conn)

    # Namespaces present in the state file (for remaining-wrap check).
    namespaces = tuple(dict.fromkeys(ns for ns, _ in affected)) or DEFAULT_NAMESPACES
    remaining = count_wrapped_accessions(conn, namespaces)
    return RepairResult(
        quoted_keys=0,
        affected_bare=len(affected),
        moved_rows=0,
        dup_loci_skipped=0,
        keys_repaired=0,
        wrapped_remaining=remaining,
        applied=True,
    )


def repair_wrapped_accessions(
    conn: sqlite3.Connection,
    *,
    namespaces: Sequence[str] = DEFAULT_NAMESPACES,
    state_dir: Path,
    apply: bool = False,
    progress_every: int = _COMMIT_EVERY,
) -> RepairResult:
    """
    Unwrap quote-wrapped symbol/alias accessions in the published DB.

    Streams matches into ``state_dir`` (no full-namespace load). Default is
    dry-run (``apply=False``): write state files and report, leave DB unchanged.
    """
    ns_tuple = tuple(namespaces)
    quoted: list[tuple[str, str]] = []
    affected_set: set[tuple[str, str]] = set()

    for ns, acc in iter_wrapped_accessions(conn, ns_tuple):
        if not is_wrapped_accession(acc):
            continue
        bare = strip_wrapping_quotes(acc)
        if not bare:
            print(
                f"repair warn: skip empty bare for {ns}:{acc!r}",
                file=sys.stderr,
            )
            continue
        quoted.append((ns, acc))
        affected_set.add((ns, bare))

    affected = sorted(affected_set)
    _write_state_files(state_dir, quoted, affected)

    print(
        f"repair: found {len(quoted)} wrapped accession(s) "
        f"→ {len(affected)} bare key(s); state → {state_dir}",
        file=sys.stderr,
        flush=True,
    )
    for ns, acc in quoted[:_SAMPLE_LIMIT]:
        print(
            f"  sample {ns}:{acc!r} → {strip_wrapping_quotes(acc)!r}",
            file=sys.stderr,
        )
    if len(quoted) > _SAMPLE_LIMIT:
        print(
            f"  … ({len(quoted) - _SAMPLE_LIMIT} more)",
            file=sys.stderr,
        )

    if not apply:
        remaining = count_wrapped_accessions(conn, ns_tuple)
        return RepairResult(
            quoted_keys=len(quoted),
            affected_bare=len(affected),
            moved_rows=0,
            dup_loci_skipped=0,
            keys_repaired=0,
            wrapped_remaining=remaining,
            applied=False,
        )

    moved_total = 0
    dups_total = 0
    for i, (ns, acc) in enumerate(quoted, start=1):
        bare = strip_wrapping_quotes(acc)
        moved, dups = _move_quoted_key(conn, ns, acc, bare)
        moved_total += moved
        dups_total += dups
        if progress_every > 0 and i % progress_every == 0:
            conn.commit()
            print(
                f"repair move: {i}/{len(quoted)} "
                f"(moved_rows={moved_total} dups={dups_total})",
                file=sys.stderr,
                flush=True,
            )
            if i % _CHECKPOINT_EVERY == 0:
                checkpoint_wal(conn)

    conn.commit()

    print(
        f"repair: recomputing xref_meta for {len(affected)} bare key(s)…",
        file=sys.stderr,
        flush=True,
    )
    _recompute_xref_meta_for_keys(conn, affected)
    rebuild_namespace_stats(conn)
    checkpoint_wal(conn)

    remaining = count_wrapped_accessions(conn, ns_tuple)
    result = RepairResult(
        quoted_keys=len(quoted),
        affected_bare=len(affected),
        moved_rows=moved_total,
        dup_loci_skipped=dups_total,
        keys_repaired=len(quoted),
        wrapped_remaining=remaining,
        applied=True,
    )
    print(
        f"repair done: keys={result.keys_repaired} "
        f"moved_rows={result.moved_rows} "
        f"dup_loci_skipped={result.dup_loci_skipped} "
        f"wrapped_remaining={result.wrapped_remaining}",
        file=sys.stderr,
        flush=True,
    )
    return result
