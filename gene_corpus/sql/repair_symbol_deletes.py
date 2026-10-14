"""
Delete-delta repair for via_level>0 symbol/alias (/hits globals).

Phase 1: two workers append delete rows to local part-*.tsv, then DELETE
from per-GFF gene_xref.
Phase 2: merge parts → patch xref_file / xref_meta (drop or decrement).

TSV line: annotation_key \\t namespace \\t accession \\t n_loci
"""
from __future__ import annotations

import os
import shutil
import sqlite3
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import IO, TextIO
from urllib.parse import quote

from gene_corpus.sql.per_gff import shard_sqlite_path
from gene_corpus.sql.schema import checkpoint_wal
from gene_corpus.sql.xref_index import ensure_schema

N_WORKERS = 2
SPILL_SUBDIR = "repair_symbol_deletes"
FINALIZED_NAME = "FINALIZED"
_SHARD_CACHE_KIB = -16_000
_SHARD_SOFT_HEAP = 128 * 1024 * 1024
_DELETE_CHUNK = 5_000
# Annotations between TSV fsync + cursor rewrite (crash redo window).
_SYNC_EVERY = 50

SpillLine = tuple[int, str, str, int]


@dataclass(frozen=True)
class ScanSummary:
    total_annotations: int
    processed: int
    skipped_resume: int
    missing_shard: int
    xref_rows_deleted: int
    tsv_lines_appended: int


@dataclass(frozen=True)
class FinalizeSummary:
    tsv_lines: int
    unique_deletes: int
    xref_file_deleted: int
    xref_meta_updated: int
    xref_meta_deleted: int


def default_spill_dir() -> Path:
    """``annotrieve-xrefs/repair_symbol_deletes/`` (repo root, not work-dir)."""
    # gene_corpus/sql/this_file → parents[2] = annotrieve-xrefs
    return Path(__file__).resolve().parents[2] / SPILL_SUBDIR


def part_tsv_path(spill_dir: Path, worker_id: int) -> Path:
    return spill_dir / f"part-{int(worker_id)}.tsv"


def cursor_path(spill_dir: Path, worker_id: int) -> Path:
    return spill_dir / f"part-{int(worker_id)}.cursor"


def finalized_path(spill_dir: Path) -> Path:
    return spill_dir / FINALIZED_NAME


def is_finalized(spill_dir: Path) -> bool:
    return finalized_path(spill_dir).is_file()


def write_finalized(spill_dir: Path) -> None:
    spill_dir.mkdir(parents=True, exist_ok=True)
    path = finalized_path(spill_dir)
    tmp = path.with_suffix(".FINALIZED.tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        fh.write("ok\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def read_cursor(spill_dir: Path, worker_id: int) -> int | None:
    path = cursor_path(spill_dir, worker_id)
    if not path.is_file():
        return None
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        return None
    return int(text)


def write_cursor(spill_dir: Path, worker_id: int, annotation_key: int) -> None:
    path = cursor_path(spill_dir, worker_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".cursor.tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        fh.write(f"{int(annotation_key)}\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def last_tsv_annotation_key(spill_dir: Path, worker_id: int) -> int | None:
    """Last annotation_key appearing in part TSV (redo hint)."""
    path = part_tsv_path(spill_dir, worker_id)
    if not path.is_file():
        return None
    last: int | None = None
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            last = int(line.split("\t", 1)[0])
    return last


def resume_from_key(spill_dir: Path, worker_id: int) -> int | None:
    """
    Key to reprocess: prefer cursor; else last TSV key.

    Callers skip jobs with key < resume_from; reprocess == resume_from.
    """
    cur = read_cursor(spill_dir, worker_id)
    if cur is not None:
        return cur
    return last_tsv_annotation_key(spill_dir, worker_id)


def list_annotations(conn: sqlite3.Connection) -> list[tuple[int, str]]:
    rows = conn.execute(
        """
        SELECT annotation_key, annotation_id
        FROM annotation
        ORDER BY annotation_key
        """
    ).fetchall()
    return [(int(r[0]), str(r[1])) for r in rows]


def _open_shard(path: Path, *, readonly: bool) -> sqlite3.Connection:
    resolved = path.resolve()
    if readonly:
        uri = f"file:{quote(resolved.as_posix())}?mode=ro"
        conn = sqlite3.connect(uri, uri=True)
    else:
        conn = sqlite3.connect(str(resolved))
    conn.execute(f"PRAGMA cache_size={_SHARD_CACHE_KIB}")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute(f"PRAGMA soft_heap_limit={_SHARD_SOFT_HEAP}")
    if readonly:
        conn.execute("PRAGMA query_only=ON")
    return conn


def _collect_delete_spill_lines(
    shard_conn: sqlite3.Connection,
    annotation_key: int,
) -> list[tuple[int, str, str, int]]:
    """
    Lines to append: via_level>0 (ns, acc) with no remaining via_level=0.

    n_loci = COUNT(DISTINCT local_id) among via_level>0 rows for that pair.
    """
    rows = shard_conn.execute(
        """
        SELECT x.namespace, x.accession, COUNT(DISTINCT x.local_id) AS n_loci
        FROM gene_xref x
        WHERE x.namespace IN ('symbol', 'alias')
          AND x.via_level > 0
          AND NOT EXISTS (
            SELECT 1 FROM gene_xref z
            WHERE z.namespace = x.namespace
              AND z.accession = x.accession
              AND z.via_level = 0
          )
        GROUP BY x.namespace, x.accession
        """
    ).fetchall()
    out: list[tuple[int, str, str, int]] = []
    for ns, acc, n_loci in rows:
        safe_acc = (
            str(acc).replace("\t", " ").replace("\n", " ").replace("\r", " ")
        )
        out.append((int(annotation_key), str(ns), safe_acc, int(n_loci)))
    return out


def _has_via_level_gt0(shard_conn: sqlite3.Connection) -> bool:
    hit = shard_conn.execute(
        """
        SELECT 1 FROM gene_xref
        WHERE namespace IN ('symbol', 'alias')
          AND via_level > 0
        LIMIT 1
        """
    ).fetchone()
    return hit is not None


def _delete_via_level_gt0(shard_conn: sqlite3.Connection) -> int:
    """DELETE all via_level>0 symbol/alias (caller already probed)."""
    shard_conn.execute(
        """
        DELETE FROM gene_xref
        WHERE namespace IN ('symbol', 'alias')
          AND via_level > 0
        """
    )
    n = int(shard_conn.execute("SELECT changes()").fetchone()[0])
    shard_conn.commit()
    return n


def _write_lines(fh: TextIO, lines: list[SpillLine]) -> None:
    for key, ns, acc, n_loci in lines:
        fh.write(f"{key}\t{ns}\t{acc}\t{n_loci}\n")


def _fsync_handle(fh: IO[str]) -> None:
    fh.flush()
    os.fsync(fh.fileno())


def _flush_cursor_batch(
    spill_dir: Path,
    worker_id: int,
    last_finished_key: int | None,
) -> None:
    """Durable-write cursor for last finished key (TSV already fsynced if needed)."""
    if last_finished_key is not None:
        write_cursor(spill_dir, worker_id, last_finished_key)


def process_annotation(
    *,
    per_gff_root: Path,
    annotation_key: int,
    annotation_id: str,
    apply: bool,
    tsv_fh: TextIO | None = None,
) -> tuple[int, bool, int]:
    """
    Returns ``(xref_rows_deleted, missing_shard, tsv_lines_written)``.

    Apply dirty path: collect → write+fsync TSV (if lines) → DELETE.
    Clean path: probe miss → no TSV/DELETE. Cursor batching is the worker's job.
    """
    shard = shard_sqlite_path(per_gff_root, annotation_id)
    if not shard.is_file():
        print(
            f"repair_symbol_deletes missing shard: {annotation_id}",
            file=sys.stderr,
            flush=True,
        )
        return 0, True, 0

    if not apply:
        shard_conn = _open_shard(shard, readonly=True)
        try:
            n = int(
                shard_conn.execute(
                    """
                    SELECT COUNT(*) FROM gene_xref
                    WHERE namespace IN ('symbol', 'alias')
                      AND via_level > 0
                    """
                ).fetchone()[0]
            )
            return n, False, 0
        finally:
            shard_conn.close()

    shard_conn = _open_shard(shard, readonly=False)
    try:
        if not _has_via_level_gt0(shard_conn):
            return 0, False, 0
        lines = _collect_delete_spill_lines(shard_conn, annotation_key)
        # Spill must be durable before DELETE (crash resume cannot re-derive).
        if lines:
            if tsv_fh is None:
                raise ValueError(
                    "tsv_fh required when spilling delete rows"
                )
            _write_lines(tsv_fh, lines)
            _fsync_handle(tsv_fh)
        deleted = _delete_via_level_gt0(shard_conn)
        return deleted, False, len(lines)
    finally:
        shard_conn.close()


def _worker_batch(
    worker_id: int,
    jobs: list[tuple[int, str]],
    per_gff_root: str,
    spill_dir: str,
    apply: bool,
) -> tuple[int, int, int, int, int]:
    """
    Returns
    ``(processed, skipped_resume, missing, deleted, tsv_lines)``.
    """
    root = Path(spill_dir)
    per_gff = Path(per_gff_root)
    if apply:
        root.mkdir(parents=True, exist_ok=True)

    resume_key = resume_from_key(root, worker_id) if apply else None
    processed = 0
    skipped = 0
    missing = 0
    deleted = 0
    tsv_lines = 0
    pending_since_sync = 0
    last_finished_key: int | None = None

    tsv_fh: TextIO | None = None
    try:
        if apply:
            tsv_path = part_tsv_path(root, worker_id)
            tsv_fh = tsv_path.open(
                "a", encoding="utf-8", buffering=1024 * 1024
            )

        for key, ann_id in jobs:
            if apply and resume_key is not None and key < resume_key:
                skipped += 1
                continue
            d, miss, n_lines = process_annotation(
                per_gff_root=per_gff,
                annotation_key=key,
                annotation_id=ann_id,
                apply=apply,
                tsv_fh=tsv_fh,
            )
            processed += 1
            deleted += d
            tsv_lines += n_lines
            if miss:
                missing += 1
                if apply:
                    _flush_cursor_batch(root, worker_id, last_finished_key)
                    pending_since_sync = 0
                    print(
                        f"repair_symbol_deletes worker-{worker_id}: "
                        f"stopping at missing shard key={key} "
                        f"({ann_id}); re-run after shard is present",
                        file=sys.stderr,
                        flush=True,
                    )
                    break
            if apply:
                last_finished_key = key
                pending_since_sync += 1
                if pending_since_sync >= _SYNC_EVERY:
                    _flush_cursor_batch(root, worker_id, last_finished_key)
                    pending_since_sync = 0
            if processed % 50 == 0 or processed + skipped == len(jobs):
                print(
                    f"repair_symbol_deletes worker-{worker_id}: "
                    f"processed={processed} skipped={skipped}/{len(jobs)} "
                    f"missing={missing} deleted_xref={deleted} "
                    f"tsv_lines={tsv_lines}",
                    file=sys.stderr,
                    flush=True,
                )

        if apply and last_finished_key is not None:
            _flush_cursor_batch(root, worker_id, last_finished_key)
    finally:
        if tsv_fh is not None:
            try:
                tsv_fh.flush()
            except OSError:
                pass
            tsv_fh.close()

    return processed, skipped, missing, deleted, tsv_lines


def partition_jobs(
    annotations: list[tuple[int, str]],
) -> list[list[tuple[int, str]]]:
    parts: list[list[tuple[int, str]]] = [[] for _ in range(N_WORKERS)]
    for i, job in enumerate(annotations):
        parts[i % N_WORKERS].append(job)
    return parts


def scan_complete(
    conn: sqlite3.Connection,
    spill_dir: Path,
) -> bool:
    """True when each worker cursor equals the last key in its partition."""
    annotations = list_annotations(conn)
    if not annotations:
        return False
    parts = partition_jobs(annotations)
    for wid, jobs in enumerate(parts):
        if not jobs:
            continue
        cur = read_cursor(spill_dir, wid)
        if cur is None:
            return False
        last_key = jobs[-1][0]
        if cur != last_key:
            return False
    return True


def run_scan(
    conn: sqlite3.Connection,
    *,
    work_dir: Path,
    spill_dir: Path | None = None,
    apply: bool = False,
) -> ScanSummary:
    annotations = list_annotations(conn)
    per_gff_root = work_dir / "per_gff"
    spill = spill_dir if spill_dir is not None else default_spill_dir()
    parts = partition_jobs(annotations)

    if apply:
        spill.mkdir(parents=True, exist_ok=True)

    if not annotations:
        return ScanSummary(
            total_annotations=0,
            processed=0,
            skipped_resume=0,
            missing_shard=0,
            xref_rows_deleted=0,
            tsv_lines_appended=0,
        )

    processed = 0
    skipped = 0
    missing = 0
    deleted = 0
    tsv_lines = 0

    if not apply:
        # Sequential dry-run (no pool).
        for wid, jobs in enumerate(parts):
            if not jobs:
                continue
            p, s, m, d, t = _worker_batch(
                wid, jobs, str(per_gff_root), str(spill), False
            )
            processed += p
            skipped += s
            missing += m
            deleted += d
            tsv_lines += t
    else:
        with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
            futs = [
                pool.submit(
                    _worker_batch,
                    wid,
                    jobs,
                    str(per_gff_root),
                    str(spill),
                    True,
                )
                for wid, jobs in enumerate(parts)
                if jobs
            ]
            for fut in as_completed(futs):
                p, s, m, d, t = fut.result()
                processed += p
                skipped += s
                missing += m
                deleted += d
                tsv_lines += t

    return ScanSummary(
        total_annotations=len(annotations),
        processed=processed,
        skipped_resume=skipped,
        missing_shard=missing,
        xref_rows_deleted=deleted,
        tsv_lines_appended=tsv_lines,
    )


def _load_delete_rows(
    spill_dir: Path,
) -> tuple[dict[tuple[int, str, str], int], int]:
    """Dedupe (key, ns, acc) → MAX(n_loci). Returns (map, raw_line_count)."""
    out: dict[tuple[int, str, str], int] = {}
    n_lines = 0
    for wid in range(N_WORKERS):
        path = part_tsv_path(spill_dir, wid)
        if not path.is_file():
            continue
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                n_lines += 1
                parts = line.split("\t")
                if len(parts) != 4:
                    raise ValueError(f"bad TSV line in {path}: {line!r}")
                key_s, ns, acc, loci_s = parts
                key = int(key_s)
                n_loci = int(loci_s)
                pair = (key, ns, acc)
                prev = out.get(pair)
                if prev is None or n_loci > prev:
                    out[pair] = n_loci
    return out, n_lines


def _batched_delete_xref_file(
    conn: sqlite3.Connection,
    deletes: dict[tuple[int, str, str], int],
) -> int:
    """DELETE xref_file rows in chunks; return rows removed."""
    batch: list[tuple[int, str, str]] = []
    deleted = 0

    def _flush() -> None:
        nonlocal batch, deleted
        if not batch:
            return
        before = conn.total_changes
        conn.executemany(
            """
            DELETE FROM xref_file
            WHERE annotation_key = ?
              AND namespace = ?
              AND accession = ?
            """,
            batch,
        )
        deleted += conn.total_changes - before
        batch = []

    for key, ns, acc in deletes:
        batch.append((key, ns, acc))
        if len(batch) >= _DELETE_CHUNK:
            _flush()
    _flush()
    return deleted


def finalize(
    conn: sqlite3.Connection,
    *,
    spill_dir: Path | None = None,
) -> FinalizeSummary:
    """Patch xref_file / xref_meta from delete TSVs."""
    ensure_schema(conn)
    spill = spill_dir if spill_dir is not None else default_spill_dir()

    if is_finalized(spill):
        print(
            "repair_symbol_deletes finalize: already FINALIZED (no-op)",
            file=sys.stderr,
            flush=True,
        )
        return FinalizeSummary(
            tsv_lines=0,
            unique_deletes=0,
            xref_file_deleted=0,
            xref_meta_updated=0,
            xref_meta_deleted=0,
        )

    if not scan_complete(conn, spill):
        raise RuntimeError(
            "repair_symbol_deletes finalize: scan incomplete "
            "(every partition cursor must cover its keys; missing shards block)"
        )

    deletes, n_lines = _load_delete_rows(spill)

    # Aggregate per (ns, acc)
    sum_loci: dict[tuple[str, str], int] = defaultdict(int)
    keys_per_acc: dict[tuple[str, str], set[int]] = defaultdict(set)
    for (key, ns, acc), n_loci in deletes.items():
        sum_loci[(ns, acc)] += n_loci
        keys_per_acc[(ns, acc)].add(key)

    file_deleted = _batched_delete_xref_file(conn, deletes)

    meta_updated = 0
    meta_deleted = 0
    for (ns, acc), loci in sum_loci.items():
        n_ann_rm = len(keys_per_acc[(ns, acc)])
        row = conn.execute(
            """
            SELECT n_annotations, n_loci FROM xref_meta
            WHERE namespace = ? AND accession = ?
            """,
            (ns, acc),
        ).fetchone()
        if row is None:
            continue
        n_ann = max(0, int(row[0]) - n_ann_rm)
        n_loci = max(0, int(row[1]) - loci)
        if n_ann == 0:
            conn.execute(
                """
                DELETE FROM xref_meta
                WHERE namespace = ? AND accession = ?
                """,
                (ns, acc),
            )
            meta_deleted += int(conn.execute("SELECT changes()").fetchone()[0])
        else:
            conn.execute(
                """
                UPDATE xref_meta
                SET n_annotations = ?, n_loci = ?
                WHERE namespace = ? AND accession = ?
                """,
                (n_ann, n_loci, ns, acc),
            )
            meta_updated += int(conn.execute("SELECT changes()").fetchone()[0])

    conn.commit()
    checkpoint_wal(conn)
    write_finalized(spill)

    print(
        f"repair_symbol_deletes finalize: tsv_lines={n_lines:,} "
        f"unique_deletes={len(deletes):,} "
        f"xref_file_deleted={file_deleted:,} "
        f"meta_updated={meta_updated:,} meta_deleted={meta_deleted:,}",
        file=sys.stderr,
        flush=True,
    )

    return FinalizeSummary(
        tsv_lines=n_lines,
        unique_deletes=len(deletes),
        xref_file_deleted=file_deleted,
        xref_meta_updated=meta_updated,
        xref_meta_deleted=meta_deleted,
    )


def cleanup_spill(spill_dir: Path | None = None) -> None:
    spill = spill_dir if spill_dir is not None else default_spill_dir()
    if spill.is_dir():
        shutil.rmtree(spill, ignore_errors=True)


__all__ = [
    "FinalizeSummary",
    "N_WORKERS",
    "ScanSummary",
    "cleanup_spill",
    "cursor_path",
    "default_spill_dir",
    "finalize",
    "finalized_path",
    "is_finalized",
    "last_tsv_annotation_key",
    "list_annotations",
    "part_tsv_path",
    "partition_jobs",
    "process_annotation",
    "read_cursor",
    "resume_from_key",
    "run_scan",
    "scan_complete",
    "write_cursor",
    "write_finalized",
]
