"""
Delete stale symbol/alias rows from gene_hit via xref_file anti-join.

Phase 1 (scan): for each of symbol/alias, stream gene_hit rows whose
(namespace, accession, annotation_key) is missing from xref_file (and
whose annotation_key is in xref_file_status) into a local part-*.tsv.

Phase 2 (finalize): PK-ordered chunked DELETE from gene_hit.

TSV line: namespace \\t accession \\t annotation_key \\t local_id
"""
from __future__ import annotations

import os
import shutil
import sqlite3
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterator, TextIO
from urllib.parse import quote

from gene_corpus.sql.schema import checkpoint_wal
from gene_corpus.sql.xref_index import ensure_schema

NAMESPACES: tuple[str, ...] = ("alias", "symbol")
SPILL_SUBDIR = "repair_gene_hit_deletes"
FINALIZED_NAME = "FINALIZED"
_DELETE_CHUNK = 20_000
_COMMIT_EVERY = 250_000
_CHECKPOINT_EVERY_COMMITS = 4
_SAMPLE_ACCESSIONS = 20
_STALE_FRACTION_ABORT = 0.5
_MMAP_BYTES = 1 * 1024 * 1024 * 1024
_RO_CACHE_KIB = -131_072  # 128 MiB

SpillRow = tuple[str, str, int, int]
# namespace -> (gene_hit_rows, xref_file_rows)
NsCounts = dict[str, tuple[int, int]]


@dataclass(frozen=True)
class NamespaceScanSummary:
    namespace: str
    gene_hit_rows: int
    xref_file_rows: int
    stale_rows: int
    distinct_groups: int
    sample_accessions: tuple[str, ...]


@dataclass(frozen=True)
class ScanSummary:
    namespaces: tuple[NamespaceScanSummary, ...]
    tsv_lines_appended: int

    @property
    def stale_rows(self) -> int:
        return sum(n.stale_rows for n in self.namespaces)


@dataclass(frozen=True)
class FinalizeSummary:
    tsv_lines: int
    deleted: int
    skipped_resume: int


class XrefFileGuardError(RuntimeError):
    """xref_file looks empty/incomplete for symbol/alias; refuse without --force."""


def default_spill_dir() -> Path:
    """``annotrieve-xrefs/repair_gene_hit_deletes/`` (repo root, not work-dir)."""
    return Path(__file__).resolve().parents[2] / SPILL_SUBDIR


def part_tsv_path(spill_dir: Path, namespace: str) -> Path:
    return spill_dir / f"part-{namespace}.tsv"


def scan_done_path(spill_dir: Path, namespace: str) -> Path:
    return spill_dir / f"part-{namespace}.done"


def finalize_cursor_path(spill_dir: Path) -> Path:
    return spill_dir / "finalize.cursor"


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


def is_scan_done(spill_dir: Path, namespace: str) -> bool:
    return scan_done_path(spill_dir, namespace).is_file()


def write_scan_done(spill_dir: Path, namespace: str) -> None:
    path = scan_done_path(spill_dir, namespace)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".done.tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        fh.write("ok\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def scan_complete(spill_dir: Path) -> bool:
    return all(is_scan_done(spill_dir, ns) for ns in NAMESPACES)


def _fsync_handle(fh: IO[str]) -> None:
    fh.flush()
    os.fsync(fh.fileno())


def _open_readonly(db_path: Path) -> sqlite3.Connection:
    resolved = db_path.resolve()
    uri = f"file:{quote(resolved.as_posix())}?mode=ro"
    conn = sqlite3.connect(uri, uri=True)
    conn.execute(f"PRAGMA cache_size={_RO_CACHE_KIB}")
    conn.execute(f"PRAGMA mmap_size={_MMAP_BYTES}")
    conn.execute("PRAGMA query_only=ON")
    return conn


def _count_table(
    conn: sqlite3.Connection,
    table: str,
    namespace: str,
) -> int:
    if table not in ("gene_hit", "xref_file"):
        raise ValueError(f"unsupported count table: {table}")
    return int(
        conn.execute(
            f"SELECT COUNT(*) FROM {table} WHERE namespace = ?",
            (namespace,),
        ).fetchone()[0]
    )


def check_xref_file_guard(
    conn: sqlite3.Connection,
    *,
    force: bool = False,
) -> NsCounts:
    """
    Count once and return ``{ns: (gene_hit_rows, xref_file_rows)}``.

    Raises ``XrefFileGuardError`` when gene_hit has symbol/alias rows but
    xref_file has none for that namespace, unless ``force``.
    """
    ensure_schema(conn)
    counts: NsCounts = {}
    problems: list[str] = []
    for ns in NAMESPACES:
        gh = _count_table(conn, "gene_hit", ns)
        xf = _count_table(conn, "xref_file", ns)
        counts[ns] = (gh, xf)
        if gh > 0 and xf == 0:
            problems.append(
                f"{ns}: gene_hit={gh:,} but xref_file=0"
            )
    if problems and not force:
        raise XrefFileGuardError(
            "xref_file looks empty for symbol/alias "
            f"({'; '.join(problems)}); pass --force to override"
        )
    return counts


def _check_stale_fraction(
    summaries: list[NamespaceScanSummary],
    *,
    force: bool,
) -> None:
    """Abort when stale_rows >= 50% of gene_hit for a namespace unless force."""
    if force:
        return
    problems: list[str] = []
    for ns in summaries:
        if ns.gene_hit_rows <= 0:
            continue
        frac = ns.stale_rows / ns.gene_hit_rows
        if frac >= _STALE_FRACTION_ABORT:
            problems.append(
                f"{ns.namespace}: stale={ns.stale_rows:,}/"
                f"{ns.gene_hit_rows:,} ({frac:.0%})"
            )
    if problems:
        raise XrefFileGuardError(
            "stale fraction looks suspiciously high "
            f"({'; '.join(problems)}); pass --force to override"
        )


def _iter_stale_rows(
    conn: sqlite3.Connection,
    namespace: str,
) -> Iterator[SpillRow]:
    """
    Yield stale gene_hit PKs for ``namespace`` in PK order.

    Only annotation_keys present in xref_file_status are considered.
    """
    cur = conn.execute(
        """
        SELECT g.accession, g.annotation_key, g.local_id
        FROM gene_hit g
        WHERE g.namespace = ?
          AND EXISTS (
            SELECT 1 FROM xref_file_status s
            WHERE s.annotation_key = g.annotation_key
          )
          AND NOT EXISTS (
            SELECT 1 FROM xref_file f
            WHERE f.namespace = g.namespace
              AND f.accession = g.accession
              AND f.annotation_key = g.annotation_key
          )
        ORDER BY g.accession, g.annotation_key, g.local_id
        """,
        (namespace,),
    )
    for acc, key, local_id in cur:
        safe_acc = (
            str(acc).replace("\t", " ").replace("\n", " ").replace("\r", " ")
        )
        yield (namespace, safe_acc, int(key), int(local_id))


def _scan_one_namespace(
    db_path: str,
    namespace: str,
    spill_dir: str,
    apply: bool,
) -> tuple[str, int, int, list[str]]:
    """
    Worker entry: anti-join scan for one namespace (no COUNT(*)).

    Returns ``(namespace, stale_rows, distinct_groups, sample_accessions)``.
    """
    path = Path(db_path)
    spill = Path(spill_dir)
    conn = _open_readonly(path)
    try:
        if apply and is_scan_done(spill, namespace):
            # Resume: already spilled; recount from TSV for summary.
            stale, groups, samples = _summarize_tsv(
                part_tsv_path(spill, namespace)
            )
            return (namespace, stale, groups, samples)

        stale = 0
        groups: set[tuple[str, int]] = set()
        samples: list[str] = []
        seen_acc: set[str] = set()

        tsv_fh: TextIO | None = None
        try:
            if apply:
                spill.mkdir(parents=True, exist_ok=True)
                tsv_path = part_tsv_path(spill, namespace)
                # Rewrite from scratch (scan is idempotent / restartable).
                tsv_fh = tsv_path.open(
                    "w", encoding="utf-8", buffering=1024 * 1024
                )

            for ns, acc, key, local_id in _iter_stale_rows(conn, namespace):
                stale += 1
                groups.add((acc, key))
                if acc not in seen_acc and len(samples) < _SAMPLE_ACCESSIONS:
                    seen_acc.add(acc)
                    samples.append(acc)
                if tsv_fh is not None:
                    tsv_fh.write(f"{ns}\t{acc}\t{key}\t{local_id}\n")
                    if stale % 50_000 == 0:
                        _fsync_handle(tsv_fh)
                        print(
                            f"repair_gene_hit_deletes scan {namespace}: "
                            f"stale={stale:,}",
                            file=sys.stderr,
                            flush=True,
                        )

            if tsv_fh is not None:
                _fsync_handle(tsv_fh)
        finally:
            if tsv_fh is not None:
                tsv_fh.close()

        # .done is written by run_scan after the stale-fraction guard passes.
        print(
            f"repair_gene_hit_deletes scan {namespace}: "
            f"stale={stale:,} groups={len(groups):,}",
            file=sys.stderr,
            flush=True,
        )
        return (namespace, stale, len(groups), samples)
    finally:
        conn.close()


def _summarize_tsv(
    path: Path,
) -> tuple[int, int, list[str]]:
    if not path.is_file():
        return 0, 0, []
    stale = 0
    groups: set[tuple[str, int]] = set()
    samples: list[str] = []
    seen_acc: set[str] = set()
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 4:
                raise ValueError(f"bad TSV line in {path}: {line!r}")
            _ns, acc, key_s, _local = parts
            stale += 1
            groups.add((acc, int(key_s)))
            if acc not in seen_acc and len(samples) < _SAMPLE_ACCESSIONS:
                seen_acc.add(acc)
                samples.append(acc)
    return stale, len(groups), samples


def run_scan(
    conn: sqlite3.Connection,
    *,
    db_path: Path,
    spill_dir: Path | None = None,
    apply: bool = False,
    force: bool = False,
) -> ScanSummary:
    """
    Anti-join scan for symbol/alias.

    Dry-run (``apply=False``): count only, no spill.
    Apply: write part-*.tsv + done markers (2 threads).
    """
    spill = spill_dir if spill_dir is not None else default_spill_dir()
    counts = check_xref_file_guard(conn, force=force)

    if apply:
        spill.mkdir(parents=True, exist_ok=True)

    by_ns: dict[str, tuple[int, int, list[str]]] = {}
    with ThreadPoolExecutor(max_workers=len(NAMESPACES)) as pool:
        futs = [
            pool.submit(
                _scan_one_namespace,
                str(db_path),
                ns,
                str(spill),
                apply,
            )
            for ns in NAMESPACES
        ]
        for fut in as_completed(futs):
            ns_name, stale, groups, samples = fut.result()
            by_ns[ns_name] = (stale, groups, samples)

    results: list[NamespaceScanSummary] = []
    tsv_lines = 0
    for ns in NAMESPACES:
        stale, groups, samples = by_ns[ns]
        gh, xf = counts[ns]
        if apply:
            tsv_lines += stale
        results.append(
            NamespaceScanSummary(
                namespace=ns,
                gene_hit_rows=gh,
                xref_file_rows=xf,
                stale_rows=stale,
                distinct_groups=groups,
                sample_accessions=tuple(samples),
            )
        )

    _check_stale_fraction(results, force=force)

    # Mark scan durable only after the fraction guard accepts the result.
    if apply:
        for ns in NAMESPACES:
            if not is_scan_done(spill, ns):
                write_scan_done(spill, ns)

    return ScanSummary(
        namespaces=tuple(results),
        tsv_lines_appended=tsv_lines,
    )


def _parse_spill_row(line: str, path: Path) -> SpillRow:
    line = line.rstrip("\n")
    parts = line.split("\t")
    if len(parts) != 4:
        raise ValueError(f"bad TSV line in {path}: {line!r}")
    ns, acc, key_s, local_s = parts
    return (ns, acc, int(key_s), int(local_s))


def _iter_spill_rows(spill_dir: Path) -> Iterator[SpillRow]:
    """Stream spill rows in PK order (namespaces sorted: alias, symbol)."""
    for ns in NAMESPACES:
        path = part_tsv_path(spill_dir, ns)
        if not path.is_file():
            continue
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                if not line.strip():
                    continue
                yield _parse_spill_row(line, path)


def read_finalize_cursor(spill_dir: Path) -> SpillRow | None:
    path = finalize_cursor_path(spill_dir)
    if not path.is_file():
        return None
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        return None
    parts = text.split("\t")
    if len(parts) != 4:
        raise ValueError(f"bad finalize cursor: {text!r}")
    ns, acc, key_s, local_s = parts
    return (ns, acc, int(key_s), int(local_s))


def write_finalize_cursor(spill_dir: Path, row: SpillRow) -> None:
    path = finalize_cursor_path(spill_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".cursor.tmp")
    ns, acc, key, local_id = row
    with tmp.open("w", encoding="utf-8") as fh:
        fh.write(f"{ns}\t{acc}\t{key}\t{local_id}\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def finalize(
    conn: sqlite3.Connection,
    *,
    spill_dir: Path | None = None,
) -> FinalizeSummary:
    """PK-ordered chunked DELETE from gene_hit using spill TSVs."""
    spill = spill_dir if spill_dir is not None else default_spill_dir()

    if is_finalized(spill):
        print(
            "repair_gene_hit_deletes finalize: already FINALIZED (no-op)",
            file=sys.stderr,
            flush=True,
        )
        return FinalizeSummary(tsv_lines=0, deleted=0, skipped_resume=0)

    if not scan_complete(spill):
        raise RuntimeError(
            "repair_gene_hit_deletes finalize: scan incomplete "
            "(need part-{alias,symbol}.done)"
        )

    resume = read_finalize_cursor(spill)
    batch: list[SpillRow] = []
    deleted = 0
    skipped = 0
    seen = 0
    since_commit = 0
    commits = 0
    last_row: SpillRow | None = None
    past_resume = resume is None

    def _flush_batch() -> None:
        nonlocal batch, deleted, since_commit, commits, last_row
        if not batch:
            return
        before = conn.total_changes
        conn.executemany(
            """
            DELETE FROM gene_hit
            WHERE namespace = ?
              AND accession = ?
              AND annotation_key = ?
              AND local_id = ?
            """,
            batch,
        )
        deleted += conn.total_changes - before
        since_commit += len(batch)
        last_row = batch[-1]
        batch = []
        if since_commit >= _COMMIT_EVERY:
            conn.commit()
            if last_row is not None:
                write_finalize_cursor(spill, last_row)
            commits += 1
            since_commit = 0
            if commits % _CHECKPOINT_EVERY_COMMITS == 0:
                checkpoint_wal(conn)
            print(
                f"repair_gene_hit_deletes finalize: "
                f"deleted={deleted:,} skipped={skipped:,}",
                file=sys.stderr,
                flush=True,
            )

    for row in _iter_spill_rows(spill):
        seen += 1
        if not past_resume:
            assert resume is not None
            if row <= resume:
                skipped += 1
                if row == resume:
                    past_resume = True
                continue
            past_resume = True
        batch.append(row)
        if len(batch) >= _DELETE_CHUNK:
            _flush_batch()

    _flush_batch()
    if since_commit > 0 or last_row is not None:
        conn.commit()
        if last_row is not None:
            write_finalize_cursor(spill, last_row)

    checkpoint_wal(conn)
    write_finalized(spill)

    print(
        f"repair_gene_hit_deletes finalize: tsv_lines={seen:,} "
        f"deleted={deleted:,} skipped_resume={skipped:,}",
        file=sys.stderr,
        flush=True,
    )
    return FinalizeSummary(
        tsv_lines=seen,
        deleted=deleted,
        skipped_resume=skipped,
    )


def cleanup_spill(spill_dir: Path | None = None) -> None:
    spill = spill_dir if spill_dir is not None else default_spill_dir()
    if spill.is_dir():
        shutil.rmtree(spill, ignore_errors=True)


__all__ = [
    "FinalizeSummary",
    "NAMESPACES",
    "NamespaceScanSummary",
    "ScanSummary",
    "XrefFileGuardError",
    "check_xref_file_guard",
    "cleanup_spill",
    "default_spill_dir",
    "finalize",
    "finalize_cursor_path",
    "finalized_path",
    "is_finalized",
    "is_scan_done",
    "part_tsv_path",
    "read_finalize_cursor",
    "run_scan",
    "scan_complete",
    "scan_done_path",
    "write_finalize_cursor",
    "write_finalized",
    "write_scan_done",
]
