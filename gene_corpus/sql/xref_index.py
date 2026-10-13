"""
Build xref_file from per-GFF tier_a_counts (annotation presence index).

Tables are created by ``ensure_schema``; merge schema stays unchanged.
``index_annotation`` is safe to call from sync after ``attach_shard``.
"""
from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

from gene_corpus.namespaces.tiers import SEARCHABLE_NAMESPACES
from gene_corpus.sql.per_gff import shard_sqlite_path
from gene_corpus.sql.schema import checkpoint_wal


def _searchable_ns_sql() -> tuple[str, list[str]]:
    ns = sorted(SEARCHABLE_NAMESPACES)
    return ",".join("?" for _ in ns), ns

_TABLES_SQL = """
CREATE TABLE IF NOT EXISTS xref_file (
    namespace TEXT NOT NULL,
    accession TEXT NOT NULL,
    annotation_key INTEGER NOT NULL,
    PRIMARY KEY (namespace, accession, annotation_key)
) WITHOUT ROWID;

CREATE TABLE IF NOT EXISTS xref_file_status (
    annotation_key INTEGER PRIMARY KEY,
    placed INTEGER NOT NULL
) WITHOUT ROWID;
"""


def ensure_schema(conn: sqlite3.Connection) -> None:
    """Create xref_file / xref_file_status if missing."""
    conn.executescript(_TABLES_SQL)
    conn.commit()


def _status_exists(conn: sqlite3.Connection, annotation_key: int) -> bool:
    row = conn.execute(
        "SELECT 1 FROM xref_file_status WHERE annotation_key = ? LIMIT 1",
        (int(annotation_key),),
    ).fetchone()
    return row is not None


def index_annotation(
    conn: sqlite3.Connection,
    per_gff_root: Path,
    annotation_key: int,
    annotation_id: str,
) -> bool:
    """
    Copy searchable rows from one shard's ``tier_a_counts`` into ``xref_file``.

    Only ``SEARCHABLE_NAMESPACES`` are written (API presence allowlist).
    Returns immediately when ``xref_file_status`` already has
    ``annotation_key``. Zero searchable rows record ``placed = 0`` and return
    False. One or more rows insert presence and record ``placed = 1``.

    TODO (return contract): False means both "already in xref_file_status"
    and "just wrote placed=0". Clarify before wiring sync (tri-state or
    require callers to read placed).
    """
    key_i = int(annotation_key)
    if _status_exists(conn, key_i):
        return False

    shard = shard_sqlite_path(per_gff_root, annotation_id)
    if not shard.is_file():
        # TODO (missing shard): FileNotFoundError does not write
        # xref_file_status, so the annotation stays pending forever and is
        # retried every run. Decide on a permanent skip or failed status later.
        raise FileNotFoundError(f"missing shard for {annotation_id}: {shard}")

    try:
        conn.execute("ATTACH DATABASE ? AS shard", (str(shard.resolve()),))
    except Exception:
        conn.rollback()
        raise

    try:
        try:
            placeholders, ns_list = _searchable_ns_sql()
            n = int(
                conn.execute(
                    f"""
                    SELECT COUNT(*) FROM shard.tier_a_counts
                    WHERE namespace IN ({placeholders})
                    """,
                    ns_list,
                ).fetchone()[0]
            )
            if n == 0:
                conn.execute(
                    "INSERT INTO xref_file_status (annotation_key, placed) "
                    "VALUES (?, 0)",
                    (key_i,),
                )
                conn.commit()
                return False

            conn.execute(
                f"""
                INSERT INTO xref_file (namespace, accession, annotation_key)
                SELECT namespace, accession, ?
                FROM shard.tier_a_counts
                WHERE namespace IN ({placeholders})
                """,
                (key_i, *ns_list),
            )
            conn.execute(
                "INSERT INTO xref_file_status (annotation_key, placed) "
                "VALUES (?, 1)",
                (key_i,),
            )
            conn.commit()
            return True
        except Exception:
            conn.rollback()
            raise
    finally:
        conn.execute("DETACH DATABASE shard")


def index_missing_annotations(
    conn: sqlite3.Connection,
    per_gff_root: Path,
    *,
    checkpoint_every: int = 50,
) -> tuple[int, int, int]:
    """
    Index every ``annotation`` row missing from ``xref_file_status``.

    Returns ``(placed, empty, errors)``.
    """
    ensure_schema(conn)
    pending = conn.execute(
        """
        SELECT a.annotation_key, a.annotation_id
        FROM annotation a
        LEFT JOIN xref_file_status s
            ON s.annotation_key = a.annotation_key
        WHERE s.annotation_key IS NULL
        ORDER BY a.annotation_key
        """
    ).fetchall()

    total = len(pending)
    progress_every = 50
    if checkpoint_every < 1:
        checkpoint_every = progress_every

    print(
        f"index: {total} annotation(s) missing from xref_file_status",
        file=sys.stderr,
        flush=True,
    )

    placed = 0
    empty = 0
    errors = 0
    newly = 0
    for i, (annotation_key, annotation_id) in enumerate(pending, start=1):
        try:
            ok = index_annotation(
                conn,
                per_gff_root,
                int(annotation_key),
                str(annotation_id),
            )
            # Pending keys only: False means we just recorded placed=0.
            if ok:
                placed += 1
            else:
                empty += 1
            newly += 1
        except Exception as exc:  # noqa: BLE001
            # TODO (missing shard): failures (e.g. FileNotFoundError) do not
            # write xref_file_status, so the key stays pending and is retried
            # every run. Decide on a permanent skip or failed status later.
            errors += 1
            print(
                f"index fail {annotation_id}: {exc}",
                file=sys.stderr,
            )
        if i % progress_every == 0 or i == total:
            print(
                f"index progress: {i}/{total} "
                f"(placed={placed} empty={empty} errors={errors})",
                file=sys.stderr,
                flush=True,
            )
        if newly > 0 and newly % checkpoint_every == 0:
            checkpoint_wal(conn)

    if newly > 0 and newly % checkpoint_every != 0:
        checkpoint_wal(conn)

    return placed, empty, errors
