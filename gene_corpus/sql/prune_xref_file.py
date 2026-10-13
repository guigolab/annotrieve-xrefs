"""
Prune ``xref_file`` to ``SEARCHABLE_NAMESPACES`` and reclaim disk.

Dry-run reports what would change; apply deletes demoted rows, recomputes
``xref_file_status.placed``, then ``VACUUM``.
"""
from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path

from gene_corpus.namespaces.tiers import SEARCHABLE_NAMESPACES
from gene_corpus.sql.schema import checkpoint_wal
from gene_corpus.sql.xref_index import ensure_schema


def _ns_placeholders() -> tuple[str, list[str]]:
    ns = sorted(SEARCHABLE_NAMESPACES)
    return ",".join("?" for _ in ns), ns


@dataclass(frozen=True)
class PruneSummary:
    """Counts for dry-run or post-apply reporting."""

    rows_total: int
    rows_to_delete: int
    namespaces_to_delete: tuple[tuple[str, int], ...]
    placed_to_zero: int
    rows_kept: int


@dataclass(frozen=True)
class PruneResult:
    """Outcome of ``apply_prune``."""

    summary: PruneSummary
    rows_deleted: int
    placed_flipped: int
    size_before: int | None
    size_after: int | None
    vacuumed: bool


def summarize_prune(conn: sqlite3.Connection) -> PruneSummary:
    """Count demoted ``xref_file`` rows and annotations that would lose presence."""
    ensure_schema(conn)
    placeholders, ns_list = _ns_placeholders()

    rows_total = int(
        conn.execute("SELECT COUNT(*) FROM xref_file").fetchone()[0]
    )
    rows_to_delete = int(
        conn.execute(
            f"""
            SELECT COUNT(*) FROM xref_file
            WHERE namespace NOT IN ({placeholders})
            """,
            ns_list,
        ).fetchone()[0]
    )
    ns_rows = conn.execute(
        f"""
        SELECT namespace, COUNT(*) AS n
        FROM xref_file
        WHERE namespace NOT IN ({placeholders})
        GROUP BY namespace
        ORDER BY n DESC, namespace
        """,
        ns_list,
    ).fetchall()
    namespaces_to_delete = tuple(
        (str(row[0]), int(row[1])) for row in ns_rows
    )

    # Keys currently placed=1 that have no searchable xref_file row remaining
    # after the delete (includes keys that only had demoted namespaces).
    placed_to_zero = int(
        conn.execute(
            f"""
            SELECT COUNT(*)
            FROM xref_file_status s
            WHERE s.placed = 1
              AND NOT EXISTS (
                SELECT 1 FROM xref_file x
                WHERE x.annotation_key = s.annotation_key
                  AND x.namespace IN ({placeholders})
              )
            """,
            ns_list,
        ).fetchone()[0]
    )
    return PruneSummary(
        rows_total=rows_total,
        rows_to_delete=rows_to_delete,
        namespaces_to_delete=namespaces_to_delete,
        placed_to_zero=placed_to_zero,
        rows_kept=rows_total - rows_to_delete,
    )


def apply_prune(
    conn: sqlite3.Connection,
    *,
    db_path: Path | None = None,
) -> PruneResult:
    """
    Delete non-searchable ``xref_file`` rows, fix ``placed``, then VACUUM.

    ``db_path`` is used only for before/after file size reporting.
    """
    summary = summarize_prune(conn)
    placeholders, ns_list = _ns_placeholders()

    size_before: int | None = None
    if db_path is not None and db_path.is_file():
        size_before = db_path.stat().st_size

    conn.execute(
        f"""
        DELETE FROM xref_file
        WHERE namespace NOT IN ({placeholders})
        """,
        ns_list,
    )
    rows_deleted = int(conn.execute("SELECT changes()").fetchone()[0])

    conn.execute(
        """
        UPDATE xref_file_status
        SET placed = 0
        WHERE placed = 1
          AND NOT EXISTS (
            SELECT 1 FROM xref_file x
            WHERE x.annotation_key = xref_file_status.annotation_key
          )
        """
    )
    placed_flipped = int(conn.execute("SELECT changes()").fetchone()[0])
    conn.commit()

    checkpoint_wal(conn)
    conn.execute("VACUUM")
    vacuumed = True

    size_after: int | None = None
    if db_path is not None and db_path.is_file():
        size_after = db_path.stat().st_size

    return PruneResult(
        summary=summary,
        rows_deleted=rows_deleted,
        placed_flipped=placed_flipped,
        size_before=size_before,
        size_after=size_after,
        vacuumed=vacuumed,
    )
