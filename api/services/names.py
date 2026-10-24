"""Gene-name prefix discovery (``GET /names``)."""
from __future__ import annotations

import sqlite3
from collections.abc import Sequence

from helpers.cursor import QueryError
from helpers.name_category import classify_name
from helpers.names import (
    MAX_PREFIX_EXPAND,
    parse_prefix_query,
    parse_taxid_list,
)


def _table_exists(conn: sqlite3.Connection, name: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1",
        (name,),
    ).fetchone()
    return row is not None


def _tax_predicate_sql(
    include: list[int],
    exclude: list[int],
) -> tuple[str, list]:
    """Return AND-clauses + params restricting ``f.annotation_key`` by lineage."""
    clauses: list[str] = []
    params: list = []
    if include:
        placeholders = ",".join("?" for _ in include)
        clauses.append(
            f"""
            f.annotation_key IN (
                SELECT annotation_key FROM annotation_lineage
                WHERE taxid IN ({placeholders})
            )
            """
        )
        params.extend(include)
    if exclude:
        placeholders = ",".join("?" for _ in exclude)
        clauses.append(
            f"""
            f.annotation_key NOT IN (
                SELECT annotation_key FROM annotation_lineage
                WHERE taxid IN ({placeholders})
            )
            """
        )
        params.extend(exclude)
    if not clauses:
        return "", []
    return " AND " + " AND ".join(clauses), params


def list_names_by_prefix(
    conn: sqlite3.Connection,
    prefix: str,
    *,
    include_taxids: str | Sequence[str] | Sequence[int] | None = None,
    exclude_taxids: str | Sequence[str] | Sequence[int] | None = None,
) -> dict:
    """
    Distinct symbol∪alias accessions matching ``prefix``, with scoped counts.

    Raises ``prefix_too_short``, ``invalid_prefix``, or ``prefix_too_broad``.
    """
    stem = parse_prefix_query(prefix)
    include = parse_taxid_list(include_taxids, field="include_taxids")
    exclude = parse_taxid_list(exclude_taxids, field="exclude_taxids")

    if not _table_exists(conn, "xref_file"):
        return {
            "prefix": stem,
            "include_taxids": include,
            "exclude_taxids": exclude,
            "total": 0,
            "results": [],
        }

    tax_sql, tax_params = _tax_predicate_sql(include, exclude)
    rows = conn.execute(
        f"""
        SELECT f.accession AS name,
               COUNT(DISTINCT f.annotation_key) AS n_annotations
        FROM xref_file f
        WHERE f.namespace IN ('symbol', 'alias')
          AND f.accession >= ?
          AND f.accession < ?
          {tax_sql}
        GROUP BY f.accession
        ORDER BY n_annotations DESC, f.accession ASC
        LIMIT ?
        """,
        (stem, stem + "\uffff", *tax_params, MAX_PREFIX_EXPAND + 1),
    ).fetchall()

    if len(rows) > MAX_PREFIX_EXPAND:
        raise QueryError(
            f"prefix expands to more than {MAX_PREFIX_EXPAND} accessions",
            code="prefix_too_broad",
            max=MAX_PREFIX_EXPAND,
            expanded=len(rows),
            prefix=stem,
        )

    results = [
        {
            "name": str(row["name"]),
            "category": classify_name(str(row["name"])),
            "n_annotations": int(row["n_annotations"]),
        }
        for row in rows
    ]
    return {
        "prefix": stem,
        "include_taxids": include,
        "exclude_taxids": exclude,
        "total": len(results),
        "results": results,
    }
