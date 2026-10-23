"""
One-time backfill: ALTER annotation for release_date, then UPDATE from
the Annotrieve annotations report TSV (selected_fields=release_date).
"""
from __future__ import annotations

import csv
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from gene_corpus.sql.schema import SCHEMA_VERSION


RELEASE_DATE_COLUMN = "release_date"
# Extended report column; defaults (incl. annotation_id) are always present.
REPORT_RELEASE_DATE_FIELDS = "release_date"


@dataclass(frozen=True)
class SeedReleaseDateStats:
    report_rows: int
    with_date: int
    matched: int
    updated: int
    missing_in_db: int
    empty_in_report: int
    column_added: bool


def annotation_has_release_date_column(conn: sqlite3.Connection) -> bool:
    cols = {
        str(row[1])
        for row in conn.execute("PRAGMA table_info(annotation)")
    }
    return RELEASE_DATE_COLUMN in cols


def ensure_release_date_column(conn: sqlite3.Connection) -> bool:
    """ADD COLUMN release_date TEXT when missing. Returns True if ALTER ran."""
    if annotation_has_release_date_column(conn):
        return False
    conn.execute(
        f"ALTER TABLE annotation ADD COLUMN {RELEASE_DATE_COLUMN} TEXT"
    )
    conn.execute(
        "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
        ("schema_version", str(SCHEMA_VERSION)),
    )
    conn.commit()
    return True


def parse_release_date_tsv(report: Path) -> tuple[dict[str, str], int, int]:
    """
    Parse report TSV → (annotation_id → release_date, total_rows, empty_dates).

    Only non-empty release_date values are included in the map.
    """
    out: dict[str, str] = {}
    total = 0
    empty = 0
    with report.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"empty or headerless TSV: {report}")
        names = list(reader.fieldnames)
        if "annotation_id" not in names:
            raise ValueError(
                f"TSV missing annotation_id column: {report} "
                f"(columns={names!r})"
            )
        if RELEASE_DATE_COLUMN not in names:
            raise ValueError(
                f"TSV missing {RELEASE_DATE_COLUMN} column: {report} "
                f"(columns={names!r}); "
                "fetch with selected_fields=release_date"
            )
        for row in reader:
            ann_id = (row.get("annotation_id") or "").strip()
            if not ann_id:
                continue
            total += 1
            raw = (row.get(RELEASE_DATE_COLUMN) or "").strip()
            if not raw:
                empty += 1
                continue
            out[ann_id] = raw
    return out, total, empty


def seed_release_dates(
    conn: sqlite3.Connection,
    id_to_date: Mapping[str, str],
    *,
    apply: bool,
    report_rows: int = 0,
    empty_in_report: int = 0,
) -> SeedReleaseDateStats:
    """
    Ensure column exists (on --apply) and UPDATE annotation.release_date
    by annotation_id. Dry-run counts matches without writing.
    """
    column_added = False
    if apply:
        column_added = ensure_release_date_column(conn)
    elif not annotation_has_release_date_column(conn):
        column_added = True  # would be added on --apply

    db_ids = {
        str(r[0])
        for r in conn.execute("SELECT annotation_id FROM annotation")
    }
    matched_ids = sorted(set(id_to_date) & db_ids)
    missing_in_db = len(set(id_to_date) - db_ids)

    updated = 0
    if apply and matched_ids:
        conn.executemany(
            f"""
            UPDATE annotation
            SET {RELEASE_DATE_COLUMN} = ?
            WHERE annotation_id = ?
            """,
            [(id_to_date[aid], aid) for aid in matched_ids],
        )
        updated = len(matched_ids)
        conn.commit()
    elif not apply:
        updated = len(matched_ids)

    return SeedReleaseDateStats(
        report_rows=report_rows or len(id_to_date),
        with_date=len(id_to_date),
        matched=len(matched_ids),
        updated=updated,
        missing_in_db=missing_in_db,
        empty_in_report=empty_in_report,
        column_added=column_added,
    )
