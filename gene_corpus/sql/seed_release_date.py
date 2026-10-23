"""
One-time backfill: ALTER annotation for release_date + assembly_name, then
UPDATE from the Annotrieve annotations report TSV
(selected_fields=release_date; assembly_name is a default report column).
"""
from __future__ import annotations

import csv
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from gene_corpus.sql.schema import SCHEMA_VERSION


RELEASE_DATE_COLUMN = "release_date"
ASSEMBLY_NAME_COLUMN = "assembly_name"
SEED_COLUMNS: tuple[str, ...] = (RELEASE_DATE_COLUMN, ASSEMBLY_NAME_COLUMN)

# Extended report column; defaults (incl. annotation_id, assembly_name) always present.
REPORT_SEED_FIELDS = "release_date"
# Back-compat alias for callers/tests that imported the old name.
REPORT_RELEASE_DATE_FIELDS = REPORT_SEED_FIELDS


@dataclass(frozen=True)
class AnnotationSeedValues:
    release_date: str | None
    assembly_name: str | None


@dataclass(frozen=True)
class SeedReleaseDateStats:
    report_rows: int
    with_date: int
    with_assembly_name: int
    matched: int
    updated: int
    missing_in_db: int
    empty_release_date: int
    empty_assembly_name: int
    columns_added: tuple[str, ...]


def annotation_columns(conn: sqlite3.Connection) -> set[str]:
    return {
        str(row[1])
        for row in conn.execute("PRAGMA table_info(annotation)")
    }


def annotation_has_release_date_column(conn: sqlite3.Connection) -> bool:
    return RELEASE_DATE_COLUMN in annotation_columns(conn)


def annotation_has_assembly_name_column(conn: sqlite3.Connection) -> bool:
    return ASSEMBLY_NAME_COLUMN in annotation_columns(conn)


def ensure_seed_columns(conn: sqlite3.Connection) -> tuple[str, ...]:
    """
    ADD COLUMN for each missing seed column. Returns names that were added.
    """
    existing = annotation_columns(conn)
    added: list[str] = []
    for col in SEED_COLUMNS:
        if col in existing:
            continue
        conn.execute(f"ALTER TABLE annotation ADD COLUMN {col} TEXT")
        added.append(col)
    if added:
        conn.execute(
            "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
            ("schema_version", str(SCHEMA_VERSION)),
        )
        conn.commit()
    return tuple(added)


def ensure_release_date_column(conn: sqlite3.Connection) -> bool:
    """Back-compat: ensure release_date only. Prefer ensure_seed_columns."""
    return RELEASE_DATE_COLUMN in ensure_seed_columns(conn)


def parse_release_date_tsv(
    report: Path,
) -> tuple[dict[str, AnnotationSeedValues], int, int, int]:
    """
    Parse report TSV → (annotation_id → seed values, total_rows,
    empty_release_date, empty_assembly_name).

    Rows without annotation_id are skipped. Empty fields become None.
    """
    out: dict[str, AnnotationSeedValues] = {}
    total = 0
    empty_date = 0
    empty_name = 0
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
        if ASSEMBLY_NAME_COLUMN not in names:
            raise ValueError(
                f"TSV missing {ASSEMBLY_NAME_COLUMN} column: {report} "
                f"(columns={names!r})"
            )
        for row in reader:
            ann_id = (row.get("annotation_id") or "").strip()
            if not ann_id:
                continue
            total += 1
            date_raw = (row.get(RELEASE_DATE_COLUMN) or "").strip()
            name_raw = (row.get(ASSEMBLY_NAME_COLUMN) or "").strip()
            if not date_raw:
                empty_date += 1
            if not name_raw:
                empty_name += 1
            out[ann_id] = AnnotationSeedValues(
                release_date=date_raw or None,
                assembly_name=name_raw or None,
            )
    return out, total, empty_date, empty_name


def seed_release_dates(
    conn: sqlite3.Connection,
    id_to_seed: Mapping[str, AnnotationSeedValues] | Mapping[str, str],
    *,
    apply: bool,
    report_rows: int = 0,
    empty_in_report: int = 0,
    empty_assembly_name: int = 0,
) -> SeedReleaseDateStats:
    """
    Ensure columns exist (on --apply) and UPDATE release_date + assembly_name
    by annotation_id. Dry-run counts matches without writing.

    ``id_to_seed`` may be ``AnnotationSeedValues`` (preferred) or a legacy
    ``annotation_id → release_date`` string map.
    """
    seeds = _normalize_seed_map(id_to_seed)

    columns_added: tuple[str, ...] = ()
    if apply:
        columns_added = ensure_seed_columns(conn)
    else:
        existing = annotation_columns(conn)
        columns_added = tuple(c for c in SEED_COLUMNS if c not in existing)

    db_ids = {
        str(r[0])
        for r in conn.execute("SELECT annotation_id FROM annotation")
    }
    matched_ids = sorted(set(seeds) & db_ids)
    missing_in_db = len(set(seeds) - db_ids)

    updated = 0
    if apply and matched_ids:
        conn.executemany(
            f"""
            UPDATE annotation
            SET {RELEASE_DATE_COLUMN} = ?,
                {ASSEMBLY_NAME_COLUMN} = ?
            WHERE annotation_id = ?
            """,
            [
                (
                    seeds[aid].release_date,
                    seeds[aid].assembly_name,
                    aid,
                )
                for aid in matched_ids
            ],
        )
        updated = len(matched_ids)
        conn.commit()
    elif not apply:
        updated = len(matched_ids)

    with_date = sum(1 for s in seeds.values() if s.release_date)
    with_name = sum(1 for s in seeds.values() if s.assembly_name)

    return SeedReleaseDateStats(
        report_rows=report_rows or len(seeds),
        with_date=with_date,
        with_assembly_name=with_name,
        matched=len(matched_ids),
        updated=updated,
        missing_in_db=missing_in_db,
        empty_release_date=empty_in_report,
        empty_assembly_name=empty_assembly_name,
        columns_added=columns_added,
    )


def _normalize_seed_map(
    id_to_seed: Mapping[str, AnnotationSeedValues] | Mapping[str, str],
) -> dict[str, AnnotationSeedValues]:
    out: dict[str, AnnotationSeedValues] = {}
    for aid, value in id_to_seed.items():
        if isinstance(value, AnnotationSeedValues):
            out[aid] = value
        else:
            raw = str(value).strip() if value is not None else ""
            out[aid] = AnnotationSeedValues(
                release_date=raw or None,
                assembly_name=None,
            )
    return out
