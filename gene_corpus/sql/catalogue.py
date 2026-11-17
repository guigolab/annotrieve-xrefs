"""
Catalogue loader: Annotrieve's annotation metadata into the manifest (decision D8).

Source: ``GET {api_base}/annotations/report?selected_fields=release_date,taxon_lineage``,
a TSV streamed line by line (default columns: annotation_id, assembly_accession,
assembly_name, organism_name, taxid, database, provider, source_url, bgzip_path,
csi_path; extended: release_date, taxon_lineage as ``;``-joined taxids).

Memory is flat: rows are read one at a time and inserted in small batches. The loader
only stores and reports; it never removes anything (reconciliation is D6, postponed).
"""
from __future__ import annotations

import csv
import sqlite3
import urllib.request
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field

DEFAULT_API_BASE = "https://genome.crg.es/annotrieve/api/v0"
REPORT_FIELDS = "release_date,taxon_lineage"
REQUIRED_COLUMNS = (
    "annotation_id", "assembly_accession", "assembly_name", "organism_name", "taxid",
    "database", "provider", "release_date", "taxon_lineage",
)
_BATCH = 2000


class CatalogueError(Exception):
    """The report is unreachable, malformed or lacks required columns."""


@dataclass
class CatalogueStats:
    rows: int = 0
    skipped_blank_id: int = 0
    duplicate_ids: int = 0
    matched_annotations: int = 0          # manifest rows found in the catalogue
    shard_only: list[str] = field(default_factory=list)       # in manifest, not in catalogue
    catalogue_only: int = 0               # in catalogue, no shard in the manifest
    database_mismatch: int = 0            # shard source_database != catalogue database
    lineage_rows: int = 0


def report_url(api_base: str = DEFAULT_API_BASE) -> str:
    return f"{api_base.rstrip('/')}/annotations/report?selected_fields={REPORT_FIELDS}"


def fetch_report_lines(url: str, *, timeout: float = 120.0) -> Iterator[str]:
    """Stream the report as decoded lines (no full-body buffering)."""
    req = urllib.request.Request(url, headers={"Accept": "text/tab-separated-values",
                                               "User-Agent": "annotrieve-xrefs-catalogue"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            for raw in resp:
                yield raw.decode("utf-8").rstrip("\r\n")
    except Exception as exc:
        raise CatalogueError(f"cannot fetch {url}: {exc}") from exc


def _lineage_taxids(taxid: int | None, lineage: str) -> set[int]:
    out: set[int] = set()
    for part in lineage.split(";"):
        part = part.strip()
        if part.isdigit():
            out.add(int(part))
    if taxid is not None:
        out.add(taxid)           # an annotation belongs to its own taxon as well
    return out


def _int_or_none(v: str) -> int | None:
    v = (v or "").strip()
    return int(v) if v.isdigit() else None


def load_catalogue(conn: sqlite3.Connection, lines: Iterable[str]) -> CatalogueStats:
    """
    Replace the ``catalogue`` table with the report rows, then apply it to the manifest:
    fill ``assembly_name`` / ``release_date`` (and ``taxid`` / ``organism_name`` when the shard
    had none), set ``in_catalogue`` and rebuild ``lineage`` for registered annotations.
    """
    reader = csv.reader(iter(lines), delimiter="\t", quoting=csv.QUOTE_NONE)
    try:
        header = next(reader)
    except StopIteration:
        raise CatalogueError("empty report") from None
    index = {name: i for i, name in enumerate(header)}
    missing = [c for c in REQUIRED_COLUMNS if c not in index]
    if missing:
        raise CatalogueError(f"report lacks columns {missing}; got {header}")

    def col(row: list[str], name: str) -> str:
        i = index[name]
        return row[i] if i < len(row) else ""

    try:
        return _load(conn, reader, index, col)
    except CatalogueError:
        conn.rollback()           # keep the previous catalogue on a bad report
        raise


def _load(conn, reader, index, col) -> CatalogueStats:
    stats = CatalogueStats()
    conn.execute("DELETE FROM catalogue")
    batch: list[tuple] = []
    seen: set[str] = set()
    for row in reader:
        if not row:
            continue
        stats.rows += 1
        aid = col(row, "annotation_id").strip()
        if not aid:
            stats.skipped_blank_id += 1
            continue
        if aid in seen:
            stats.duplicate_ids += 1
            continue
        seen.add(aid)
        batch.append((
            aid, col(row, "assembly_accession"), col(row, "assembly_name"),
            col(row, "organism_name"), _int_or_none(col(row, "taxid")),
            col(row, "database"), col(row, "provider"),
            col(row, "source_url") if "source_url" in index else "",
            col(row, "bgzip_path") if "bgzip_path" in index else "",
            col(row, "release_date"), col(row, "taxon_lineage"),
        ))
        if len(batch) >= _BATCH:
            conn.executemany("INSERT INTO catalogue VALUES (?,?,?,?,?,?,?,?,?,?,?)", batch)
            batch.clear()
    if batch:
        conn.executemany("INSERT INTO catalogue VALUES (?,?,?,?,?,?,?,?,?,?,?)", batch)
    del seen
    if stats.rows == 0:
        raise CatalogueError("report has a header but no rows")

    apply_catalogue(conn, stats)
    conn.commit()
    return stats


def apply_catalogue(conn: sqlite3.Connection, stats: CatalogueStats | None = None) -> CatalogueStats:
    """
    (Re)apply the stored ``catalogue`` table to the manifest. Idempotent; run it again
    after registering new shards so they pick up their metadata and lineage.
    Lineage rows are written for every registered annotation present in the catalogue,
    whatever its status; the presence build filters on status.
    """
    stats = stats or CatalogueStats()
    conn.execute("UPDATE annotation SET in_catalogue = 0")
    conn.execute(
        "UPDATE annotation SET in_catalogue = 1 WHERE annotation_id IN "
        "(SELECT annotation_id FROM catalogue)"
    )
    conn.execute(
        """
        UPDATE annotation SET
            assembly_name = (SELECT c.assembly_name FROM catalogue c WHERE c.annotation_id = annotation.annotation_id),
            release_date  = (SELECT c.release_date  FROM catalogue c WHERE c.annotation_id = annotation.annotation_id),
            taxid = COALESCE(taxid, (SELECT c.taxid FROM catalogue c WHERE c.annotation_id = annotation.annotation_id)),
            organism_name = COALESCE(NULLIF(organism_name, ''),
                            (SELECT c.organism_name FROM catalogue c WHERE c.annotation_id = annotation.annotation_id))
        WHERE in_catalogue = 1
        """
    )
    stats.matched_annotations = conn.execute(
        "SELECT COUNT(*) FROM annotation WHERE in_catalogue = 1").fetchone()[0]
    stats.shard_only = [r[0] for r in conn.execute(
        "SELECT annotation_id FROM annotation WHERE in_catalogue = 0 ORDER BY ord LIMIT 1000")]
    stats.catalogue_only = conn.execute(
        "SELECT COUNT(*) FROM catalogue c WHERE NOT EXISTS "
        "(SELECT 1 FROM annotation a WHERE a.annotation_id = c.annotation_id)").fetchone()[0]
    stats.database_mismatch = conn.execute(
        "SELECT COUNT(*) FROM annotation a JOIN catalogue c USING (annotation_id) "
        "WHERE lower(trim(a.source_database)) != lower(trim(c.database))").fetchone()[0]

    conn.execute("DELETE FROM lineage")
    batch: list[tuple[int, int]] = []
    cur = conn.execute(
        "SELECT a.ord, a.taxid, c.taxon_lineage FROM annotation a "
        "JOIN catalogue c USING (annotation_id)"
    )
    for ord_, taxid, lineage in cur:
        for t in _lineage_taxids(taxid, lineage or ""):
            batch.append((t, ord_))
        if len(batch) >= 20000:
            conn.executemany("INSERT OR IGNORE INTO lineage VALUES (?, ?)", batch)
            stats.lineage_rows += len(batch)
            batch.clear()
    if batch:
        conn.executemany("INSERT OR IGNORE INTO lineage VALUES (?, ?)", batch)
        stats.lineage_rows += len(batch)
    return stats
