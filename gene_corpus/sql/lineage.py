"""Fill annotation_lineage from an Annotrieve flattened taxonomy TSV."""
from __future__ import annotations

import csv
import hashlib
import sqlite3
import sys
from collections.abc import Mapping
from pathlib import Path

_MAX_LINEAGE_DEPTH = 128
TAXONOMY_SHA256_META_KEY = "taxonomy_sha256"


def load_parent_map(taxonomy_tsv: Path) -> dict[int, int | None]:
    """
    Load taxid → parent_taxid from a flattened-tree TSV.

    Expects header columns ``taxid`` and ``parent_taxid`` (Annotrieve export).
    Empty / missing parent becomes None.
    """
    parents: dict[int, int | None] = {}
    with taxonomy_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"empty taxonomy TSV: {taxonomy_tsv}")
        fields = {name.strip().lower(): name for name in reader.fieldnames}
        if "taxid" not in fields or "parent_taxid" not in fields:
            raise ValueError(
                f"taxonomy TSV needs taxid and parent_taxid columns: {taxonomy_tsv}"
            )
        tax_col = fields["taxid"]
        parent_col = fields["parent_taxid"]
        for row in reader:
            raw_tax = (row.get(tax_col) or "").strip()
            if not raw_tax or not raw_tax.isdigit():
                continue
            taxid = int(raw_tax)
            raw_parent = (row.get(parent_col) or "").strip()
            if raw_parent and raw_parent.isdigit():
                parents[taxid] = int(raw_parent)
            else:
                parents[taxid] = None
    return parents


def taxonomy_file_sha256(taxonomy_tsv: Path) -> str:
    """SHA-256 hex digest of the taxonomy TSV file bytes."""
    digest = hashlib.sha256()
    with taxonomy_tsv.open("rb") as handle:
        while True:
            chunk = handle.read(1 << 20)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def get_meta_value(conn: sqlite3.Connection, key: str) -> str | None:
    row = conn.execute(
        "SELECT value FROM meta WHERE key = ?", (key,)
    ).fetchone()
    return None if row is None else str(row[0])


def set_meta_value(conn: sqlite3.Connection, key: str, value: str) -> None:
    conn.execute(
        "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
        (key, value),
    )
    conn.commit()


def lineage_taxids(
    species_taxid: int,
    parents: Mapping[int, int | None],
    *,
    max_depth: int = _MAX_LINEAGE_DEPTH,
) -> list[int]:
    """
    Species taxid plus ancestors, walking ``parents``.

    Always includes ``species_taxid``. Stops on missing parent, self-loop,
    or depth cap.
    """
    out: list[int] = []
    seen: set[int] = set()
    current: int | None = int(species_taxid)
    depth = 0
    while current is not None and current not in seen and depth < max_depth:
        out.append(current)
        seen.add(current)
        if current not in parents:
            break
        current = parents[current]
        depth += 1
    return out


def _lineage_pairs_for_taxid(
    annotation_key: int,
    species_taxid: int,
    parents: Mapping[int, int | None],
    *,
    warn: bool = True,
) -> list[tuple[int, int]]:
    """Build (annotation_key, taxid) rows for one species taxid."""
    taxid_i = int(species_taxid)
    key_i = int(annotation_key)
    if taxid_i not in parents:
        if warn:
            print(
                f"merge warn: taxid {taxid_i} not in taxonomy TSV "
                f"(annotation_key={key_i}); species-only lineage",
                file=sys.stderr,
            )
        return [(key_i, taxid_i)]
    return [(key_i, t) for t in lineage_taxids(taxid_i, parents)]


def fill_lineage_for_annotation(
    conn: sqlite3.Connection,
    annotation_key: int,
    species_taxid: int,
    parents: Mapping[int, int | None],
    *,
    commit: bool = True,
) -> int:
    """
    Insert lineage rows for one annotation_key (no delete of other keys).

    Returns number of rows inserted for this key.
    When ``commit`` is False, leaves the transaction open for the caller.
    """
    batch = _lineage_pairs_for_taxid(
        annotation_key, species_taxid, parents, warn=True
    )
    if batch:
        conn.executemany(
            "INSERT OR IGNORE INTO annotation_lineage "
            "(annotation_key, taxid) VALUES (?, ?)",
            batch,
        )
    if commit:
        conn.commit()
    return len(batch)


def fill_annotation_lineage(
    conn: sqlite3.Connection,
    parents: Mapping[int, int | None],
) -> int:
    """
    Insert lineage rows for every annotation.

    Warns when a species taxid is absent from the parent map (still inserts
    the species row alone). Returns number of lineage rows written.
    """
    conn.execute("DELETE FROM annotation_lineage")
    ann_rows = conn.execute(
        "SELECT annotation_key, taxid FROM annotation ORDER BY annotation_key"
    ).fetchall()
    batch: list[tuple[int, int]] = []
    missing = 0
    for annotation_key, taxid in ann_rows:
        taxid_i = int(taxid)
        if taxid_i not in parents:
            missing += 1
        batch.extend(
            _lineage_pairs_for_taxid(
                int(annotation_key), taxid_i, parents, warn=taxid_i not in parents
            )
        )

    if batch:
        conn.executemany(
            "INSERT OR IGNORE INTO annotation_lineage "
            "(annotation_key, taxid) VALUES (?, ?)",
            batch,
        )
    conn.commit()
    if missing:
        print(
            f"merge lineage: {missing} annotation(s) missing from taxonomy map",
            file=sys.stderr,
            flush=True,
        )
    n = int(conn.execute("SELECT COUNT(*) FROM annotation_lineage").fetchone()[0])
    return n


def refresh_lineage_if_taxonomy_changed(
    conn: sqlite3.Connection,
    taxonomy_tsv: Path,
    parents: Mapping[int, int | None],
) -> tuple[bool, int]:
    """
    Rebuild ``annotation_lineage`` when the TSV content hash changed.

    Stores ``taxonomy_sha256`` in ``meta``. Returns
    ``(refreshed, n_lineage_rows)``.
    """
    digest = taxonomy_file_sha256(taxonomy_tsv)
    previous = get_meta_value(conn, TAXONOMY_SHA256_META_KEY)
    if previous == digest:
        n = int(
            conn.execute("SELECT COUNT(*) FROM annotation_lineage").fetchone()[0]
        )
        return False, n
    n = fill_annotation_lineage(conn, parents)
    set_meta_value(conn, TAXONOMY_SHA256_META_KEY, digest)
    return True, n
