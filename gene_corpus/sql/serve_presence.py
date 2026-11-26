"""
Shared pieces of the served-index builds: the ``presence.sqlite`` schema, ``BuildError`` and the
atomic ``publish_version`` switch. ``presence_from_sidecars`` and ``loci_from_sidecars`` build the
files; this module holds what they share.

``presence.sqlite`` tables:

    meta         key/value (data_version, schema_version, counts, partial flag)
    annotation   one row per annotation, keyed by ``ord``
    lineage      (taxid, ord)
    name_dict    casefolded gene name -> name_id, n_sym, n_alias_only
    posting      name_id -> sym / alias posting blobs (see ``postings.py``)

Semantics: a *match* is a symbol hit. ``n_alias_only`` is the number of annotations where the
name occurs only as an alias, which the API reports as a suggestion.
"""
from __future__ import annotations

import json
import os
import sqlite3
from pathlib import Path

SCHEMA_VERSION = 1
PAGE_SIZE = 16384

# Must match api/helpers/source_codes.py (1/2/3, else 0).
_SOURCE_CODES = {"ensembl": 1, "genbank": 2, "refseq": 3}

_OUT_SCHEMA = """
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);

CREATE TABLE annotation (
    ord INTEGER PRIMARY KEY,
    annotation_id TEXT NOT NULL UNIQUE,
    taxid INTEGER NOT NULL,
    organism_name TEXT,
    assembly_accession TEXT,
    assembly_name TEXT,
    source INTEGER NOT NULL,
    release_date TEXT,
    gff_path TEXT
);

CREATE TABLE lineage (
    taxid INTEGER NOT NULL,
    ord INTEGER NOT NULL,
    PRIMARY KEY (taxid, ord)
) WITHOUT ROWID;

-- Annotations whose GFF carries at least one gene name (the ones a name search can reach).
CREATE TABLE named (
    ord INTEGER PRIMARY KEY
);

CREATE TABLE name_dict (
    name TEXT PRIMARY KEY,
    name_id INTEGER NOT NULL UNIQUE,
    n_sym INTEGER NOT NULL,
    n_alias_only INTEGER NOT NULL
) WITHOUT ROWID;

CREATE TABLE posting (
    name_id INTEGER PRIMARY KEY,
    sym BLOB,
    alias BLOB
);
"""


class BuildError(Exception):
    """Build or validation failure; the output must not be published."""


def _source_code(value: str | None) -> int:
    return _SOURCE_CODES.get(str(value or "").strip().casefold(), 0)


def _open_output(path: Path) -> sqlite3.Connection:
    if path.exists():
        path.unlink()
    conn = sqlite3.connect(str(path))
    conn.execute(f"PRAGMA page_size={PAGE_SIZE}")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA cache_size=-65536")
    conn.executescript(_OUT_SCHEMA)
    return conn


def _alias_only_count(sym: list[int], alias: list[int]) -> int:
    if not alias:
        return 0
    if not sym:
        return len(alias)
    return len(set(alias).difference(sym))


def publish_version(out_dir: Path, built_file: Path, stats: dict) -> Path:
    """
    Move ``built_file`` into ``out_dir/v<data_version>/presence.sqlite`` and switch
    the ``out_dir/current`` symlink to that directory in one atomic step.
    Older version directories are left in place (the API may still have them open).
    """
    version_dir = out_dir / f"v{stats['data_version']}"
    version_dir.mkdir(parents=True, exist_ok=True)
    final = version_dir / "presence.sqlite"
    os.replace(built_file, final)
    (version_dir / "build_report.json").write_text(json.dumps(stats, indent=2))
    tmp_link = out_dir / ".current.tmp"
    if tmp_link.is_symlink() or tmp_link.exists():
        tmp_link.unlink()
    tmp_link.symlink_to(version_dir.name)  # relative target
    os.replace(tmp_link, out_dir / "current")
    return final
