"""
Build ``presence.sqlite``: the compact, read-only name index served by the API.

Source of truth for this build is the global ``gene_corpus.sqlite`` (tables
``annotation``, ``annotation_lineage`` and ``xref_file``). Output tables:

    meta         key/value (data_version, schema_version, counts, partial flag)
    annotation   one row per annotation, keyed by ``ord`` (= annotation_key)
    lineage      (taxid, ord) copy of annotation_lineage
    name_dict    casefolded gene name -> name_id, n_sym, n_alias_only
    posting      name_id -> sym / alias posting blobs (see ``postings.py``)

Semantics (see API-DESIGN.md): a *match* is a symbol hit. ``n_alias_only`` is
the number of annotations where the name occurs only as an alias, which the API
reports as a suggestion. Postings are filtered to annotations present in
``annotation`` (the old /hits path ignored orphan keys the same way).

Memory is constant: both ``xref_file`` streams are read once, in primary-key
order, and merged by accession.
"""
from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import sys
import time
from collections.abc import Callable, Iterator
from itertools import groupby
from operator import itemgetter
from pathlib import Path

from gene_corpus.sql.postings import decode_postings, encode_postings

SCHEMA_VERSION = 1
PAGE_SIZE = 16384
INSERT_BATCH = 20_000

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


def open_source(path: Path, cache_mib: int = 64) -> sqlite3.Connection:
    conn = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
    conn.execute(f"PRAGMA cache_size=-{cache_mib * 1024}")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute("PRAGMA mmap_size=0")
    conn.execute("PRAGMA query_only=ON")
    return conn


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


def _name_groups(
    src: sqlite3.Connection,
    namespace: str,
    lo: str | None,
    hi: str | None,
    valid: set[int],
    counters: dict,
) -> Iterator[tuple[str, list[int]]]:
    """Yield ``(accession, sorted valid ordinals)`` in accession order."""
    if lo is None:
        cur = src.execute(
            "SELECT accession, annotation_key FROM xref_file WHERE namespace = ? "
            "ORDER BY accession, annotation_key",
            (namespace,),
        )
    else:
        cur = src.execute(
            "SELECT accession, annotation_key FROM xref_file WHERE namespace = ? "
            "AND accession >= ? AND accession < ? ORDER BY accession, annotation_key",
            (namespace, lo, hi),
        )
    prev: str | None = None
    for accession, rows in groupby(cur, key=itemgetter(0)):
        if prev is not None and not accession > prev:
            raise BuildError(
                f"xref_file stream for {namespace!r} is not in accession order "
                f"({prev!r} then {accession!r}); refusing to build"
            )
        prev = accession
        ords = []
        for _acc, key in rows:
            counters["rows_" + namespace] += 1
            if key in valid:
                ords.append(key)
            else:
                counters["orphan_rows"] += 1
        yield accession, ords


def _merge_names(
    sym_stream: Iterator[tuple[str, list[int]]],
    alias_stream: Iterator[tuple[str, list[int]]],
) -> Iterator[tuple[str, list[int], list[int]]]:
    """Merge two accession-ordered streams into ``(name, sym_ords, alias_ords)``."""
    end = (None, None)
    s = next(sym_stream, end)
    a = next(alias_stream, end)
    while s[0] is not None or a[0] is not None:
        if a[0] is None or (s[0] is not None and s[0] < a[0]):
            yield s[0], s[1], []
            s = next(sym_stream, end)
        elif s[0] is None or a[0] < s[0]:
            yield a[0], [], a[1]
            a = next(alias_stream, end)
        else:
            yield s[0], s[1], a[1]
            s = next(sym_stream, end)
            a = next(alias_stream, end)


def _alias_only_count(sym: list[int], alias: list[int]) -> int:
    if not alias:
        return 0
    if not sym:
        return len(alias)
    return len(set(alias).difference(sym))


def build_presence(
    src_db: Path,
    out_file: Path,
    *,
    name_prefix: str | None = None,
    limit_names: int | None = None,
    verify_sample: int = 200,
    log: Callable[[str], None] | None = None,
) -> dict:
    """
    Build ``out_file`` from ``src_db`` and return a stats dict.

    ``name_prefix`` and ``limit_names`` produce a *partial* build (for trials);
    it is flagged in ``meta`` and must not be published as current.
    Raises ``BuildError`` when validation fails.
    """
    if log is None:
        log = lambda m: print(m, file=sys.stderr, flush=True)  # noqa: E731
    t0 = time.time()
    partial = bool(name_prefix or limit_names)
    lo = hi = None
    if name_prefix:
        lo = name_prefix.casefold()
        hi = lo + "￿"

    src = open_source(src_db)
    out = _open_output(out_file)
    try:
        # ---- annotations -------------------------------------------------
        ann_rows = src.execute(
            "SELECT annotation_key, annotation_id, taxid, organism_name, assembly_accession, "
            "assembly_name, source_database, release_date, gff_path "
            "FROM annotation ORDER BY annotation_key"
        ).fetchall()
        if not ann_rows:
            raise BuildError("source has no annotations")
        valid = {int(r[0]) for r in ann_rows}
        max_ord = max(valid)
        out.executemany(
            "INSERT INTO annotation (ord, annotation_id, taxid, organism_name, "
            "assembly_accession, assembly_name, source, release_date, gff_path) "
            "VALUES (?,?,?,?,?,?,?,?,?)",
            [(int(r[0]), r[1], int(r[2]), r[3], r[4], r[5], _source_code(r[6]), r[7], r[8])
             for r in ann_rows],
        )
        lineage = [
            (int(t), int(k))
            for t, k in src.execute("SELECT taxid, annotation_key FROM annotation_lineage")
            if int(k) in valid
        ]
        out.executemany("INSERT OR IGNORE INTO lineage (taxid, ord) VALUES (?, ?)", lineage)
        log(f"annotations: {len(ann_rows)}, lineage rows: {len(lineage)}, max ord: {max_ord}")

        # ---- names -------------------------------------------------------
        counters = {"rows_symbol": 0, "rows_alias": 0, "orphan_rows": 0}
        n_names = n_alias_only_names = 0
        blob_bytes = varint_blobs = bitset_blobs = 0
        dict_batch: list[tuple] = []
        post_batch: list[tuple] = []
        last_log = time.time()

        def flush() -> None:
            if dict_batch:
                out.executemany(
                    "INSERT INTO name_dict (name, name_id, n_sym, n_alias_only) VALUES (?,?,?,?)",
                    dict_batch)
                out.executemany(
                    "INSERT INTO posting (name_id, sym, alias) VALUES (?,?,?)", post_batch)
                dict_batch.clear()
                post_batch.clear()

        merged = _merge_names(
            _name_groups(src, "symbol", lo, hi, valid, counters),
            _name_groups(src, "alias", lo, hi, valid, counters),
        )
        for name, sym, alias in merged:
            if not sym and not alias:
                continue  # only orphan annotations
            n_names += 1
            sym_blob = encode_postings(sym, max_ord) if sym else None
            alias_blob = encode_postings(alias, max_ord) if alias else None
            only = _alias_only_count(sym, alias)
            if not sym:
                n_alias_only_names += 1
            for blob in (sym_blob, alias_blob):
                if blob is not None:
                    blob_bytes += len(blob)
                    if blob[0] == 0x01:
                        varint_blobs += 1
                    else:
                        bitset_blobs += 1
            dict_batch.append((name, n_names, len(sym), only))
            post_batch.append((n_names, sym_blob, alias_blob))
            if len(dict_batch) >= INSERT_BATCH:
                flush()
            if time.time() - last_log > 30:
                last_log = time.time()
                log(f"names: {n_names:,} written (at {name!r})")
            if limit_names and n_names >= limit_names:
                log(f"stopping early at --limit-names={limit_names}")
                break
        flush()
        out.commit()

        # ---- meta --------------------------------------------------------
        st = src_db.stat()
        version_src = (
            f"{st.st_mtime_ns}:{st.st_size}:{len(ann_rows)}:{n_names}:"
            f"{counters['rows_symbol']}:{counters['rows_alias']}:{lo}:{limit_names}"
        )
        data_version = hashlib.sha256(version_src.encode()).hexdigest()[:16]
        meta = {
            "schema_version": str(SCHEMA_VERSION),
            "data_version": data_version,
            "built_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "source_db": str(src_db),
            "source_mtime_ns": str(st.st_mtime_ns),
            "n_annotations": str(len(ann_rows)),
            "n_names": str(n_names),
            "max_ord": str(max_ord),
            "partial": "1" if partial else "0",
        }
        out.executemany("INSERT INTO meta (key, value) VALUES (?, ?)", list(meta.items()))
        out.commit()

        # ---- validation --------------------------------------------------
        checked = _verify(src, out, valid, max_ord, verify_sample, lo, hi, log)
    finally:
        out.close()
        src.close()

    stats = {
        "data_version": data_version,
        "partial": partial,
        "n_annotations": len(ann_rows),
        "n_names": n_names,
        "names_alias_only": n_alias_only_names,
        "rows_symbol": counters["rows_symbol"],
        "rows_alias": counters["rows_alias"],
        "orphan_rows_dropped": counters["orphan_rows"],
        "varint_blobs": varint_blobs,
        "bitset_blobs": bitset_blobs,
        "blob_bytes": blob_bytes,
        "output_bytes": out_file.stat().st_size,
        "verified_names": checked,
        "seconds": round(time.time() - t0, 1),
    }
    log("build done: " + json.dumps(stats))
    return stats


def _verify(src, out, valid, max_ord, sample, lo, hi, log) -> int:
    """Recompute postings for evenly spaced names straight from xref_file and compare."""
    n_names = out.execute("SELECT COUNT(*) FROM name_dict").fetchone()[0]
    if sample <= 0 or n_names == 0:
        return 0
    step = max(1, n_names // sample)
    picks = out.execute(
        "SELECT name, name_id, n_sym, n_alias_only FROM name_dict "
        "WHERE (name_id - 1) % ? = 0 OR name_id = ? ORDER BY name_id",
        (step, n_names),
    ).fetchall()
    for name, name_id, n_sym, n_alias_only in picks:
        want_sym = [int(k) for (k,) in src.execute(
            "SELECT annotation_key FROM xref_file WHERE namespace='symbol' AND accession=? "
            "ORDER BY annotation_key", (name,)) if int(k) in valid]
        want_alias = [int(k) for (k,) in src.execute(
            "SELECT annotation_key FROM xref_file WHERE namespace='alias' AND accession=? "
            "ORDER BY annotation_key", (name,)) if int(k) in valid]
        sym_blob, alias_blob = out.execute(
            "SELECT sym, alias FROM posting WHERE name_id = ?", (name_id,)).fetchone()
        got_sym = decode_postings(sym_blob) if sym_blob else []
        got_alias = decode_postings(alias_blob) if alias_blob else []
        if got_sym != want_sym or got_alias != want_alias:
            raise BuildError(f"verification failed for name {name!r}: postings differ from xref_file")
        if n_sym != len(want_sym) or n_alias_only != _alias_only_count(want_sym, want_alias):
            raise BuildError(f"verification failed for name {name!r}: counts differ")
    # global cross-checks
    total_sym = out.execute("SELECT COALESCE(SUM(n_sym), 0) FROM name_dict").fetchone()[0]
    ann = out.execute("SELECT COUNT(*) FROM annotation").fetchone()[0]
    if ann != len(valid):
        raise BuildError("annotation count mismatch")
    log(f"verified {len(picks)} names; sum(n_sym)={total_sym:,}")
    return len(picks)


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
