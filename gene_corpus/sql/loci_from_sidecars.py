"""
Build ``loci.sqlite`` from per-annotation sidecars (PLAN-build-tables.md, /loci phase).

One row per (name, gene): the coordinates sit in the row, rows are clustered by name, so
a request reads one short key range per name. ``name_id`` is the id used by the matching
``presence.sqlite`` (same merge order); the build checks, name by name and in lockstep,
that the two agree and refuses otherwise. The output is written next to that presence file
(same ``v<data_version>/`` directory) and carries the same ``data_version``.

    locus(name_id, ord, local_id, seqid_id, start, end, strand, feature_type_id,
          biotype_id, via, gene_id, primary_name)  PRIMARY KEY (name_id, ord, local_id)

``via`` is 0 for a symbol match and 1 for an alias-only match of that gene (a gene that is
both keeps 0). ``gene_id`` is the source gene id, ``primary_name`` the gene's own name in
its original case. Coordinates are 1-based inclusive, as in the GFF.

Algorithm, same two levels as the presence build:
1. Annotations are grouped into batches whose named genes fit a budget (``gene_budget``).
   For each annotation the ``#genes`` section is loaded (a dict, or an on-disk SQLite table
   for very large annotations), then its sorted ``#names`` rows are joined to it; the batch
   is merged into one gzip run.
2. The runs are merged again; each name's rows are sorted by (ord, local_id) and inserted.

Memory: about ``gene_budget`` gene tuples per batch (default 300,000, roughly 100 MiB) plus
one reader per run; very large annotations use a scratch SQLite file under ``tmp_dir``
(on the NFS work area, never local disk).
"""
from __future__ import annotations

import gzip
import heapq
import json
import os
import resource
import shutil
import sqlite3
import sys
import tempfile
import time
from collections.abc import Callable, Iterator
from itertools import groupby
from operator import itemgetter
from pathlib import Path

from gene_corpus.sql.serve_presence import SCHEMA_VERSION, BuildError
from gene_corpus.sql.sidecar import (
    SHARD_NAME,
    SidecarError,
    iter_genes,
    iter_names,
    read_header,
    sidecar_path_for,
)

LOCI_FILE = "loci.sqlite"
INSERT_BATCH = 20_000
DICT_GENES_LIMIT = 150_000       # above this an annotation's genes go to a scratch SQLite table
_PAGE_SIZE = 16384

_SCHEMA = """
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE locus (
    name_id INTEGER NOT NULL,
    ord INTEGER NOT NULL,
    local_id INTEGER NOT NULL,
    seqid_id INTEGER NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    strand INTEGER NOT NULL,
    feature_type_id INTEGER NOT NULL,
    biotype_id INTEGER,
    via INTEGER NOT NULL,
    gene_id TEXT,
    primary_name TEXT,
    PRIMARY KEY (name_id, ord, local_id)
) WITHOUT ROWID;
CREATE TABLE seqid (ord INTEGER NOT NULL, seqid_id INTEGER NOT NULL, seqid TEXT NOT NULL,
                    PRIMARY KEY (ord, seqid_id)) WITHOUT ROWID;
CREATE TABLE dict_feature_type (id INTEGER PRIMARY KEY, name TEXT UNIQUE NOT NULL);
CREATE TABLE dict_biotype (id INTEGER PRIMARY KEY, name TEXT UNIQUE NOT NULL);
CREATE TABLE name_loci (name_id INTEGER PRIMARY KEY, n_loci INTEGER NOT NULL);
"""


def _peak_rss_mib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _none(v: str) -> str | None:
    return v if v != "" else None


class _Dicts:
    """Global small dictionaries (feature types, biotypes), ids assigned on first use."""

    def __init__(self) -> None:
        self.ft: dict[str, int] = {}
        self.bt: dict[str, int] = {}

    def ft_id(self, v: str) -> int:
        return self.ft.setdefault(v, len(self.ft) + 1)

    def bt_id(self, v: str | None) -> int | None:
        return None if v is None else self.bt.setdefault(v, len(self.bt) + 1)


class _GeneStore:
    """local_id -> (seqid_id, start, end, strand, ft_id, bt_id, gene_id, primary_name)."""

    def __init__(self, path: Path, ord_: int, n_genes: int, dicts: _Dicts,
                 out: sqlite3.Connection, tmp_dir: Path) -> None:
        self._dict: dict[int, tuple] | None = None
        self._db: sqlite3.Connection | None = None
        self._scratch: Path | None = None
        seqids: dict[str, int] = {}
        rows: list[tuple] = []
        if n_genes <= DICT_GENES_LIMIT:
            self._dict = {}
        else:
            tmp_dir.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(dir=tmp_dir, suffix=".genes.sqlite")
            os.close(fd)
            self._scratch = Path(name)
            self._db = sqlite3.connect(name)
            self._db.execute("PRAGMA journal_mode=OFF")
            self._db.execute("PRAGMA synchronous=OFF")
            self._db.execute("PRAGMA cache_size=-8192")
            self._db.execute(
                "CREATE TABLE g (local_id INTEGER PRIMARY KEY, seqid_id INTEGER, start INTEGER, "
                "end INTEGER, strand INTEGER, ft INTEGER, bt INTEGER, gene_id TEXT, pname TEXT)")
        batch: list[tuple] = []
        n = 0
        for local_id, gene_id, seqid, start, end, strand, ftype, biotype, pname in iter_genes(path):
            sid = seqids.get(seqid)
            if sid is None:
                sid = seqids[seqid] = len(seqids) + 1
                rows.append((ord_, sid, seqid))
            rec = (sid, start, end, strand, dicts.ft_id(ftype), dicts.bt_id(_none(biotype)),
                   _none(gene_id), _none(pname))
            n += 1
            if self._dict is not None:
                self._dict[local_id] = rec
            else:
                batch.append((local_id, *rec))
                if len(batch) >= 20_000:
                    self._db.executemany("INSERT INTO g VALUES (?,?,?,?,?,?,?,?,?)", batch)
                    batch.clear()
        if batch:
            self._db.executemany("INSERT INTO g VALUES (?,?,?,?,?,?,?,?,?)", batch)
        out.executemany("INSERT INTO seqid VALUES (?,?,?)", rows)
        self.n = n

    def get(self, local_id: int) -> tuple | None:
        if self._dict is not None:
            return self._dict.get(local_id)
        row = self._db.execute(
            "SELECT seqid_id, start, end, strand, ft, bt, gene_id, pname FROM g WHERE local_id=?",
            (local_id,)).fetchone()
        return tuple(row) if row else None

    def close(self) -> None:
        self._dict = None
        if self._db is not None:
            self._db.close()
        if self._scratch is not None:
            self._scratch.unlink(missing_ok=True)


def _annotation_loci(path: Path, ord_: int, store: _GeneStore, counts: dict[int, int]
                     ) -> Iterator[tuple]:
    n = 0
    prev = None
    try:
        for name, kind, local_id in iter_names(path):
            key = (name, kind)
            if prev is not None and key < prev:
                raise BuildError(f"{path}: names not sorted ({prev!r} then {key!r})")
            prev = key
            g = store.get(local_id)
            if g is None:
                raise BuildError(f"{path}: name {name!r} refers to unknown gene {local_id}")
            n += 1
            yield (name, kind, ord_, local_id, *g)
    except SidecarError as exc:
        raise BuildError(f"{path}: {exc}") from exc
    counts[ord_] = n


def _fmt(v) -> str:
    return "" if v is None else str(v)


def _write_run(run_path: Path, merged: Iterator[tuple]) -> int:
    lines = 0
    with gzip.open(run_path, "wt", encoding="utf-8", compresslevel=1, newline="\n") as fh:
        for row in merged:
            fh.write("\t".join(map(_fmt, row)))
            fh.write("\n")
            lines += 1
    return lines


def _read_run(run_path: Path) -> Iterator[tuple]:
    with gzip.open(run_path, "rt", encoding="utf-8", newline="\n") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            yield (f[0], f[1], int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]),
                   int(f[8]), int(f[9]) if f[9] != "" else None, _none(f[10]), _none(f[11]))


_KEY = itemgetter(0, 1, 2, 3)


def build_loci(
    presence_file: Path,
    manifest_path: Path,
    per_gff: Path,
    out_file: Path,
    *,
    tmp_dir: Path,
    gene_budget: int = 300_000,
    verify_sample: int = 30,
    keep_runs: bool = False,
    log: Callable[[str], None] | None = None,
) -> dict:
    if log is None:
        log = lambda m: print(f"[{time.strftime('%H:%M:%S')}] {m}", file=sys.stderr, flush=True)  # noqa: E731
    t0 = time.time()
    pres = sqlite3.connect(f"file:{presence_file}?mode=ro", uri=True)
    pmeta = dict(pres.execute("SELECT key, value FROM meta"))
    if pmeta.get("partial") == "1":
        raise BuildError("presence.sqlite is a partial build; loci need a full one")
    p_ords = {int(r[0]) for r in pres.execute("SELECT ord FROM annotation")}

    manifest = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True)
    rows = manifest.execute(
        "SELECT ord, annotation_id, n_genes_named, n_name_rows, shard_fingerprint "
        "FROM annotation WHERE status='indexed' AND in_catalogue=1 ORDER BY ord").fetchall()
    if {int(r[0]) for r in rows} != p_ords:
        raise BuildError("the manifest's annotation set differs from presence.sqlite "
                         "(rebuild presence first, or restore the manifest)")
    named = [r for r in rows if (r[3] or 0) > 0]
    if len(named) != int(pmeta["n_annotations_with_names"]):
        raise BuildError("annotations with names differ from presence.sqlite")
    log(f"annotations with names: {len(named):,}; presence data_version {pmeta['data_version']}")

    tmp_dir.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix="loci-runs-", dir=tmp_dir))
    if out_file.exists():
        out_file.unlink()
    out = sqlite3.connect(str(out_file))
    out.execute(f"PRAGMA page_size={_PAGE_SIZE}")
    out.execute("PRAGMA journal_mode=OFF")
    out.execute("PRAGMA synchronous=OFF")
    out.execute("PRAGMA cache_size=-65536")
    out.executescript(_SCHEMA)

    dicts = _Dicts()
    counts: dict[int, int] = {}
    runs: list[Path] = []
    paths: dict[int, Path] = {}
    try:
        # ---- level 1 -----------------------------------------------------
        i = 0
        while i < len(named):
            budget = 0
            batch = []
            while i < len(named) and (not batch or budget + (named[i][2] or 0) <= gene_budget):
                batch.append(named[i])
                budget += named[i][2] or 0
                i += 1
            stores: list[_GeneStore] = []
            streams = []
            try:
                for ord_, aid, n_genes, n_rows, _fp in batch:
                    path = sidecar_path_for(per_gff / aid / SHARD_NAME)
                    header = read_header(path)
                    if int(header.get("n_name_rows", -1)) != n_rows:
                        raise BuildError(f"{aid}: sidecar differs from the manifest (re-run `sidecars`)")
                    paths[ord_] = path
                    try:
                        store = _GeneStore(path, ord_, int(n_genes or 0), dicts, out, tmp_dir)
                    except SidecarError as exc:
                        raise BuildError(f"{aid}: {exc}") from exc
                    stores.append(store)
                    streams.append(_annotation_loci(path, ord_, store, counts))
                run_path = run_dir / f"run-{len(runs):05d}.tsv.gz"
                _write_run(run_path, heapq.merge(*streams, key=_KEY))
                runs.append(run_path)
            finally:
                for s in stores:
                    s.close()
            out.commit()
            log(f"level 1: {i:,}/{len(named):,} annotations, run {len(runs)}, "
                f"peak_rss={_peak_rss_mib():.0f}MiB")
        for r in named:
            if counts.get(r[0]) != r[3]:
                raise BuildError(f"{r[1]}: joined {counts.get(r[0])} name rows, manifest says {r[3]}")

        # ---- level 2 -----------------------------------------------------
        merged = heapq.merge(*[_read_run(p) for p in runs], key=_KEY)
        names_cur = pres.execute("SELECT name, n_sym, n_alias_only FROM name_dict ORDER BY name_id")
        n_names = n_loci = 0
        locus_batch: list[tuple] = []
        nl_batch: list[tuple] = []
        last_log = time.time()

        def flush() -> None:
            if locus_batch:
                out.executemany("INSERT INTO locus VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", locus_batch)
                out.executemany("INSERT INTO name_loci VALUES (?,?)", nl_batch)
                locus_batch.clear()
                nl_batch.clear()

        for name, grp in groupby(merged, key=itemgetter(0)):
            n_names += 1
            want = names_cur.fetchone()
            if want is None or want[0] != name:
                raise BuildError(f"name #{n_names}: loci has {name!r}, presence has "
                                 f"{None if want is None else want[0]!r}")
            by_gene: dict[tuple[int, int], list] = {}
            sym_ords: set[int] = set()
            alias_ords: set[int] = set()
            for _n, kind, ord_, local_id, sid, start, end, strand, ft, bt, gid, pname in grp:
                key = (ord_, local_id)
                if kind == "s":
                    sym_ords.add(ord_)
                    by_gene[key] = [sid, start, end, strand, ft, bt, 0, gid, pname]
                else:
                    alias_ords.add(ord_)
                    if key not in by_gene:
                        by_gene[key] = [sid, start, end, strand, ft, bt, 1, gid, pname]
            if len(sym_ords) != want[1] or len(alias_ords - sym_ords) != want[2]:
                raise BuildError(f"{name!r}: loci counts sym={len(sym_ords)} "
                                 f"alias_only={len(alias_ords - sym_ords)} differ from presence "
                                 f"({want[1]}, {want[2]})")
            for (ord_, local_id) in sorted(by_gene):
                locus_batch.append((n_names, ord_, local_id, *by_gene[(ord_, local_id)]))
            n_loci += len(by_gene)
            nl_batch.append((n_names, len(by_gene)))
            if len(locus_batch) >= INSERT_BATCH:
                flush()
            if time.time() - last_log > 30:
                last_log = time.time()
                log(f"level 2: {n_names:,} names, {n_loci:,} loci (at {name!r}) "
                    f"peak_rss={_peak_rss_mib():.0f}MiB")
        flush()
        if names_cur.fetchone() is not None or n_names != int(pmeta["n_names"]):
            raise BuildError(f"name counts differ: loci {n_names:,}, presence {pmeta['n_names']}")

        out.executemany("INSERT INTO dict_feature_type VALUES (?,?)",
                        [(v, k) for k, v in dicts.ft.items()])
        out.executemany("INSERT INTO dict_biotype VALUES (?,?)",
                        [(v, k) for k, v in dicts.bt.items()])
        meta = {
            "schema_version": str(SCHEMA_VERSION),
            "data_version": pmeta["data_version"],
            "built_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "n_names": str(n_names),
            "n_loci": str(n_loci),
            "n_annotations_with_names": str(len(named)),
        }
        out.executemany("INSERT INTO meta VALUES (?,?)", list(meta.items()))
        out.commit()

        checked = _verify(out, pres, named, paths, dicts, verify_sample, log)
    finally:
        out.close()
        pres.close()
        manifest.close()
        if not keep_runs:
            shutil.rmtree(run_dir, ignore_errors=True)

    stats = {
        "data_version": pmeta["data_version"], "n_names": n_names, "n_loci": n_loci,
        "runs": len(runs), "output_bytes": out_file.stat().st_size,
        "verified_annotations": checked, "peak_rss_mib": round(_peak_rss_mib(), 1),
        "seconds": round(time.time() - t0, 1),
    }
    log("loci build done: " + json.dumps(stats))
    return stats


def _verify(out: sqlite3.Connection, pres: sqlite3.Connection, named: list[tuple],
            paths: dict[int, Path], dicts: _Dicts, sample: int, log) -> int:
    """For evenly spaced annotations, every name row must be in loci.sqlite with its coordinates."""
    if sample <= 0 or not named:
        return 0
    inv_ft = {v: k for k, v in dicts.ft.items()}
    step = max(1, len(named) // sample)
    picks = named[::step]
    for ord_, aid, *_ in picks:
        genes = {g[0]: g for g in iter_genes(paths[ord_])}
        seqids = dict(out.execute("SELECT seqid_id, seqid FROM seqid WHERE ord=?", (ord_,)))
        for name, kind, local_id in iter_names(paths[ord_]):
            nid = pres.execute("SELECT name_id FROM name_dict WHERE name=?", (name,)).fetchone()
            row = None if nid is None else out.execute(
                "SELECT seqid_id, start, end, strand, feature_type_id, via, gene_id FROM locus "
                "WHERE name_id=? AND ord=? AND local_id=?", (nid[0], ord_, local_id)).fetchone()
            if row is None:
                raise BuildError(f"verification: {aid} gene {local_id} ({name}) missing from loci")
            g = genes[local_id]
            want = (g[2], g[3], g[4], g[5], g[6], g[1] or None)
            got = (seqids[row[0]], row[1], row[2], row[3], inv_ft[row[4]], row[6])
            if got != want or (kind == "s" and row[5] != 0):
                raise BuildError(f"verification: {aid} gene {local_id} ({name}): {got} != {want}")
    log(f"verified {len(picks)} annotations")
    return len(picks)
