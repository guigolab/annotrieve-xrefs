"""
Build ``presence.sqlite`` from per-annotation sidecars.

Input: ``manifest.sqlite`` (which annotations, which ordinals) plus the ``#names``
section of each sidecar. Annotations included: ``status='indexed'`` and
``in_catalogue=1`` (stale shards that the live catalogue no longer lists are left out).
The output has the same schema as the legacy build, so the API reader is unchanged.

Algorithm, two-level external merge (flat memory, scratch on the NFS work area):

1. Annotations are taken in ordinal order in batches of ``batch`` (default 256). The
   sorted ``#names`` streams of a batch are merged with ``heapq.merge`` and written to one
   gzip *run* file with one line per ``(name, kind)``: ``name TAB kind TAB ord,ord,...``.
2. The runs are merged again (``heapq.merge``); a name's ordinals arrive ascending because
   batches cover increasing ordinal ranges. Each name becomes one ``name_dict`` row and one
   ``posting`` row.

Memory: one open gzip reader per sidecar in a batch (about 256 x ~50 KiB) and one per run
(about 70), plus one name's ordinal lists. Expected peak well under 200 MiB.

Validation (no legacy tables involved): per-annotation row counts read from each sidecar
equal the manifest, sidecar and manifest fingerprints agree, sorted order is enforced,
(annotation, name, kind) pair totals equal the sum of the posting sizes, and the postings
of a sample of annotations are checked name by name.
"""
from __future__ import annotations

import gzip
import hashlib
import heapq
import json
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

from gene_corpus.sql.postings import decode_postings, encode_postings
from gene_corpus.sql.serve_presence import (
    SCHEMA_VERSION,
    BuildError,
    _alias_only_count,
    _open_output,
    _source_code,
)
from gene_corpus.sql.sidecar import SHARD_NAME, SidecarError, iter_names, read_header, sidecar_path_for

INSERT_BATCH = 20_000


def _peak_rss_mib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _annotation_names(path: Path, ord_: int, counts: dict[int, int],
                      prefix: str | None) -> Iterator[tuple[str, str, int]]:
    """Yield ``(name, kind, ord)`` in sidecar order; enforce sorting; record the row count."""
    prev: tuple[str, str] | None = None
    n = 0
    try:
        for name, kind, _local in iter_names(path):
            key = (name, kind)
            if prev is not None and key < prev:
                raise BuildError(f"{path}: names not sorted ({prev!r} then {key!r})")
            prev = key
            n += 1
            if prefix is None or name.startswith(prefix):
                yield name, kind, ord_
    except SidecarError as exc:
        raise BuildError(f"{path}: {exc}") from exc
    counts[ord_] = n


def _write_run(run_path: Path, merged: Iterator[tuple[str, str, int]]) -> tuple[int, int, int]:
    """Write one run; return (lines, sym_pairs, alias_pairs)."""
    lines = sym_pairs = alias_pairs = 0
    with gzip.open(run_path, "wt", encoding="utf-8", compresslevel=1, newline="\n") as fh:
        for (name, kind), grp in groupby(merged, key=itemgetter(0, 1)):
            ords: list[int] = []
            for _n, _k, o in grp:
                if not ords or ords[-1] != o:
                    ords.append(o)
            fh.write(f"{name}\t{kind}\t{','.join(map(str, ords))}\n")
            lines += 1
            if kind == "s":
                sym_pairs += len(ords)
            else:
                alias_pairs += len(ords)
    return lines, sym_pairs, alias_pairs


def _read_run(run_path: Path, idx: int) -> Iterator[tuple[str, str, int, list[int]]]:
    with gzip.open(run_path, "rt", encoding="utf-8", newline="\n") as fh:
        for line in fh:
            name, kind, ords = line.rstrip("\n").split("\t")
            yield name, kind, idx, [int(x) for x in ords.split(",")]


def _name_groups(runs: list[Path]) -> Iterator[tuple[str, list[int], list[int]]]:
    merged = heapq.merge(*[_read_run(p, i) for i, p in enumerate(runs)],
                         key=lambda t: (t[0], t[1]))
    for name, grp in groupby(merged, key=itemgetter(0)):
        sym: list[int] = []
        alias: list[int] = []
        for _n, kind, _i, ords in grp:
            (sym if kind == "s" else alias).extend(ords)
        yield name, sym, alias


def _select_annotations(manifest: sqlite3.Connection, limit: int | None) -> list[tuple]:
    sql = ("SELECT ord, annotation_id, taxid, organism_name, assembly_accession, assembly_name, "
           "source_database, release_date, gff_path, shard_fingerprint, sidecar_fingerprint, "
           "n_name_rows FROM annotation WHERE status='indexed' AND in_catalogue=1 ORDER BY ord")
    if limit:
        sql += f" LIMIT {int(limit)}"
    return manifest.execute(sql).fetchall()   # one small tuple per annotation (~17k), not gene data


def build_from_sidecars(
    manifest_path: Path,
    per_gff: Path,
    out_file: Path,
    *,
    tmp_dir: Path,
    batch: int = 256,
    name_prefix: str | None = None,
    limit_annotations: int | None = None,
    verify_annotations: int = 30,
    keep_runs: bool = False,
    log: Callable[[str], None] | None = None,
) -> dict:
    if log is None:
        log = lambda m: print(f"[{time.strftime('%H:%M:%S')}] {m}", file=sys.stderr, flush=True)  # noqa: E731
    t0 = time.time()
    prefix = name_prefix.casefold() if name_prefix else None
    partial = bool(prefix or limit_annotations)

    manifest = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True)
    rows = _select_annotations(manifest, limit_annotations)
    if not rows:
        raise BuildError("no indexed, in-catalogue annotations in the manifest")
    for r in rows:
        if r[2] is None:
            raise BuildError(f"annotation {r[1]} has no taxid (run `catalogue`)")
    max_ord = max(r[0] for r in rows)
    named = [r for r in rows if (r[11] or 0) > 0]
    log(f"annotations in build: {len(rows):,} ({len(named):,} with names), max ord {max_ord}")

    # ---- level 1: runs ---------------------------------------------------
    tmp_dir.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix="presence-runs-", dir=tmp_dir))
    runs: list[Path] = []
    counts: dict[int, int] = {}
    pairs_sym = pairs_alias = run_lines = 0
    paths: dict[int, Path] = {}
    sha = hashlib.sha256()
    try:
        for i in range(0, len(named), batch):
            chunk = named[i:i + batch]
            streams = []
            for r in chunk:
                ord_, aid = r[0], r[1]
                path = sidecar_path_for(per_gff / aid / SHARD_NAME)
                try:
                    header = read_header(path)
                except (OSError, SidecarError) as exc:
                    raise BuildError(f"{aid}: cannot read sidecar header: {exc}") from exc
                if header.get("shard_fingerprint") != r[9]:
                    raise BuildError(f"{aid}: sidecar is stale vs manifest (re-run `sidecars`)")
                if int(header.get("n_name_rows", -1)) != r[11]:
                    raise BuildError(f"{aid}: n_name_rows differs between sidecar and manifest")
                paths[ord_] = path
                sha.update(f"{ord_}:{r[10]};".encode())
                streams.append(_annotation_names(path, ord_, counts, prefix))
            run_path = run_dir / f"run-{len(runs):05d}.tsv.gz"
            n, s, a = _write_run(run_path, heapq.merge(*streams))
            runs.append(run_path)
            run_lines += n
            pairs_sym += s
            pairs_alias += a
            log(f"level 1: {min(i + batch, len(named)):,}/{len(named):,} annotations, "
                f"run {len(runs)}, peak_rss={_peak_rss_mib():.0f}MiB")
        for r in named:
            if counts.get(r[0]) != r[11]:
                raise BuildError(f"{r[1]}: read {counts.get(r[0])} name rows, manifest says {r[11]}")

        # ---- level 2: merge runs into presence.sqlite --------------------
        out = _open_output(out_file)
        try:
            out.executemany(
                "INSERT INTO annotation (ord, annotation_id, taxid, organism_name, "
                "assembly_accession, assembly_name, source, release_date, gff_path) "
                "VALUES (?,?,?,?,?,?,?,?,?)",
                [(r[0], r[1], int(r[2]), r[3], r[4], r[5], _source_code(r[6]), r[7], r[8])
                 for r in rows])
            lineage_rows = 0
            batch_l: list[tuple[int, int]] = []
            wanted = {r[0] for r in rows}
            for taxid, ord_ in manifest.execute("SELECT taxid, ord FROM lineage"):
                if ord_ in wanted:
                    batch_l.append((taxid, ord_))
                    if len(batch_l) >= 50_000:
                        out.executemany("INSERT OR IGNORE INTO lineage (taxid, ord) VALUES (?,?)", batch_l)
                        lineage_rows += len(batch_l)
                        batch_l.clear()
            if batch_l:
                out.executemany("INSERT OR IGNORE INTO lineage (taxid, ord) VALUES (?,?)", batch_l)
                lineage_rows += len(batch_l)
            del wanted

            n_names = n_alias_only_names = blob_bytes = varint_blobs = bitset_blobs = 0
            got_sym = got_alias = 0
            dict_batch: list[tuple] = []
            post_batch: list[tuple] = []
            last_log = time.time()

            def flush() -> None:
                if dict_batch:
                    out.executemany("INSERT INTO name_dict (name, name_id, n_sym, n_alias_only) "
                                    "VALUES (?,?,?,?)", dict_batch)
                    out.executemany("INSERT INTO posting (name_id, sym, alias) VALUES (?,?,?)",
                                    post_batch)
                    dict_batch.clear()
                    post_batch.clear()

            try:
                for name, sym, alias in _name_groups(runs):
                    n_names += 1
                    got_sym += len(sym)
                    got_alias += len(alias)
                    sym_blob = encode_postings(sym, max_ord) if sym else None
                    alias_blob = encode_postings(alias, max_ord) if alias else None
                    if not sym:
                        n_alias_only_names += 1
                    for blob in (sym_blob, alias_blob):
                        if blob is not None:
                            blob_bytes += len(blob)
                            if blob[0] == 0x01:
                                varint_blobs += 1
                            else:
                                bitset_blobs += 1
                    dict_batch.append((name, n_names, len(sym), _alias_only_count(sym, alias)))
                    post_batch.append((n_names, sym_blob, alias_blob))
                    if len(dict_batch) >= INSERT_BATCH:
                        flush()
                    if time.time() - last_log > 30:
                        last_log = time.time()
                        log(f"level 2: {n_names:,} names (at {name!r}) peak_rss={_peak_rss_mib():.0f}MiB")
            except ValueError as exc:
                raise BuildError(f"posting encode failed: {exc}") from exc
            flush()
            out.commit()

            if got_sym != pairs_sym or got_alias != pairs_alias:
                raise BuildError(
                    f"pair totals differ: runs sym={pairs_sym:,} alias={pairs_alias:,}, "
                    f"postings sym={got_sym:,} alias={got_alias:,}")

            sha.update(f"{n_names}:{pairs_sym}:{pairs_alias}:{prefix}:{limit_annotations}".encode())
            data_version = sha.hexdigest()[:16]
            meta = {
                "schema_version": str(SCHEMA_VERSION),
                "data_version": data_version,
                "built_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
                "source": "sidecars",
                "manifest": str(manifest_path),
                "n_annotations": str(len(rows)),
                "n_annotations_with_names": str(len(named)),
                "n_names": str(n_names),
                "max_ord": str(max_ord),
                "partial": "1" if partial else "0",
            }
            out.executemany("INSERT INTO meta (key, value) VALUES (?, ?)", list(meta.items()))
            out.commit()

            checked_ann, checked_pairs = _verify(out, named, paths, verify_annotations, prefix, log)
        finally:
            out.close()
    finally:
        manifest.close()
        if not keep_runs:
            shutil.rmtree(run_dir, ignore_errors=True)

    stats = {
        "data_version": data_version,
        "partial": partial,
        "source": "sidecars",
        "n_annotations": len(rows),
        "n_annotations_with_names": len(named),
        "n_names": n_names,
        "names_alias_only": n_alias_only_names,
        "pairs_sym": pairs_sym,
        "pairs_alias": pairs_alias,
        "lineage_rows": lineage_rows,
        "runs": len(runs),
        "varint_blobs": varint_blobs,
        "bitset_blobs": bitset_blobs,
        "blob_bytes": blob_bytes,
        "output_bytes": out_file.stat().st_size,
        "verified_annotations": checked_ann,
        "verified_pairs": checked_pairs,
        "peak_rss_mib": round(_peak_rss_mib(), 1),
        "seconds": round(time.time() - t0, 1),
    }
    log("build done: " + json.dumps(stats))
    return stats


def _verify(out: sqlite3.Connection, named: list[tuple], paths: dict[int, Path],
            sample: int, prefix: str | None, log) -> tuple[int, int]:
    """Check every (name, kind) of evenly spaced annotations against the postings."""
    if sample <= 0 or not named:
        return 0, 0
    step = max(1, len(named) // sample)
    picks = named[::step]
    pairs = 0
    for r in picks:
        ord_ = r[0]
        for name, kind, _local in iter_names(paths[ord_]):
            if prefix is not None and not name.startswith(prefix):
                continue
            row = out.execute(
                "SELECT p.sym, p.alias FROM name_dict d JOIN posting p USING (name_id) "
                "WHERE d.name = ?", (name,)).fetchone()
            if row is None:
                raise BuildError(f"verification: name {name!r} of {r[1]} missing")
            blob = row[0] if kind == "s" else row[1]
            if not blob or ord_ not in decode_postings(blob):
                raise BuildError(f"verification: {name!r} ({kind}) lacks ordinal {ord_} ({r[1]})")
            pairs += 1
    log(f"verified {len(picks)} annotations, {pairs:,} name rows")
    return len(picks), pairs
