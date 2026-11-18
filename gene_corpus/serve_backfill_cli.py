"""
Sub-commands for the sidecar / manifest side of ``serve_build`` (PLAN-build-tables.md M3).

    python -m gene_corpus.serve_build catalogue --work-dir DIR [--api-base URL | --tsv-file F]
    python -m gene_corpus.serve_build sidecars  --work-dir DIR [--workers 2] [--force] [--sample N]
    python -m gene_corpus.serve_build status    --work-dir DIR

Defaults: shards in ``<work-dir>/per_gff``, manifest at ``<work-dir>/serve/manifest.sqlite``.
Sidecars are written next to the shards (``per_gff/<id>/sidecar.tsv.gz``). Memory stays
small: workers hold one shard's streaming state (about 50 MiB each in M2), the manifest
process holds nothing but the current result.
"""
from __future__ import annotations

import argparse
import os
import random
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from gene_corpus.sql.catalogue import (
    DEFAULT_API_BASE,
    CatalogueError,
    apply_catalogue,
    fetch_report_lines,
    load_catalogue,
    report_url,
)
from gene_corpus.sql.manifest import (
    STATUS_FAILED,
    open_manifest,
    record_result,
    work_shard,
)

COMMANDS = ("catalogue", "sidecars", "status", "build", "missing", "loci")


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def _common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--work-dir", type=Path, required=True, help="corpus root (contains per_gff/)")
    p.add_argument("--per-gff", type=Path, default=None, help="default: <work-dir>/per_gff")
    p.add_argument("--manifest", type=Path, default=None,
                   help="default: <work-dir>/serve/manifest.sqlite")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="serve_build")
    sub = ap.add_subparsers(dest="command", required=True)

    c = sub.add_parser("catalogue", help="load Annotrieve annotation metadata into the manifest")
    _common(c)
    c.add_argument("--api-base", default=DEFAULT_API_BASE)
    c.add_argument("--url", default=None, help="full report URL (overrides --api-base)")
    c.add_argument("--tsv-file", type=Path, default=None, help="read a saved report instead of fetching")
    c.add_argument("--timeout", type=float, default=120.0)

    s = sub.add_parser("sidecars", help="write missing/stale sidecars and register shards")
    _common(s)
    s.add_argument("--workers", type=int, default=2)
    s.add_argument("--force", action="store_true", help="rewrite every sidecar")
    s.add_argument("--ids", type=Path, default=None, help="file with one shard directory name per line")
    s.add_argument("--sample", type=int, default=None, help="only N random shards (testing)")
    s.add_argument("--limit", type=int, default=None, help="only the first N shards (testing)")
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--max-worker-mib", type=int, default=768,
                   help="per-worker address-space cap; 0 disables (default 768)")

    b = sub.add_parser("build", help="build presence.sqlite from the manifest and the sidecars")
    _common(b)
    b.add_argument("--out-dir", type=Path, default=None, help="default: <work-dir>/serve")
    b.add_argument("--tmp-dir", type=Path, default=None,
                   help="scratch for merge runs; default <work-dir>/serve/tmp (NFS, not local disk)")
    b.add_argument("--batch", type=int, default=256, help="sidecars merged per run (default 256)")
    b.add_argument("--name-prefix", default=None, help="trial build: only names starting with this")
    b.add_argument("--limit-annotations", type=int, default=None, help="trial build: first N annotations")
    b.add_argument("--verify-annotations", type=int, default=30)
    b.add_argument("--keep-runs", action="store_true")

    lo = sub.add_parser("loci", help="build loci.sqlite next to the current presence.sqlite")
    _common(lo)
    lo.add_argument("--out-dir", type=Path, default=None, help="default: <work-dir>/serve")
    lo.add_argument("--tmp-dir", type=Path, default=None, help="default: <out-dir>/tmp (NFS)")
    lo.add_argument("--gene-budget", type=int, default=300_000,
                    help="named genes held in memory per batch (default 300000)")
    lo.add_argument("--verify-annotations", type=int, default=30)
    lo.add_argument("--keep-runs", action="store_true")

    m = sub.add_parser("missing", help="list catalogue annotations that have no shard yet")
    _common(m)
    m.add_argument("--out", type=Path, default=None,
                   help="write harvestable annotation_ids here (default: <work-dir>/serve/missing_ids.txt)")
    m.add_argument("--limit", type=int, default=None, help="only the first N ids (pilot)")

    t = sub.add_parser("status", help="show manifest counts")
    _common(t)
    return ap


def _paths(args) -> tuple[Path, Path]:
    per_gff = args.per_gff or (args.work_dir / "per_gff")
    manifest = args.manifest or (args.work_dir / "serve" / "manifest.sqlite")
    manifest.parent.mkdir(parents=True, exist_ok=True)
    return per_gff, manifest


def _cmd_catalogue(args) -> int:
    _per_gff, manifest = _paths(args)
    conn = open_manifest(manifest)
    try:
        if args.tsv_file:
            with open(args.tsv_file, encoding="utf-8") as fh:
                stats = load_catalogue(conn, (line.rstrip("\r\n") for line in fh))
            source = str(args.tsv_file)
        else:
            source = args.url or report_url(args.api_base)
            _log(f"fetching {source}")
            stats = load_catalogue(conn, fetch_report_lines(source, timeout=args.timeout))
    except CatalogueError as exc:
        print(f"CATALOGUE FAILED: {exc}", file=sys.stderr)
        return 1
    n_manifest = conn.execute("SELECT COUNT(*) FROM annotation").fetchone()[0]
    print(f"catalogue: {stats.rows:,} rows from {source}")
    print(f"  duplicates={stats.duplicate_ids} blank_ids={stats.skipped_blank_id}")
    print(f"  manifest annotations={n_manifest:,} matched={stats.matched_annotations:,} "
          f"shard_only={n_manifest - stats.matched_annotations:,} catalogue_only={stats.catalogue_only:,}")
    print(f"  database_mismatch={stats.database_mismatch} lineage_rows={stats.lineage_rows:,}")
    if n_manifest == 0:
        print("  (manifest is empty: run `sidecars` too, then `catalogue` again or let `sidecars` re-apply)")
    if stats.shard_only:
        print(f"  first shard_only ids: {', '.join(stats.shard_only[:5])}")
    return 0


def _list_dirs(per_gff: Path, args) -> list[str]:
    if args.ids:
        dirs = [ln.strip() for ln in args.ids.read_text().splitlines() if ln.strip()]
    else:
        with os.scandir(per_gff) as it:
            dirs = sorted(e.name for e in it if e.is_dir(follow_symlinks=False))
    if args.sample is not None:
        dirs = sorted(random.Random(args.seed).sample(dirs, min(args.sample, len(dirs))))
    if args.limit is not None:
        dirs = dirs[: args.limit]
    return dirs


def _init_worker(max_mib: int) -> None:
    """
    Hard memory cap for a worker (address space). A shard that would exceed it raises
    MemoryError inside ``work_shard`` and is recorded as ``failed``, instead of the kernel
    killing a process (or a neighbouring service) when RAM runs out.
    """
    if max_mib > 0:
        try:
            import resource

            limit = max_mib * 1024 * 1024
            resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
        except (ImportError, ValueError, OSError):
            pass


def _work(task: tuple[str, str, bool]) -> dict:
    import resource

    per_gff, annotation_dir, force = task
    res = work_shard(per_gff, annotation_dir, force)
    res["maxrss_mib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    return res


def _cmd_sidecars(args) -> int:
    per_gff, manifest = _paths(args)
    if not per_gff.is_dir():
        print(f"per_gff not found: {per_gff}", file=sys.stderr)
        return 2
    dirs = _list_dirs(per_gff, args)
    conn = open_manifest(manifest)
    _log(f"{len(dirs):,} shards, workers={args.workers}, manifest={manifest}")
    tasks = [(str(per_gff), d, args.force) for d in dirs]

    counts: Counter = Counter()
    written = 0
    bytes_written = 0
    failures: list[tuple[str, str]] = []
    peak_rss = 0.0
    t0 = time.perf_counter()
    pending_commit = 0
    interrupted = False

    def handle(res: dict, i: int) -> None:
        nonlocal written, bytes_written, pending_commit, peak_rss
        peak_rss = max(peak_rss, res.get("maxrss_mib", 0.0))
        status = record_result(conn, res)
        counts[status] += 1
        if res.get("written"):
            written += 1
            bytes_written += res["info"]["size_bytes"]
        if status == STATUS_FAILED:
            failures.append((res["dir"], res.get("error") or ""))
        pending_commit += 1
        if pending_commit >= 200:
            conn.commit()
            pending_commit = 0
        if i % 500 == 0 or i == len(tasks):
            el = time.perf_counter() - t0
            rate = i / el if el else 0
            eta = (len(tasks) - i) / rate if rate else 0
            _log(f"  {i:,}/{len(tasks):,}  {rate:.1f} shards/s  eta {eta / 60:.1f} min  "
                 f"written={written:,} failed={counts[STATUS_FAILED]} peak_rss={peak_rss:.0f}MiB")

    try:
        if args.workers <= 1:
            _init_worker(args.max_worker_mib)
            for i, task in enumerate(tasks, 1):
                handle(_work(task), i)
        else:
            with ProcessPoolExecutor(max_workers=args.workers, initializer=_init_worker,
                                     initargs=(args.max_worker_mib,)) as ex:
                # ordered results: ordinals are assigned in sorted shard order, deterministically
                for i, res in enumerate(ex.map(_work, tasks, chunksize=4), 1):
                    handle(res, i)
    except KeyboardInterrupt:
        interrupted = True
        _log("interrupted: saving progress (a re-run continues where this stopped)")
    conn.commit()

    if not interrupted and conn.execute("SELECT COUNT(*) FROM catalogue").fetchone()[0]:
        st = apply_catalogue(conn)
        conn.commit()
        _log(f"catalogue re-applied: matched={st.matched_annotations:,} "
             f"shard_only={conn.execute('SELECT COUNT(*) FROM annotation WHERE in_catalogue=0').fetchone()[0]:,}")

    el = time.perf_counter() - t0
    print(f"done in {el / 60:.1f} min: {dict(counts)}; sidecars written {written:,} "
          f"({bytes_written / 2**20:.1f} MiB); peak worker RSS {peak_rss:.0f} MiB")
    for d, err in failures[:20]:
        print(f"  FAILED {d}: {err}")
    if len(failures) > 20:
        print(f"  ... and {len(failures) - 20} more (see `status`)")
    return 130 if interrupted else (1 if failures else 0)


def _cmd_loci(args) -> int:
    import json

    from gene_corpus.sql.loci_from_sidecars import LOCI_FILE, build_loci
    from gene_corpus.sql.serve_presence import BuildError

    per_gff, manifest = _paths(args)
    out_dir = args.out_dir or (args.work_dir / "serve")
    current = out_dir / "current"
    if not current.exists():
        print(f"no published presence index at {current} (run `build` first)", file=sys.stderr)
        return 2
    version_dir = current.resolve()
    building = out_dir / ".loci.building.sqlite"
    try:
        stats = build_loci(
            version_dir / "presence.sqlite", manifest, per_gff, building,
            tmp_dir=args.tmp_dir or (out_dir / "tmp"), gene_budget=args.gene_budget,
            verify_sample=args.verify_annotations, keep_runs=args.keep_runs)
    except BuildError as exc:
        print(f"LOCI BUILD FAILED: {exc}", file=sys.stderr)
        building.unlink(missing_ok=True)
        return 1
    os.replace(building, version_dir / LOCI_FILE)
    (version_dir / "loci_report.json").write_text(json.dumps(stats, indent=2))
    print(f"published {version_dir / LOCI_FILE}")
    return 0


def _cmd_missing(args) -> int:
    from gene_corpus.sql.manifest import STATUS_SKIPPED_COMMUNITY, classify_source

    _per_gff, manifest = _paths(args)
    if not manifest.exists():
        print(f"no manifest at {manifest}", file=sys.stderr)
        return 2
    conn = open_manifest(manifest)
    if not conn.execute("SELECT COUNT(*) FROM catalogue").fetchone()[0]:
        print("catalogue is empty: run `catalogue` first", file=sys.stderr)
        return 2
    by_db: Counter = Counter()
    ids: list[str] = []
    n_community = n_unknown = 0
    for aid, database in conn.execute(
        "SELECT c.annotation_id, c.database FROM catalogue c WHERE NOT EXISTS "
        "(SELECT 1 FROM annotation a WHERE a.annotation_id = c.annotation_id) "
        "ORDER BY c.annotation_id"
    ):
        by_db[database or "(empty)"] += 1
        skip = classify_source(database)
        if skip is None:
            ids.append(aid)
        elif skip == STATUS_SKIPPED_COMMUNITY:
            n_community += 1
        else:
            n_unknown += 1
    total = sum(by_db.values())
    print(f"catalogue annotations without a shard: {total:,}")
    for db, n in by_db.most_common():
        print(f"  {db:<24} {n:>7,}")
    print(f"  -> harvestable (Ensembl/GenBank/RefSeq): {len(ids):,}; community (skipped): "
          f"{n_community:,}; unknown source (skipped): {n_unknown:,}")
    if args.limit:
        ids = ids[: args.limit]
    out = args.out or (args.work_dir / "serve" / "missing_ids.txt")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(f"{i}\n" for i in ids), encoding="utf-8")
    print(f"wrote {len(ids):,} ids to {out}")
    return 0


def _cmd_status(args) -> int:
    _per_gff, manifest = _paths(args)
    if not manifest.exists():
        print(f"no manifest at {manifest}", file=sys.stderr)
        return 2
    conn = open_manifest(manifest)
    q = lambda sql: conn.execute(sql).fetchall()  # noqa: E731
    print(f"manifest: {manifest}")
    for status, n in q("SELECT status, COUNT(*) FROM annotation GROUP BY status ORDER BY 2 DESC"):
        print(f"  {status:<20} {n:>8,}")
    (tot,) = q("SELECT COUNT(*) FROM annotation")[0]
    (in_cat,) = q("SELECT COUNT(*) FROM annotation WHERE in_catalogue=1")[0]
    (cat,) = q("SELECT COUNT(*) FROM catalogue")[0]
    (nonames,) = q("SELECT COUNT(*) FROM annotation WHERE status='indexed' AND n_name_rows=0")[0]
    print(f"  annotations={tot:,} in_catalogue={in_cat:,} catalogue_rows={cat:,} "
          f"indexed_without_names={nonames:,}")
    for aid, err in q("SELECT annotation_id, error FROM annotation WHERE status='failed' LIMIT 10"):
        print(f"  failed {aid}: {err}")
    return 0


def _cmd_build(args) -> int:
    import json

    from gene_corpus.sql.presence_from_sidecars import build_from_sidecars
    from gene_corpus.sql.serve_presence import BuildError, publish_version

    per_gff, manifest = _paths(args)
    if not manifest.exists():
        print(f"no manifest at {manifest}", file=sys.stderr)
        return 2
    out_dir = args.out_dir or (args.work_dir / "serve")
    tmp_dir = args.tmp_dir or (out_dir / "tmp")
    out_dir.mkdir(parents=True, exist_ok=True)
    building = out_dir / ".presence.building.sqlite"
    try:
        stats = build_from_sidecars(
            manifest, per_gff, building, tmp_dir=tmp_dir, batch=args.batch,
            name_prefix=args.name_prefix, limit_annotations=args.limit_annotations,
            verify_annotations=args.verify_annotations, keep_runs=args.keep_runs)
    except BuildError as exc:
        print(f"BUILD FAILED: {exc}", file=sys.stderr)
        building.unlink(missing_ok=True)
        return 1
    if stats["partial"]:
        part = out_dir / f"partial-{stats['data_version']}"
        part.mkdir(parents=True, exist_ok=True)
        os.replace(building, part / "presence.sqlite")
        (part / "build_report.json").write_text(json.dumps(stats, indent=2))
        print(f"partial build written to {part / 'presence.sqlite'} (current not changed)")
    else:
        final = publish_version(out_dir, building, stats)
        print(f"published {final}; current -> {(out_dir / 'current').resolve().name}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return {"catalogue": _cmd_catalogue, "sidecars": _cmd_sidecars, "status": _cmd_status,
            "build": _cmd_build, "missing": _cmd_missing,
            "loci": _cmd_loci}[
        args.command](args)
