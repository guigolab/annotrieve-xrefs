#!/usr/bin/env python3
"""
Phase 0 measurements for the annotrieve-xrefs serve-index plan.

Read-only and standard-library only. Designed to run in the background against
the NFS-mounted corpus, with bounded memory (small SQLite cache, streaming
scans). It writes ``report.json`` and ``report.md`` after every step, so a
partial run still leaves useful results.

    nohup python3 tools/phase0_measure.py \
        --corpus-dir /data/annotrieve/gene_corpus \
        --out-dir /tmp/phase0 \
        > /tmp/phase0.log 2>&1 &

Optional live-API load test (run against the container directly, NOT through
traefik, whose 3 req/s per IP limit would dominate the result):

    ... --http-base http://127.0.0.1:5000 --api-pid <gunicorn master pid>

Steps (each can be skipped; see --help):
  env       mount type, RAM, cgroup limit, SQLite version
  files     DB / WAL sizes, shard count and total shard size
  nfs       sequential throughput and random-read latency on the DB file
  counts    annotation / lineage counts, sqlite_stat1 estimates
  names     streaming scan of symbol+alias rows: distinct names, histogram,
            estimated size of the proposed postings table, sample name sets
  queries   timings of the current /hits and /names SQL on sampled inputs
  http      optional concurrent load test with peak RSS of the API process tree
  dbstat    per-table / per-index sizes (reads the whole file; slowest, last)

Run it on a published snapshot, or when no sync/merge is writing the DB: the
file is opened with immutable=1, exactly like the API does.
"""
from __future__ import annotations

import argparse
import heapq
import json
import math
import os
import platform
import random
import signal
import sqlite3
import statistics
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

DB_NAME = "gene_corpus.sqlite"
SQLITE_MAX_VARS = 900  # same chunk size as api/services/hits.py
QUERY_TIMEOUT_S = 300.0

REPORT: dict = {}
OUT_DIR: Path | None = None
T0 = time.time()


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[{time.time() - T0:8.1f}s] {msg}", file=sys.stderr, flush=True)


def human(n: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(n) < 1024 or unit == "TiB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{int(n)} B"
        n /= 1024
    return str(n)


def chunks(items: list, size: int):
    for i in range(0, len(items), size):
        yield items[i : i + size]


def pct(values: list[float], p: float) -> float:
    if not values:
        return float("nan")
    s = sorted(values)
    k = min(len(s) - 1, max(0, int(math.ceil(p / 100 * len(s))) - 1))
    return s[k]


def timed(fn, *args, **kwargs):
    t = time.perf_counter()
    out = fn(*args, **kwargs)
    return out, time.perf_counter() - t


def open_db(path: Path, cache_kib: int, tmp_dir: str | None) -> sqlite3.Connection:
    uri = f"file:{path}?mode=ro&immutable=1"
    conn = sqlite3.connect(uri, uri=True, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    conn.execute(f"PRAGMA cache_size=-{cache_kib}")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute("PRAGMA mmap_size=0")
    conn.execute("PRAGMA query_only=ON")
    if tmp_dir:
        os.environ["SQLITE_TMPDIR"] = tmp_dir
    return conn


def with_timeout(conn: sqlite3.Connection, seconds: float):
    """Abort the running statement after ``seconds`` (progress handler)."""
    deadline = time.monotonic() + seconds
    conn.set_progress_handler(lambda: 1 if time.monotonic() > deadline else 0, 100_000)


def clear_timeout(conn: sqlite3.Connection) -> None:
    conn.set_progress_handler(None, 0)


def table_exists(conn: sqlite3.Connection, name: str) -> bool:
    return (
        conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=? LIMIT 1",
            (name,),
        ).fetchone()
        is not None
    )


# --------------------------------------------------------------------------
# report writing
# --------------------------------------------------------------------------
def save_report() -> None:
    if OUT_DIR is None:
        return
    REPORT["elapsed_s"] = round(time.time() - T0, 1)
    (OUT_DIR / "report.json").write_text(json.dumps(REPORT, indent=2, default=str))
    (OUT_DIR / "report.md").write_text(render_md())


def md_table(headers: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def render_md() -> str:
    L: list[str] = ["# Phase 0 measurements", ""]
    L.append(f"Generated {datetime.now().isoformat(timespec='seconds')} "
             f"(elapsed {REPORT.get('elapsed_s', 0)} s)")
    L.append("")
    env = REPORT.get("env")
    if env:
        L += ["## Environment", ""]
        L += [f"- `{k}`: {v}" for k, v in env.items()]
        L.append("")
    files = REPORT.get("files")
    if files:
        L += ["## Files", ""]
        L += [f"- {k}: {v}" for k, v in files.items()]
        L.append("")
    nfs = REPORT.get("nfs")
    if nfs:
        L += ["## NFS read behaviour", ""]
        L += [f"- {k}: {v}" for k, v in nfs.items()]
        L.append("")
    counts = REPORT.get("counts")
    if counts:
        L += ["## Row counts", ""]
        L += [f"- {k}: {v}" for k, v in counts.items() if k != "stat1"]
        if counts.get("stat1"):
            L += ["", "`sqlite_stat1` estimates (first number = rows):", ""]
            L.append(md_table(["table", "index", "stat"], counts["stat1"]))
        L.append("")
    names = REPORT.get("names")
    if names:
        L += ["## Symbol / alias scan", ""]
        for k in ("rows_symbol", "rows_alias", "distinct_names", "names_in_both",
                  "avg_name_bytes", "annotations", "scan_seconds"):
            if k in names:
                L.append(f"- {k}: {names[k]}")
        L += ["", "Rows-per-name histogram (symbol rows + alias rows):", ""]
        L.append(md_table(["bucket", "names", "share"], names["histogram"]))
        est = names["estimate"]
        L += ["", "Estimated size of the proposed `name` table (rough):", ""]
        L += [f"- {k}: {v}" for k, v in est.items()]
        L += ["", "Top names by row count: " + ", ".join(
            f"{n}({c})" for n, c in names["top_names"][:20]), ""]
    q = REPORT.get("queries")
    if q:
        L += ["## Current SQL timings (seconds)", "",
              "`first` is the first run (cold-ish: the OS/NFS cache may already hold "
              "pages from earlier steps); `warm_med` is the median of the repeats.", ""]
        L.append(md_table(
            ["scenario", "names", "postings", "hits", "first", "warm_med", "parts (first run)"],
            q["hits_rows"]))
        L.append("")
        L.append(md_table(
            ["prefix/scope", "results", "first", "warm_med", "note"], q["names_rows"]))
        L.append("")
    h = REPORT.get("http")
    if h:
        L += ["## HTTP load test", ""]
        L += [f"- {k}: {v}" for k, v in h.items()]
        L.append("")
    d = REPORT.get("dbstat")
    if d:
        L += ["## dbstat (bytes by table/index)", ""]
        if isinstance(d, dict) and "error" in d:
            L.append(f"- {d['error']}")
        else:
            L.append(md_table(["name", "size", "pages"], d))
        L.append("")
    return "\n".join(L)


# --------------------------------------------------------------------------
# step: env
# --------------------------------------------------------------------------
def mount_info(path: Path) -> dict:
    best = None
    try:
        real = os.path.realpath(path)
        with open("/proc/mounts") as f:
            for line in f:
                parts = line.split()
                if len(parts) < 4:
                    continue
                mnt = parts[1].replace("\\040", " ")
                if real == mnt or real.startswith(mnt.rstrip("/") + "/"):
                    if best is None or len(mnt) > len(best[1]):
                        best = (parts[0], mnt, parts[2], parts[3])
    except OSError:
        pass
    if not best:
        return {"mount": "unknown"}
    return {"device": best[0], "mountpoint": best[1], "fstype": best[2], "options": best[3]}


def step_env(corpus: Path) -> None:
    env: dict = {
        "python": platform.python_version(),
        "sqlite": sqlite3.sqlite_version,
        "cpus": os.cpu_count(),
    }
    try:
        for line in open("/proc/meminfo"):
            if line.startswith(("MemTotal", "MemAvailable")):
                k, v = line.split(":")
                env[k] = human(int(v.split()[0]) * 1024)
    except OSError:
        pass
    for cg in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            raw = open(cg).read().strip()
            env["cgroup_memory_limit"] = raw if raw == "max" else human(int(raw))
            break
        except (OSError, ValueError):
            continue
    env["corpus_mount"] = mount_info(corpus)
    REPORT["env"] = env
    save_report()


# --------------------------------------------------------------------------
# step: files
# --------------------------------------------------------------------------
def step_files(corpus: Path, skip_shards: bool) -> None:
    out: dict = {}
    for name in (DB_NAME, DB_NAME + "-wal", DB_NAME + "-shm", "gene_corpus.sqlite-journal"):
        p = corpus / name
        if p.exists():
            out[name] = human(p.stat().st_size)
    per_gff = corpus / "per_gff"
    if per_gff.is_dir() and not skip_shards:
        n = total = missing = 0
        sizes: list[int] = []
        t = time.perf_counter()
        with os.scandir(per_gff) as it:
            for i, d in enumerate(it):
                if not d.is_dir():
                    continue
                g = Path(d.path) / "genes.sqlite"
                try:
                    s = g.stat().st_size
                except OSError:
                    missing += 1
                    continue
                n += 1
                total += s
                sizes.append(s)
                if i and i % 5000 == 0:
                    log(f"files: stat'd {i} shard dirs")
        out["shards"] = n
        out["shards_missing_genes_sqlite"] = missing
        out["shards_total"] = human(total)
        if sizes:
            out["shard_median"] = human(statistics.median(sizes))
            out["shard_p95"] = human(pct(sizes, 95))
            out["shard_max"] = human(max(sizes))
        out["shard_scan_seconds"] = round(time.perf_counter() - t, 1)
    REPORT["files"] = out
    save_report()


# --------------------------------------------------------------------------
# step: nfs
# --------------------------------------------------------------------------
def step_nfs(db_path: Path, seq_mb: int, rand_reads: int) -> None:
    out: dict = {}
    size = db_path.stat().st_size
    fd = os.open(db_path, os.O_RDONLY)
    try:
        # sequential read from the middle of the file (likely cold)
        want = min(seq_mb * 1024 * 1024, size)
        off = max(0, size // 2 - want // 2)
        buf = 1024 * 1024
        read = 0
        t = time.perf_counter()
        while read < want:
            data = os.pread(fd, buf, off + read)
            if not data:
                break
            read += len(data)
        dt = time.perf_counter() - t
        out["seq_read"] = f"{human(read)} in {dt:.2f}s = {human(read / dt if dt else 0)}/s"

        # random page reads (16 KiB) across the file
        rng = random.Random(1)
        lat: list[float] = []
        page = 16384
        for _ in range(rand_reads):
            o = rng.randrange(0, max(1, size - page)) // 4096 * 4096
            t = time.perf_counter()
            os.pread(fd, page, o)
            lat.append((time.perf_counter() - t) * 1000)
        out["random_16KiB_reads"] = rand_reads
        out["random_read_ms_p50"] = round(pct(lat, 50), 2)
        out["random_read_ms_p95"] = round(pct(lat, 95), 2)
        out["random_read_ms_max"] = round(max(lat), 2)
    finally:
        os.close(fd)
    REPORT["nfs"] = out
    save_report()


# --------------------------------------------------------------------------
# step: counts
# --------------------------------------------------------------------------
def step_counts(conn: sqlite3.Connection) -> int:
    out: dict = {}
    n_ann = 0
    if table_exists(conn, "annotation"):
        n_ann = conn.execute("SELECT COUNT(*) FROM annotation").fetchone()[0]
        out["annotation"] = n_ann
        rows = conn.execute(
            "SELECT source_database, COUNT(*) FROM annotation GROUP BY 1"
        ).fetchall()
        out["annotation_by_source"] = {r[0]: r[1] for r in rows}
    if table_exists(conn, "annotation_lineage"):
        out["annotation_lineage_rows"] = conn.execute(
            "SELECT COUNT(*) FROM annotation_lineage"
        ).fetchone()[0]
        out["lineage_distinct_taxids"] = conn.execute(
            "SELECT COUNT(DISTINCT taxid) FROM annotation_lineage"
        ).fetchone()[0]
    if table_exists(conn, "sqlite_stat1"):
        out["stat1"] = [
            [r[0], r[1], r[2]] for r in conn.execute("SELECT tbl, idx, stat FROM sqlite_stat1")
        ]
    else:
        out["stat1"] = []
        out["note"] = "no sqlite_stat1 (ANALYZE was not run); row counts come from the name scan"
    tables = [r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")]
    out["tables_present"] = ", ".join(tables)
    REPORT["counts"] = out
    save_report()
    return int(n_ann)


# --------------------------------------------------------------------------
# step: names (streaming scan of symbol + alias)
# --------------------------------------------------------------------------
def posting_bytes(n: int, total_ann: int) -> int:
    """Rough size of one posting blob (varint deltas vs bitset, whichever is smaller)."""
    if n <= 0:
        return 0
    avg_gap = total_ann / n
    per = 1 if avg_gap < 64 else (2 if avg_gap < 8000 else 3)
    varint = n * per + 1
    bitset = math.ceil(total_ann / 8) + 1
    return min(varint, bitset)


BUCKETS = [(1, 1), (2, 2), (3, 5), (6, 10), (11, 100), (101, 1000), (1001, 10**9)]


def stream_ns(conn: sqlite3.Connection, ns: str):
    cur = conn.execute(
        "SELECT accession, COUNT(*) FROM xref_file WHERE namespace = ? "
        "GROUP BY accession ORDER BY accession",
        (ns,),
    )
    for acc, n in cur:
        yield acc, int(n)


def merge_streams(a, b):
    """Merge two (accession, n) iterators ordered by accession. Yields (acc, n_sym, n_alias)."""
    sentinel = (None, 0)
    x = next(a, sentinel)
    y = next(b, sentinel)
    while x[0] is not None or y[0] is not None:
        if y[0] is None or (x[0] is not None and x[0] < y[0]):
            yield x[0], x[1], 0
            x = next(a, sentinel)
        elif x[0] is None or y[0] < x[0]:
            yield y[0], 0, y[1]
            y = next(b, sentinel)
        else:
            yield x[0], x[1], y[1]
            x = next(a, sentinel)
            y = next(b, sentinel)


def step_names(conn: sqlite3.Connection, n_ann: int, sample_size: int, top_k: int,
               seed: int, max_names: int | None) -> dict:
    if not table_exists(conn, "xref_file"):
        REPORT["names"] = {"error": "xref_file missing"}
        save_report()
        return {}
    rng = random.Random(seed)
    reservoir: list[tuple[str, int]] = []
    top_heap: list[tuple[int, str]] = []
    hist = [0] * len(BUCKETS)
    distinct = both = 0
    rows_sym = rows_alias = 0
    name_bytes = 0
    est_blob = est_keys = 0
    t = time.perf_counter()
    last_log = t

    for acc, n_s, n_a in merge_streams(stream_ns(conn, "symbol"), stream_ns(conn, "alias")):
        distinct += 1
        rows_sym += n_s
        rows_alias += n_a
        if n_s and n_a:
            both += 1
        tot = n_s + n_a
        nb = len(acc.encode("utf-8"))
        name_bytes += nb
        est_keys += nb + 12  # key + cell/pointer overhead (rough)
        est_blob += posting_bytes(n_s, n_ann) + posting_bytes(n_a, n_ann) + 2
        for i, (lo, hi) in enumerate(BUCKETS):
            if lo <= tot <= hi:
                hist[i] += 1
                break
        # reservoir sample (algorithm R) and top-k
        if len(reservoir) < sample_size:
            reservoir.append((acc, tot))
        else:
            j = rng.randrange(distinct)
            if j < sample_size:
                reservoir[j] = (acc, tot)
        if len(top_heap) < top_k:
            heapq.heappush(top_heap, (tot, acc))
        elif tot > top_heap[0][0]:
            heapq.heapreplace(top_heap, (tot, acc))

        now = time.perf_counter()
        if now - last_log > 30:
            last_log = now
            log(f"names: {distinct:,} distinct so far, last={acc!r}")
        if max_names and distinct >= max_names:
            log(f"names: stopping early at --max-names={max_names}")
            break

    dt = time.perf_counter() - t
    top = sorted(top_heap, reverse=True)
    total_rows = rows_sym + rows_alias
    hist_rows = [
        [f"{lo}" if lo == hi else (f"{lo}+" if hi >= 10**9 else f"{lo}-{hi}"),
         f"{c:,}", f"{(c / distinct * 100 if distinct else 0):.1f}%"]
        for (lo, hi), c in zip(BUCKETS, hist)
    ]
    est_total = est_keys + est_blob
    out = {
        "annotations": n_ann,
        "rows_symbol": f"{rows_sym:,}",
        "rows_alias": f"{rows_alias:,}",
        "distinct_names": f"{distinct:,}",
        "names_in_both": f"{both:,}",
        "avg_name_bytes": round(name_bytes / distinct, 1) if distinct else 0,
        "scan_seconds": round(dt, 1),
        "histogram": hist_rows,
        "top_names": [(a, c) for c, a in top],
        "estimate": {
            "name_keys_plus_overhead": human(est_keys),
            "posting_blobs": human(est_blob),
            "estimated_name_table_total": human(est_total),
            "note": "excludes B-tree page slack (add roughly 15-30%); varint/bitset sizing "
                    "is approximate, real build will differ",
            "legacy_rows_symbol_alias": f"{total_rows:,}",
        },
        "early_stop": bool(max_names and distinct >= max_names),
    }
    REPORT["names"] = out
    save_report()
    return {"reservoir": reservoir, "top": [(a, c) for c, a in top]}


# --------------------------------------------------------------------------
# step: queries (copies of the current /hits and /names SQL)
# --------------------------------------------------------------------------
def q_postings(conn, names: list[str]):
    out = []
    for chunk in chunks(names, SQLITE_MAX_VARS):
        ph = ",".join("?" for _ in chunk)
        rows = conn.execute(
            f"SELECT namespace, accession, annotation_key FROM xref_file "
            f"WHERE accession IN ({ph}) AND namespace IN ('symbol','alias')",
            chunk,
        ).fetchall()
        out.extend((r[0], r[1], r[2]) for r in rows)
    return out


def q_filter(conn, keys: set[int], include: list[int], exclude: list[int]) -> set[int]:
    if not keys or (not include and not exclude):
        return set(keys)
    key_list = list(keys)
    allowed = set(keys)
    if include:
        allowed = set()
        inc_ph = ",".join("?" for _ in include)
        step = max(1, SQLITE_MAX_VARS - len(include))
        for chunk in chunks(key_list, step):
            ph = ",".join("?" for _ in chunk)
            for r in conn.execute(
                f"SELECT DISTINCT annotation_key FROM annotation_lineage "
                f"WHERE annotation_key IN ({ph}) AND taxid IN ({inc_ph})",
                [*chunk, *include],
            ):
                allowed.add(r[0])
    if exclude and allowed:
        exc_ph = ",".join("?" for _ in exclude)
        step = max(1, SQLITE_MAX_VARS - len(exclude))
        banned = set()
        for chunk in chunks(list(allowed), step):
            ph = ",".join("?" for _ in chunk)
            for r in conn.execute(
                f"SELECT DISTINCT annotation_key FROM annotation_lineage "
                f"WHERE annotation_key IN ({ph}) AND taxid IN ({exc_ph})",
                [*chunk, *exclude],
            ):
                banned.add(r[0])
        allowed -= banned
    return allowed


def q_annotations(conn, keys: list[int]) -> int:
    n = 0
    for chunk in chunks(keys, SQLITE_MAX_VARS):
        ph = ",".join("?" for _ in chunk)
        n += len(conn.execute(
            f"SELECT annotation_key, annotation_id, taxid, source_database, release_date, "
            f"assembly_name FROM annotation WHERE annotation_key IN ({ph})",
            chunk,
        ).fetchall())
    return n


def run_hits(conn, names, include, exclude) -> dict:
    parts: dict[str, float] = {}
    with_timeout(conn, QUERY_TIMEOUT_S)
    try:
        post, parts["postings"] = timed(q_postings, conn, names)
        keys = {p[2] for p in post}
        allowed, parts["tax_filter"] = timed(q_filter, conn, keys, include, exclude)
        cand = [k for k in keys if k in allowed]
        n_rows, parts["annotations"] = timed(q_annotations, conn, cand)
    finally:
        clear_timeout(conn)
    return {"postings": len(post), "hits": n_rows, "parts": parts}


def run_names(conn, prefix: str, include, exclude) -> dict:
    clauses, params = [], []
    if include:
        ph = ",".join("?" for _ in include)
        clauses.append(f"f.annotation_key IN (SELECT annotation_key FROM annotation_lineage "
                       f"WHERE taxid IN ({ph}))")
        params += include
    if exclude:
        ph = ",".join("?" for _ in exclude)
        clauses.append(f"f.annotation_key NOT IN (SELECT annotation_key FROM annotation_lineage "
                       f"WHERE taxid IN ({ph}))")
        params += exclude
    tax_sql = (" AND " + " AND ".join(clauses)) if clauses else ""
    with_timeout(conn, QUERY_TIMEOUT_S)
    try:
        rows = conn.execute(
            f"""SELECT f.accession, COUNT(DISTINCT f.annotation_key) AS n
                FROM xref_file f
                WHERE f.namespace IN ('symbol','alias') AND f.accession >= ? AND f.accession < ?
                {tax_sql}
                GROUP BY f.accession ORDER BY n DESC, f.accession ASC LIMIT 501""",
            (prefix, prefix + "￿", *params),
        ).fetchall()
    finally:
        clear_timeout(conn)
    return {"results": len(rows)}


def bench(fn, reps: int):
    times, res = [], None
    for _ in range(reps):
        t = time.perf_counter()
        try:
            res = fn()
            err = None
        except sqlite3.OperationalError as exc:
            res, err = None, str(exc)
        times.append(time.perf_counter() - t)
        if err:
            return times, res, err
    return times, res, None


def pick_scopes(conn) -> list[tuple[str, list[int], list[int]]]:
    scopes: list[tuple[str, list[int], list[int]]] = [("no scope", [], [])]
    if not table_exists(conn, "annotation_lineage"):
        return scopes
    rows = conn.execute(
        "SELECT taxid, COUNT(*) c FROM annotation_lineage GROUP BY taxid ORDER BY c DESC"
    ).fetchall()
    if not rows:
        return scopes
    targets = [("include big clade", None), ("include ~2000", 2000), ("include ~200", 200),
               ("include ~20", 20)]
    used = set()
    for label, target in targets:
        if target is None:
            tid, c = rows[min(1, len(rows) - 1)]  # skip the root-most taxon
        else:
            tid, c = min(rows, key=lambda r: abs(r[1] - target))
        if tid in used:
            continue
        used.add(tid)
        scopes.append((f"{label} (taxid {tid}, {c} ann)", [int(tid)], []))
    big = rows[min(1, len(rows) - 1)][0]
    scopes.append((f"exclude big clade (taxid {big})", [], [int(big)]))
    return scopes


def step_queries(conn, sample: dict, reps: int, seed: int) -> None:
    rng = random.Random(seed)
    reservoir = [a for a, _ in sample.get("reservoir", [])]
    top = [a for a, _ in sample.get("top", [])]
    pool_common = ["tp53", "brca1", "actb", "gapdh", "myc", "sox2", "il6", "egfr", "hox a1",
                   "shh", "pax6", "foxp2", "cdk1", "mapk1", "tnf", "apoe"]
    common = [n for n in pool_common if n in set(reservoir) | set(top)]
    if not reservoir:
        log("queries: no sampled names; skipping")
        return

    name_sets: list[tuple[str, list[str]]] = []
    if common:
        name_sets.append((f"known common ({len(common)})", common))
    name_sets.append(("top-frequency 10", top[:10]))
    name_sets.append(("random 1", rng.sample(reservoir, 1)))
    name_sets.append(("random 10", rng.sample(reservoir, min(10, len(reservoir)))))
    name_sets.append(("random 50", rng.sample(reservoir, min(50, len(reservoir)))))
    name_sets.append(("random 200", rng.sample(reservoir, min(200, len(reservoir)))))
    name_sets.append(("top-frequency 200 (worst case)", top[:200]))

    scopes = pick_scopes(conn)
    hits_rows = []
    for (set_label, names), (scope_label, inc, exc) in (
        (a, b) for a in name_sets for b in scopes
    ):
        if not names:
            continue
        # keep the matrix manageable: full scope sweep only for 3 representative sets
        if scope_label != "no scope" and not (
            set_label.startswith(("random 10", "random 200", "top-frequency 10"))
        ):
            continue
        log(f"queries: /hits [{set_label}] x [{scope_label}]")
        times, res, err = bench(lambda: run_hits(conn, names, inc, exc), reps)
        if err:
            hits_rows.append([f"{set_label} / {scope_label}", len(names), "-", "-",
                              f"{times[0]:.3f}", "-", f"ERROR {err}"])
            continue
        warm = times[1:] or times[:1]
        hits_rows.append([
            f"{set_label} / {scope_label}", len(names), res["postings"], res["hits"],
            f"{times[0]:.3f}", f"{statistics.median(warm):.3f}",
            ", ".join(f"{k}={v:.3f}" for k, v in res["parts"].items()),
        ])
        save_partial_queries(hits_rows, None)

    names_rows = []
    prefixes = ["sel", "brc", "cdk", "hox", "zn", "il", "mir", "tp5"]
    prefixes += [r[:3] for r in rng.sample(reservoir, min(6, len(reservoir))) if len(r) >= 3]
    seen = set()
    for pfx in prefixes:
        if pfx in seen:
            continue
        seen.add(pfx)
        for scope_label, inc, exc in (scopes[:1] + scopes[1:2]):
            log(f"queries: /names prefix={pfx!r} [{scope_label}]")
            times, res, err = bench(lambda: run_names(conn, pfx, inc, exc), reps)
            if err:
                names_rows.append([f"{pfx} / {scope_label}", "-", f"{times[0]:.3f}", "-",
                                   f"ERROR {err}"])
                continue
            warm = times[1:] or times[:1]
            note = "prefix_too_broad (>500)" if res["results"] > 500 else ""
            names_rows.append([f"{pfx} / {scope_label}", res["results"],
                               f"{times[0]:.3f}", f"{statistics.median(warm):.3f}", note])
            save_partial_queries(hits_rows, names_rows)
    save_partial_queries(hits_rows, names_rows)


def save_partial_queries(hits_rows, names_rows) -> None:
    q = REPORT.setdefault("queries", {"hits_rows": [], "names_rows": []})
    q["hits_rows"] = hits_rows
    if names_rows is not None:
        q["names_rows"] = names_rows
    save_report()


# --------------------------------------------------------------------------
# step: http load test + RSS of the API process tree
# --------------------------------------------------------------------------
def tree_rss_kib(root_pid: int) -> int:
    children: dict[int, list[int]] = {}
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        try:
            with open(f"/proc/{d}/stat") as f:
                ppid = int(f.read().rsplit(")", 1)[1].split()[1])
            children.setdefault(ppid, []).append(int(d))
        except (OSError, ValueError, IndexError):
            continue
    stack, total = [root_pid], 0
    while stack:
        pid = stack.pop()
        stack.extend(children.get(pid, []))
        try:
            with open(f"/proc/{pid}/status") as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        total += int(line.split()[1])
                        break
        except OSError:
            pass
    return total


def step_http(base: str, sample: dict, concurrency: int, requests: int, pid: int | None,
              seed: int) -> None:
    rng = random.Random(seed)
    pool = [a for a, _ in sample.get("reservoir", [])]
    if not pool:
        log("http: no sampled names; skipping")
        return
    base = base.rstrip("/")
    peak = {"kib": 0}
    stop = threading.Event()

    def sampler():
        while not stop.is_set():
            if pid:
                peak["kib"] = max(peak["kib"], tree_rss_kib(pid))
            time.sleep(0.5)

    th = threading.Thread(target=sampler, daemon=True)
    base_rss = tree_rss_kib(pid) if pid else 0
    th.start()

    def one(_):
        names = rng.sample(pool, min(10, len(pool)))
        url = f"{base}/hits?names={urllib.parse.quote(','.join(names))}"
        t = time.perf_counter()
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                r.read()
                code = r.status
        except urllib.error.HTTPError as e:
            code = e.code
        except Exception:
            code = -1
        return code, (time.perf_counter() - t) * 1000

    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=concurrency) as ex:
        results = list(ex.map(one, range(requests)))
    wall = time.perf_counter() - t0
    stop.set()
    th.join(timeout=2)
    lat = [ms for code, ms in results if code == 200]
    codes: dict[int, int] = {}
    for code, _ in results:
        codes[code] = codes.get(code, 0) + 1
    out = {
        "base": base,
        "requests": requests,
        "concurrency": concurrency,
        "wall_seconds": round(wall, 1),
        "status_counts": codes,
        "latency_ms_p50": round(pct(lat, 50), 1) if lat else None,
        "latency_ms_p95": round(pct(lat, 95), 1) if lat else None,
        "latency_ms_max": round(max(lat), 1) if lat else None,
    }
    if pid:
        out["api_rss_before"] = human(base_rss * 1024)
        out["api_rss_peak"] = human(peak["kib"] * 1024)
    REPORT["http"] = out
    save_report()


# --------------------------------------------------------------------------
# step: dbstat
# --------------------------------------------------------------------------
def step_dbstat(conn: sqlite3.Connection, timeout_s: float) -> None:
    try:
        with_timeout(conn, timeout_s)
        rows = conn.execute(
            "SELECT name, SUM(pgsize) AS bytes, COUNT(*) AS pages FROM dbstat "
            "GROUP BY name ORDER BY bytes DESC LIMIT 25"
        ).fetchall()
        REPORT["dbstat"] = [[r[0], human(r[1]), f"{r[2]:,}"] for r in rows]
    except sqlite3.Error as exc:
        REPORT["dbstat"] = {"error": f"dbstat unavailable or timed out: {exc}"}
    finally:
        clear_timeout(conn)
    save_report()


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main() -> int:
    global OUT_DIR
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus-dir", type=Path, default=Path("/data/annotrieve/gene_corpus"))
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="default: ./phase0_report_<timestamp>")
    ap.add_argument("--cache-mib", type=int, default=64, help="SQLite page cache for this script")
    ap.add_argument("--tmp-dir", type=str, default=None,
                    help="local directory for SQLite temp files (avoid NFS)")
    ap.add_argument("--skip", default="", help="comma list of steps to skip: "
                    "env,files,nfs,counts,names,queries,http,dbstat")
    ap.add_argument("--skip-shards", action="store_true", help="do not stat per_gff shards")
    ap.add_argument("--seq-mb", type=int, default=256, help="sequential read size for NFS test")
    ap.add_argument("--rand-reads", type=int, default=300, help="random 16 KiB reads for NFS test")
    ap.add_argument("--max-names", type=int, default=None,
                    help="stop the name scan after N distinct names (quick dry run)")
    ap.add_argument("--sample-size", type=int, default=2000, help="reservoir of sampled names")
    ap.add_argument("--top-k", type=int, default=300, help="most frequent names to keep")
    ap.add_argument("--reps", type=int, default=4, help="repeats per query scenario")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--http-base", default=None,
                    help="e.g. http://127.0.0.1:5000 (the API container, not traefik)")
    ap.add_argument("--api-pid", type=int, default=None,
                    help="PID of the gunicorn master, to sample peak RSS of its process tree")
    ap.add_argument("--http-concurrency", type=int, default=20)
    ap.add_argument("--http-requests", type=int, default=400)
    ap.add_argument("--dbstat-timeout", type=float, default=7200.0)
    args = ap.parse_args()

    skip = {s.strip() for s in args.skip.split(",") if s.strip()}
    OUT_DIR = args.out_dir or Path(f"phase0_report_{datetime.now():%Y%m%d_%H%M%S}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log(f"writing report to {OUT_DIR}")

    db_path = args.corpus_dir / DB_NAME
    if not db_path.is_file():
        log(f"ERROR: {db_path} not found")
        return 2

    def on_signal(signum, _frame):
        log(f"signal {signum}: saving partial report and exiting")
        save_report()
        sys.exit(130)

    signal.signal(signal.SIGINT, on_signal)
    signal.signal(signal.SIGTERM, on_signal)

    REPORT["args"] = {k: str(v) for k, v in vars(args).items()}
    conn = open_db(db_path, args.cache_mib * 1024, args.tmp_dir)

    if "env" not in skip:
        log("step env"); step_env(args.corpus_dir)
    if "files" not in skip:
        log("step files"); step_files(args.corpus_dir, args.skip_shards)
    if "nfs" not in skip:
        log("step nfs"); step_nfs(db_path, args.seq_mb, args.rand_reads)
    n_ann = 0
    if "counts" not in skip:
        log("step counts"); n_ann = step_counts(conn)
    sample: dict = {}
    if "names" not in skip:
        log("step names (streaming scan; this is the long one)")
        if not n_ann and table_exists(conn, "annotation"):
            n_ann = conn.execute("SELECT COUNT(*) FROM annotation").fetchone()[0]
        sample = step_names(conn, n_ann or 19000, args.sample_size, args.top_k, args.seed,
                            args.max_names)
    if "queries" not in skip and sample:
        log("step queries"); step_queries(conn, sample, args.reps, args.seed)
    if "http" not in skip and args.http_base:
        log("step http")
        step_http(args.http_base, sample, args.http_concurrency, args.http_requests,
                  args.api_pid, args.seed)
    if "dbstat" not in skip:
        log("step dbstat (reads the whole file)"); step_dbstat(conn, args.dbstat_timeout)

    save_report()
    log(f"done. See {OUT_DIR / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
