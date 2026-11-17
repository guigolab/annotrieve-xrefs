#!/usr/bin/env python3
"""
Phase 0b checks for the loci design (read-only, standard library only).

Answers three questions before anything is built:

1. TIMING   How long does it take to read all loci of one gene name straight
            from ``gene_hit`` (already clustered by name, with coordinates in
            the row)? If this is fast on a cold cache, ``/loci`` may not need
            a separate ``loci.sqlite``. Runs FIRST so the cache is as cold as
            we can get without root.
2. CASE     Do ``gene_hit.primary_name`` values differ in case from the
            casefolded ``accession``? (They do for symbols, so the served
            table needs a display-name column.)
3. CONSIST  Do ``xref_file`` and ``gene_hit`` agree on which (name, annotation)
            pairs exist for symbol and alias? Streams both tables in primary-key
            order and merges them in constant memory.

    nohup python3 tools/phase0b_check.py \
        --corpus-dir /data/annotrieve/gene_corpus \
        --out-dir /tmp/phase0b > /tmp/phase0b.log 2>&1 &

Section 3 compares 30 random windows by default (``--consist-sample``); the full scan
(``--consist-sample 0``) took hours cold in testing. Skip a section with
``--skip timing,case,consist``. The report (report.md / report.json) is
rewritten after every section; SIGINT / SIGTERM save a partial report.

Run on a published snapshot, or when no sync/merge is writing the DB.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import signal
import sqlite3
import statistics
import sys
import time
from datetime import datetime
from pathlib import Path

DB_NAME = "gene_corpus.sqlite"
QUERY_TIMEOUT_S = 600.0
DEFAULT_NAMES = ["tp53", "brca1", "actb", "rps2", "sec13", "rps3", "nop58", "cct2", "sod1"]

REPORT: dict = {}
OUT_DIR: Path | None = None
T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:8.1f}s] {msg}", file=sys.stderr, flush=True)


def human(n: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(n) < 1024 or unit == "TiB":
            return f"{int(n)} B" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return str(n)


def pct(values: list[float], p: float) -> float:
    if not values:
        return float("nan")
    s = sorted(values)
    k = min(len(s) - 1, max(0, int(-(-p / 100 * len(s) // 1)) - 1))
    return s[k]


def open_db(path: Path, cache_mib: int) -> sqlite3.Connection:
    conn = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
    conn.execute(f"PRAGMA cache_size=-{cache_mib * 1024}")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute("PRAGMA mmap_size=0")
    conn.execute("PRAGMA query_only=ON")
    return conn


def set_timeout(conn: sqlite3.Connection, seconds: float) -> None:
    deadline = time.monotonic() + seconds
    conn.set_progress_handler(lambda: 1 if time.monotonic() > deadline else 0, 100_000)


def clear_timeout(conn: sqlite3.Connection) -> None:
    conn.set_progress_handler(None, 0)


def md_table(headers: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------
def render_md() -> str:
    L = ["# Phase 0b checks", "",
         f"Generated {datetime.now().isoformat(timespec='seconds')} "
         f"(elapsed {REPORT.get('elapsed_s', 0)} s)", ""]
    t = REPORT.get("timing")
    if t:
        L += ["## 1. Cold range read of one name in gene_hit", "",
              "`first` is the first read of that name; `repeat` is the second (cached). "
              "Rows carry coordinates in-row, so this is roughly what `/loci` would cost "
              "without a separate file. Names come first from the fixed list, then from "
              "random positions in the table (probably cold).", ""]
        L.append(md_table(
            ["name", "kind", "rows", "annotations", "payload", "first (s)", "repeat (s)",
             "rows/s first"], t["rows"]))
        L += ["", f"- random-name read, median first-read: {t['random_median_s']} s "
                  f"(p95 {t['random_p95_s']} s) over {t['random_n']} names",
              f"- average row payload (text fields plus about 24 bytes of integers): "
              f"{t['avg_payload_bytes']} bytes", ""]
    c = REPORT.get("case")
    if c:
        L += ["## 2. primary_name vs accession (case)", ""]
        L += [f"- {k}: {v}" for k, v in c.items()]
        L.append("")
    k = REPORT.get("consist")
    if k:
        L += ["## 3. xref_file vs gene_hit consistency", ""]
        for ns, r in k.items():
            L.append(f"### namespace `{ns}`")
            L.append("")
            L += [f"- {a}: {b}" for a, b in r["summary"].items()]
            if r.get("by_annotation"):
                L += ["", "Annotations with the most mismatches:", ""]
                L.append(md_table(["annotation_key", "only in xref_file", "only in gene_hit"],
                                  r["by_annotation"]))
            if r.get("samples"):
                L += ["", "Samples (name, annotation_key):", ""]
                L += [f"- {s}" for s in r["samples"]]
            L.append("")
        L += ["Reading: a few mismatches may be harmless leftovers; a mismatch share above "
              "about 0.1% means `/hits` and `/loci` would disagree and the build needs a "
              "single source or a repair first.", ""]
    return "\n".join(L)


def save_report() -> None:
    if OUT_DIR is None:
        return
    REPORT["elapsed_s"] = round(time.time() - T0, 1)
    (OUT_DIR / "report.json").write_text(json.dumps(REPORT, indent=2, default=str))
    (OUT_DIR / "report.md").write_text(render_md())


# --------------------------------------------------------------------------
# 1. timing
# --------------------------------------------------------------------------
LOCI_SQL = (
    "SELECT annotation_key, local_id, seqid, start, end, strand, feature_type, biotype, "
    "primary_name FROM gene_hit WHERE namespace = ? AND accession = ?"
)


def read_name(conn: sqlite3.Connection, ns: str, name: str) -> tuple[list, float]:
    set_timeout(conn, QUERY_TIMEOUT_S)
    try:
        t = time.perf_counter()
        rows = conn.execute(LOCI_SQL, (ns, name)).fetchall()
        return rows, time.perf_counter() - t
    finally:
        clear_timeout(conn)


def payload_bytes(row) -> int:
    return 24 + sum(len(str(row[i] or "")) for i in (2, 6, 7, 8))


def random_names(conn: sqlite3.Connection, n: int, seed: int) -> list[str]:
    rng = random.Random(seed)
    alphabet = "abcdefghijklmnopqrstuvwxyz0123456789"
    out: list[str] = []
    tries = 0
    while len(out) < n and tries < n * 10:
        tries += 1
        stem = "".join(rng.choice(alphabet) for _ in range(rng.choice((2, 3, 4))))
        row = conn.execute(
            "SELECT accession FROM gene_hit WHERE namespace='symbol' AND accession >= ? LIMIT 1",
            (stem,),
        ).fetchone()
        if row and row[0] not in out:
            out.append(row[0])
    return out


def step_timing(conn: sqlite3.Connection, names: list[str], n_random: int, seed: int) -> None:
    rows_out: list[list] = []
    payloads: list[int] = []
    rand_first: list[float] = []

    def one(name: str, kind: str) -> None:
        rows, t1 = read_name(conn, "symbol", name)
        _, t2 = read_name(conn, "symbol", name)
        anns = len({r[0] for r in rows})
        pay = [payload_bytes(r) for r in rows[:2000]]
        payloads.extend(pay)
        size = sum(payload_bytes(r) for r in rows) if len(rows) <= 2000 else int(
            statistics.mean(pay) * len(rows)) if pay else 0
        rows_out.append([name, kind, f"{len(rows):,}", f"{anns:,}", human(size),
                         f"{t1:.3f}", f"{t2:.3f}", f"{len(rows) / t1:,.0f}" if t1 else "-"])
        if kind == "random":
            rand_first.append(t1)
        log(f"timing: {kind} {name!r}: {len(rows):,} rows first={t1:.3f}s repeat={t2:.3f}s")
        REPORT["timing"] = {
            "rows": rows_out,
            "random_median_s": round(statistics.median(rand_first), 3) if rand_first else "-",
            "random_p95_s": round(pct(rand_first, 95), 3) if rand_first else "-",
            "random_n": len(rand_first),
            "avg_payload_bytes": round(statistics.mean(payloads), 1) if payloads else 0,
        }
        save_report()

    # random (likely cold) names first, then the fixed list
    for name in random_names(conn, n_random, seed):
        one(name, "random")
    for name in names:
        one(name, "listed")


# --------------------------------------------------------------------------
# 2. case
# --------------------------------------------------------------------------
def step_case(conn: sqlite3.Connection, sample_names: int, seed: int) -> None:
    out: dict = {}
    for ns in ("symbol", "alias"):
        tot = same_exact = same_fold = differs = 0
        examples: list[str] = []
        names = random_names(conn, sample_names, seed) if ns == "symbol" else []
        if ns == "alias":
            rng = random.Random(seed)
            for _ in range(sample_names * 5):
                stem = "".join(rng.choice("abcdefghijklmnopqrstuvwxyz") for _ in range(3))
                r = conn.execute(
                    "SELECT accession FROM gene_hit WHERE namespace='alias' AND accession >= ? "
                    "LIMIT 1", (stem,)).fetchone()
                if r and r[0] not in names:
                    names.append(r[0])
                if len(names) >= sample_names:
                    break
        for name in names:
            for acc, pn in conn.execute(
                "SELECT accession, primary_name FROM gene_hit WHERE namespace=? AND accession=? "
                "LIMIT 500", (ns, name)):
                tot += 1
                if pn is None:
                    differs += 1
                elif pn == acc:
                    same_exact += 1
                elif pn.casefold() == acc.casefold():
                    same_fold += 1
                    if len(examples) < 8:
                        examples.append(f"{acc} -> {pn}")
                else:
                    differs += 1
                    if len(examples) < 8:
                        examples.append(f"{acc} -> {pn}  (different name)")
        out[f"{ns}: rows sampled"] = tot
        if tot:
            out[f"{ns}: primary_name identical to accession"] = f"{same_exact / tot:.1%}"
            out[f"{ns}: same name, different case"] = f"{same_fold / tot:.1%}"
            out[f"{ns}: primary_name is another name / NULL"] = f"{differs / tot:.1%}"
            out[f"{ns}: examples"] = "; ".join(examples)
    out["conclusion"] = ("If 'different case' is large for symbols, the served table needs a "
                         "display-name column (display_id). If 'another name' is large for "
                         "aliases, the alias row's own symbol must be kept as well.")
    REPORT["case"] = out
    save_report()


# --------------------------------------------------------------------------
# 3. consistency
# --------------------------------------------------------------------------
def stream_xref(conn, ns):
    cur = conn.execute(
        "SELECT accession, annotation_key FROM xref_file WHERE namespace = ? "
        "ORDER BY accession, annotation_key", (ns,))
    yield from cur


def stream_hit(conn, ns):
    cur = conn.execute(
        "SELECT accession, annotation_key, COUNT(*) FROM gene_hit WHERE namespace = ? "
        "GROUP BY accession, annotation_key ORDER BY accession, annotation_key", (ns,))
    for acc, ak, n in cur:
        yield (acc, ak), n


def step_consist(conn: sqlite3.Connection, namespaces: list[str], max_pairs: int | None) -> None:
    out: dict = {}
    for ns in namespaces:
        log(f"consist: namespace {ns}")
        a = ((acc, ak) for acc, ak in stream_xref(conn, ns))
        b = stream_hit(conn, ns)
        only_x = only_h = both = n_x = n_h = multi = loci_rows = 0
        by_ann_x: dict[int, int] = {}
        by_ann_h: dict[int, int] = {}
        samples: list[str] = []
        last = time.perf_counter()
        sentinel = None
        x = next(a, sentinel)
        h = next(b, sentinel)
        while x is not None or h is not None:
            if max_pairs and (n_x + n_h) >= max_pairs:
                log("consist: stopping early (--max-pairs)")
                break
            if h is None or (x is not None and x < h[0]):
                only_x += 1
                n_x += 1
                by_ann_x[x[1]] = by_ann_x.get(x[1], 0) + 1
                if len(samples) < 25:
                    samples.append(f"only in xref_file: {x[0]!r}, {x[1]}")
                x = next(a, sentinel)
            elif x is None or h[0] < x:
                only_h += 1
                n_h += 1
                by_ann_h[h[0][1]] = by_ann_h.get(h[0][1], 0) + 1
                if len(samples) < 25:
                    samples.append(f"only in gene_hit: {h[0][0]!r}, {h[0][1]}")
                h = next(b, sentinel)
            else:
                both += 1
                n_x += 1
                n_h += 1
                loci_rows += h[1]
                if h[1] > 1:
                    multi += 1
                x = next(a, sentinel)
                h = next(b, sentinel)
            now = time.perf_counter()
            if now - last > 30:
                last = now
                log(f"consist[{ns}]: both={both:,} only_xref={only_x:,} only_hit={only_h:,}")
        total = both + only_x + only_h
        worst = sorted(set(by_ann_x) | set(by_ann_h),
                       key=lambda k: -(by_ann_x.get(k, 0) + by_ann_h.get(k, 0)))[:15]
        out[ns] = {
            "summary": {
                "pairs in both": f"{both:,}",
                "pairs only in xref_file": f"{only_x:,}",
                "pairs only in gene_hit": f"{only_h:,}",
                "mismatch share": f"{((only_x + only_h) / total if total else 0):.4%}",
                "gene_hit rows for matching pairs": f"{loci_rows:,}",
                "pairs with more than one locus": f"{multi:,}",
                "early_stop": bool(max_pairs and (n_x + n_h) >= max_pairs),
            },
            "by_annotation": [[k, by_ann_x.get(k, 0), by_ann_h.get(k, 0)] for k in worst
                              if by_ann_x.get(k, 0) or by_ann_h.get(k, 0)],
            "samples": samples,
        }
        REPORT["consist"] = out
        save_report()


def step_consist_sample(conn: sqlite3.Connection, namespaces: list[str], n_ranges: int,
                        range_pairs: int, seed: int) -> None:
    """Compare xref_file and gene_hit on random name windows (cheap, cold-read friendly)."""
    out: dict = {}
    rng = random.Random(seed)
    alphabet = "abcdefghijklmnopqrstuvwxyz0123456789"
    for ns in namespaces:
        only_x = only_h = both = windows = 0
        by_x: dict[int, int] = {}
        by_h: dict[int, int] = {}
        samples: list[str] = []
        for i in range(n_ranges):
            stem = "".join(rng.choice(alphabet) for _ in range(rng.choice((2, 3, 4))))
            set_timeout(conn, QUERY_TIMEOUT_S)
            try:
                rows = conn.execute(
                    "SELECT accession, annotation_key FROM xref_file WHERE namespace=? "
                    "AND accession >= ? ORDER BY accession, annotation_key LIMIT ?",
                    (ns, stem, range_pairs)).fetchall()
                if len(rows) < 2:
                    continue
                lo, hi = rows[0][0], rows[-1][0]
                xs = {(a, k) for a, k in rows if a != hi}
                hs = {(a, k) for a, k, _n in conn.execute(
                    "SELECT accession, annotation_key, COUNT(*) FROM gene_hit WHERE namespace=? "
                    "AND accession >= ? AND accession < ? GROUP BY accession, annotation_key",
                    (ns, lo, hi))}
            finally:
                clear_timeout(conn)
            windows += 1
            ox, oh = xs - hs, hs - xs
            both += len(xs & hs)
            only_x += len(ox)
            only_h += len(oh)
            for a, k in ox:
                by_x[k] = by_x.get(k, 0) + 1
                if len(samples) < 25:
                    samples.append(f"only in xref_file: {a!r}, {k}")
            for a, k in oh:
                by_h[k] = by_h.get(k, 0) + 1
                if len(samples) < 25:
                    samples.append(f"only in gene_hit: {a!r}, {k}")
            log(f"consist-sample[{ns}] window {i + 1}/{n_ranges} from {lo!r}: "
                f"both={len(xs & hs):,} only_xref={len(ox):,} only_hit={len(oh):,}")
        total = both + only_x + only_h
        worst = sorted(set(by_x) | set(by_h), key=lambda k: -(by_x.get(k, 0) + by_h.get(k, 0)))[:15]
        out[ns] = {
            "summary": {
                "windows sampled": windows,
                "pairs per window": range_pairs,
                "pairs in both": f"{both:,}",
                "pairs only in xref_file": f"{only_x:,}",
                "pairs only in gene_hit": f"{only_h:,}",
                "mismatch share": f"{((only_x + only_h) / total if total else 0):.4%}",
            },
            "by_annotation": [[k, by_x.get(k, 0), by_h.get(k, 0)] for k in worst],
            "samples": samples,
        }
        REPORT["consist"] = out
        save_report()


# --------------------------------------------------------------------------
def main() -> int:
    global OUT_DIR
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus-dir", type=Path, default=Path("/data/annotrieve/gene_corpus"))
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--cache-mib", type=int, default=64)
    ap.add_argument("--tmp-dir", default=None, help="local dir for SQLite temp files")
    ap.add_argument("--skip", default="", help="comma list: timing,case,consist")
    ap.add_argument("--names", default=",".join(DEFAULT_NAMES),
                    help="comma list of names for the timed reads")
    ap.add_argument("--random-names", type=int, default=20)
    ap.add_argument("--case-sample", type=int, default=40)
    ap.add_argument("--namespaces", default="symbol,alias")
    ap.add_argument("--max-pairs", type=int, default=None,
                    help="stop each consistency scan after N pairs (quick dry run)")
    ap.add_argument("--consist-sample", type=int, default=30,
                    help="number of random name windows to compare (default 30); "
                         "0 = full scan, which takes hours cold")
    ap.add_argument("--range-pairs", type=int, default=20000, help="pairs per sampled window")
    ap.add_argument("--evict-cache", action="store_true",
                    help="posix_fadvise(DONTNEED) on the DB file before the timing section "
                         "(advisory, no root needed)")
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()

    skip = {s.strip() for s in args.skip.split(",") if s.strip()}
    OUT_DIR = args.out_dir or Path(f"phase0b_report_{datetime.now():%Y%m%d_%H%M%S}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.tmp_dir:
        os.environ["SQLITE_TMPDIR"] = args.tmp_dir

    db = args.corpus_dir / DB_NAME
    if not db.is_file():
        log(f"ERROR: {db} not found")
        return 2

    def on_signal(signum, _f):
        log(f"signal {signum}: saving partial report")
        save_report()
        sys.exit(130)

    signal.signal(signal.SIGINT, on_signal)
    signal.signal(signal.SIGTERM, on_signal)

    REPORT["args"] = {k: str(v) for k, v in vars(args).items()}
    conn = open_db(db, args.cache_mib)
    names = [n.strip().casefold() for n in args.names.split(",") if n.strip()]

    if args.evict_cache:
        try:
            fd = os.open(db, os.O_RDONLY)
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
            os.close(fd)
            log("evicted the DB file from the page cache (advisory)")
        except (OSError, AttributeError) as exc:
            log(f"evict-cache failed: {exc}")

    if "timing" not in skip:
        log("section 1: timed range reads from gene_hit (run first, cache coldest)")
        step_timing(conn, names, args.random_names, args.seed)
    if "case" not in skip:
        log("section 2: primary_name vs accession")
        step_case(conn, args.case_sample, args.seed)
    if "consist" not in skip:
        nss = [n.strip() for n in args.namespaces.split(",") if n.strip()]
        if args.consist_sample > 0:
            log("section 3: sampled xref_file vs gene_hit consistency")
            step_consist_sample(conn, nss, args.consist_sample, args.range_pairs, args.seed)
        else:
            log("section 3: FULL xref_file vs gene_hit consistency (hours when cold)")
            step_consist(conn, nss, args.max_pairs)

    save_report()
    log(f"done. See {OUT_DIR / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
