#!/usr/bin/env python3
"""
M2 gate: measure sidecar writing on a random sample of real shards.

Read-only on the corpus. Sidecars are written to --out-dir (NOT next to the shards),
so nothing under per_gff/ is modified. Put --out-dir and --report-dir on the NFS
work area, not on local disk.

    python3 tools/m2_measure_sidecars.py \
        --per-gff /data/annotrieve/gene_corpus/per_gff \
        --out-dir /data/annotrieve/gene_corpus/m2_sidecars \
        --report-dir /data/annotrieve/gene_corpus/m2_report \
        --sample 200 --largest 5 --workers 4 --seed 1

What it reports (see PLAN-build-tables.md section 8, M2):
  * time per shard and throughput; peak memory of the writer
  * sidecar bytes (names + genes) and the ratio to shard bytes
  * projected full backfill time and total sidecar size
  * D4: how many names changed under casefold (are shards already casefolded?)
  * D9: the source_database / source_provider values actually present
  * validation 1: sidecar names/genes equal a fresh query of the shard

Run from the repo root (annotrieve-xrefs) so `gene_corpus` is importable.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import resource
import sqlite3
import statistics
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gene_corpus.sql.manifest import classify_source  # noqa: E402
from gene_corpus.sql.sidecar import (  # noqa: E402
    SHARD_NAME,
    iter_genes,
    iter_names,
    write_sidecar,
)


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def evict(path: Path) -> None:
    """Ask the kernel to drop cached pages of the shard (advisory, esp. on NFS)."""
    try:
        fd = os.open(path, os.O_RDONLY)
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)
    except (OSError, AttributeError):
        pass


def list_shards(per_gff: Path) -> list[tuple[str, int]]:
    out = []
    with os.scandir(per_gff) as it:
        for entry in it:
            if entry.is_dir(follow_symlinks=False):
                p = Path(entry.path) / SHARD_NAME
                try:
                    out.append((entry.name, p.stat().st_size))
                except OSError:
                    continue
    return sorted(out)


FULL_VERIFY_MAX_ROWS = 1_000_000
SPOT_CHECK = 2000


def verify_sidecar(shard: Path, out: Path, n_skipped: int) -> dict:
    """
    Validation 1 (sidecar vs a fresh query of the shard), bounded in memory.

    Small shards: full set comparison. Big shards: gene rows are compared by streaming,
    and ``SPOT_CHECK`` random name rows must be present in the sidecar.
    """
    conn = sqlite3.connect(f"file:{shard}?mode=ro", uri=True)
    try:
        n_src = conn.execute(
            "SELECT COUNT(*) FROM gene_xref WHERE namespace IN ('symbol','alias')"
        ).fetchone()[0]
        full = n_src <= FULL_VERIFY_MAX_ROWS
        # genes: streamed side by side
        genes_ok = True
        # the sidecar holds only genes that have a symbol/alias row
        gcur = conn.execute(
            "SELECT local_id, start, end, strand FROM gene WHERE local_id IN "
            "(SELECT local_id FROM gene_xref WHERE namespace IN ('symbol','alias') "
            "AND accession IS NOT NULL) ORDER BY local_id"
        )
        n_got = 0
        for fresh, got in zip(gcur, iter_genes(out)):
            n_got += 1
            if tuple(fresh) != (got[0], got[3], got[4], got[5]):
                genes_ok = False
                break
        else:
            n_fresh = conn.execute(
                "SELECT COUNT(DISTINCT local_id) FROM gene_xref WHERE namespace IN "
                "('symbol','alias') AND accession IS NOT NULL").fetchone()[0]
            n_side = sum(1 for _ in iter_genes(out))
            # unsafe names can drop a gene whose only names were skipped
            genes_ok = n_side == n_fresh if n_skipped == 0 else n_side <= n_fresh
        if full:
            fresh_names = {
                (a.casefold(), "s" if ns == "symbol" else "a", int(lid))
                for a, ns, lid in conn.execute(
                    "SELECT accession, namespace, local_id FROM gene_xref "
                    "WHERE namespace IN ('symbol','alias') AND accession IS NOT NULL"
                )
            }
            got_names = set(iter_names(out))
            missing = len(fresh_names - got_names)
            extra = len(got_names - fresh_names)
            ok = genes_ok and extra == 0 and missing <= n_skipped
            return {"verify_ok": ok, "verify_mode": "full",
                    "verify_detail": None if ok else
                    {"missing": missing, "extra": extra, "skipped": n_skipped, "genes_ok": genes_ok}}
        # spot check: random rows by rowid
        mx = conn.execute("SELECT MAX(rowid) FROM gene_xref").fetchone()[0] or 0
        rng = random.Random(7)
        ids = [rng.randint(1, mx) for _ in range(SPOT_CHECK * 3)]
        targets = set()
        for i in range(0, len(ids), 500):
            chunk = ids[i:i + 500]
            q = ",".join("?" for _ in chunk)
            for a, ns, lid in conn.execute(
                f"SELECT accession, namespace, local_id FROM gene_xref WHERE rowid IN ({q}) "
                "AND namespace IN ('symbol','alias') AND accession IS NOT NULL", chunk):
                nm = a.casefold()
                if nm == nm.translate(str.maketrans({"\t": " ", "\n": " ", "\r": " ", "\0": ""})):
                    targets.add((nm, "s" if ns == "symbol" else "a", int(lid)))
            if len(targets) >= SPOT_CHECK:
                break
        found = 0
        remaining = set(targets)
        for row in iter_names(out):
            if row in remaining:
                remaining.discard(row)
                found += 1
                if not remaining:
                    break
        ok = genes_ok and not remaining
        return {"verify_ok": ok, "verify_mode": f"spot({len(targets)})",
                "verify_detail": None if ok else
                {"not_found": len(remaining), "genes_ok": genes_ok}}
    finally:
        conn.close()


def measure_one(args: tuple[str, str, str, bool, bool]) -> dict:
    annotation_dir, per_gff, out_dir, verify, do_evict = args
    shard = Path(per_gff) / annotation_dir / SHARD_NAME
    out = Path(out_dir) / f"{annotation_dir}.sidecar.tsv.gz"
    rec: dict = {"dir": annotation_dir, "shard_bytes": shard.stat().st_size}
    try:
        if do_evict:
            evict(shard)
        t0 = time.perf_counter()
        info = write_sidecar(shard, out, tmp_dir=Path(out_dir) / "_tmp")
        rec["seconds"] = time.perf_counter() - t0
        rec.update(
            n_genes=info.n_genes,
            n_genes_named=info.n_genes_named,
            n_name_rows=info.n_name_rows,
            n_distinct_names=info.n_distinct_names,
            n_casefold_changed=info.n_casefold_changed,
            n_skipped_names=info.n_skipped_names,
            sidecar_bytes=info.size_bytes,
            source_database=info.meta.get("source_database", ""),
            source_provider=info.meta.get("source_provider", ""),
        )
        if verify:
            rec.update(verify_sidecar(shard, out, info.n_skipped_names))
    except Exception as exc:  # report, do not stop the sample
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["maxrss_mib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    return rec


def run_isolated(w: tuple) -> dict:
    """Run one shard in its own process so a killed worker (OOM) names the culprit."""
    try:
        with ProcessPoolExecutor(max_workers=1) as ex:
            return ex.submit(measure_one, w).result()
    except BrokenProcessPool:
        return {"dir": w[0], "shard_bytes": Path(w[1], w[0], SHARD_NAME).stat().st_size,
                "error": "worker killed (most likely out of memory)"}


def fmt_bytes(n: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(n) < 1024 or unit == "TiB":
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TiB"


def fmt_secs(s: float) -> str:
    return f"{s / 3600:.1f} h" if s >= 3600 else (f"{s / 60:.1f} min" if s >= 60 else f"{s:.1f} s")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-gff", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True, help="sidecars go here (NFS work area)")
    ap.add_argument("--report-dir", type=Path, required=True)
    ap.add_argument("--sample", type=int, default=200)
    ap.add_argument("--largest", type=int, default=5, help="also measure the N largest shards")
    ap.add_argument("--workers", type=int, default=1, help="parallel processes")
    ap.add_argument("--assume-workers", type=int, default=8,
                    help="worker count used for the projected backfill time")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-verify", action="store_true")
    ap.add_argument("--resume", action="store_true",
                    help="continue from report-dir/results.jsonl, skipping measured shards")
    ap.add_argument("--no-evict", action="store_true", help="do not drop shard pages before reading")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    log("listing shards (stat on every shard file)")
    t0 = time.perf_counter()
    shards = list_shards(args.per_gff)
    if not shards:
        print("no shards found", file=sys.stderr)
        return 2
    total_bytes = sum(b for _d, b in shards)
    log(f"{len(shards)} shards, {fmt_bytes(total_bytes)} in {time.perf_counter() - t0:.1f}s")

    rng = random.Random(args.seed)
    by_size = sorted(shards, key=lambda x: -x[1])
    picked = {d for d, _b in by_size[: args.largest]}
    rest = [d for d, _b in shards if d not in picked]
    picked |= set(rng.sample(rest, min(args.sample, len(rest))))
    work = [(d, str(args.per_gff), str(args.out_dir), not args.no_verify, not args.no_evict)
            for d in sorted(picked)]
    log(f"measuring {len(work)} shards ({args.largest} largest + random), workers={args.workers}")

    jsonl = args.report_dir / "results.jsonl"
    results: list[dict] = []
    done: set[str] = set()
    if args.resume and jsonl.exists():
        for line in jsonl.read_text().splitlines():
            r = json.loads(line)
            results.append(r)
            done.add(r["dir"])
        log(f"resuming: {len(done)} shards already measured")
    todo = [w for w in work if w[0] not in done]
    out_fh = open(jsonl, "a" if args.resume else "w")

    def keep(rec: dict) -> None:
        results.append(rec)
        out_fh.write(json.dumps(rec) + "\n")
        out_fh.flush()

    wall0 = time.perf_counter()
    if args.workers <= 1:
        for i, w in enumerate(todo, 1):
            keep(measure_one(w))
            if i % 20 == 0 or i == len(todo):
                log(f"  {i}/{len(todo)}")
    else:
        finished: set[str] = set()
        try:
            with ProcessPoolExecutor(max_workers=args.workers) as ex:
                futs = {ex.submit(measure_one, w): w for w in todo}
                for i, f in enumerate(as_completed(futs), 1):
                    keep(f.result())
                    finished.add(futs[f][0])
                    if i % 20 == 0 or i == len(todo):
                        log(f"  {i}/{len(todo)}")
        except BrokenProcessPool:
            left = [w for w in todo if w[0] not in finished]
            log(f"a worker was killed (likely out of memory); re-running {len(left)} "
                "unfinished shards one at a time to find the culprit")
            for w in left:
                keep(run_isolated(w))
    out_fh.close()
    wall = time.perf_counter() - wall0

    ok = [r for r in results if "error" not in r]
    errors = [r for r in results if "error" in r]
    if not ok:
        print("every shard failed; see report", file=sys.stderr)

    # ---- aggregates -----------------------------------------------------
    sample_shard_bytes = sum(r["shard_bytes"] for r in ok)
    sample_side_bytes = sum(r["sidecar_bytes"] for r in ok)
    sample_secs = sum(r["seconds"] for r in ok)
    n_ok = max(len(ok), 1)
    ratio_side = sample_side_bytes / sample_shard_bytes if sample_shard_bytes else 0.0
    secs_per_byte = sample_secs / sample_shard_bytes if sample_shard_bytes else 0.0
    secs = [r["seconds"] for r in ok] or [0.0]

    proj_total_side = ratio_side * total_bytes
    # random-only estimate is less biased than size-weighted when the sample is mixed;
    # report both the per-shard mean and the per-byte rate.
    proj_secs_bytes = secs_per_byte * total_bytes
    proj_secs_mean = statistics.mean(secs) * len(shards)

    changed_shards = sum(1 for r in ok if r["n_casefold_changed"] > 0)
    changed_rows = sum(r["n_casefold_changed"] for r in ok)
    name_rows = sum(r["n_name_rows"] for r in ok)
    skipped = sum(r["n_skipped_names"] for r in ok)
    bad_verify = [r for r in ok if r.get("verify_ok") is False]

    src_db = Counter(r["source_database"] for r in ok)
    src_prov = Counter(r["source_provider"] for r in ok)
    classes = Counter(classify_source(k) or "indexed" for k in src_db for _ in range(src_db[k]))

    report = {
        "n_shards_total": len(shards),
        "total_shard_bytes": total_bytes,
        "n_measured": len(results),
        "n_errors": len(errors),
        "wall_seconds": wall,
        "workers": args.workers,
        "sample_shard_bytes": sample_shard_bytes,
        "sample_sidecar_bytes": sample_side_bytes,
        "sidecar_to_shard_ratio": ratio_side,
        "seconds_per_shard_mean": statistics.mean(secs),
        "seconds_per_shard_p50": statistics.median(secs),
        "seconds_per_shard_p95": sorted(secs)[int(0.95 * (len(secs) - 1))],
        "seconds_per_shard_max": max(secs),
        "projected_total_sidecar_bytes": proj_total_side,
        "projected_single_process_seconds_by_bytes": proj_secs_bytes,
        "projected_single_process_seconds_by_mean": proj_secs_mean,
        "assume_workers": args.assume_workers,
        "projected_wall_seconds_by_bytes": proj_secs_bytes / args.assume_workers,
        "projected_wall_seconds_by_mean": proj_secs_mean / args.assume_workers,
        "peak_maxrss_mib": max((r["maxrss_mib"] for r in results), default=0.0),
        "name_rows": name_rows,
        "casefold_changed_rows": changed_rows,
        "casefold_changed_shards": changed_shards,
        "skipped_names": skipped,
        "verify_failures": len(bad_verify),
        "source_database": dict(src_db),
        "source_provider": dict(src_prov),
        "source_classification": dict(classes),
        "slowest": sorted(ok, key=lambda r: -r["seconds"])[:5],
        "largest": sorted(ok, key=lambda r: -r["shard_bytes"])[:5],
        "errors": errors[:20],
        "verify_failed_examples": bad_verify[:10],
        "evicted_cache": not args.no_evict,
    }
    (args.report_dir / "m2_report.json").write_text(json.dumps(report, indent=2))

    # ---- markdown -------------------------------------------------------
    d4 = ("shards ARE already casefolded (no value changed)" if changed_rows == 0
          else f"{changed_rows:,} of {name_rows:,} name rows ({changed_rows / max(name_rows, 1):.2%}) "
               f"changed under casefold, in {changed_shards} of {len(ok)} shards: shards are NOT all casefolded")
    lines = [
        "# M2 report: sidecar measurement",
        "",
        f"- shards on disk: **{len(shards):,}**, total **{fmt_bytes(total_bytes)}**",
        f"- measured: {len(results)} ({len(errors)} errors), workers={args.workers}, wall {fmt_secs(wall)}"
        f", page cache eviction {'on (advisory)' if not args.no_evict else 'off'}",
        f"- peak writer memory (max RSS of any process): **{report['peak_maxrss_mib']:.0f} MiB**",
        "",
        "## Time",
        "| metric | value |",
        "|---|---|",
        f"| per shard mean / p50 / p95 / max | {fmt_secs(report['seconds_per_shard_mean'])} / "
        f"{fmt_secs(report['seconds_per_shard_p50'])} / {fmt_secs(report['seconds_per_shard_p95'])} / "
        f"{fmt_secs(report['seconds_per_shard_max'])} |",
        f"| projected single-process backfill (by bytes) | {fmt_secs(proj_secs_bytes)} |",
        f"| projected single-process backfill (by mean x count) | {fmt_secs(proj_secs_mean)} |",
        f"| projected wall with {args.assume_workers} workers (by bytes / by mean) | "
        f"{fmt_secs(report['projected_wall_seconds_by_bytes'])} / {fmt_secs(report['projected_wall_seconds_by_mean'])} |",
        "",
        "## Size",
        "| metric | value |",
        "|---|---|",
        f"| sample shards / sample sidecars | {fmt_bytes(sample_shard_bytes)} / {fmt_bytes(sample_side_bytes)} |",
        f"| sidecar / shard ratio | {ratio_side:.4f} |",
        f"| projected total sidecar size (all shards) | **{fmt_bytes(proj_total_side)}** |",
        f"| mean sidecar | {fmt_bytes(sample_side_bytes / n_ok)} |",
        "",
        "## D4: are shard names casefolded?",
        d4,
        f"skipped unsafe names (tab/newline/empty): {skipped:,}",
        "",
        "## D9: source_database / source_provider seen",
        "| source_database | shards | status the manifest would give |",
        "|---|---|---|",
        *[f"| `{k}` | {v} | {classify_source(k) or 'indexed'} |" for k, v in src_db.most_common()],
        "",
        "| source_provider | shards |",
        "|---|---|",
        *[f"| `{k}` | {v} |" for k, v in src_prov.most_common(30)],
        "",
        "## Validation (sidecar vs fresh shard query)",
        f"- failures: **{len(bad_verify)}** of {sum(1 for r in ok if 'verify_ok' in r)} verified",
        *[f"  - `{r['dir']}`: {r.get('verify_detail')}" for r in bad_verify[:10]],
        "",
        "## Slowest shards",
        "| dir | shard | seconds | genes (named) | name rows |",
        "|---|---|---|---|---|",
        *[f"| `{r['dir']}` | {fmt_bytes(r['shard_bytes'])} | {r['seconds']:.1f} | {r['n_genes']:,} ({r['n_genes_named']:,}) | {r['n_name_rows']:,} |"
          for r in report["slowest"]],
        "",
        "## Errors",
        *( [f"- `{r['dir']}`: {r['error']}" for r in errors[:20]] or ["none"] ),
        "",
        "## Reading this",
        "- Cold timing: eviction is advisory on NFS; a repeated run on the same shards can be warmer. "
        "Use a different `--seed` for a fresh sample.",
        "- Gate: stop and re-decide if the projected wall time or total size is unacceptable.",
    ]
    (args.report_dir / "m2_report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    log(f"report written to {args.report_dir}")
    return 1 if (errors or bad_verify) else 0


if __name__ == "__main__":
    raise SystemExit(main())
