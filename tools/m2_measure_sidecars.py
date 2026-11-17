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


def fresh_sets(shard: Path) -> tuple[set, list]:
    """Independent query used only for validation (does not reuse writer code)."""
    conn = sqlite3.connect(f"file:{shard}?mode=ro", uri=True)
    try:
        names = {
            (a.casefold(), "s" if ns == "symbol" else "a", int(lid))
            for a, ns, lid in conn.execute(
                "SELECT accession, namespace, local_id FROM gene_xref "
                "WHERE namespace IN ('symbol','alias') AND accession IS NOT NULL"
            )
        }
        genes = conn.execute(
            "SELECT local_id, start, end, strand FROM gene ORDER BY local_id"
        ).fetchall()
    finally:
        conn.close()
    return names, genes


def measure_one(args: tuple[str, str, str, bool, bool]) -> dict:
    annotation_dir, per_gff, out_dir, verify, do_evict = args
    shard = Path(per_gff) / annotation_dir / SHARD_NAME
    out = Path(out_dir) / f"{annotation_dir}.sidecar.tsv.gz"
    rec: dict = {"dir": annotation_dir, "shard_bytes": shard.stat().st_size}
    try:
        if do_evict:
            evict(shard)
        t0 = time.perf_counter()
        info = write_sidecar(shard, out)
        rec["seconds"] = time.perf_counter() - t0
        rec.update(
            n_genes=info.n_genes,
            n_name_rows=info.n_name_rows,
            n_distinct_names=info.n_distinct_names,
            n_casefold_changed=info.n_casefold_changed,
            n_skipped_names=info.n_skipped_names,
            sidecar_bytes=info.size_bytes,
            source_database=info.meta.get("source_database", ""),
            source_provider=info.meta.get("source_provider", ""),
        )
        if verify:
            names, genes = fresh_sets(shard)
            got_names = set(iter_names(out))
            got_genes = [(g[0], g[3], g[4], g[5]) for g in iter_genes(out)]
            # sidecar skips unsafe names; compare allowing exactly that many missing
            missing = len(names - got_names)
            extra = len(got_names - names)
            rec["verify_ok"] = (
                extra == 0
                and missing <= info.n_skipped_names
                and got_genes == [tuple(g) for g in genes]
            )
            if not rec["verify_ok"]:
                rec["verify_detail"] = {
                    "missing": missing, "extra": extra,
                    "skipped": info.n_skipped_names,
                    "genes_equal": got_genes == [tuple(g) for g in genes],
                }
    except Exception as exc:  # report, do not stop the sample
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["maxrss_mib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    return rec


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

    results: list[dict] = []
    wall0 = time.perf_counter()
    if args.workers <= 1:
        for i, w in enumerate(work, 1):
            results.append(measure_one(w))
            if i % 20 == 0 or i == len(work):
                log(f"  {i}/{len(work)}")
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(measure_one, w) for w in work]
            for i, f in enumerate(as_completed(futs), 1):
                results.append(f.result())
                if i % 20 == 0 or i == len(work):
                    log(f"  {i}/{len(work)}")
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
        "| dir | shard | seconds | genes | name rows |",
        "|---|---|---|---|---|",
        *[f"| `{r['dir']}` | {fmt_bytes(r['shard_bytes'])} | {r['seconds']:.1f} | {r['n_genes']:,} | {r['n_name_rows']:,} |"
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
