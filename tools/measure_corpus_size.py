#!/usr/bin/env python3
"""
Measure how big the xrefs corpus is, to size a transcript/exon table.

Read-only, standard library only, streaming. It never loads a table into memory (peak RSS is
reported per step and must stay well under 2 GB). Results are rewritten to ``report.json`` and
``report.md`` after every step, so a partial run still leaves something useful.

    python3 tools/measure_corpus_size.py \
        --presence-dir "$PRESENCE_DIR" --corpus-dir "$GENE_CORPUS_DIR" \
        --out-dir /path/on/nfs/measure --shard-sample 400 --gff-sample 20

Steps (``--steps`` picks a subset, comma separated; default is all):
  presence  annotations (by source), taxa, distinct names, symbols vs alias-only, posting sizes
  loci      indexed loci: rows, distinct named genes, rows per gene, per-annotation spread, top names
  shards    per-annotation gene tables (sampled; ``--all-shards`` reads every one): genes, named
            genes, feature types, xref namespaces; extrapolated to the whole corpus
  gff       streams a sample of GFFs and counts genes, transcripts, exons and CDS lines, then
            estimates the transcript/exon table (rows and bytes) for several scopes

Scratch space: SQLite sorts spill to ``--tmp-dir`` (set it to the NFS work area, not local disk).
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import random
import resource
import sqlite3
import statistics
import sys
import time
from collections import Counter
from pathlib import Path

PRESENCE_FILE = "presence.sqlite"
LOCI_FILE = "loci.sqlite"
REF_FILE = "loci_ref.sqlite"
SOURCES = {1: "ensembl", 2: "genbank", 3: "refseq", 0: "other"}

GENE_TYPES = {"gene", "pseudogene", "ncRNA_gene", "transposable_element_gene", "tRNA_gene",
              "rRNA_gene", "snRNA_gene", "snoRNA_gene", "miRNA_gene", "lnc_RNA_gene"}
TRANSCRIPT_TYPES = {
    "mRNA", "transcript", "lnc_RNA", "lncRNA", "ncRNA", "rRNA", "tRNA", "snRNA", "snoRNA",
    "miRNA", "scRNA", "antisense_RNA", "primary_transcript", "pseudogenic_transcript",
    "unconfirmed_transcript", "guide_RNA", "RNase_P_RNA", "SRP_RNA", "telomerase_RNA", "Y_RNA",
    "vault_RNA", "V_gene_segment", "C_gene_segment", "J_gene_segment", "D_gene_segment",
}

# Size model for the proposed transcript table (bytes). Documented guesses, shown in the report.
ROW_OVERHEAD_B = 48            # key, start/end, strand, biotype, counts, flags (WITHOUT ROWID row)
BYTES_PER_EXON_LOW = 3         # two small delta varints
BYTES_PER_EXON_HIGH = 6        # larger introns, plus a CDS interval per exon


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def peak_rss_mb() -> float:
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)   # KiB on Linux


def open_ro(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
    conn.execute("PRAGMA cache_size=-16384")      # 16 MiB
    conn.execute("PRAGMA mmap_size=0")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute("PRAGMA query_only=ON")
    return conn


def scalar(conn: sqlite3.Connection, sql: str, *args):
    row = conn.execute(sql, args).fetchone()
    return row[0] if row else None


def has_table(conn: sqlite3.Connection, name: str) -> bool:
    return bool(scalar(conn, "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name=?", name))


def size_mb(path: Path) -> float:
    return round(path.stat().st_size / 1_048_576, 1)


def quantiles(values: list[int]) -> dict:
    if not values:
        return {}
    values = sorted(values)
    pick = lambda q: values[min(len(values) - 1, int(q * (len(values) - 1)))]
    return {"n": len(values), "min": values[0], "p25": pick(.25), "median": pick(.5),
            "p75": pick(.75), "p95": pick(.95), "max": values[-1],
            "mean": round(statistics.fmean(values), 1)}


class Report:
    def __init__(self, out_dir: Path):
        self.out_dir = out_dir
        self.data: dict = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "steps": {}}
        out_dir.mkdir(parents=True, exist_ok=True)

    def set(self, step: str, result: dict) -> None:
        result["peak_rss_mb"] = peak_rss_mb()
        self.data["steps"][step] = result
        self.flush()

    def flush(self) -> None:
        (self.out_dir / "report.json").write_text(json.dumps(self.data, indent=2, default=str))
        lines = [f"# Corpus size measurement ({self.data['generated']})", ""]
        for step, res in self.data["steps"].items():
            lines += [f"## {step}", ""]
            lines += render(res)
            lines.append("")
        (self.out_dir / "report.md").write_text("\n".join(lines))


def render(value, indent: int = 0) -> list[str]:
    pad = "  " * indent
    out: list[str] = []
    if isinstance(value, dict):
        for k, v in value.items():
            if isinstance(v, (dict, list)) and v:
                out.append(f"{pad}- **{k}**:")
                out += render(v, indent + 1)
            else:
                out.append(f"{pad}- **{k}**: {v}")
    elif isinstance(value, list):
        for item in value:
            if isinstance(item, dict):
                out.append(pad + "- " + ", ".join(f"{k}={v}" for k, v in item.items()))
            else:
                out.append(f"{pad}- {item}")
    return out


# ---------------------------------------------------------------- presence

def step_presence(pdir: Path) -> dict:
    path = pdir / PRESENCE_FILE
    if not path.exists():
        return {"skipped": f"{path} not found"}
    conn = open_ro(path)
    res: dict = {"file_mb": size_mb(path)}
    res["meta"] = dict(conn.execute("SELECT key, value FROM meta"))
    res["annotations"] = scalar(conn, "SELECT COUNT(*) FROM annotation")
    res["annotations_by_source"] = {
        SOURCES.get(s, str(s)): n
        for s, n in conn.execute("SELECT source, COUNT(*) FROM annotation GROUP BY source")}
    res["distinct_taxa"] = scalar(conn, "SELECT COUNT(DISTINCT taxid) FROM annotation")
    res["lineage_rows"] = scalar(conn, "SELECT COUNT(*) FROM lineage")
    res["distinct_names"] = scalar(conn, "SELECT COUNT(*) FROM name_dict")
    res["names_with_symbol"] = scalar(conn, "SELECT COUNT(*) FROM name_dict WHERE n_sym > 0")
    res["names_alias_only"] = scalar(conn, "SELECT COUNT(*) FROM name_dict WHERE n_sym = 0")
    res["sum_annotations_as_symbol"] = scalar(conn, "SELECT SUM(n_sym) FROM name_dict")
    res["sum_annotations_as_alias_only"] = scalar(conn, "SELECT SUM(n_alias_only) FROM name_dict")
    res["posting_bytes_symbol"] = scalar(conn, "SELECT SUM(LENGTH(sym)) FROM posting")
    res["posting_bytes_alias"] = scalar(conn, "SELECT SUM(LENGTH(alias)) FROM posting")
    buckets = conn.execute(
        "SELECT CASE WHEN n_sym=0 THEN 'a: 0 (alias only)' WHEN n_sym=1 THEN 'b: 1 annotation' "
        "WHEN n_sym<5 THEN 'c: 2-4' WHEN n_sym<10 THEN 'd: 5-9' WHEN n_sym<100 THEN 'e: 10-99' "
        "WHEN n_sym<1000 THEN 'f: 100-999' WHEN n_sym<10000 THEN 'g: 1000-9999' "
        "ELSE 'h: 10000+' END b, COUNT(*) FROM name_dict GROUP BY b ORDER BY b")
    res["names_by_annotation_count"] = {b: n for b, n in buckets}
    conn.close()
    return res


# ---------------------------------------------------------------- loci

def step_loci(pdir: Path, top: int, deep: bool) -> dict:
    lpath, rpath, ppath = pdir / LOCI_FILE, pdir / REF_FILE, pdir / PRESENCE_FILE
    res: dict = {}
    if lpath.exists():
        conn = open_ro(lpath)
        res["loci_file_mb"] = size_mb(lpath)
        res["meta"] = dict(conn.execute("SELECT key, value FROM meta"))
        res["locus_rows"] = scalar(conn, "SELECT COUNT(*) FROM locus")
        res["rows_by_via"] = {str(v): n for v, n in
                              conn.execute("SELECT via, COUNT(*) FROM locus GROUP BY via")}
        res["seqid_rows"] = scalar(conn, "SELECT COUNT(*) FROM seqid")
        res["name_loci_rows"] = scalar(conn, "SELECT COUNT(*) FROM name_loci")
        res["sum_n_loci"] = scalar(conn, "SELECT SUM(n_loci) FROM name_loci")
        buckets = conn.execute(
            "SELECT CASE WHEN n_loci=1 THEN 'a: 1' WHEN n_loci<10 THEN 'b: 2-9' "
            "WHEN n_loci<100 THEN 'c: 10-99' WHEN n_loci<1000 THEN 'd: 100-999' "
            "WHEN n_loci<10000 THEN 'e: 1000-9999' ELSE 'f: 10000+' END b, COUNT(*) "
            "FROM name_loci GROUP BY b ORDER BY b")
        res["names_by_locus_count"] = {b: n for b, n in buckets}
        top_rows = conn.execute(
            "SELECT name_id, n_loci FROM name_loci ORDER BY n_loci DESC LIMIT ?", (top,)).fetchall()
        names: dict[int, str] = {}
        if ppath.exists():
            pconn = open_ro(ppath)
            for nid, _ in top_rows:
                row = pconn.execute("SELECT name FROM name_dict WHERE name_id=?", (nid,)).fetchone()
                names[nid] = row[0] if row else str(nid)
            pconn.close()
        res["top_names_by_loci"] = [{"name": names.get(nid, str(nid)), "n_loci": n}
                                    for nid, n in top_rows]
        if deep:
            log("loci: deep scan (feature types, biotypes)")
            res["rows_by_feature_type"] = {
                n: c for n, c in conn.execute(
                    "SELECT d.name, COUNT(*) FROM locus l JOIN dict_feature_type d "
                    "ON d.id = l.feature_type_id GROUP BY d.name ORDER BY 2 DESC LIMIT 30")}
            res["rows_by_biotype"] = {
                n: c for n, c in conn.execute(
                    "SELECT d.name, COUNT(*) FROM locus l JOIN dict_biotype d "
                    "ON d.id = l.biotype_id GROUP BY d.name ORDER BY 2 DESC LIMIT 30")}
        conn.close()
    else:
        res["loci_file"] = "missing"
    if rpath.exists():
        conn = open_ro(rpath)
        res["ref_file_mb"] = size_mb(rpath)
        res["named_genes_indexed"] = scalar(conn, "SELECT COUNT(*) FROM ref")
        per_ord = [n for _, n in conn.execute("SELECT ord, COUNT(*) FROM ref GROUP BY ord")]
        res["named_genes_per_annotation"] = quantiles(per_ord)
        res["annotations_with_named_genes"] = len(per_ord)
        if res.get("locus_rows") and res["named_genes_indexed"]:
            res["locus_rows_per_named_gene"] = round(
                res["locus_rows"] / res["named_genes_indexed"], 2)
        conn.close()
    else:
        res["ref_file"] = "missing"
    return res


# ---------------------------------------------------------------- shards

def annotation_rows(pdir: Path) -> list[tuple[int, str, int, str | None]]:
    """(ord, annotation_id, source, gff_path) for every annotation (a few tens of thousands)."""
    conn = open_ro(pdir / PRESENCE_FILE)
    rows = conn.execute("SELECT ord, annotation_id, source, gff_path FROM annotation").fetchall()
    conn.close()
    return rows


def pick_sample(rows: list, n: int, seed: int) -> list:
    """Sample proportionally per source (rows are tuples with the source at index 2)."""
    if n <= 0 or n >= len(rows):
        return list(rows)
    rng = random.Random(seed)
    by_src: dict[int, list] = {}
    for r in rows:
        by_src.setdefault(r[2], []).append(r)
    picked: list = []
    for src, group in by_src.items():
        k = max(1, round(n * len(group) / len(rows)))
        picked += rng.sample(group, min(k, len(group)))
    return picked


def step_shards(pdir: Path, cdir: Path | None, n: int, all_shards: bool, seed: int) -> dict:
    if cdir is None or not (cdir / "per_gff").is_dir():
        return {"skipped": "no corpus dir with per_gff/ (use --corpus-dir)"}
    rows = annotation_rows(pdir)
    sample = rows if all_shards else pick_sample(rows, n, seed)
    log(f"shards: reading {len(sample)} of {len(rows)} annotation shards")
    per_src: dict[str, dict[str, list[int]]] = {}
    feature_types: Counter = Counter()
    biotypes: Counter = Counter()
    namespaces: Counter = Counter()
    xref_rows = genes_total = named_total = missing = 0
    t0 = time.time()
    for i, (_ord, aid, src, _gff) in enumerate(sample, 1):
        shard = cdir / "per_gff" / aid / "genes.sqlite"
        if not shard.exists():
            missing += 1
            continue
        conn = open_ro(shard)
        n_genes = scalar(conn, "SELECT COUNT(*) FROM gene")
        n_named = scalar(conn, "SELECT COUNT(*) FROM gene WHERE primary_name IS NOT NULL "
                               "AND primary_name <> ''")
        s = per_src.setdefault(SOURCES.get(src, "other"), {"genes": [], "named": []})
        s["genes"].append(n_genes)
        s["named"].append(n_named)
        genes_total += n_genes
        named_total += n_named
        for name, c in conn.execute("SELECT feature_type, COUNT(*) FROM gene GROUP BY feature_type"):
            feature_types[name] += c
        for name, c in conn.execute("SELECT COALESCE(biotype,'(none)'), COUNT(*) FROM gene "
                                    "GROUP BY 1 ORDER BY 2 DESC LIMIT 15"):
            biotypes[name] += c
        for ns, c in conn.execute("SELECT namespace, COUNT(*) FROM gene_xref GROUP BY namespace"):
            namespaces[ns] += c
            xref_rows += c
        conn.close()
        if i % 100 == 0:
            log(f"shards: {i}/{len(sample)} ({time.time() - t0:.0f}s, rss {peak_rss_mb()} MB)")
    src_counts = Counter(SOURCES.get(r[2], "other") for r in rows)
    res: dict = {"annotations_total": len(rows), "shards_read": len(sample) - missing,
                 "shards_missing": missing, "sampled": not all_shards}
    est_genes = est_named = 0.0
    by_source = {}
    for src, d in per_src.items():
        mg, mn = statistics.fmean(d["genes"]), statistics.fmean(d["named"])
        by_source[src] = {"annotations": src_counts[src], "sampled": len(d["genes"]),
                          "genes_per_annotation": quantiles(d["genes"]),
                          "named_per_annotation_mean": round(mn, 1),
                          "named_fraction": round(mn / mg, 3) if mg else None}
        est_genes += mg * src_counts[src]
        est_named += mn * src_counts[src]
    res["by_source"] = by_source
    res["estimated_total_genes"] = round(est_genes)
    res["estimated_named_genes"] = round(est_named)
    res["feature_types_in_sample"] = dict(feature_types.most_common(15))
    res["top_biotypes_in_sample"] = dict(biotypes.most_common(15))
    res["xref_rows_in_sample"] = xref_rows
    res["xref_rows_per_gene_in_sample"] = round(xref_rows / genes_total, 2) if genes_total else None
    res["xref_namespaces_in_sample"] = dict(namespaces.most_common(15))
    return res


# ---------------------------------------------------------------- gff

def resolve_gff(gff_root: Path | None, raw: str | None) -> Path | None:
    if not raw:
        return None
    p = Path(raw)
    if p.is_absolute() and p.exists():
        return p
    if gff_root is not None:
        for cand in (gff_root / raw, gff_root / p.name):
            if cand.exists():
                return cand
    return p if p.exists() else None


def count_gff(path: Path, timeout_s: float) -> tuple[Counter, bool]:
    """Count column-3 feature types, streaming. Returns (counts, finished)."""
    counts: Counter = Counter()
    opener = gzip.open if path.suffix in (".gz", ".bgz") else open
    deadline = time.time() + timeout_s
    with opener(path, "rt", errors="replace") as fh:
        for n, line in enumerate(fh):
            if line.startswith("#"):
                continue
            cols = line.split("\t", 3)
            if len(cols) >= 3:
                counts[cols[2]] += 1
            if n % 200_000 == 0 and time.time() > deadline:
                return counts, False
    return counts, True


def table_size(rows: float, exons_per_tx: float) -> dict:
    lo = rows * (ROW_OVERHEAD_B + exons_per_tx * BYTES_PER_EXON_LOW) / 1e9
    hi = rows * (ROW_OVERHEAD_B + exons_per_tx * BYTES_PER_EXON_HIGH) / 1e9
    return {"rows": round(rows), "size_gb_low": round(lo, 2), "size_gb_high": round(hi, 2)}


def step_gff(pdir: Path, cdir: Path | None, gff_root: Path | None, n: int, seed: int,
             timeout_s: float, prior: dict) -> dict:
    rows = annotation_rows(pdir)
    sample = pick_sample([r for r in rows if r[3]], n, seed)
    log(f"gff: streaming {len(sample)} GFF files (timeout {timeout_s:.0f}s each)")
    totals: Counter = Counter()
    per_file = []
    for _ord, aid, src, raw in sample:
        path = resolve_gff(gff_root, raw)
        if path is None:
            per_file.append({"annotation_id": aid, "status": "gff not found", "path": raw})
            continue
        t0 = time.time()
        counts, done = count_gff(path, timeout_s)
        genes = sum(c for t, c in counts.items() if t in GENE_TYPES)
        tx = sum(c for t, c in counts.items() if t in TRANSCRIPT_TYPES)
        info = {"annotation_id": aid, "source": SOURCES.get(src, "other"),
                "status": "ok" if done else "partial (timeout)", "genes": genes,
                "transcripts": tx, "mRNA": counts.get("mRNA", 0), "exons": counts.get("exon", 0),
                "cds": counts.get("CDS", 0), "seconds": round(time.time() - t0, 1)}
        per_file.append(info)
        if done:
            totals.update({"genes": genes, "transcripts": tx, "mRNA": counts.get("mRNA", 0),
                           "exons": counts.get("exon", 0), "cds": counts.get("CDS", 0),
                           "files": 1})
    res: dict = {"files_ok": totals["files"], "files_total": len(sample),
                 "totals_in_sample": dict(totals), "per_file": per_file[:40]}
    if not totals["genes"]:
        res["note"] = "no complete GFF read; set --gff-root and/or raise --gff-timeout"
        return res
    tx_per_gene = totals["transcripts"] / totals["genes"]
    exons_per_tx = (totals["exons"] / totals["transcripts"]) if totals["transcripts"] else 0
    res["ratios"] = {"transcripts_per_gene": round(tx_per_gene, 2),
                     "mRNA_per_gene": round(totals["mRNA"] / totals["genes"], 2),
                     "exons_per_transcript": round(exons_per_tx, 2),
                     "cds_lines_per_mRNA": round(totals["cds"] / totals["mRNA"], 2)
                     if totals["mRNA"] else None}
    named = (prior.get("loci", {}) or {}).get("named_genes_indexed")
    total_genes = (prior.get("shards", {}) or {}).get("estimated_total_genes")
    est: dict = {"assumptions": {
        "row_overhead_bytes": ROW_OVERHEAD_B,
        "bytes_per_exon": [BYTES_PER_EXON_LOW, BYTES_PER_EXON_HIGH],
        "note": "transcripts per gene measured on the sample are assumed to hold for the named "
                "genes; real sizes need a build on a larger sample"}}
    if named:
        est["named_genes_all_transcripts"] = table_size(named * tx_per_gene, exons_per_tx)
        est["named_genes_one_transcript_each"] = table_size(named, exons_per_tx)
    if total_genes:
        est["all_genes_all_transcripts"] = table_size(total_genes * tx_per_gene, exons_per_tx)
    if not named and not total_genes:
        est["note2"] = "run the loci and shards steps too, so counts can be scaled"
    res["estimated_transcript_table"] = est
    return res


# ---------------------------------------------------------------- main

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--presence-dir", default=os.getenv("PRESENCE_DIR", ""))
    ap.add_argument("--corpus-dir", default=os.getenv("GENE_CORPUS_DIR", ""))
    ap.add_argument("--gff-root", default="", help="base dir for relative gff_path values")
    ap.add_argument("--out-dir", default="measure-out")
    ap.add_argument("--tmp-dir", default="", help="SQLite temp files (use the NFS work area)")
    ap.add_argument("--steps", default="presence,loci,shards,gff")
    ap.add_argument("--shard-sample", type=int, default=400)
    ap.add_argument("--all-shards", action="store_true", help="read every shard (slow on NFS)")
    ap.add_argument("--gff-sample", type=int, default=20)
    ap.add_argument("--gff-timeout", type=float, default=180.0, help="seconds per GFF")
    ap.add_argument("--top-names", type=int, default=25)
    ap.add_argument("--deep", action="store_true", help="also scan locus rows by type and biotype")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    if args.tmp_dir:
        os.environ["SQLITE_TMPDIR"] = args.tmp_dir
    if not args.presence_dir:
        ap.error("--presence-dir (or PRESENCE_DIR) is required")
    pdir = Path(args.presence_dir)
    if (pdir / "current").exists():
        pdir = (pdir / "current").resolve()
    cdir = Path(args.corpus_dir) if args.corpus_dir else None
    gff_root = Path(args.gff_root) if args.gff_root else None
    steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    report = Report(Path(args.out_dir))
    report.data["inputs"] = {"presence_dir": str(pdir), "corpus_dir": str(cdir),
                             "shard_sample": args.shard_sample, "gff_sample": args.gff_sample}
    done: dict = {}
    t0 = time.time()
    for step in steps:
        log(f"step {step}")
        if step == "presence":
            result = step_presence(pdir)
        elif step == "loci":
            result = step_loci(pdir, args.top_names, args.deep)
        elif step == "shards":
            result = step_shards(pdir, cdir, args.shard_sample, args.all_shards, args.seed)
        elif step == "gff":
            result = step_gff(pdir, cdir, gff_root, args.gff_sample, args.seed,
                              args.gff_timeout, done)
        else:
            log(f"unknown step {step!r}")
            continue
        done[step] = result
        report.set(step, result)
    log(f"done in {time.time() - t0:.0f}s, peak RSS {peak_rss_mb()} MB; "
        f"report in {report.out_dir}/report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
