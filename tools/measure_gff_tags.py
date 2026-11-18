"""
Measure which `tag=` values real GFFs carry, to validate the canonical-transcript list.

    python3 tools/measure_gff_tags.py --root /data/annotrieve/annotations --sample 300

Walks --root for *.gff / *.gff3 (optionally .gz), takes a deterministic sample (every k-th file of
the sorted list), and STREAMS each file line by line (memory stays small: one counter per gene
parent of the current file). For each file class (RefSeq GCF_, GenBank GCA_, other = Ensembl or
unknown) it reports:

  - the most frequent tag values on non-structural, non-gene features (transcript lines);
  - how many files carry any canonical-looking tag, and which ones;
  - for genes with a canonical-flagged transcript: how many have exactly 1 flag vs several;
  - for tags that LOOK canonical but are not in transcript_select.CANONICAL_TAGS.

Prints peak RSS at the end. Read-only.
"""
from __future__ import annotations

import argparse
import gzip
import resource
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "api")]
from helpers.transcript_select import CANONICAL_TAGS, normalise_tag  # noqa: E402

STRUCTURAL = {"exon", "cds", "utr", "five_prime_utr", "three_prime_utr", "5utr", "3utr", "intron"}
GENE_TYPES = {"gene", "pseudogene"}
LOOKS_CANONICAL = ("canon", "select", "mane", "principal", "primary", "ccds", "appris")


def classify(path: Path) -> str:
    name = path.name
    if "GCF_" in name:
        return "refseq"
    if "GCA_" in name:
        return "genbank"
    return "other"


def open_text(path: Path):
    if path.suffix in (".gz", ".bgz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def find_files(root: Path) -> list[Path]:
    out: list[Path] = []
    for p in root.rglob("*"):
        n = p.name.casefold()
        if ".gff" in n and not n.endswith((".csi", ".tbi", ".gzi", ".tsv", ".json", ".md5")):
            out.append(p)
    return sorted(out)


def scan_file(path: Path, tags: Counter, per_gene: Counter) -> tuple[int, int]:
    """Return (transcript lines, tagged transcript lines); fill tag counts and flags per gene."""
    n_tx = n_tagged = 0
    for line in open_text(path):
        if not line or line[0] == "#":
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        ftype = cols[2].casefold()
        if ftype in STRUCTURAL or ftype.endswith("_utr") or ftype in GENE_TYPES:
            continue
        attrs = cols[8]
        parent = ""
        raw_tag = ""
        for part in attrs.split(";"):
            key, _, value = part.partition("=")
            k = key.strip().casefold()
            if k == "parent":
                parent = value.split(",")[0]
            elif k in ("tag", "tags"):
                raw_tag = unquote(value)
        if not parent:
            continue
        n_tx += 1
        if raw_tag:
            n_tagged += 1
            vals = [normalise_tag(v) for v in raw_tag.split(",") if v.strip()]
            for v in vals:
                tags[v] += 1
            if any(v in CANONICAL_TAGS for v in vals):
                per_gene[parent] += 1
    return n_tx, n_tagged


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--sample", type=int, default=300, help="number of files (deterministic)")
    ap.add_argument("--top", type=int, default=25)
    a = ap.parse_args()

    files = find_files(a.root)
    if not files:
        seen = []
        for i, p in enumerate(a.root.rglob("*") if a.root.is_dir() else []):
            if i >= 8:
                break
            seen.append(str(p))
        sys.exit(f"no GFF files under {a.root} (exists={a.root.is_dir()}); first entries: {seen}")
    step = max(1, len(files) // max(1, a.sample))
    sample = files[::step][: a.sample]
    print(f"{len(files)} GFF files found, scanning {len(sample)} (every {step}th)")

    cls_tags: dict[str, Counter] = defaultdict(Counter)
    cls_files: Counter = Counter()
    cls_files_canon: Counter = Counter()
    cls_tx: Counter = Counter()
    cls_tagged: Counter = Counter()
    cls_genes_one: Counter = Counter()
    cls_genes_many: Counter = Counter()
    no_tag_examples: dict[str, list[str]] = defaultdict(list)

    for i, path in enumerate(sample, 1):
        c = classify(path)
        tags: Counter = Counter()
        per_gene: Counter = Counter()
        try:
            n_tx, n_tagged = scan_file(path, tags, per_gene)
        except (OSError, EOFError, UnicodeError) as exc:
            print(f"  skip {path.name}: {type(exc).__name__}: {exc}")
            continue
        cls_files[c] += 1
        cls_tx[c] += n_tx
        cls_tagged[c] += n_tagged
        cls_tags[c].update(tags)
        if per_gene:
            cls_files_canon[c] += 1
            cls_genes_one[c] += sum(1 for v in per_gene.values() if v == 1)
            cls_genes_many[c] += sum(1 for v in per_gene.values() if v > 1)
        elif len(no_tag_examples[c]) < 5:
            no_tag_examples[c].append(path.name)
        if i % 50 == 0:
            print(f"  {i}/{len(sample)} files", flush=True)

    for c in sorted(cls_files):
        print(f"\n=== {c}: {cls_files[c]} files ===")
        print(f"transcript-like lines {cls_tx[c]:,}; with a tag attribute {cls_tagged[c]:,}")
        print(f"files with a recognised canonical tag: {cls_files_canon[c]} / {cls_files[c]}")
        print(f"genes with exactly 1 flagged transcript: {cls_genes_one[c]:,}; with several: "
              f"{cls_genes_many[c]:,}")
        print(f"top tag values (recognised canonical marked *):")
        for tag, n in cls_tags[c].most_common(a.top):
            mark = "*" if tag in CANONICAL_TAGS else " "
            print(f"  {mark} {tag:<32} {n:>12,}")
        unknown = [t for t in cls_tags[c]
                   if t not in CANONICAL_TAGS and any(w in t for w in LOOKS_CANONICAL)]
        if unknown:
            print("  tags that look canonical but are NOT recognised:", ", ".join(sorted(unknown)))
        if no_tag_examples[c]:
            print("  examples without a recognised canonical tag:", ", ".join(no_tag_examples[c]))

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"\npeak RSS MiB {peak:.0f}")


if __name__ == "__main__":
    main()
