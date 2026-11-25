#!/usr/bin/env python3
"""Free disk on the server: remove what production does not use. Dry run unless --apply.

    python3 tools/free_disk.py --work-dir /data/annotrieve/gene_corpus \\
        --results-dir /data/annotrieve-xrefs/sequences [--keep-versions 2] [--apply]

What it looks at (each line of the report is one item with its size):

  versions    serve/v<data_version>/ directories that are not ``current`` and not among the newest
              ``--keep-versions`` (default 2: the live one and one to roll back to). Stop the API
              or wait for a restart first if it still has an old version open (it reopens on its own
              when ``current`` changes, so a version older than two builds is safe).
  partial     serve/partial-*/ trial builds, and half-written build files (.*.building*, *.part, *.tmp).
  downloads   download folders under --results-dir that the API itself would purge: finished and
              older than 2 hours, failed and older than 1 hour, or without a status file and older
              than 1 hour. Running or queued ones are never touched.
  shard-temp  per_gff/*/ leftovers of an interrupted harvest (genes.sqlite.tmp*, *-wal, *-shm).
  legacy-db   the old global gene_corpus.sqlite (and its -wal/-shm): no route reads it any more.
              Only with --include-legacy-db, because it takes hours to rebuild.
  orphans     per_gff/<id>/ shards of annotations that are not in manifest.sqlite or are marked
              removed. Only with --orphan-shards.
  docker      dangling Docker images and build cache. Only with --docker (and --apply).

Nothing is followed through symlinks and nothing outside --work-dir and --results-dir is touched.
Memory use is constant: directories are walked one at a time.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

HOUR = 3600.0


@dataclass
class Item:
    kind: str
    path: Path
    size: int
    why: str


def human(n: int) -> str:
    x = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if x < 1024 or unit == "TB":
            return f"{x:,.1f} {unit}" if unit != "B" else f"{int(x)} B"
        x /= 1024
    return f"{n} B"


def size_of(path: Path) -> int:
    """Bytes under ``path`` (a file or a directory); symlinks are counted as themselves."""
    try:
        st = path.lstat()
    except OSError:
        return 0
    if not path.is_symlink() and path.is_dir():
        total = 0
        for root, dirs, files in os.walk(path, followlinks=False):
            for name in files:
                try:
                    total += os.lstat(os.path.join(root, name)).st_size
                except OSError:
                    pass
        return total
    return st.st_size


def _age_s(path: Path, now: float) -> float:
    try:
        return now - path.lstat().st_mtime
    except OSError:
        return 0.0


# --- the categories ---------------------------------------------------------------------------------

def find_versions(work: Path, keep: int) -> list[Item]:
    serve = work / "serve"
    if not serve.is_dir():
        return []
    current = (serve / "current").resolve() if (serve / "current").exists() else None
    versions = sorted((d for d in serve.iterdir() if d.is_dir() and not d.is_symlink()
                       and d.name.startswith("v")), key=lambda d: d.stat().st_mtime, reverse=True)
    keep_set = {v for v in versions[: max(keep, 1)]}
    if current is not None:
        keep_set.add(current)
    items = []
    for d in versions:
        if d in keep_set or (current is not None and d.resolve() == current):
            continue
        items.append(Item("versions", d, size_of(d), "old published version"))
    return items


def find_partial(work: Path) -> list[Item]:
    serve = work / "serve"
    if not serve.is_dir():
        return []
    items = []
    for p in serve.iterdir():
        name = p.name
        if p.is_symlink():
            continue
        if p.is_dir() and name.startswith("partial-"):
            items.append(Item("partial", p, size_of(p), "trial build"))
        elif p.is_file() and (name.startswith(".") and "building" in name
                              or name.endswith((".part", ".tmp"))):
            items.append(Item("partial", p, size_of(p), "half-written build file"))
    return items


def find_downloads(results: Path, now: float, keep_done_s: float, keep_failed_s: float,
                   stale_run_s: float) -> list[Item]:
    if not results.is_dir():
        return []
    items = []
    for d in results.iterdir():
        if d.is_symlink() or not d.is_dir() or d.name.startswith("."):
            continue
        status_file = d / "status.json"
        try:
            st = json.loads(status_file.read_text())
        except (OSError, ValueError):
            st = None
        if st is None:
            if _age_s(d, now) > HOUR:
                items.append(Item("downloads", d, size_of(d), "no status file"))
            continue
        state, upd = st.get("state"), float(st.get("updated", 0))
        if state == "done" and now - upd > keep_done_s:
            items.append(Item("downloads", d, size_of(d), "finished, past keep time"))
        elif state == "failed" and now - upd > keep_failed_s:
            items.append(Item("downloads", d, size_of(d), "failed, past keep time"))
        elif state == "running" and now - upd > stale_run_s:
            items.append(Item("downloads", d, size_of(d), "running status that stopped updating"))
    scratch = results / "tmp"
    if scratch.is_dir() and not scratch.is_symlink():
        for p in scratch.iterdir():
            if p.name.startswith(".nfs"):
                continue    # NFS silly-rename of a file still open elsewhere: cannot be removed
            if _age_s(p, now) > HOUR:
                items.append(Item("downloads", p, size_of(p), "old scratch file"))
    return items


def find_shard_temp(work: Path) -> list[Item]:
    root = work / "per_gff"
    if not root.is_dir():
        return []
    items = []
    for d in root.iterdir():
        if d.is_symlink() or not d.is_dir():
            continue
        for p in d.iterdir():
            n = p.name
            if p.is_file() and (n.startswith("genes.sqlite.tmp") or n in ("genes.sqlite-wal", "genes.sqlite-shm")
                                or n.endswith((".part", ".tmp"))):
                items.append(Item("shard-temp", p, size_of(p), "interrupted harvest leftover"))
    return items


def find_legacy_db(work: Path) -> list[Item]:
    items = []
    for name in ("gene_corpus.sqlite", "gene_corpus.sqlite-wal", "gene_corpus.sqlite-shm",
                 "gene_corpus.sqlite-journal"):
        p = work / name
        if p.is_file() and not p.is_symlink():
            items.append(Item("legacy-db", p, size_of(p), "old global corpus, unused by the API"))
    return items


def find_orphans(work: Path) -> list[Item]:
    root, manifest = work / "per_gff", work / "serve" / "manifest.sqlite"
    if not root.is_dir() or not manifest.is_file():
        return []
    conn = sqlite3.connect(f"file:{manifest}?mode=ro", uri=True)
    try:
        known = {}
        for aid, status in conn.execute("SELECT annotation_id, status FROM annotation"):
            known[aid] = status
    finally:
        conn.close()
    items = []
    for d in root.iterdir():
        if d.is_symlink() or not d.is_dir():
            continue
        if d.name not in known:
            items.append(Item("orphans", d, size_of(d), "not in manifest.sqlite"))
        elif known[d.name] == "removed":
            items.append(Item("orphans", d, size_of(d), "annotation removed from the catalogue"))
    return items


# --- running ------------------------------------------------------------------------------------------

def remove(item: Item) -> None:
    p = item.path
    if p.is_symlink() or not p.is_dir():
        p.unlink(missing_ok=True)
    else:
        shutil.rmtree(p)


def run_docker() -> int:
    freed = 0
    for cmd in (["docker", "image", "prune", "-f"], ["docker", "builder", "prune", "-f"]):
        r = subprocess.run(cmd, capture_output=True, text=True)
        print(f"$ {' '.join(cmd)}\n{r.stdout.strip() or r.stderr.strip()}")
        freed += r.returncode != 0
    return freed


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work-dir", type=Path, required=True,
                    help="corpus root: contains per_gff/ and serve/ (GENE_CORPUS_HOST_DIR)")
    ap.add_argument("--results-dir", type=Path, default=None,
                    help="downloads folder (SEQUENCES_RESULTS_DIR on the host)")
    ap.add_argument("--keep-versions", type=int, default=2,
                    help="published versions to keep, current included (default 2)")
    ap.add_argument("--include-legacy-db", action="store_true",
                    help="also remove the old global gene_corpus.sqlite")
    ap.add_argument("--orphan-shards", action="store_true",
                    help="also remove shards of annotations that are not in the manifest")
    ap.add_argument("--docker", action="store_true", help="also prune dangling Docker images and build cache")
    ap.add_argument("--apply", action="store_true", help="delete (default: only list)")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    work: Path = args.work_dir
    if not work.is_dir():
        print(f"--work-dir is not a directory: {work}", file=sys.stderr)
        return 2
    now = time.time()
    items: list[Item] = []
    items += find_versions(work, args.keep_versions)
    items += find_partial(work)
    if args.results_dir is not None:
        items += find_downloads(args.results_dir, now, keep_done_s=2 * HOUR, keep_failed_s=HOUR,
                                stale_run_s=2 * HOUR)
    items += find_shard_temp(work)
    skipped_notes = []
    legacy = find_legacy_db(work)
    if args.include_legacy_db:
        items += legacy
    elif legacy:
        skipped_notes.append(f"legacy-db: {human(sum(i.size for i in legacy))} not removed "
                             "(add --include-legacy-db)")
    orphans = find_orphans(work)
    if args.orphan_shards:
        items += orphans
    elif orphans:
        skipped_notes.append(f"orphans: {len(orphans)} shard(s), {human(sum(i.size for i in orphans))} "
                             "not removed (add --orphan-shards)")

    total = sum(i.size for i in items)
    for i in items:
        print(f"{i.kind:11} {human(i.size):>10}  {i.path}  ({i.why})")
    for note in skipped_notes:
        print(f"note: {note}")
    print(f"{len(items)} item(s), {human(total)} {'to remove' if not args.apply else 'selected'}")

    if not args.apply:
        print("dry run: nothing was deleted (add --apply)")
        if args.docker:
            print("docker: would prune dangling images and build cache")
        return 0
    freed, errors = 0, 0
    for i in items:
        try:
            remove(i)
            freed += i.size
        except OSError as exc:
            errors += 1
            print(f"could not remove {i.path}: {exc}", file=sys.stderr)
    if args.docker:
        errors += run_docker()
    print(f"freed {human(freed)}" + (f", {errors} error(s)" if errors else ""))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
