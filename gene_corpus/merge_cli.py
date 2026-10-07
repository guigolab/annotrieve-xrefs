"""
Merge CLI: scan per_gff shards → global gene_corpus.sqlite.

    python -m gene_corpus.merge \\
        --work-dir /data/gene_corpus_build \\
        --taxonomy-tsv /path/to/flattened-tree.tsv

    # Resume an interrupted merge (NFS-safe; does not wipe the DB):
    python -m gene_corpus.merge \\
        --work-dir /data/gene_corpus_build \\
        --taxonomy-tsv /path/to/flattened-tree.tsv \\
        --resume --resume-after-key 14050
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path

from gene_corpus.sql.lineage import load_parent_map
from gene_corpus.sql.merge import merge_shards_into_db
from gene_corpus.sql.per_gff import (
    discover_shard_paths,
    meta_complete,
    read_shard_meta,
)
from gene_corpus.sql.schema import connect_for_build, init_schema, seed_annotations


def _seed_rows_from_shards(per_gff_root: Path) -> list[dict]:
    """
    Discover shards, read meta, assign dense annotation_key by sorted id.

    Skips shards with incomplete meta (warns on stderr).
    """
    paths = discover_shard_paths(per_gff_root)
    total = len(paths)
    progress_every = 50
    print(f"merge seed: scanning {total} shard(s)", file=sys.stderr, flush=True)

    by_id: dict[str, dict] = {}
    skipped_meta = 0
    done = 0
    for shard in paths:
        conn = sqlite3.connect(f"file:{shard.resolve()}?mode=ro", uri=True)
        try:
            meta = read_shard_meta(conn)
        finally:
            conn.close()
        if not meta_complete(meta):
            skipped_meta += 1
            print(
                f"merge warn: incomplete meta in {shard}",
                file=sys.stderr,
            )
        elif not meta.get("annotation_id"):
            skipped_meta += 1
            print(
                f"merge warn: empty annotation_id in {shard}",
                file=sys.stderr,
            )
        else:
            ann_id = meta["annotation_id"]
            by_id[ann_id] = {
                "annotation_id": ann_id,
                "taxid": int(meta["taxid"]),
                "assembly_accession": meta["assembly_accession"],
                "organism_name": meta["organism_name"] or None,
                "source_database": meta["source_database"],
                "source_provider": meta["source_provider"] or None,
                "profile_id": meta["profile_id"],
                "gff_path": meta["gff_path"] or None,
            }
        done += 1
        if done % progress_every == 0 or done == total:
            print(
                f"merge seed progress: {done}/{total} "
                f"(ready={len(by_id)} skipped_meta={skipped_meta})",
                file=sys.stderr,
                flush=True,
            )

    rows: list[dict] = []
    for i, ann_id in enumerate(sorted(by_id), start=1):
        row = dict(by_id[ann_id])
        row["annotation_key"] = i
        rows.append(row)
    return rows


def _validate_annotation_seed(
    conn: sqlite3.Connection, seed_rows: list[dict]
) -> None:
    """Ensure DB annotation keys match a fresh shard rescan (no key shift)."""
    existing = conn.execute(
        "SELECT annotation_key, annotation_id FROM annotation "
        "ORDER BY annotation_key"
    ).fetchall()
    expected = [
        (int(r["annotation_key"]), str(r["annotation_id"])) for r in seed_rows
    ]
    got = [(int(k), str(aid)) for k, aid in existing]
    if got != expected:
        raise ValueError(
            "resume annotation mismatch: DB keys/ids do not match shard rescan "
            f"(db={len(got)} shards_ready={len(expected)}). "
            "Refusing to resume (re-run without --resume only if you intend "
            "to wipe and rebuild)."
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Merge per_gff/*/genes.sqlite into global gene_corpus.sqlite "
            "(gene_hit + xref_meta + annotation_lineage). "
            "Annotation keys assigned at merge time."
        )
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("tmp/gene_corpus_build"),
        help="directory containing per_gff/ (writes gene_corpus.sqlite)",
    )
    parser.add_argument(
        "--taxonomy-tsv",
        type=Path,
        required=True,
        help=(
            "Annotrieve flattened taxonomy TSV (taxid, parent_taxid columns); "
            "used to fill annotation_lineage"
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "continue an interrupted merge without wiping gene_corpus.sqlite; "
            "skips annotation_keys already in merge_done"
        ),
    )
    parser.add_argument(
        "--resume-after-key",
        type=int,
        default=None,
        metavar="N",
        help=(
            "with --resume: treat annotation_keys 1..N as done (safe floor "
            "from merge.log, e.g. 14050); avoids scanning gene_hit on NFS"
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    work_dir: Path = args.work_dir
    per_gff_root = work_dir / "per_gff"
    db_path = work_dir / "gene_corpus.sqlite"
    taxonomy_tsv: Path = args.taxonomy_tsv
    resume: bool = bool(args.resume)
    resume_after_key: int | None = args.resume_after_key

    if resume_after_key is not None and not resume:
        print(
            "--resume-after-key requires --resume",
            file=sys.stderr,
        )
        return 1
    if resume_after_key is not None and resume_after_key < 0:
        print("--resume-after-key must be >= 0", file=sys.stderr)
        return 1

    if not per_gff_root.is_dir():
        print(f"per_gff directory missing: {per_gff_root}", file=sys.stderr)
        return 1
    if not taxonomy_tsv.is_file():
        print(f"taxonomy TSV missing: {taxonomy_tsv}", file=sys.stderr)
        return 1

    try:
        parents = load_parent_map(taxonomy_tsv)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    print(
        f"taxonomy: loaded {len(parents)} taxid→parent edges from {taxonomy_tsv}",
        file=sys.stderr,
        flush=True,
    )

    seed_rows = _seed_rows_from_shards(per_gff_root)
    if not seed_rows:
        print("no shards with complete meta to merge", file=sys.stderr)
        return 1

    if resume:
        if not db_path.is_file():
            print(
                f"resume failed: gene_corpus.sqlite missing: {db_path}",
                file=sys.stderr,
            )
            return 1
        conn = connect_for_build(db_path)
        init_schema(conn)
        try:
            _validate_annotation_seed(conn, seed_rows)
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            conn.close()
            return 1
        print(
            f"resume: validated {len(seed_rows)} annotations → {db_path}",
            file=sys.stderr,
            flush=True,
        )
    else:
        for path in (
            db_path,
            Path(str(db_path) + "-wal"),
            Path(str(db_path) + "-shm"),
        ):
            if path.exists():
                path.unlink()

        conn = connect_for_build(db_path)
        init_schema(conn)
        seed_annotations(conn, seed_rows)
        print(
            f"seeded {len(seed_rows)} annotations from shards → {db_path}",
            file=sys.stderr,
        )

    n_hit, n_meta, n_lineage = merge_shards_into_db(
        conn,
        per_gff_root,
        taxonomy_parents=parents,
        resume=resume,
        resume_after_key=resume_after_key,
    )
    conn.close()
    print(
        f"merged gene_hit={n_hit} xref_meta={n_meta} "
        f"annotation_lineage={n_lineage} → {db_path}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
