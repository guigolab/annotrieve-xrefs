"""Resolve published gene_corpus paths (global DB + per-GFF shards)."""
from __future__ import annotations

from pathlib import Path

from settings import gene_corpus_dir

GENE_CORPUS_SQLITE = "gene_corpus.sqlite"
PER_GFF_DIRNAME = "per_gff"
SHARD_GENES_SQLITE = "genes.sqlite"
MAX_ANNOTATION_ID_CHARS = 128


def gene_corpus_db_path(base: Path | None = None) -> Path | None:
    """Return ``{dir}/gene_corpus.sqlite`` if the parent dir is configured."""
    root = base if base is not None else gene_corpus_dir()
    if root is None:
        return None
    return root / GENE_CORPUS_SQLITE


def sanitize_annotation_id(annotation_id: str) -> str:
    """Replace path-unsafe characters; IDs today are hex-like."""
    text = (
        (annotation_id or "")
        .replace("\0", "")
        .replace("/", "_")
        .replace("\\", "_")
        .replace("..", "_")
    )
    return text.strip()


def per_gff_root(base: Path | None = None) -> Path | None:
    """Return ``{GENE_CORPUS_DIR}/per_gff`` if configured."""
    root = base if base is not None else gene_corpus_dir()
    if root is None:
        return None
    return root / PER_GFF_DIRNAME


def shard_genes_path(annotation_id: str, base: Path | None = None) -> Path | None:
    """
    Return ``{per_gff}/{sanitized_id}/genes.sqlite`` if corpus dir is set.

    Returns None when the id is empty/too long or would escape ``per_gff/``.
    """
    root = per_gff_root(base)
    if root is None:
        return None
    cleaned = sanitize_annotation_id(annotation_id)
    if not cleaned or len(cleaned) > MAX_ANNOTATION_ID_CHARS:
        return None
    path = (root / cleaned / SHARD_GENES_SQLITE).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError:
        return None
    return path
