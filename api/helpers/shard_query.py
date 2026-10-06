"""Read-only queries against one per-GFF genes.sqlite shard."""
from __future__ import annotations

import sqlite3
from pathlib import Path

from helpers.paths import per_gff_root, shard_genes_path
from helpers.session import GeneCorpusUnavailable

_SHARD_CACHE_SIZE = -16384  # 16 MiB
# Pathological gene with huge attr dumps; one locus card stays bounded.
MAX_XREFS_PER_GENE = 2000


class ShardUnavailable(Exception):
    """Raised when a per-GFF genes.sqlite is missing."""


def open_shard(
    annotation_id: str,
    *,
    base: Path | None = None,
) -> sqlite3.Connection:
    """
    Open one annotation's genes.sqlite read-only.

    Caller must close the connection. Not process-cached.
    """
    if per_gff_root(base) is None:
        raise GeneCorpusUnavailable("Gene corpus is not published")
    path = shard_genes_path(annotation_id, base=base)
    if path is None or not path.is_file():
        raise ShardUnavailable(
            f"Gene shard not found for annotation {annotation_id}"
        )
    uri = f"file:{path.resolve()}?mode=ro&immutable=1"
    conn = sqlite3.connect(uri, uri=True)
    conn.row_factory = sqlite3.Row
    conn.execute(f"PRAGMA cache_size={_SHARD_CACHE_SIZE}")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute("PRAGMA mmap_size=0")
    conn.execute("PRAGMA query_only=ON")
    return conn


def get_gene_with_xrefs(conn: sqlite3.Connection, local_id: int) -> dict | None:
    """Return one gene row plus gene_xref chips (capped), or None if missing."""
    if local_id < 0:
        return None
    row = conn.execute(
        """
        SELECT
            local_id, source_gene_id, feature_type, seqid, start, end,
            strand, biotype, primary_name, prose, prose_kind
        FROM gene
        WHERE local_id = ?
        """,
        (int(local_id),),
    ).fetchone()
    if row is None:
        return None
    xrefs = conn.execute(
        """
        SELECT namespace, accession, display, origin, via_level
        FROM gene_xref
        WHERE local_id = ?
        ORDER BY namespace, accession
        LIMIT ?
        """,
        (int(local_id), MAX_XREFS_PER_GENE),
    ).fetchall()
    return {
        "local_id": int(row["local_id"]),
        "source_gene_id": row["source_gene_id"],
        "feature_type": row["feature_type"],
        "seqid": row["seqid"],
        "start": int(row["start"]),
        "end": int(row["end"]),
        "strand": int(row["strand"]),
        "biotype": row["biotype"],
        "primary_name": row["primary_name"],
        "prose": row["prose"],
        "prose_kind": row["prose_kind"],
        "xrefs": [
            {
                "namespace": x["namespace"],
                "accession": x["accession"],
                "display": x["display"],
                "origin": x["origin"],
                "via_level": int(x["via_level"]),
            }
            for x in xrefs
        ],
    }
