"""HTTP-facing DB open helpers shared by routers."""
from __future__ import annotations

import sqlite3

from helpers.errors import http_error
from helpers.paths import MAX_ANNOTATION_ID_CHARS
from services.annotations import GeneCorpusUnavailable, ShardUnavailable, open_shard


def annotation_shard(annotation_id: str) -> sqlite3.Connection:
    """Open one annotation genes.sqlite, or raise 404/503."""
    if (
        not annotation_id
        or len(annotation_id) > MAX_ANNOTATION_ID_CHARS
        or "\0" in annotation_id
    ):
        raise http_error(404, "shard_not_found", "Gene shard not found")
    try:
        return open_shard(annotation_id)
    except GeneCorpusUnavailable as exc:
        raise http_error(
            503,
            "gene_corpus_unavailable",
            str(exc),
        ) from exc
    except ShardUnavailable as exc:
        raise http_error(404, "shard_not_found", str(exc)) from exc
