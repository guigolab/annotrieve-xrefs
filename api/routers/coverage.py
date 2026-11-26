"""How much of each clade a name search can reach (``/coverage``)."""
from __future__ import annotations

from fastapi import APIRouter, Query, Request

from helpers.cursor import QueryError
from helpers.errors import http_for_query_error
from helpers.names import parse_taxid_list
from helpers.presence import presence_reader
from helpers.strict import reject_unknown_params

router = APIRouter()

#: Taxids answered in one call (one overview page of clades is far fewer).
MAX_COVERAGE_TAXIDS = 500
_CHUNK = 200


def _has_named_table(conn) -> bool:
    row = conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='named'").fetchone()
    return row is not None


@router.get("/coverage", summary="Annotations and organisms per clade, and how many carry gene names")
def get_coverage(
    request: Request,
    taxids: str = Query(..., description=f"Comma-separated NCBI taxids (at most {MAX_COVERAGE_TAXIDS})"),
):
    """
    For each taxid: annotations in the index under that clade, those whose GFF carries at least
    one gene name, and the same two counts in distinct organisms. Older index builds do not
    record which annotations carry names: then ``available`` is false and ``results`` is empty.
    """
    reject_unknown_params(request, ("taxids",))
    try:
        ids = parse_taxid_list(taxids, field="taxids", max_taxids=MAX_COVERAGE_TAXIDS)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc
    reader = presence_reader()
    conn = reader.conn
    if not _has_named_table(conn):
        return {"available": False, "data_version": reader.data_version, "results": []}
    found: dict[int, dict] = {}
    for i in range(0, len(ids), _CHUNK):
        chunk = ids[i:i + _CHUNK]
        marks = ",".join("?" * len(chunk))
        for taxid, n_ann, n_named, n_org, n_org_named in conn.execute(
            f"""
            SELECT l.taxid,
                   COUNT(*),
                   COUNT(n.ord),
                   COUNT(DISTINCT a.taxid),
                   COUNT(DISTINCT CASE WHEN n.ord IS NOT NULL THEN a.taxid END)
            FROM lineage l
            JOIN annotation a ON a.ord = l.ord
            LEFT JOIN named n ON n.ord = l.ord
            WHERE l.taxid IN ({marks})
            GROUP BY l.taxid
            """,
            chunk,
        ):
            found[taxid] = {
                "taxid": taxid,
                "n_annotations": n_ann,
                "n_annotations_with_names": n_named,
                "n_organisms": n_org,
                "n_organisms_with_names": n_org_named,
            }
    results = [
        found.get(t, {"taxid": t, "n_annotations": 0, "n_annotations_with_names": 0,
                      "n_organisms": 0, "n_organisms_with_names": 0})
        for t in ids
    ]
    return {"available": True, "data_version": reader.data_version, "results": results}
