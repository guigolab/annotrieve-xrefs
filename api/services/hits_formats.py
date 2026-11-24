"""TSV writer for /hits (annotation-level table with ``#`` provenance lines and a completeness trailer)."""
from __future__ import annotations

import time
from collections.abc import Iterator

TSV_COLUMNS = [
    "annotation_id", "organism_name", "taxid", "assembly_accession", "assembly_name", "source",
    "release_date", "n_matched", "n_terms", "matched_terms", "alias_only_terms", "missing_terms",
    "download_url",
]


_ROWS_PER_CHUNK = 500


def _cell(v) -> str:
    if v is None:
        return ""
    return str(v).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def stream_hits_tsv(result: dict) -> Iterator[str]:
    """``result`` is a ``shape=records`` answer of ``list_name_hits_presence``."""
    terms = [t["input"] for t in result["terms"]]

    def names(indices) -> str:
        return ",".join(terms[i] for i in indices)

    yield (f"#format=tsv data_version={result['data_version']} "
           f"generated={time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")
    yield (f"#query names={','.join(terms)} alias={result['alias']} "
           f"include_taxids={','.join(map(str, result['include_taxids']))} "
           f"exclude_taxids={','.join(map(str, result['exclude_taxids']))} "
           f"sources={','.join(result['source_filter'])} min_matched={result['min_matched'] or ''}\n")
    yield "\t".join(TSV_COLUMNS) + "\n"
    n = 0
    batch: list[str] = []  # the response is pulled from a worker thread: many rows per hop
    try:
        for r in result["results"]:
            vals = dict(r)
            vals["matched_terms"] = names(r["matched"])
            vals["alias_only_terms"] = names(r["alias_only"])
            vals["missing_terms"] = names(r["missing"])
            batch.append("\t".join(_cell(vals.get(c)) for c in TSV_COLUMNS) + "\n")
            n += 1
            if len(batch) >= _ROWS_PER_CHUNK:
                yield "".join(batch)
                batch.clear()
        if batch:
            yield "".join(batch)
    except Exception as exc:  # noqa: BLE001 - the status line is already sent
        if batch:
            yield "".join(batch)
        yield f"#error {type(exc).__name__}: {_cell(exc)}\n"
        return
    yield f"#end rows={n} complete=true\n"
