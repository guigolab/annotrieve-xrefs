"""Streaming writers for /loci: tsv, bed, ndjson (each ends with a completeness trailer)."""
from __future__ import annotations

import json
import time
from collections.abc import Iterable, Iterator

from services.presence_loci import LociQuery

TSV_COLUMNS = [
    "annotation_id", "organism_name", "assembly_accession", "source", "symbol", "terms",
    "matched_names", "matched_via", "seqid", "start", "end", "strand", "feature_type",
    "biotype", "locus_ref", "gene_id", "gff_url",
]
_BATCH = 500


def _cell(v) -> str:
    if v is None:
        return ""
    return str(v).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _provenance(q: LociQuery, fmt: str) -> list[str]:
    return [
        f"#format={fmt} data_version={q.data_version} generated={time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
        f"#query names={','.join(t['input'] for t in q.terms)} alias={q.mode} "
        f"include_taxids={','.join(map(str, q.include))} exclude_taxids={','.join(map(str, q.exclude))} "
        f"case_sensitive={str(q.case_sensitive).lower()} flank_bp={q.flank_bp}",
    ]


def _limited(rows: Iterable[dict], limit: int, state: dict) -> Iterator[dict]:
    for row in rows:
        if state["rows"] >= limit:
            state["truncated"] = True
            return
        state["rows"] += 1
        yield row


def _trailer(state: dict) -> str:
    return f"#end rows={state['rows']} complete={str(not state['truncated']).lower()}\n"


def _tsv_row(r: dict) -> str:
    vals = dict(r)
    vals["terms"] = ",".join(map(str, r["terms"]))
    vals["matched_names"] = ",".join(r["matched_names"])
    return "\t".join(_cell(vals.get(c)) for c in TSV_COLUMNS) + "\n"


def _guard(inner: Iterator[str]) -> Iterator[str]:
    """Report a mid-stream failure in-band (the HTTP status is already sent)."""
    try:
        yield from inner
    except Exception as exc:  # noqa: BLE001
        yield f"#error {type(exc).__name__}: {_cell(exc)}\n"


def stream_tsv(q: LociQuery, rows: Iterable[dict], limit: int) -> Iterator[str]:
    return _guard(_stream_tsv(q, rows, limit))


def stream_bed(q: LociQuery, rows: Iterable[dict], limit: int) -> Iterator[str]:
    return _guard(_stream_bed(q, rows, limit))


def _stream_tsv(q: LociQuery, rows: Iterable[dict], limit: int) -> Iterator[str]:
    state = {"rows": 0, "truncated": False}
    yield "".join(line + "\n" for line in _provenance(q, "tsv"))
    yield "\t".join(TSV_COLUMNS) + "\n"
    buf: list[str] = []
    for r in _limited(rows, limit, state):
        buf.append(_tsv_row(r))
        if len(buf) >= _BATCH:
            yield "".join(buf)
            buf.clear()
    if buf:
        yield "".join(buf)
    yield _trailer(state)


def _stream_bed(q: LociQuery, rows: Iterable[dict], limit: int) -> Iterator[str]:
    state = {"rows": 0, "truncated": False}
    yield "".join(line + "\n" for line in _provenance(q, "bed"))
    buf: list[str] = []
    for r in _limited(rows, limit, state):
        name = f"{r['symbol']}|{r['annotation_id']}"
        buf.append("\t".join([_cell(r["seqid"]), str(max(0, r["start"] - 1)), str(r["end"]),
                              _cell(name), "0", r["strand"]]) + "\n")
        if len(buf) >= _BATCH:
            yield "".join(buf)
            buf.clear()
    if buf:
        yield "".join(buf)
    yield _trailer(state)


def stream_ndjson(q: LociQuery, rows: Iterable[dict], limit: int) -> Iterator[str]:
    state = {"rows": 0, "truncated": False}
    buf: list[str] = []
    try:
        for r in _limited(rows, limit, state):
            buf.append(json.dumps(r, ensure_ascii=False) + "\n")
            if len(buf) >= _BATCH:
                yield "".join(buf)
                buf.clear()
        if buf:
            yield "".join(buf)
        yield json.dumps({"end": {"rows": state["rows"], "complete": not state["truncated"]}}) + "\n"
    except Exception as exc:  # report the failure in-band; the HTTP status is already sent
        yield "".join(buf)
        yield json.dumps({"error": {"type": type(exc).__name__, "message": str(exc)}}) + "\n"
