"""Parse and normalise the body of ``POST /downloads`` (and ``/transcripts``).

Exactly one source of genes:

- ``query``: the fields of ``/loci`` (names, include_taxids, exclude_taxids, alias,
  annotation_ids, case_sensitive);
- ``loci``: a list of ``"annotation_id:local_id"`` strings or objects
  ``{"locus_ref": ..., "transcript_id": ...}``; a pinned ``transcript_id`` skips the rule;
- ``items``: the legacy form ``{"annotation_id", "local_id", "transcript_id"}``, same meaning.

``result_id`` is a hash of the normalised request and the data version, so equal requests share
one cached result.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Any

from helpers.cursor import QueryError
from helpers.names import MAX_SQLITE_INT
from helpers.gffread_runner import SequenceAction
from helpers.pagination import MAX_SEQUENCES_BUFFER_BP, MAX_SEQUENCES_ID_CHARS
from helpers.transcript_select import Rule, RuleError, parse_rule
from helpers.pagination import MAX_NAMES_POST
from services.terms import parse_terms
from settings import settings

QUERY_FIELDS = ("names", "terms", "include_taxids", "exclude_taxids", "alias", "annotation_ids",
                "case_sensitive")
MAX_LOCI_ENTRIES = 5000
ON_ERROR = ("skip", "fail")
ALIAS_MODES = ("off", "suggest", "include")


@dataclass(frozen=True)
class LociEntry:
    annotation_id: str
    local_id: int
    transcript_id: str | None = None

    @property
    def locus_ref(self) -> str:
        return f"{self.annotation_id}:{self.local_id}"


@dataclass(frozen=True)
class SeqRequest:
    action: SequenceAction
    rule: Rule
    rule_spec: Any                       # what the client sent (normalised), for the hash
    buffer_bp: int
    include_exon_coords: bool
    on_error: str
    max_genes: int
    query: dict | None = None
    loci: tuple[LociEntry, ...] | None = None

    @property
    def can_give_several_per_gene(self) -> bool:
        """True when a gene can yield more than one sequence (the rule keeps several transcripts)."""
        if self.rule.keep == 1:
            return False
        return self.loci is None or any(e.transcript_id is None for e in self.loci)

    def normalised(self) -> dict:
        return {
            "action": self.action.value,
            "transcript": self.rule_spec,
            "buffer_bp": self.buffer_bp,
            "exon_coords": self.include_exon_coords,
            "on_error": self.on_error,
            "max_genes": self.max_genes,
            "query": self.query,
            "loci": None if self.loci is None else [
                [e.annotation_id, e.local_id, e.transcript_id] for e in self.loci],
        }

    def result_id(self, data_version: str) -> str:
        blob = json.dumps([self.normalised(), data_version], sort_keys=True, separators=(",", ":"))
        return hashlib.sha1(blob.encode()).hexdigest()[:20]


def _bad(message: str, code: str, **extra) -> QueryError:
    return QueryError(message, code=code, **extra)


def _text(value: Any, field: str) -> str:
    text = str(value if value is not None else "").strip()
    if not text:
        raise _bad(f"{field} is required", "validation_error", field=field)
    if len(text) > MAX_SEQUENCES_ID_CHARS or any(ord(ch) < 32 for ch in text):
        raise _bad(f"{field} is too long or has control characters", "validation_error", field=field)
    return text


def parse_ref(ref: Any, index: int) -> tuple[str, int]:
    """``"annotation_id:local_id"`` -> (annotation_id, local_id)."""
    text = _text(ref, "locus_ref")
    annotation_id, sep, local = text.rpartition(":")
    if not sep or not annotation_id:
        raise _bad("locus_ref must look like 'annotation_id:local_id'", "invalid_locus_ref",
                   index=index, value=text)
    if not (local.isascii() and local.isdigit()):
        raise _bad("locus_ref must end with an integer local_id", "invalid_locus_ref",
                   index=index, value=text)
    local_id = int(local)
    if local_id > MAX_SQLITE_INT:
        raise _bad("local_id is too large", "invalid_locus_ref", index=index, value=text)
    return annotation_id, local_id


def _entries(raw: Any, field: str) -> tuple[LociEntry, ...]:
    if not isinstance(raw, list) or not raw:
        raise _bad(f"{field} must be a non-empty list", "items_required", field=field)
    if len(raw) > MAX_LOCI_ENTRIES:
        raise _bad(f"at most {MAX_LOCI_ENTRIES} entries allowed", "too_many_items",
                   max_items=MAX_LOCI_ENTRIES, got=len(raw))
    out: list[LociEntry] = []
    seen: set[tuple[str, int, str | None]] = set()
    for i, item in enumerate(raw):
        pinned = None
        if isinstance(item, str):
            annotation_id, local_id = parse_ref(item, i)
        elif isinstance(item, dict):
            if field == "items":
                annotation_id = _text(item.get("annotation_id"), "annotation_id")
                raw_id = item.get("local_id")
                if isinstance(raw_id, bool) or not (
                        isinstance(raw_id, int)
                        or (isinstance(raw_id, str) and raw_id.isascii() and raw_id.isdigit())):
                    raise _bad("local_id must be an integer >= 0", "validation_error", index=i)
                local_id = int(raw_id)
                if not 0 <= local_id <= MAX_SQLITE_INT:
                    raise _bad("local_id must be an integer >= 0", "validation_error", index=i)
            else:
                annotation_id, local_id = parse_ref(item.get("locus_ref"), i)
            if item.get("transcript_id") not in (None, ""):
                pinned = _text(item.get("transcript_id"), "transcript_id")
        else:
            raise _bad(f"{field} entries must be strings or objects", "validation_error", index=i)
        key = (annotation_id, local_id, pinned)
        if key in seen:
            continue
        seen.add(key)
        out.append(LociEntry(annotation_id, local_id, pinned))
    return tuple(out)


def _int_list(raw: Any, field: str) -> list[int]:
    if raw is None:
        return []
    if not isinstance(raw, list):
        raise _bad(f"query.{field} must be a list of integers", "validation_error", field=field)
    out = []
    for v in raw:
        text = v.strip() if isinstance(v, str) else v
        if isinstance(v, bool) or not (
                isinstance(v, int) or (isinstance(text, str) and text.isascii() and text.isdigit())):
            raise _bad(f"query.{field} must contain integers", "invalid_taxid", field=field, value=v)
        if not 1 <= int(text) <= MAX_SQLITE_INT:
            raise _bad(f"query.{field} must contain positive integers", "invalid_taxid",
                       field=field, value=v)
        out.append(int(text))
    return sorted(set(out))


def _query(raw: Any) -> dict:
    if not isinstance(raw, dict):
        raise _bad("query must be an object", "validation_error", field="query")
    unknown = sorted(set(raw) - set(QUERY_FIELDS))
    if unknown:
        raise _bad(f"unknown query fields {unknown}", "unknown_query_field",
                   allowed=list(QUERY_FIELDS), got=unknown)
    names = raw.get("names")
    if isinstance(names, str):
        names = [n for n in names.split(",") if n.strip()]
    terms_raw = raw.get("terms")
    if names is None and terms_raw is not None:
        names = []
    if not isinstance(names, list) or not (names or terms_raw):
        raise _bad("query.names must be a non-empty list", "names_required")
    if not all(isinstance(n, str) for n in names):
        raise _bad("query.names must contain strings", "validation_error", field="names")
    case = raw.get("case_sensitive", False)
    if not isinstance(case, bool):
        raise _bad("query.case_sensitive must be true or false", "validation_error",
                   field="case_sensitive")
    norm = [n.strip() if case else n.strip().casefold() for n in names if n.strip()]
    if not norm and not terms_raw:
        raise _bad("query.names must be a non-empty list", "names_required")
    terms = []
    if terms_raw is not None:
        # Validated here so a bad term fails the request, not the job; stored canonically so the
        # same request always hashes to the same result_id.
        terms = sorted((t.as_dict() for t in parse_terms(None, terms_raw, max_terms=MAX_NAMES_POST)),
                       key=lambda t: json.dumps(t, sort_keys=True))
    alias = raw.get("alias") if raw.get("alias") is not None else "suggest"
    if not isinstance(alias, str) or alias.lower() not in ALIAS_MODES:
        raise _bad(f"query.alias must be one of {', '.join(ALIAS_MODES)}", "invalid_alias_mode",
                   field="alias", value=alias)
    annotation_ids = raw.get("annotation_ids")
    if annotation_ids is not None and (
            not isinstance(annotation_ids, list) or not all(isinstance(a, str) for a in annotation_ids)):
        raise _bad("query.annotation_ids must be a list of strings", "validation_error",
                   field="annotation_ids")
    return {
        "names": sorted(set(norm)),
        "terms": terms,
        "include_taxids": _int_list(raw.get("include_taxids"), "include_taxids"),
        "exclude_taxids": _int_list(raw.get("exclude_taxids"), "exclude_taxids"),
        "alias": alias.lower(),
        "annotation_ids": sorted(set(annotation_ids or [])),
        "case_sensitive": case,
    }


def parse_request(body: dict) -> SeqRequest:
    if not isinstance(body, dict):
        raise _bad("request body must be a JSON object", "validation_error")
    sources = [k for k in ("query", "loci", "items") if body.get(k) is not None]
    if len(sources) != 1:
        raise _bad("give exactly one of query, loci or items", "exactly_one_source", got=sources)
    source = sources[0]

    action_raw = body.get("action")
    try:
        action = SequenceAction(action_raw)
    except ValueError:
        raise _bad("action must be cds, transcript or protein", "validation_error",
                   field="action", value=action_raw) from None

    try:
        buffer_bp = int(body.get("buffer_bp") or 0)
    except (TypeError, ValueError):
        raise _bad("buffer_bp must be an integer", "validation_error", field="buffer_bp") from None
    if not 0 <= buffer_bp <= MAX_SEQUENCES_BUFFER_BP:
        raise _bad(f"buffer_bp must be between 0 and {MAX_SEQUENCES_BUFFER_BP}",
                   "validation_error", field="buffer_bp", max_buffer_bp=MAX_SEQUENCES_BUFFER_BP)
    if buffer_bp and action is not SequenceAction.transcript:
        raise _bad("buffer_bp is only supported with action=transcript", "buffer_not_supported",
                   action=action.value, buffer_bp=buffer_bp)

    spec = body.get("transcript")
    try:
        rule = parse_rule(spec, action.value)
    except RuleError as exc:
        raise _bad(str(exc), exc.code, **exc.extra) from exc

    on_error = body.get("on_error") or ("skip" if source == "query" else "fail")
    if on_error not in ON_ERROR:
        raise _bad("on_error must be skip or fail", "validation_error", field="on_error")

    hard = min(int(settings.SEQUENCES_MAX_GENES), int(settings.SEQUENCES_MAX_SEQUENCES))
    try:
        max_genes = hard if body.get("max_genes") is None else int(body["max_genes"])
    except (TypeError, ValueError):
        raise _bad("max_genes must be an integer", "validation_error", field="max_genes") from None
    if max_genes < 1:
        raise _bad("max_genes must be >= 1", "validation_error", field="max_genes")
    max_genes = min(max_genes, hard)

    return SeqRequest(
        action=action,
        rule=rule,
        rule_spec=spec if spec is not None else "longest",
        buffer_bp=buffer_bp,
        include_exon_coords=bool(body.get("include_exon_coords", False)),
        on_error=on_error,
        max_genes=max_genes,
        query=_query(body["query"]) if source == "query" else None,
        loci=_entries(body[source], source) if source != "query" else None,
    )
