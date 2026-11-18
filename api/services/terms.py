"""Search terms shared by /hits, /loci and /sequences: exact names and prefix terms (``seleno*``).

A *term* is what the user typed. It resolves to one or more index names (``name_dict``): an exact
name gives one, a prefix gives every symbol that starts with it. A term matches an annotation or a
gene if any of its names does, so coverage (``min_matched``) is counted over terms, not names.

Prefix terms may carry options (POST only): ``exclude`` (names to drop), ``categories`` (keep these
name categories, as ``GET /names`` reports them) and ``split`` (one term per expanded name, so
coverage is counted per name). Limits: at least 3 characters, at most 500 names per prefix term and
2,000 expanded names per request; beyond that the request is refused with ``prefix_too_broad``.
"""
from __future__ import annotations

import sqlite3
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from helpers.cursor import QueryError
from helpers.name_category import NAME_CATEGORIES, classify_name
from helpers.names import normalize_name_token, parse_exact_name
from helpers.pagination import MAX_NAMES_GET, MAX_PREFIX_EXPAND, MIN_PREFIX_LEN

MAX_EXPANDED_NAMES = 2000
_SCAN_LIMIT = 5000                    # rows read from name_dict before a prefix is called too broad
_TERM_FIELDS = ("q", "exclude", "categories", "split")
_IN_CHUNK = 900


@dataclass(frozen=True)
class RawTerm:
    q: str
    exclude: tuple[str, ...] = ()
    categories: tuple[str, ...] = ()
    split: bool = False

    @property
    def is_prefix(self) -> bool:
        return "*" in self.q

    def key(self) -> tuple:
        return (normalize_name_token(self.q), self.exclude, self.categories, self.split)

    def as_dict(self) -> dict | str:
        """Canonical form (also used to hash a request): a bare string unless options are set."""
        if not (self.exclude or self.categories or self.split):
            return self.q
        out: dict[str, Any] = {"q": self.q}
        if self.exclude:
            out["exclude"] = list(self.exclude)
        if self.categories:
            out["categories"] = list(self.categories)
        if self.split:
            out["split"] = True
        return out


def _from_item(item: Any) -> RawTerm:
    if isinstance(item, str):
        return RawTerm(item.strip())
    if not isinstance(item, dict):
        raise QueryError("each term must be a string or an object", code="invalid_terms")
    unknown = sorted(set(item) - set(_TERM_FIELDS))
    if unknown:
        raise QueryError(f"unknown term field(s): {', '.join(unknown)}", code="unknown_term_field",
                         unknown=unknown, allowed=list(_TERM_FIELDS))
    q = item.get("q")
    if not isinstance(q, str) or not q.strip():
        raise QueryError("term.q must be a non-empty string", code="invalid_terms", field="q")
    exclude = item.get("exclude") or []
    if not isinstance(exclude, list) or not all(isinstance(x, str) for x in exclude):
        raise QueryError("term.exclude must be a list of names", code="invalid_terms",
                         field="exclude")
    categories = item.get("categories") or []
    if not isinstance(categories, list) or not all(isinstance(x, str) for x in categories):
        raise QueryError("term.categories must be a list of category names", code="invalid_terms",
                         field="categories")
    bad = sorted({c.strip().casefold() for c in categories} - set(NAME_CATEGORIES))
    if bad:
        raise QueryError(f"unknown categories: {', '.join(bad)}", code="invalid_category",
                         unknown=bad, allowed=list(NAME_CATEGORIES))
    split = item.get("split", False)
    if not isinstance(split, bool):
        raise QueryError("term.split must be true or false", code="invalid_terms", field="split")
    term = RawTerm(
        q.strip(),
        exclude=tuple(sorted({normalize_name_token(x) for x in exclude if x.strip()})),
        categories=tuple(sorted({c.strip().casefold() for c in categories})),
        split=split)
    if (term.exclude or term.categories or term.split) and not term.is_prefix:
        raise QueryError("exclude, categories and split apply to prefix terms (name*) only",
                         code="option_needs_prefix", term=term.q)
    return term


def parse_terms(names: str | Sequence[str] | None = None, terms: Sequence[Any] | None = None,
                *, max_terms: int = MAX_NAMES_GET) -> list[RawTerm]:
    """``names`` (CSV string or list of strings) and/or ``terms`` (strings or option objects)."""
    items: list[Any] = []
    if names is not None:
        pieces = (names,) if isinstance(names, str) else names
        for piece in pieces:
            for part in str(piece).split(","):
                if part.strip():
                    items.append(part.strip())
    if terms is not None:
        if not isinstance(terms, list):
            raise QueryError("terms must be a list", code="invalid_terms")
        items.extend(terms)
    out: list[RawTerm] = []
    seen: set[tuple] = set()
    for item in items:
        term = _from_item(item)
        if term.key() in seen:
            continue
        seen.add(term.key())
        out.append(term)
    if not out:
        raise QueryError("at least one name is required", code="names_required")
    if len(out) > max_terms:
        raise QueryError(f"names accepts at most {max_terms} values", code="too_many_names",
                         max=max_terms, received=len(out))
    return out


def _stem(term: RawTerm) -> tuple[str, str]:
    """(typed stem, folded stem) of a prefix term; ``*`` must be the single trailing character."""
    text = term.q
    if text.count("*") != 1 or not text.endswith("*"):
        raise QueryError("a wildcard is only allowed once, as the last character (seleno*)",
                         code="invalid_prefix", term=text)
    typed = text[:-1].strip()
    if len(normalize_name_token(typed)) < MIN_PREFIX_LEN:
        raise QueryError(f"a prefix needs at least {MIN_PREFIX_LEN} characters before the *",
                         code="prefix_too_short", min_length=MIN_PREFIX_LEN, term=text)
    return typed, normalize_name_token(typed)


def _lookup_exact(conn: sqlite3.Connection, names: list[str]) -> dict[str, int]:
    found: dict[str, int] = {}
    for i in range(0, len(names), _IN_CHUNK):
        chunk = names[i:i + _IN_CHUNK]
        ph = ",".join("?" for _ in chunk)
        for row in conn.execute(f"SELECT name, name_id FROM name_dict WHERE name IN ({ph})", chunk):
            found[row[0]] = int(row[1])
    return found


def _expand_prefix(conn: sqlite3.Connection, term: RawTerm, stem: str, mode: str
                   ) -> list[tuple[str, int]]:
    rows = conn.execute(
        "SELECT name, name_id, n_sym FROM name_dict WHERE name >= ? AND name < ? "
        "ORDER BY name LIMIT ?", (stem, stem + "￿", _SCAN_LIMIT + 1)).fetchall()
    scanned_all = len(rows) <= _SCAN_LIMIT
    keep = []
    for name, name_id, n_sym in rows[:_SCAN_LIMIT]:
        if mode != "include" and not n_sym:          # alias-only names stay hidden, as in /names
            continue
        if name in term.exclude:
            continue
        if term.categories and classify_name(name) not in term.categories:
            continue
        keep.append((name, int(name_id)))
    if len(keep) > MAX_PREFIX_EXPAND or not scanned_all:
        raise QueryError(
            f"{term.q!r} expands to more than {MAX_PREFIX_EXPAND} names; use a longer prefix, "
            "`categories` or `exclude`",
            code="prefix_too_broad", max=MAX_PREFIX_EXPAND, n_names=len(keep),
            exact_count=scanned_all, term=term.q)
    return keep


def resolve_terms(conn: sqlite3.Connection, raws: Sequence[RawTerm], mode: str) -> list[dict]:
    """Resolve terms against ``name_dict``. Every returned dict has:

    ``index``, ``input`` (as typed), ``kind`` (exact | prefix), ``name`` (folded name or stem),
    ``typed`` (the user's spelling, for case-sensitive matching), ``names`` and ``name_ids``
    (what the term expands to; empty when nothing matches), ``n_names`` and, for ``split`` children,
    ``from_prefix``.
    """
    exact_names = []
    for raw in raws:
        if not raw.is_prefix:
            exact_names.append(parse_exact_name(raw.q)[1])
    ids = _lookup_exact(conn, exact_names)

    out: list[dict] = []
    total = 0
    for raw in raws:
        if not raw.is_prefix:
            text, name = parse_exact_name(raw.q)
            nid = ids.get(name)
            out.append({"input": text, "kind": "exact", "name": name, "typed": text,
                        "names": [name] if nid is not None else [],
                        "name_ids": [nid] if nid is not None else [],
                        "n_names": 1 if nid is not None else 0, "from_prefix": None})
            total += 1
            continue
        typed, stem = _stem(raw)
        found = _expand_prefix(conn, raw, stem, mode)
        total += len(found)
        if raw.split and found:
            for name, nid in found:
                out.append({"input": name, "kind": "exact", "name": name, "typed": typed,
                            "names": [name], "name_ids": [nid], "n_names": 1,
                            "from_prefix": raw.q})
        else:
            out.append({"input": raw.q, "kind": "prefix", "name": stem, "typed": typed,
                        "names": [n for n, _ in found], "name_ids": [i for _, i in found],
                        "n_names": len(found), "from_prefix": None})
        if total > MAX_EXPANDED_NAMES:
            raise QueryError(
                f"the terms expand to more than {MAX_EXPANDED_NAMES} names in total",
                code="prefix_too_broad", max=MAX_EXPANDED_NAMES, n_names=total)
    for i, term in enumerate(out):
        term["index"] = i
    return out


def term_echo(term: dict) -> dict:
    """The expansion as a response reports it (no internal ids)."""
    echo = {"index": term["index"], "input": term["input"], "kind": term["kind"],
            "name": term["name"], "names": term["names"], "n_names": term["n_names"],
            "found": bool(term["name_ids"])}
    if term["kind"] == "prefix":
        echo["truncated"] = False
    if term.get("from_prefix"):
        echo["from_prefix"] = term["from_prefix"]
    return echo
