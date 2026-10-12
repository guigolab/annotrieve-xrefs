"""Intersect ontology terms with xref_meta and write label tables + FTS."""
from __future__ import annotations

import sqlite3
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from gene_corpus.labels.providers import PROVIDERS, LabelProvider, LabelTerm
from gene_corpus.labels.providers.pfam import normalize_pfam_accession
from gene_corpus.sql.schema import connect_for_build, init_schema
from gene_corpus.stream.emit_util import normalize_go_token

_INSERT_BATCH = 1000


@dataclass(frozen=True)
class LoadStats:
    namespace: str
    meta_accessions: int
    labels_written: int
    aliases_written: int
    fts_written: int
    meta_missing_label: int


def _normalize_meta_accession(namespace: str, accession: str) -> str:
    if namespace == "go":
        go = normalize_go_token(accession)
        return go if go is not None else accession.strip()
    if namespace == "pfam":
        # Labels are keyed unversioned (PF00001). The labels API aggregates
        # versioned xref_meta rows under that key (see _labels_meta_join).
        pfam = normalize_pfam_accession(accession)
        return pfam if pfam is not None else accession.strip()
    return accession.strip()


def load_meta_accessions(conn: sqlite3.Connection, namespace: str) -> set[str]:
    rows = conn.execute(
        "SELECT accession FROM xref_meta WHERE namespace = ?",
        (namespace,),
    )
    return {
        _normalize_meta_accession(namespace, str(acc))
        for (acc,) in rows
        if acc is not None and str(acc).strip()
    }


def _write_batches(
    conn: sqlite3.Connection,
    labels: list[tuple[str, str, str]],
    aliases: list[tuple[str, str, str]],
    fts_rows: list[tuple[str, str, str]],
) -> None:
    if labels:
        conn.executemany(
            "INSERT INTO xref_label (namespace, accession, label) "
            "VALUES (?, ?, ?)",
            labels,
        )
    if aliases:
        conn.executemany(
            "INSERT INTO xref_label_alias (namespace, accession, alias) "
            "VALUES (?, ?, ?)",
            aliases,
        )
    if fts_rows:
        conn.executemany(
            "INSERT INTO xref_label_fts (namespace, accession, text) "
            "VALUES (?, ?, ?)",
            fts_rows,
        )
    conn.commit()
    labels.clear()
    aliases.clear()
    fts_rows.clear()


def replace_namespace_labels(
    conn: sqlite3.Connection,
    namespace: str,
    terms: Iterable[LabelTerm],
    meta_accessions: set[str],
) -> LoadStats:
    """
    Replace all label/alias/FTS rows for *namespace* from *terms* ∩ meta.

    Deletes existing rows for the namespace first (idempotent re-run).
    """
    conn.execute(
        "DELETE FROM xref_label_alias WHERE namespace = ?", (namespace,)
    )
    conn.execute("DELETE FROM xref_label WHERE namespace = ?", (namespace,))
    conn.execute("DELETE FROM xref_label_fts WHERE namespace = ?", (namespace,))
    conn.commit()

    covered: set[str] = set()
    label_buf: list[tuple[str, str, str]] = []
    alias_buf: list[tuple[str, str, str]] = []
    fts_buf: list[tuple[str, str, str]] = []
    labels_written = 0
    aliases_written = 0
    fts_written = 0

    for term in terms:
        candidates = {term.accession, *term.alt_ids}
        write_ids = sorted(candidates & meta_accessions)
        if not write_ids:
            continue
        for acc in write_ids:
            label_buf.append((namespace, acc, term.label))
            labels_written += 1
            covered.add(acc)
            fts_buf.append((namespace, acc, term.label))
            fts_written += 1
            for alias in term.aliases:
                alias_buf.append((namespace, acc, alias))
                aliases_written += 1
                fts_buf.append((namespace, acc, alias))
                fts_written += 1
        if (
            len(label_buf) >= _INSERT_BATCH
            or len(alias_buf) >= _INSERT_BATCH
            or len(fts_buf) >= _INSERT_BATCH
        ):
            _write_batches(conn, label_buf, alias_buf, fts_buf)

    _write_batches(conn, label_buf, alias_buf, fts_buf)

    return LoadStats(
        namespace=namespace,
        meta_accessions=len(meta_accessions),
        labels_written=labels_written,
        aliases_written=aliases_written,
        fts_written=fts_written,
        meta_missing_label=len(meta_accessions - covered),
    )


def load_namespace_labels(
    work_dir: Path,
    namespace: str,
    *,
    cache_dir: Path | None = None,
    force_download: bool = False,
    provider: LabelProvider | None = None,
    source_path: Path | None = None,
) -> LoadStats:
    """
    Load labels for one namespace into ``{work_dir}/gene_corpus.sqlite``.

    When *source_path* is set (tests), skip download and parse that file.
    Also rebuilds ``xref_label_fts`` (trigram) for the namespace.
    """
    if provider is None:
        try:
            provider = PROVIDERS[namespace]
        except KeyError as exc:
            raise ValueError(
                f"no label provider registered for namespace: {namespace}"
            ) from exc

    db_path = work_dir / "gene_corpus.sqlite"
    if not db_path.is_file():
        raise FileNotFoundError(f"gene_corpus.sqlite missing: {db_path}")

    cache = cache_dir if cache_dir is not None else work_dir / "label_cache"
    if source_path is not None:
        path = source_path
    else:
        path = provider.download(cache, force=force_download)

    conn = connect_for_build(db_path)
    try:
        init_schema(conn)
        meta = load_meta_accessions(conn, namespace)
        print(
            f"labels {namespace}: {len(meta)} accession(s) in xref_meta",
            file=sys.stderr,
            flush=True,
        )
        stats = replace_namespace_labels(
            conn, namespace, provider.iter_terms(path), meta
        )
    finally:
        conn.close()

    print(
        f"labels {namespace}: wrote labels={stats.labels_written} "
        f"aliases={stats.aliases_written} fts={stats.fts_written} "
        f"meta_missing_label={stats.meta_missing_label}",
        file=sys.stderr,
        flush=True,
    )
    return stats
