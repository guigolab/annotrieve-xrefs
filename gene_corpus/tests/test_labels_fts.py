"""xref_label_fts trigram: substring MATCH after labels load."""
from __future__ import annotations

import sqlite3
import tempfile
from pathlib import Path

import pytest

from gene_corpus.labels.load import load_namespace_labels
from gene_corpus.sql.schema import connect_for_build, init_schema

# Minimal OBO with a mid-string token for trigram substring search.
_OBO_APOPTOTIC = """\
format-version: 1.2
ontology: go

[Term]
id: GO:1904709
name: negative regulation of granulosa cell apoptotic process
namespace: biological_process
synonym: "down-regulation of granulosa cell apoptosis" EXACT []
is_a: GO:0006915 ! apoptotic process

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function
"""


def _seed_meta(conn: sqlite3.Connection, accessions: list[str]) -> None:
    conn.executemany(
        "INSERT INTO xref_meta (namespace, accession, n_annotations, n_loci) "
        "VALUES ('go', ?, 1, 1)",
        [(a,) for a in accessions],
    )
    conn.commit()


def _require_trigram(conn: sqlite3.Connection) -> None:
    try:
        conn.execute(
            "CREATE VIRTUAL TABLE IF NOT EXISTS _trigram_probe "
            "USING fts5(x, tokenize='trigram')"
        )
        conn.execute("DROP TABLE IF EXISTS _trigram_probe")
    except sqlite3.OperationalError as exc:
        pytest.skip(f"SQLite FTS5 trigram unavailable (need 3.34+): {exc}")


def test_fts_substring_match_and_reload_replaces() -> None:
    with tempfile.TemporaryDirectory() as d:
        work = Path(d)
        db = work / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        _require_trigram(conn)
        _seed_meta(conn, ["GO:1904709"])
        conn.close()

        obo = work / "fixture.obo"
        obo.write_text(_OBO_APOPTOTIC, encoding="utf-8")

        stats = load_namespace_labels(
            work, "go", source_path=obo, force_download=False
        )
        assert stats.labels_written == 1
        assert stats.aliases_written == 1
        assert stats.fts_written == 2  # label + alias

        conn = sqlite3.connect(str(db))
        mid = [
            row[0]
            for row in conn.execute(
                "SELECT accession FROM xref_label_fts "
                "WHERE namespace = ? AND text MATCH ?",
                ("go", "apoptotic"),
            )
        ]
        assert mid == ["GO:1904709"]

        prefix = [
            row[0]
            for row in conn.execute(
                "SELECT accession FROM xref_label_fts "
                "WHERE namespace = ? AND text MATCH ?",
                ("go", "negative"),
            )
        ]
        assert "GO:1904709" in prefix

        alias_hit = [
            row[0]
            for row in conn.execute(
                "SELECT accession FROM xref_label_fts "
                "WHERE namespace = ? AND text MATCH ?",
                ("go", "apoptosis"),
            )
        ]
        assert alias_hit == ["GO:1904709"]

        n_before = conn.execute(
            "SELECT COUNT(*) FROM xref_label_fts WHERE namespace = 'go'"
        ).fetchone()[0]
        assert n_before == 2
        conn.close()

        # Re-load must replace, not duplicate.
        stats2 = load_namespace_labels(
            work, "go", source_path=obo, force_download=False
        )
        assert stats2.fts_written == 2

        conn = sqlite3.connect(str(db))
        n_after = conn.execute(
            "SELECT COUNT(*) FROM xref_label_fts WHERE namespace = 'go'"
        ).fetchone()[0]
        assert n_after == 2
        conn.close()


def test_init_schema_creates_xref_label_fts() -> None:
    with tempfile.TemporaryDirectory() as d:
        conn = connect_for_build(Path(d) / "db.sqlite")
        init_schema(conn)
        _require_trigram(conn)
        row = conn.execute(
            "SELECT sql FROM sqlite_master WHERE name = 'xref_label_fts'"
        ).fetchone()
        assert row is not None
        assert "fts5" in row[0].lower()
        assert "trigram" in row[0].lower()
        conn.close()
