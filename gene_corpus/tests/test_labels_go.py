"""GO OBO parse + xref_label load against xref_meta."""
from __future__ import annotations

import sqlite3
import tempfile
from pathlib import Path

from gene_corpus.labels.load import load_namespace_labels
from gene_corpus.labels.providers.go import (
    MAX_ALIASES_PER_TERM,
    GoLabelProvider,
    _parse_term_stanza,
    _select_aliases,
)
from gene_corpus.sql.schema import SCHEMA_VERSION, connect_for_build, init_schema

_OBO_FIXTURE = """\
format-version: 1.2
ontology: go

[Term]
id: GO:0000010
name: heptaprenyl diphosphate synthase activity
namespace: molecular_function
alt_id: GO:0036422
def: "Catalysis of a reaction." [PMID:1]
synonym: "HepPP synthase activity" RELATED [EC:2.5.1.30]
synonym: "all-trans-heptaprenyl-diphosphate synthase activity" EXACT [EC:2.5.1.30]
is_a: GO:0120531 ! prenyl diphosphate synthase activity

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process
is_obsolete: true

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Typedef]
id: part_of
name: part of
"""


def _seed_meta(conn: sqlite3.Connection, accessions: list[str]) -> None:
    conn.executemany(
        "INSERT INTO xref_meta (namespace, accession, n_annotations, n_loci) "
        "VALUES ('go', ?, 1, 1)",
        [(a,) for a in accessions],
    )
    conn.commit()


def test_schema_version_stays_one() -> None:
    assert SCHEMA_VERSION == 1


def test_parse_term_stanza_maps_fields() -> None:
    term = _parse_term_stanza(
        [
            "id: GO:0000010",
            "name: heptaprenyl diphosphate synthase activity",
            "namespace: molecular_function",
            "alt_id: GO:0036422",
            'synonym: "HepPP synthase activity" RELATED [EC:2.5.1.30]',
            'synonym: "all-trans-heptaprenyl-diphosphate synthase activity" '
            "EXACT [EC:2.5.1.30]",
            "def: \"Catalysis.\" [PMID:1]",
            "is_a: GO:0120531 ! prenyl diphosphate synthase activity",
        ]
    )
    assert term is not None
    assert term.accession == "GO:0000010"
    assert term.label == "heptaprenyl diphosphate synthase activity"
    assert term.alt_ids == ("GO:0036422",)
    assert "HepPP synthase activity" in term.aliases
    assert (
        "all-trans-heptaprenyl-diphosphate synthase activity" in term.aliases
    )


def test_parse_keeps_obsolete() -> None:
    term = _parse_term_stanza(
        [
            "id: GO:0008150",
            "name: biological_process",
            "is_obsolete: true",
        ]
    )
    assert term is not None
    assert term.accession == "GO:0008150"
    assert term.label == "biological_process"


def test_select_aliases_prefers_exact_related_and_caps() -> None:
    scoped = [
        ("narrow one", "NARROW"),
        ("exact one", "EXACT"),
        ("related one", "RELATED"),
        ("exact two", "EXACT"),
        ("narrow two", "NARROW"),
        ("related two", "RELATED"),
        ("exact three", "EXACT"),
        ("broad one", "BROAD"),
    ]
    picked = _select_aliases(scoped, limit=5)
    assert len(picked) == MAX_ALIASES_PER_TERM
    assert picked == [
        "exact one",
        "exact two",
        "exact three",
        "related one",
        "related two",
    ]
    assert "narrow one" not in picked
    assert "broad one" not in picked


def test_parse_caps_aliases_preferring_exact() -> None:
    lines = [
        "id: GO:1904709",
        "name: negative regulation of granulosa cell apoptotic process",
        'synonym: "down regulation of granulosa cell apoptosis" NARROW []',
        'synonym: "down regulation of granulosa cell apoptotic process" EXACT []',
        'synonym: "inhibition of granulosa cell apoptosis" NARROW []',
        'synonym: "negative regulation of granulosa cell of ovary '
        'apoptotic process" EXACT []',
        'synonym: "down-regulation of granulosa cell apoptotic process" EXACT []',
        'synonym: "downregulation of granulosa cell apoptotic process" EXACT []',
        'synonym: "HepPP-like short" RELATED []',
        'synonym: "another exact" EXACT []',
        'synonym: "yet another exact" EXACT []',
    ]
    term = _parse_term_stanza(lines)
    assert term is not None
    assert len(term.aliases) == MAX_ALIASES_PER_TERM
    for alias in term.aliases:
        assert "inhibition" not in alias
        assert alias != "down regulation of granulosa cell apoptosis"


def test_iter_terms_streams_fixture(tmp_path: Path) -> None:
    obo = tmp_path / "go-basic.obo"
    obo.write_text(_OBO_FIXTURE, encoding="utf-8")
    terms = list(GoLabelProvider().iter_terms(obo))
    by_id = {t.accession: t for t in terms}
    assert "GO:0000010" in by_id
    assert "GO:0003674" in by_id
    assert "GO:0008150" in by_id  # obsolete kept for corpus labels
    assert by_id["GO:0008150"].label == "biological_process"


def test_load_writes_labels_and_aliases_for_meta_ids() -> None:
    with tempfile.TemporaryDirectory() as d:
        work = Path(d)
        db = work / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        _seed_meta(conn, ["GO:0000010", "GO:0036422", "GO:0008150"])
        conn.close()

        obo = work / "fixture.obo"
        obo.write_text(_OBO_FIXTURE, encoding="utf-8")

        stats = load_namespace_labels(
            work, "go", source_path=obo, force_download=False
        )
        assert stats.labels_written == 3
        assert stats.aliases_written == 4  # 2 synonyms × 2 accessions
        assert stats.meta_missing_label == 0

        conn = sqlite3.connect(str(db))
        labels = {
            (ns, acc, lab)
            for ns, acc, lab in conn.execute(
                "SELECT namespace, accession, label FROM xref_label"
            )
        }
        assert (
            "go",
            "GO:0000010",
            "heptaprenyl diphosphate synthase activity",
        ) in labels
        assert (
            "go",
            "GO:0036422",
            "heptaprenyl diphosphate synthase activity",
        ) in labels
        assert ("go", "GO:0008150", "biological_process") in labels
        # Absent from meta must not appear.
        accessions = {acc for _, acc, _ in labels}
        assert "GO:0003674" not in accessions

        aliases = {
            alias
            for (_, _, alias) in conn.execute(
                "SELECT namespace, accession, alias FROM xref_label_alias "
                "WHERE accession = 'GO:0000010'"
            )
        }
        assert "HepPP synthase activity" in aliases
        ver = conn.execute(
            "SELECT value FROM meta WHERE key='schema_version'"
        ).fetchone()[0]
        assert ver == "1"
        conn.close()


def test_load_skips_term_not_in_meta() -> None:
    with tempfile.TemporaryDirectory() as d:
        work = Path(d)
        db = work / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        _seed_meta(conn, ["GO:0009999"])
        conn.close()

        obo = work / "fixture.obo"
        obo.write_text(_OBO_FIXTURE, encoding="utf-8")
        stats = load_namespace_labels(work, "go", source_path=obo)
        assert stats.labels_written == 0
        assert stats.meta_missing_label == 1
