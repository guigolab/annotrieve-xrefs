"""InterPro label parse + load against xref_meta."""
from __future__ import annotations

import sqlite3
import tempfile
from pathlib import Path

from gene_corpus.labels.load import load_namespace_labels
from gene_corpus.labels.providers.interpro import (
    InterproLabelProvider,
    _iter_entry_list,
    _load_short_names,
)
from gene_corpus.sql.schema import connect_for_build, init_schema


def _write_interpro_fixture(cache: Path) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    (cache / "entry.list").write_text(
        "ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\n"
        "IPR000008\tDomain\tC2 domain\n"
        "IPR000126\tActive_site\t"
        "Serine proteases, V8 family, serine active site\n"
        "IPR000001\tDomain\tKringle\n",
        encoding="utf-8",
    )
    (cache / "short_names.dat").write_text(
        "IPR000008\tC2_dom\n"
        "IPR000126\tV8_ser_AS\n"
        "IPR000001\tKringle\n",
        encoding="utf-8",
    )


def _seed_meta(conn: sqlite3.Connection, accessions: list[str]) -> None:
    conn.executemany(
        "INSERT INTO xref_meta (namespace, accession, n_annotations, n_loci) "
        "VALUES ('interpro', ?, 1, 1)",
        [(a,) for a in accessions],
    )
    conn.commit()


def test_iter_entry_list_skips_type_uses_short_alias(tmp_path: Path) -> None:
    cache = tmp_path / "interpro"
    _write_interpro_fixture(cache)
    shorts = _load_short_names(cache / "short_names.dat")
    by_id = {
        t.accession: t
        for t in _iter_entry_list(cache / "entry.list", short_names=shorts)
    }
    assert by_id["IPR000008"].label == "C2 domain"
    assert by_id["IPR000008"].aliases == ("C2_dom",)
    assert "Domain" not in by_id["IPR000008"].aliases
    assert by_id["IPR000126"].label.startswith("Serine proteases")
    assert by_id["IPR000126"].aliases == ("V8_ser_AS",)
    # short name equals label → no alias
    assert by_id["IPR000001"].label == "Kringle"
    assert by_id["IPR000001"].aliases == ()


def test_iter_terms_directory(tmp_path: Path) -> None:
    cache = tmp_path / "interpro"
    _write_interpro_fixture(cache)
    terms = list(InterproLabelProvider().iter_terms(cache))
    assert {t.accession for t in terms} == {
        "IPR000008",
        "IPR000126",
        "IPR000001",
    }


def test_load_writes_labels_for_meta_ids() -> None:
    with tempfile.TemporaryDirectory() as d:
        work = Path(d)
        db = work / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        _seed_meta(conn, ["IPR000008", "IPR000126", "IPR999999"])
        conn.close()

        cache = work / "fixture_interpro"
        _write_interpro_fixture(cache)
        stats = load_namespace_labels(
            work, "interpro", source_path=cache, force_download=False
        )
        assert stats.labels_written == 2
        assert stats.meta_missing_label == 1
        assert stats.fts_written == stats.labels_written + stats.aliases_written

        conn = sqlite3.connect(str(db))
        labels = {
            (acc, lab)
            for acc, lab in conn.execute(
                "SELECT accession, label FROM xref_label "
                "WHERE namespace='interpro'"
            )
        }
        assert ("IPR000008", "C2 domain") in labels
        aliases = {
            alias
            for (alias,) in conn.execute(
                "SELECT alias FROM xref_label_alias "
                "WHERE namespace='interpro' AND accession='IPR000008'"
            )
        }
        assert aliases == {"C2_dom"}
        assert "Domain" not in aliases
        conn.close()
