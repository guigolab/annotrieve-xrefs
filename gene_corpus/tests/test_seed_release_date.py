"""One-time annotation.release_date seed from report TSV."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

from gene_corpus.catalog.report import download_annotation_report
from gene_corpus.seed_release_date_cli import main as seed_main
from gene_corpus.sql.schema import connect_for_build, init_schema, seed_annotations
from gene_corpus.sql.seed_release_date import (
    REPORT_RELEASE_DATE_FIELDS,
    annotation_has_release_date_column,
    ensure_release_date_column,
    parse_release_date_tsv,
    seed_release_dates,
)


class _FakeResponse:
    def __init__(self, body: bytes):
        self._body = body

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def read(self, size: int = -1) -> bytes:
        if size < 0:
            data, self._body = self._body, b""
            return data
        data, self._body = self._body[:size], self._body[size:]
        return data


def _legacy_db(path: Path) -> None:
    """DB whose annotation table lacks release_date (pre-schema bump)."""
    conn = connect_for_build(path)
    conn.executescript(
        """
        CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE annotation (
            annotation_key INTEGER PRIMARY KEY,
            annotation_id TEXT UNIQUE NOT NULL,
            taxid INTEGER NOT NULL,
            assembly_accession TEXT NOT NULL,
            organism_name TEXT,
            source_database TEXT NOT NULL,
            source_provider TEXT,
            profile_id TEXT NOT NULL,
            gff_path TEXT
        );
        """
    )
    seed_annotations(
        conn,
        [
            {
                "annotation_key": 1,
                "annotation_id": "aaa",
                "taxid": 9606,
                "assembly_accession": "GCA_1",
                "organism_name": "human",
                "source_database": "RefSeq",
                "source_provider": None,
                "profile_id": "refseq",
                "gff_path": "/a.gff.gz",
            },
            {
                "annotation_key": 2,
                "annotation_id": "bbb",
                "taxid": 10090,
                "assembly_accession": "GCA_2",
                "organism_name": "mouse",
                "source_database": "Ensembl",
                "source_provider": None,
                "profile_id": "ensembl",
                "gff_path": "/b.gff.gz",
            },
        ],
    )
    conn.close()


def test_ensure_column_and_seed_apply() -> None:
    with tempfile.TemporaryDirectory() as d:
        db = Path(d) / "gene_corpus.sqlite"
        _legacy_db(db)
        conn = connect_for_build(db)
        assert not annotation_has_release_date_column(conn)
        assert ensure_release_date_column(conn) is True
        assert annotation_has_release_date_column(conn)
        assert ensure_release_date_column(conn) is False

        stats = seed_release_dates(
            conn,
            {"aaa": "2024-01-02T00:00:00", "ccc": "2020-01-01"},
            apply=True,
            report_rows=2,
            empty_in_report=0,
        )
        assert stats.matched == 1
        assert stats.updated == 1
        assert stats.missing_in_db == 1
        row = conn.execute(
            "SELECT release_date FROM annotation WHERE annotation_id = 'aaa'"
        ).fetchone()
        assert row[0] == "2024-01-02T00:00:00"
        row_b = conn.execute(
            "SELECT release_date FROM annotation WHERE annotation_id = 'bbb'"
        ).fetchone()
        assert row_b[0] is None
        conn.close()


def test_parse_release_date_tsv() -> None:
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "r.tsv"
        path.write_text(
            "annotation_id\trelease_date\tdatabase\n"
            "aaa\t2024-01-02T00:00:00\tRefSeq\n"
            "bbb\t\tEnsembl\n"
            "ccc\t2020-06-01\tGenBank\n",
            encoding="utf-8",
        )
        mapping, total, empty = parse_release_date_tsv(path)
        assert total == 3
        assert empty == 1
        assert mapping == {
            "aaa": "2024-01-02T00:00:00",
            "ccc": "2020-06-01",
        }


def test_download_selected_fields_override() -> None:
    captured: dict = {}

    def opener(request, timeout=None):
        captured["body"] = json.loads(request.data.decode("utf-8"))
        return _FakeResponse(
            b"annotation_id\trelease_date\tdatabase\n"
            b"aaa\t2024-01-02\tRefSeq\n"
        )

    with tempfile.TemporaryDirectory() as d:
        dest = Path(d) / "r.tsv"
        download_annotation_report(
            "https://example.test/api/v0",
            dest,
            selected_fields=REPORT_RELEASE_DATE_FIELDS,
            opener=opener,
        )
        assert captured["body"]["selected_fields"] == REPORT_RELEASE_DATE_FIELDS


def test_cli_apply_from_local_report() -> None:
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        db = root / "gene_corpus.sqlite"
        _legacy_db(db)
        report = root / "report.tsv"
        report.write_text(
            "annotation_id\trelease_date\n"
            "aaa\t2024-05-01T12:00:00\n"
            "bbb\t2023-01-01\n",
            encoding="utf-8",
        )
        assert seed_main(["--work-dir", str(root), "--report", str(report)]) == 0
        conn = connect_for_build(db)
        assert not annotation_has_release_date_column(conn)
        conn.close()

        assert (
            seed_main(
                [
                    "--work-dir",
                    str(root),
                    "--report",
                    str(report),
                    "--apply",
                ]
            )
            == 0
        )
        conn = connect_for_build(db)
        rows = dict(
            conn.execute(
                "SELECT annotation_id, release_date FROM annotation"
            ).fetchall()
        )
        assert rows["aaa"] == "2024-05-01T12:00:00"
        assert rows["bbb"] == "2023-01-01"
        ver = conn.execute(
            "SELECT value FROM meta WHERE key = 'schema_version'"
        ).fetchone()
        assert ver is not None
        conn.close()


def test_init_schema_includes_release_date() -> None:
    with tempfile.TemporaryDirectory() as d:
        db = Path(d) / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        assert annotation_has_release_date_column(conn)
        conn.close()
