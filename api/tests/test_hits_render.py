"""The JSON joined from stored pieces must equal the plain Python answer, for every option."""
from __future__ import annotations

import json

import pytest
from fastapi.encoders import jsonable_encoder

import helpers.presence as presence_mod
from services.presence_hits import list_name_hits_presence, render_hits_json
from tests.test_presence_api import presence_client  # noqa: F401

CASES = [
    dict(names="tp53"),
    dict(names="tp53,brca1,p53"),
    dict(names="TP53,selenow", shape="records"),
    dict(names="tp53,brca1,p53,nothere", shape="records", alias="include"),
    dict(names="tp53,brca1", alias="off"),
    dict(names="seleno*,selenow", shape="records"),
    dict(names="tp53,brca1", min_matched=2),
    dict(names="tp53,brca1", min_matched=2, shape="records"),
    dict(names="tp53", sources="RefSeq"),
    dict(names="tp53,p53", sources="Ensembl,RefSeq", shape="records"),
    dict(names="tp53", include_taxids="40674"),
    dict(names="tp53", include_taxids="40674", exclude_taxids="10090", shape="records"),
    dict(names="nothere"),
    dict(names="nothere", shape="records"),
]


@pytest.mark.parametrize("kw", CASES, ids=lambda kw: str(kw))
def test_render_equals_plain_answer(presence_client, kw):  # noqa: F811
    reader = presence_mod.presence_reader()
    plain = jsonable_encoder(list_name_hits_presence(reader, **kw))
    assert json.loads(render_hits_json(reader, **kw)) == plain


def test_answer_order_and_unicode(presence_client):  # noqa: F811
    reader = presence_mod.presence_reader()
    body = json.loads(render_hits_json(reader, "tp53,brca1"))
    scores = [len(r[3]) for r in body["rows"]]
    assert scores == sorted(scores, reverse=True)
    raw = render_hits_json(reader, "tp53").decode()
    assert ", " not in raw.split('"rows"')[1]  # compact separators


def test_count_loci_uses_totals_and_falls_back_for_older_indexes():
    import sqlite3
    from services.presence_hits import _count_loci

    def make(with_sym_column: bool):
        c = sqlite3.connect(":memory:")
        c.execute("CREATE TABLE locus (name_id, ord, local_id, via)")
        c.executemany("INSERT INTO locus VALUES (7,?,?,?)", [(1, 1, 0), (2, 1, 0), (2, 2, 1)])
        if with_sym_column:
            c.execute("CREATE TABLE name_loci (name_id, n_loci, n_sym_loci)")
            c.execute("INSERT INTO name_loci VALUES (7, 3, 2)")
        else:
            c.execute("CREATE TABLE name_loci (name_id, n_loci)")
            c.execute("INSERT INTO name_loci VALUES (7, 3)")
        return c

    for new in (True, False):
        c = make(new)
        assert _count_loci(c, 7, "suggest", None) == 2
        assert _count_loci(c, 7, "include", None) == 3
        if new:   # prove the totals are what is read: tamper with the rows
            c.execute("DELETE FROM locus")
            assert _count_loci(c, 7, "suggest", None) == 2
