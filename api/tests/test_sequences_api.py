"""/sequences, /sequences/estimate, /sequences/plan, jobs and cached results (new index only)."""
from __future__ import annotations

import sqlite3
import time
from pathlib import Path
from unittest.mock import patch

import httpx
import pysam
import pytest
from fastapi.testclient import TestClient

from gene_corpus.serve_build_cli import main as serve_build_main
from gene_corpus.tests.test_presence_from_sidecars import corpus
from helpers.gffread_runner import SequenceAction
from helpers.ncbi_eutils import NcbiEutilsClient, close_ncbi_client, reset_ncbi_client

# Gene 1 (chr1:100-150): tx1a exon 100-150, CDS 110-150 (spliced 51, CDS 41, 1 exon);
#                        tx1b exons 100-120 + 130-150, CDS 100-120 + 130-140 (spliced 42, CDS 32, 2 exons).
GENE1 = """\
chr1	s	gene	100	150	.	+	.	ID=gene-1
chr1	s	mRNA	100	150	.	+	.	ID=tx1a;Parent=gene-1
chr1	s	exon	100	150	.	+	.	Parent=tx1a
chr1	s	CDS	110	150	.	+	0	Parent=tx1a
chr1	s	mRNA	100	150	.	+	.	ID=tx1b;Parent=gene-1
chr1	s	exon	100	120	.	+	.	Parent=tx1b
chr1	s	exon	130	150	.	+	.	Parent=tx1b
chr1	s	CDS	100	120	.	+	0	Parent=tx1b
chr1	s	CDS	130	140	.	+	0	Parent=tx1b
"""
# Gene 2 (indexed as chr1:200-250) is coding and its transcript reaches 300, beyond the gene.
GENE2_CODING = """\
chr1	s	gene	200	250	.	+	.	ID=gene-2
chr1	s	mRNA	200	300	.	+	.	ID=tx2;Parent=gene-2
chr1	s	exon	200	250	.	+	.	Parent=tx2
chr1	s	exon	270	300	.	+	.	Parent=tx2
chr1	s	CDS	210	250	.	+	0	Parent=tx2
chr1	s	CDS	270	290	.	+	0	Parent=tx2
"""
GENE2_NONCODING = """\
chr1	s	gene	200	250	.	+	.	ID=gene-2
chr1	s	lnc_RNA	200	250	.	+	.	ID=tx2;Parent=gene-2
chr1	s	exon	200	250	.	+	.	Parent=tx2
"""


def write_tabix(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted((ln.split("\t") for ln in text.splitlines() if ln.strip()),
                  key=lambda c: (c[0], int(c[3])))
    plain = path.with_suffix(".gff")
    plain.write_text("\n".join("\t".join(c) for c in rows) + "\n")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="gff", force=True)
    plain.unlink()


def fake_gffread(*, gff_path, fasta_path, out_path, action, buffer_bp=0,
                 include_exon_coords=False, timeout_s=None):
    fake_gffread.calls += 1
    ids = []
    for line in Path(gff_path).read_text().splitlines():
        cols = line.split("\t")
        if len(cols) >= 9 and cols[2] in {"mRNA", "transcript", "lnc_RNA"}:
            ids.append(dict(p.split("=", 1) for p in cols[8].split(";") if "=" in p)["ID"])
    with open(out_path, "w") as fh:
        for tid in ids:
            fh.write(f">{tid}\nACGTACGT\n")


fake_gffread.calls = 0


def mock_ncbi() -> NcbiEutilsClient:
    def handler(request: httpx.Request) -> httpx.Response:
        p = dict(request.url.params)
        n = int(p["seq_stop"]) - int(p["seq_start"]) + 1
        return httpx.Response(200, text=">chr1\n" + ("ACGT" * (n // 4 + 1))[:n] + "\n")
    return NcbiEutilsClient(rate_per_sec=1000, cache_mb=0, transport=httpx.MockTransport(handler))


@pytest.fixture()
def seq_client(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import settings as settings_mod
    from helpers import presence as presence_mod
    from main import create_app

    work = corpus(tmp_path)
    assert serve_build_main(["build", "--work-dir", str(work)]) == 0
    assert serve_build_main(["loci", "--work-dir", str(work)]) == 0
    version = (work / "serve" / "current").resolve()
    ann_root = tmp_path / "annotations"
    conn = sqlite3.connect(version / "presence.sqlite")
    for aid, in conn.execute("SELECT annotation_id FROM annotation").fetchall():
        gene2 = GENE2_NONCODING if aid == "ann1" else GENE2_CODING
        write_tabix(ann_root / "9606" / "GCF_000001405.40" / f"{aid}.gff.gz", GENE1 + gene2)
        conn.execute("UPDATE annotation SET gff_path=? WHERE annotation_id=?",
                     (f"/gff/{aid}.gff.gz", aid))
    conn.commit()
    conn.close()

    s = settings_mod.settings
    saved = {k: getattr(s, k) for k in ("PRESENCE_DIR", "ANNOTATIONS_DIR", "SEQUENCES_RESULTS_DIR",
                                        "SEQUENCES_SYNC_MAX_GENES", "SEQUENCES_BATCH_GENES",
                                        "NCBI_WORKERS", "SEQUENCES_DEADLINE_S")}
    s.PRESENCE_DIR = str(work / "serve")
    s.ANNOTATIONS_DIR = str(ann_root)
    s.SEQUENCES_RESULTS_DIR = str(tmp_path / "results")
    s.SEQUENCES_SYNC_MAX_GENES = 50
    s.SEQUENCES_BATCH_GENES = 3
    s.NCBI_WORKERS = 1
    s.SEQUENCES_DEADLINE_S = 60.0
    presence_mod.close_presence()
    close_ncbi_client()
    reset_ncbi_client(mock_ncbi())
    fake_gffread.calls = 0
    with patch("services.seq_export.run_gffread", side_effect=fake_gffread):
        with TestClient(create_app()) as client:
            yield client, work
    close_ncbi_client()
    presence_mod.close_presence()
    for k, v in saved.items():
        setattr(s, k, v)


def headers(fasta: str) -> list[str]:
    return [ln[1:] for ln in fasta.splitlines() if ln.startswith(">")]


TP53 = {"query": {"names": ["tp53"], "alias": "off"}, "action": "cds"}


def test_query_sync_fasta_headers_manifest_and_cache(seq_client):
    client, _ = seq_client
    r = client.post("/sequences", json=TP53)
    assert r.status_code == 200 and r.headers["content-type"].startswith("text/x-fasta")
    manifest = client.get(f"/sequences/results/{r.headers['X-Result-Id']}/manifest.tsv").text
    assert r.headers["X-Sequence-Count"] == "5" and r.headers["X-Skipped-Count"] == "0", manifest
    hs = headers(r.text)
    assert len(hs) == 5 and all(h.startswith("tp53|Homo_sapiens|GCF_000001405.40|ann") for h in hs)
    assert all(h.split("|")[4].split()[0] == "tx1a" for h in hs)          # longest CDS (41 > 32)
    assert "chr1:100-150(+)" in hs[0]
    man = client.get(r.headers["X-Manifest-Url"].replace("http://testserver", ""))
    lines = man.text.strip().splitlines()
    assert lines[0].split("\t")[:5] == ["locus_ref", "annotation_id", "symbol", "status", "transcript_id"]
    assert len(lines) == 6 and all(ln.split("\t")[3] == "ok" for ln in lines[1:])
    calls = fake_gffread.calls
    again = client.post("/sequences", json=TP53)                         # same request -> cached
    assert again.headers["X-Result-Id"] == r.headers["X-Result-Id"] and again.text == r.text
    assert fake_gffread.calls == calls


def test_transcript_rules(seq_client):
    client, _ = seq_client
    base = {"query": {"names": ["tp53"], "alias": "off"}, "action": "cds"}
    most = client.post("/sequences", json={**base, "transcript": "most_exons"})
    assert {h.split("|")[4].split()[0] for h in headers(most.text)} == {"tx1b"}
    allt = client.post("/sequences", json={**base, "transcript": "all"})
    assert len(headers(allt.text)) == 10 and allt.headers["X-Sequence-Count"] == "10"
    custom = client.post("/sequences", json={**base, "transcript": {
        "rank_by": ["cds_len:asc"], "require": ["coding"], "keep": 1}})
    assert {h.split("|")[4].split()[0] for h in headers(custom.text)} == {"tx1b"}


def test_loci_mode_pinned_transcript_and_fail_by_default(seq_client):
    client, _ = seq_client
    aid = client.post("/sequences/plan", json=TP53).json()["genes"][0]["annotation_id"]
    ok = client.post("/sequences", json={"action": "transcript",
                                         "loci": [{"locus_ref": f"{aid}:1", "transcript_id": "tx1b"}]})
    assert ok.status_code == 200 and headers(ok.text)[0].split("|")[4].split()[0] == "tx1b"
    bad = client.post("/sequences", json={"action": "cds", "loci": [f"{aid}:1", "nope:1"]})
    assert bad.status_code == 404
    assert bad.json()["detail"]["code"] == "annotation_not_found"
    skip = client.post("/sequences", json={"action": "cds", "on_error": "skip",
                                           "loci": [f"{aid}:1", "nope:1", f"{aid}:999"]})
    assert skip.status_code == 200 and skip.headers["X-Skipped-Count"] == "2"
    man = client.get(f"/sequences/results/{skip.headers['X-Result-Id']}/manifest.tsv").text
    reasons = [ln.split("\t")[7] for ln in man.strip().splitlines()[1:]]
    assert sorted(reasons) == ["", "annotation_not_found", "locus_not_found"]


def test_noncoding_gene_is_skipped_with_reason_and_fails_when_asked(seq_client):
    client, _ = seq_client
    q = {"query": {"names": ["uniq1"], "alias": "off"}, "action": "cds"}           # ann1 gene 2
    r = client.post("/sequences", json=q)
    assert r.status_code == 200 and r.headers["X-Sequence-Count"] == "0"
    man = client.get(f"/sequences/results/{r.headers['X-Result-Id']}/manifest.tsv").text.splitlines()
    assert man[1].split("\t")[3:8:4] == ["skipped", "no_cds"]
    strict = client.post("/sequences", json={**q, "on_error": "fail"})
    assert strict.status_code == 422 and strict.json()["detail"]["code"] == "no_cds"
    tx = client.post("/sequences", json={**q, "action": "transcript"})
    assert tx.headers["X-Sequence-Count"] == "1"                                   # lncRNA is fine


def test_transcript_beyond_gene_bounds_is_fetched_whole(seq_client):
    client, _ = seq_client
    plan = client.post("/sequences/plan", json={
        "query": {"names": ["uniq0"], "alias": "off"}, "action": "cds"}).json()
    tx = plan["genes"][0]["transcripts"][0]
    assert tx["transcript_id"] == "tx2" and tx["spliced_len"] == 51 + 31 and tx["n_exons"] == 2
    assert tx["cds_len"] == 41 + 21 and plan["genes"][0]["end"] == 250


def test_plan_lists_candidates_and_a_postable_selection(seq_client):
    client, _ = seq_client
    plan = client.post("/sequences/plan", json={**TP53, "transcript": "most_exons"}).json()
    assert plan["n_genes"] == 5 and plan["n_skipped"] == 0 and len(plan["loci"]) == 5
    g = plan["genes"][0]
    assert [t["transcript_id"] for t in g["transcripts"]] == ["tx1a", "tx1b"]
    assert [t["selected"] for t in g["transcripts"]] == [False, True]
    assert g["transcripts"][1]["n_exons"] == 2 and g["rule"] == "most_exons"
    again = client.post("/sequences", json={"action": "cds", "loci": plan["loci"]})
    assert again.status_code == 200 and again.headers["X-Sequence-Count"] == "5"
    assert {h.split("|")[4].split()[0] for h in headers(again.text)} == {"tx1b"}


def test_estimate_and_limits(seq_client):
    client, _ = seq_client
    est = client.post("/sequences/estimate", json=TP53).json()
    assert est["n_genes"] == 5 and est["mode"] == "sync" and not est["cached"]
    assert est["est_bp_upper_bound"] == 5 * 51 and est["terms"][0]["found"] is True
    client.post("/sequences", json=TP53)
    assert client.post("/sequences/estimate", json=TP53).json()["cached"] is True
    big = client.post("/sequences", json={**TP53, "max_genes": 2})
    assert big.status_code == 413 and big.json()["detail"]["code"] == "too_many_genes"
    none = client.post("/sequences", json={"query": {"names": ["zzzzz"]}, "action": "cds"})
    assert none.status_code == 404 and none.json()["detail"]["code"] == "no_genes_matched"


def test_large_request_runs_as_job_then_serves_result(seq_client):
    import settings as settings_mod
    client, _ = seq_client
    settings_mod.settings.SEQUENCES_SYNC_MAX_GENES = 2
    r = client.post("/sequences", json=TP53)
    assert r.status_code == 202
    job = r.json()
    assert job["job_id"] == job["result_id"] and r.headers["Location"].endswith(job["job_id"])
    status = {}
    for _ in range(100):
        status = client.get(f"/sequence-jobs/{job['job_id']}").json()
        if status["state"] in ("done", "failed"):
            break
        time.sleep(0.05)
    assert status["state"] == "done" and status["n_sequences"] == 5, status
    res = client.get(f"/sequences/results/{job['job_id']}")
    assert res.status_code == 200 and len(headers(res.text)) == 5
    assert client.post("/sequences", json=TP53).status_code == 200            # now cached
    assert client.get("/sequence-jobs/" + "0" * 20).status_code == 404
    assert client.get("/sequences/results/not-an-id").status_code == 404


@pytest.mark.parametrize("body,code", [
    ({"action": "cds"}, "exactly_one_source"),
    ({"action": "cds", "query": {"names": ["a"]}, "loci": ["x:1"]}, "exactly_one_source"),
    ({"query": {"names": ["a"]}}, "validation_error"),
    ({"action": "cds", "query": {"names": ["a"], "bogus": 1}}, "unknown_query_field"),
    ({"action": "cds", "query": {"names": []}}, "names_required"),
    ({"action": "cds", "query": {"names": ["a"]}, "buffer_bp": 5}, "buffer_not_supported"),
    ({"action": "cds", "query": {"names": ["a"]}, "transcript": "canonical"}, "unknown_transcript_preset"),
    ({"action": "cds", "query": {"names": ["a"]}, "transcript": {"rank_by": ["canonical"]}}, "unknown_rank_key"),
    ({"action": "cds", "loci": ["no-colon"]}, "invalid_locus_ref"),
    ({"action": "cds", "loci": ["a:1"], "on_error": "maybe"}, "validation_error"),
])
def test_invalid_requests_have_coded_errors(seq_client, body, code):
    client, _ = seq_client
    r = client.post("/sequences", json=body)
    assert r.status_code == 422, r.text
    assert r.json()["detail"]["code"] == code


def test_missing_ref_index_is_a_clear_503_for_loci_mode_only(seq_client):
    client, work = seq_client
    (work / "serve" / "current").resolve().joinpath("loci_ref.sqlite").unlink()
    r = client.post("/sequences", json={"action": "cds", "loci": ["ann0:1"]})
    assert r.status_code == 503 and r.json()["detail"]["code"] == "locus_ref_index_unavailable"
    assert client.post("/sequences", json=TP53).status_code == 200


def test_links_use_the_public_base_url_behind_a_proxy(seq_client):
    import settings as settings_mod
    client, _ = seq_client
    settings_mod.settings.PUBLIC_BASE_URL = "https://example.org/annotrieve-xrefs/api/v1/"
    try:
        r = client.post("/sequences", json=TP53)
        assert r.headers["X-Manifest-Url"] == (
            f"https://example.org/annotrieve-xrefs/api/v1/sequences/results/{r.headers['X-Result-Id']}/manifest.tsv")
        job = client.get(f"/sequence-jobs/{r.headers['X-Result-Id']}")
        assert job.headers["cache-control"] == "no-store"
        assert job.json()["status_url"].startswith("https://example.org/annotrieve-xrefs/api/v1/sequence-jobs/")
    finally:
        settings_mod.settings.PUBLIC_BASE_URL = ""
