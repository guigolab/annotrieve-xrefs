"""/downloads (queue, quota, gzip FASTA), /transcripts (new index only)."""
from __future__ import annotations

import gzip
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
                                        "SEQUENCES_MAX_SEQUENCES", "SEQUENCES_BATCH_GENES",
                                        "NCBI_WORKERS", "SEQUENCES_DEADLINE_S", "DOWNLOADS_PER_DAY",
                                        "DOWNLOADS_POLL_S", "DOWNLOADS_LOCK_FILE",
                                        "SEQUENCES_MAX_GENES")}
    s.PRESENCE_DIR = str(work / "serve")
    s.ANNOTATIONS_DIR = str(ann_root)
    s.SEQUENCES_RESULTS_DIR = str(tmp_path / "results")
    s.SEQUENCES_MAX_SEQUENCES = 200
    s.SEQUENCES_BATCH_GENES = 3
    s.NCBI_WORKERS = 1
    s.SEQUENCES_DEADLINE_S = 60.0
    s.SEQUENCES_MAX_GENES = 200
    s.DOWNLOADS_PER_DAY = 0                      # the quota has its own tests
    s.DOWNLOADS_POLL_S = 0.02
    s.DOWNLOADS_LOCK_FILE = str(tmp_path / "downloads.lock")
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


class Reply:
    """A finished download seen as one reply: the FASTA text and its headers."""

    def __init__(self, resp, text: str) -> None:
        self.status_code, self.headers, self.text = 200, resp.headers, text
        self.id = resp.headers["X-Download-Id"]


def dl(client, body, headers=None, wait: float = 30.0):
    """POST /downloads, wait for the queue, fetch the FASTA. A refusal or failure is returned as is."""
    r = client.post("/downloads", json=body, headers=headers or {})
    if r.status_code not in (200, 202):
        return r
    st = r.json()
    end = time.time() + wait
    while st["state"] in ("queued", "running") and time.time() < end:
        time.sleep(0.02)
        st = client.get(f"/downloads/{st['id']}").json()
    f = client.get(f"/downloads/{st['id']}/fasta")
    if f.status_code != 200:
        return f
    assert f.headers["content-type"] == "application/gzip"
    return Reply(f, gzip.decompress(f.content).decode())


def report(client, r) -> str:
    return client.get(f"/downloads/{r.id}/report").text


def estimate(client, body):
    return client.post("/downloads?dry_run=true", json=body)


TP53 = {"query": {"names": ["tp53"], "alias": "off"}, "action": "cds"}


def test_query_sync_fasta_headers_manifest_and_cache(seq_client):
    client, _ = seq_client
    r = dl(client, TP53)
    assert r.status_code == 200 and r.headers["content-type"] == "application/gzip"
    manifest = report(client, r)
    assert r.headers["X-Sequence-Count"] == "5" and r.headers["X-Skipped-Count"] == "0", manifest
    hs = headers(r.text)
    assert len(hs) == 5 and all(h.startswith("tp53|Homo_sapiens|GCF_000001405.40|ann") for h in hs)
    assert all(h.split("|")[4].split()[0] == "tx1a" for h in hs)          # longest CDS (41 > 32)
    assert "chr1:100-150(+)" in hs[0]
    lines = manifest.strip().splitlines()
    assert lines[0].split("\t")[:5] == ["locus_ref", "annotation_id", "symbol", "status", "transcript_id"]
    assert len(lines) == 6 and all(ln.split("\t")[3] == "ok" for ln in lines[1:])
    calls = fake_gffread.calls
    again = dl(client, TP53)                         # same request -> cached
    assert again.headers["X-Download-Id"] == r.headers["X-Download-Id"] and again.text == r.text
    assert fake_gffread.calls == calls


def test_transcript_rules(seq_client):
    client, _ = seq_client
    base = {"query": {"names": ["tp53"], "alias": "off"}, "action": "cds"}
    most = dl(client, {**base, "transcript": "most_exons"})
    assert {h.split("|")[4].split()[0] for h in headers(most.text)} == {"tx1b"}
    allt = dl(client, {**base, "transcript": "all"})
    assert len(headers(allt.text)) == 10 and allt.headers["X-Sequence-Count"] == "10"
    custom = dl(client, {**base, "transcript": {
        "rank_by": ["cds_len:asc"], "require": ["coding"], "keep": 1}})
    assert {h.split("|")[4].split()[0] for h in headers(custom.text)} == {"tx1b"}


def test_loci_mode_pinned_transcript_and_fail_by_default(seq_client):
    client, _ = seq_client
    aid = client.post("/transcripts", json=TP53).json()["genes"][0]["annotation_id"]
    ok = dl(client, {"action": "transcript",
                                         "loci": [{"locus_ref": f"{aid}:1", "transcript_id": "tx1b"}]})
    assert ok.status_code == 200 and headers(ok.text)[0].split("|")[4].split()[0] == "tx1b"
    bad = dl(client, {"action": "cds", "loci": [f"{aid}:1", "nope:1"]})
    assert bad.status_code == 404
    assert bad.json()["detail"]["code"] == "annotation_not_found"
    skip = dl(client, {"action": "cds", "on_error": "skip",
                                           "loci": [f"{aid}:1", "nope:1", f"{aid}:999"]})
    assert skip.status_code == 200 and skip.headers["X-Skipped-Count"] == "2"
    man = report(client, skip)
    reasons = [ln.split("\t")[7] for ln in man.strip().splitlines()[1:]]
    assert sorted(reasons) == ["", "annotation_not_found", "locus_not_found"]


def test_noncoding_gene_is_skipped_with_reason_and_fails_when_asked(seq_client):
    client, _ = seq_client
    q = {"query": {"names": ["uniq1"], "alias": "off"}, "action": "cds"}           # ann1 gene 2
    r = dl(client, q)
    assert r.status_code == 200 and r.headers["X-Sequence-Count"] == "0"
    man = report(client, r).splitlines()
    assert man[1].split("\t")[3:8:4] == ["skipped", "no_cds"]
    strict = dl(client, {**q, "on_error": "fail"})
    assert strict.status_code == 422 and strict.json()["detail"]["code"] == "no_cds"
    tx = dl(client, {**q, "action": "transcript"})
    assert tx.headers["X-Sequence-Count"] == "1"                                   # lncRNA is fine


def test_transcript_beyond_gene_bounds_is_fetched_whole(seq_client):
    client, _ = seq_client
    plan = client.post("/transcripts", json={
        "query": {"names": ["uniq0"], "alias": "off"}, "action": "cds"}).json()
    tx = plan["genes"][0]["transcripts"][0]
    assert tx["transcript_id"] == "tx2" and tx["spliced_len"] == 51 + 31 and tx["n_exons"] == 2
    assert tx["cds_len"] == 41 + 21 and plan["genes"][0]["end"] == 250


def test_plan_lists_candidates_and_a_postable_selection(seq_client):
    client, _ = seq_client
    plan = client.post("/transcripts", json={**TP53, "transcript": "most_exons"}).json()
    assert plan["n_genes"] == 5 and plan["n_skipped"] == 0 and len(plan["loci"]) == 5
    g = plan["genes"][0]
    assert [t["transcript_id"] for t in g["transcripts"]] == ["tx1a", "tx1b"]
    assert [t["selected"] for t in g["transcripts"]] == [False, True]
    assert g["transcripts"][1]["n_exons"] == 2 and g["rule"] == "most_exons"
    again = dl(client, {"action": "cds", "loci": plan["loci"]})
    assert again.status_code == 200 and again.headers["X-Sequence-Count"] == "5"
    assert {h.split("|")[4].split()[0] for h in headers(again.text)} == {"tx1b"}


def test_estimate_and_limits(seq_client):
    client, _ = seq_client
    est = estimate(client, TP53).json()
    assert est["n_genes"] == 5 and est["max_sequences"] == 200 and not est["cached"]
    assert est["est_bp_upper_bound"] == 5 * 51 and est["terms"][0]["found"] is True
    assert client.get("/downloads/quota").json()["used"] == 0          # a dry run uses no quota
    dl(client, TP53)
    assert estimate(client, TP53).json()["cached"] is True
    big = dl(client, {**TP53, "max_genes": 2})
    assert big.status_code == 413 and big.json()["detail"]["code"] == "too_many_genes"
    none = dl(client, {"query": {"names": ["zzzzz"]}, "action": "cds"})
    assert none.status_code == 404 and none.json()["detail"]["code"] == "no_genes_matched"


def test_a_download_is_queued_polled_and_fetched(seq_client):
    client, _ = seq_client
    from services.seq_runner import stop_runner
    stop_runner()                                          # nothing runs until we say so
    r = client.post("/downloads", json=TP53)
    assert r.status_code == 202
    st = r.json()
    assert r.headers["Location"].endswith(f"/downloads/{st['id']}") and r.headers["Retry-After"] == "5"
    assert st["state"] == "queued" and st["queue_position"] == 1 and st["output"] == "fasta"
    assert st["progress"] == {"genes_done": 0, "genes_total": 5} and "fasta" not in st["links"]
    assert client.get(f"/downloads/{st['id']}/fasta").status_code == 409
    again = client.post("/downloads", json=TP53)                    # equal request: same download
    assert again.status_code == 202 and again.json()["id"] == st["id"]
    from services.seq_runner import run_pending_once
    assert run_pending_once() == 1
    done = client.get(f"/downloads/{st['id']}").json()
    assert done["state"] == "done" and done["n_sequences"] == 5 and done["queue_position"] is None
    assert set(done["links"]) == {"self", "fasta", "report"} and done["size_bytes"] > 0
    assert "Retry-After" not in client.get(f"/downloads/{st['id']}").headers
    f = client.get(done["links"]["fasta"].replace("http://testserver", ""))
    assert f.status_code == 200 and f.headers["x-sequence-count"] == "5"
    assert f.content[:2] == b"\x1f\x8b" and len(headers(gzip.decompress(f.content).decode())) == 5
    assert client.get("/downloads/" + "0" * 20).status_code == 404
    assert client.get("/downloads/not-an-id").status_code == 404
    stored = client.post("/downloads", json=TP53)                   # finished: 200, no new work
    assert stored.status_code == 200 and stored.json()["state"] == "done"


def test_the_file_is_removed_after_the_first_complete_download(seq_client):
    import settings as settings_mod
    from services.seq_results import purge_expired
    client, _ = seq_client
    r = dl(client, TP53)
    purge_expired()
    assert client.get(f"/downloads/{r.id}").status_code == 200       # still inside the grace time
    assert client.get(f"/downloads/{r.id}").json()["expires_at"] > time.time()
    settings_mod.settings.DOWNLOAD_KEEP_AFTER_FETCH_S = 0
    try:
        time.sleep(0.01)
        purge_expired()
    finally:
        settings_mod.settings.DOWNLOAD_KEEP_AFTER_FETCH_S = 600
    assert client.get(f"/downloads/{r.id}").status_code == 404
    assert client.get(f"/downloads/{r.id}/fasta").status_code == 404


def test_an_unfetched_download_waits_longer_and_a_half_read_one_is_not_counted(seq_client):
    import settings as settings_mod
    from services.seq_results import purge_expired, read_status
    client, _ = seq_client
    r = client.post("/downloads", json=TP53).json()
    for _ in range(300):
        if client.get(f"/downloads/{r['id']}").json()["state"] == "done":
            break
        time.sleep(0.02)
    settings_mod.settings.DOWNLOAD_KEEP_AFTER_FETCH_S = 0
    try:
        purge_expired()                                  # never fetched: kept (2 h rule)
        assert client.get(f"/downloads/{r['id']}").status_code == 200
        client.get(f"/downloads/{r['id']}/fasta", headers={"Range": "bytes=0-5"})
        resp = client.get(f"/downloads/{r['id']}/fasta")
        assert resp.status_code == 200
    finally:
        settings_mod.settings.DOWNLOAD_KEEP_AFTER_FETCH_S = 600
    assert read_status(r["id"])["downloaded_at"]


def test_the_queue_runs_one_download_at_a_time_in_order(seq_client, limits):
    from services.seq_results import read_status
    from services.seq_runner import run_pending_once, stop_runner
    client, _ = seq_client
    stop_runner()
    limits.SEQUENCES_IP_MAX_ACTIVE = 0
    ids = [client.post("/downloads", json=_variant(i)).json()["id"] for i in range(3)]
    pos = [client.get(f"/downloads/{i}").json()["queue_position"] for i in ids]
    assert pos == [1, 2, 3]
    assert run_pending_once() == 3
    started = [read_status(i)["started"] for i in ids]
    assert started == sorted(started) and all(read_status(i)["state"] == "done" for i in ids)


def test_only_one_runner_works_at_a_time(seq_client):
    import fcntl
    import settings as settings_mod
    from services.seq_runner import run_pending_once, stop_runner
    client, _ = seq_client
    stop_runner()
    st = client.post("/downloads", json=TP53).json()
    with open(settings_mod.settings.DOWNLOADS_LOCK_FILE, "a+") as other:
        fcntl.flock(other, fcntl.LOCK_EX)                    # another worker is the runner
        assert run_pending_once() == 0
        assert client.get(f"/downloads/{st['id']}").json()["state"] == "queued"
    assert run_pending_once() == 1


def test_a_download_can_be_cancelled_only_by_its_client(seq_client):
    from services.seq_runner import stop_runner
    client, _ = seq_client
    stop_runner()
    st = client.post("/downloads", json=TP53).json()
    other = client.delete(f"/downloads/{st['id']}", headers={"X-Forwarded-For": "203.0.113.5"})
    assert other.status_code == 403 and other.json()["detail"]["code"] == "not_your_download"
    assert client.delete(f"/downloads/{st['id']}").status_code == 204
    assert client.get(f"/downloads/{st['id']}").status_code == 404


def test_more_sequences_than_the_limit_are_refused_before_any_fetch(seq_client, monkeypatch):
    import settings as settings_mod
    from services import seq_runner
    client, _ = seq_client
    settings_mod.settings.SEQUENCES_MAX_SEQUENCES = 3
    # 5 genes, one sequence each: the gene count already exceeds the limit
    r = dl(client, TP53)
    assert r.status_code == 413 and r.json()["detail"]["code"] == "too_many_genes"
    # several transcripts per gene: counted from the GFF slices, still before queueing
    settings_mod.settings.SEQUENCES_MAX_SEQUENCES = 200
    many = {**TP53, "transcript": "all"}
    n = dl(client, many)
    assert n.status_code == 200
    got = int(n.headers["X-Sequence-Count"])
    assert got == 10                                   # 5 genes, two transcripts each in the fixture
    settings_mod.settings.SEQUENCES_MAX_SEQUENCES = 9
    monkeypatch.setattr(seq_runner, "run_export",
                        lambda *a, **k: (_ for _ in ()).throw(AssertionError("fetched")))
    r2 = client.post("/downloads", json={**many, "max_genes": 5})
    assert r2.status_code == 413 and r2.json()["detail"]["code"] == "too_many_sequences"
    assert r2.json()["detail"]["n_sequences_at_least"] == 10


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
    r = dl(client, body)
    assert r.status_code == 422, r.text
    assert r.json()["detail"]["code"] == code


def test_missing_ref_index_is_a_clear_503_for_loci_mode_only(seq_client):
    client, work = seq_client
    (work / "serve" / "current").resolve().joinpath("loci_ref.sqlite").unlink()
    r = dl(client, {"action": "cds", "loci": ["ann0:1"]})
    assert r.status_code == 503 and r.json()["detail"]["code"] == "locus_ref_index_unavailable"
    assert dl(client, TP53).status_code == 200


def test_links_use_the_public_base_url_behind_a_proxy(seq_client):
    import settings as settings_mod
    client, _ = seq_client
    settings_mod.settings.PUBLIC_BASE_URL = "https://example.org/annotrieve-xrefs/api/v1/"
    try:
        r = client.post("/downloads", json=TP53)
        base = "https://example.org/annotrieve-xrefs/api/v1/downloads/"
        assert r.headers["Location"] == base + r.json()["id"]
        assert r.json()["links"]["self"] == base + r.json()["id"]
    finally:
        settings_mod.settings.PUBLIC_BASE_URL = ""


def test_tabix_cache_keeps_at_most_two_files_open_and_releases_them(tmp_path):
    from helpers.annotation_files import AnnotationFiles
    from helpers.tabix_fetch import FetchError, TabixCache

    cache = TabixCache(max_handles=2)
    files = []
    for i in range(5):
        gff = tmp_path / f"a{i}" / "x.gff.gz"
        write_tabix(gff, GENE1)
        files.append(AnnotationFiles(gff, Path(str(gff) + ".tbi"), 9606, "GCF_1", gff.parent, ""))
    for f in files:
        assert any("gene-1" in ln for ln in cache.fetch(f, seqid="chr1", start=100, end=150))
        assert len(cache) <= 2
    with pytest.raises(FetchError) as exc:
        cache.fetch(files[0], seqid="no-such-contig", start=1, end=10)
    assert exc.value.code == "seqid_not_resolved"
    cache.close()
    assert len(cache) == 0


def test_oversized_gff_index_is_skipped_with_a_reason(seq_client):
    import settings as settings_mod
    client, _ = seq_client
    settings_mod.settings.SEQUENCES_MAX_INDEX_MB = 1e-9
    try:
        r = dl(client, TP53)
        assert r.status_code == 200 and r.headers["X-Sequence-Count"] == "0"
        assert r.headers["X-Skipped-Count"] == "5"
        man = report(client, r)
        assert {ln.split("\t")[7] for ln in man.strip().splitlines()[1:]} == {"index_too_large"}
    finally:
        settings_mod.settings.SEQUENCES_MAX_INDEX_MB = 0.0


def test_contig_unknown_to_ncbi_is_skipped_with_its_own_reason(seq_client):
    client, _ = seq_client

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, text="Error: F a i l e d  t o  u n d e r s t a n d  i d : chr1")

    reset_ncbi_client(NcbiEutilsClient(rate_per_sec=1000, cache_mb=0,
                                       transport=httpx.MockTransport(handler)))
    r = dl(client, TP53)
    assert r.status_code == 200 and r.headers["X-Sequence-Count"] == "0"
    man = report(client, r)
    assert {ln.split("\t")[7] for ln in man.strip().splitlines()[1:]} == {"ncbi_sequence_not_found"}


def test_request_order_is_kept_when_genes_are_resolved_by_annotation(seq_client):
    client, _ = seq_client
    plan = client.post("/transcripts", json=TP53).json()
    refs = [g["locus_ref"] for g in plan["genes"]][::-1]           # reverse: annotation order differs
    r = dl(client, {"action": "cds", "loci": refs})
    got = [h.split("|")[3] for h in headers(r.text)]
    assert got == [ref.split(":")[0] for ref in refs]


FUZZ_BODIES = [
    {"action": "cds", "query": {"names": ["tp53"], "include_taxids": ["abc"]}},
    {"action": "cds", "query": {"names": ["tp53"], "include_taxids": 9606}},
    {"action": "cds", "query": {"names": ["tp53"], "exclude_taxids": [None]}},
    {"action": "cds", "query": {"names": ["tp53"], "annotation_ids": 5}},
    {"action": "cds", "query": {"names": [1, None, {}]}},
    {"action": "cds", "query": {"names": ["tp53"], "alias": "bogus"}},
    {"action": "cds", "query": {"names": ["tp53"], "alias": 7}},
    {"action": "cds", "query": {"names": ["tp53"], "case_sensitive": "false"}},
    {"action": "cds", "query": {"names": "tp53,"}},
    {"action": "cds", "query": {"names": []}},
    {"action": "cds", "query": "tp53"},
    {"action": "cds", "query": {"names": ["tp53"]}, "max_genes": 0},
    {"action": "cds", "query": {"names": ["tp53"]}, "max_genes": -3},
    {"action": "cds", "query": {"names": ["tp53"]}, "buffer_bp": -1},
    {"action": "transcript", "query": {"names": ["tp53"]}, "buffer_bp": 10**9},
    {"action": "cds", "query": {"names": ["tp53"]}, "transcript": 5},
    {"action": "cds", "query": {"names": ["tp53"]}, "transcript": {"keep": 0}},
    {"action": "cds", "query": {"names": ["tp53"]}, "transcript": {"keep": "x"}},
    {"action": "cds", "query": {"names": ["tp53"]}, "transcript": {"rank_by": "cds_len:desc"}},
    {"action": "cds", "query": {"names": ["tp53"]}, "transcript": {"rank_by": [3]}},
    {"action": "cds", "query": {"names": ["tp53"]}, "transcript": {"require": "coding"}},
    {"action": "cds", "query": {"names": ["tp53"]}, "transcript": {"nonsense": 1}},
    {"action": "cds", "loci": [5]},
    {"action": "cds", "loci": [[]]},
    {"action": "cds", "loci": [{"locus_ref": 7}]},
    {"action": "cds", "loci": [{"locus_ref": "a:1", "transcript_id": {"x": 1}}]},
    {"action": "cds", "loci": "ann0:1"},
    {"action": "cds", "items": [{"annotation_id": "a", "local_id": "x"}]},
    {"action": "cds", "items": [{"annotation_id": "a"}]},
    {"action": "cds", "loci": ["a:1"], "typo_field": True},
    {"action": ["cds"], "loci": ["a:1"]},
    {"action": "cds", "loci": ["a:99999999999999999999"]},
    {"action": "cds", "loci": ["ann1:99999999999999999999"]},
    {"action": "cds", "loci": ["ann1:\u00b2"]},
    {"action": "cds", "items": [{"annotation_id": "ann1", "local_id": 10**30}]},
    {"action": "cds", "query": {"names": ["tp53"], "include_taxids": ["99999999999999999999"]}},
    {"action": "cds", "query": {"names": ["tp53"], "include_taxids": ["\u00b2"]}},
    {"action": "cds", "query": {"names": ["tp53"]}, "transcript": {"rank_by": 5}},
    {"action": "cds", "query": {"names": ["tp53"]}, "transcript": {"require": 5}},
    [],
    "text",
]


@pytest.mark.parametrize("body", FUZZ_BODIES)
@pytest.mark.parametrize("path", ["/downloads", "/downloads?dry_run=true", "/transcripts"])
def test_malformed_bodies_never_cause_a_server_error(seq_client, path, body):
    client, _ = seq_client
    r = client.post(path, json=body)
    assert r.status_code < 500, (path, body, r.status_code, r.text[:300])
    if r.status_code >= 400:
        assert "code" in r.json()["detail"], r.text[:300]


def test_cors_exposes_the_headers_a_browser_client_needs(seq_client):
    client, _ = seq_client
    r = client.options("/downloads", headers={"Origin": "https://x.org",
                                              "Access-Control-Request-Method": "POST"})
    assert r.status_code == 200
    r = client.post("/downloads", json=TP53, headers={"Origin": "https://x.org"})
    exposed = {h.strip().lower() for h in r.headers["access-control-expose-headers"].split(",")}
    assert {"x-download-id", "x-skipped-count", "x-gene-count", "x-data-version",
            "x-sequence-count", "location", "retry-after", "content-disposition"} <= exposed
    d = client.options("/downloads/" + "0" * 20, headers={"Origin": "https://x.org",
                                                         "Access-Control-Request-Method": "DELETE"})
    assert d.status_code == 200


def test_a_crashing_export_ends_failed_and_can_be_asked_again(seq_client, monkeypatch):
    from services import seq_runner
    client, _ = seq_client

    def boom(*a, **k):
        raise RuntimeError("unexpected")

    monkeypatch.setattr(seq_runner, "run_export", boom)
    r = dl(client, TP53)
    assert r.status_code == 500 and r.json()["detail"]["code"] == "download_failed"
    monkeypatch.undo()
    again = dl(client, TP53)
    assert again.status_code == 200 and again.headers["X-Sequence-Count"] == "5"


def test_a_file_removed_while_listed_is_404_not_500(seq_client):
    import shutil
    client, _ = seq_client
    r = dl(client, TP53)
    from services.seq_results import result_dir
    (result_dir(r.id) / "result.fa.gz").unlink()
    assert client.get(f"/downloads/{r.id}/fasta").status_code == 404
    assert client.get(f"/downloads/{r.id}/report").status_code == 404
    shutil.rmtree(result_dir(r.id), ignore_errors=True)


# ---- quota, queue and storage limits -----------------------------------------------------------

@pytest.fixture()
def limits():
    import settings as settings_mod
    st = settings_mod.settings
    names = ("SEQUENCES_MAX_RESULTS", "SEQUENCES_IP_MAX_ACTIVE", "DOWNLOADS_PER_DAY",
             "DOWNLOADS_MAX_QUEUED", "TRUSTED_PROXY_HOPS")
    saved = {n: getattr(st, n) for n in names}
    yield st
    for n, v in saved.items():
        setattr(st, n, v)


PRESETS = ["longest", "longest_cds", "most_exons", "shortest_cds", "fewest_exons"]


def _variant(i):
    return {**TP53, "transcript": PRESETS[i]}


def test_downloads_per_day_are_capped_per_ip_but_equal_requests_are_free(seq_client, limits):
    client, _ = seq_client
    limits.DOWNLOADS_PER_DAY = 2
    assert dl(client, _variant(0)).status_code == 200
    assert client.get("/downloads/quota").json() | {"resets_at": None} == {
        "limit": 2, "used": 1, "remaining": 1, "window_hours": 24, "resets_at": None,
        "active": 0, "max_active": 2}
    assert dl(client, _variant(1)).status_code == 200
    assert dl(client, _variant(0)).status_code == 200                        # equal: free
    r = client.post("/downloads", json=_variant(2))
    assert r.status_code == 429 and r.json()["detail"]["code"] == "daily_quota_exceeded"
    assert 0 < int(r.headers["Retry-After"]) <= 24 * 3600
    q = client.get("/downloads/quota").json()
    assert q["remaining"] == 0 and q["resets_at"] > time.time()
    other = dl(client, _variant(2), headers={"X-Forwarded-For": "203.0.113.9"})
    assert other.status_code == 200                                          # another client
    assert client.get("/downloads/quota", headers={"X-Forwarded-For": "203.0.113.9"}).json()["used"] == 1


def test_the_quota_window_is_rolling_24_hours(seq_client, limits):
    import json
    from helpers.client_ip import owner_key
    from services.seq_results import _quota_file
    client, _ = seq_client
    limits.DOWNLOADS_PER_DAY = 1
    f = _quota_file(owner_key("testclient"))
    f.write_text(json.dumps([time.time() - 24 * 3600 - 5]))        # used just over a day ago
    assert client.get("/downloads/quota").json()["used"] == 0
    assert dl(client, _variant(0)).status_code == 200
    assert client.post("/downloads", json=_variant(1)).status_code == 429


def test_downloads_waiting_or_running_per_ip_are_capped(seq_client, limits):
    from services.seq_runner import stop_runner
    client, _ = seq_client
    stop_runner()
    limits.SEQUENCES_IP_MAX_ACTIVE = 1
    assert client.post("/downloads", json=_variant(0)).status_code == 202
    r = client.post("/downloads", json=_variant(1))
    assert r.status_code == 429 and r.json()["detail"]["code"] == "too_many_active_downloads"
    ok = client.post("/downloads", json=_variant(1), headers={"X-Forwarded-For": "198.51.100.7"})
    assert ok.status_code == 202


def test_the_client_ip_cannot_be_chosen_with_x_forwarded_for(seq_client, limits):
    client, _ = seq_client
    limits.DOWNLOADS_PER_DAY = 1
    h = {"X-Forwarded-For": "1.1.1.1, 192.0.2.50"}            # proxy appended 192.0.2.50
    assert dl(client, _variant(0), headers=h).status_code == 200
    spoof = {"X-Forwarded-For": "9.9.9.9, 192.0.2.50"}         # same real client, new fake prefix
    assert client.post("/downloads", json=_variant(1), headers=spoof).status_code == 429


def test_a_full_queue_refuses_new_downloads(seq_client, limits):
    from services.seq_runner import stop_runner
    client, _ = seq_client
    stop_runner()
    limits.SEQUENCES_IP_MAX_ACTIVE = 0
    limits.DOWNLOADS_MAX_QUEUED = 2
    for i in range(2):
        assert client.post("/downloads", json=_variant(i)).status_code == 202
    r = client.post("/downloads", json=_variant(2))
    assert r.status_code == 503 and r.json()["detail"]["code"] == "downloads_queue_full"
    assert "Retry-After" in r.headers


def test_a_full_store_refuses_and_never_removes_finished_downloads(seq_client, limits):
    from services.seq_results import results_root
    client, _ = seq_client
    limits.SEQUENCES_IP_MAX_ACTIVE = 0
    limits.SEQUENCES_MAX_RESULTS = 2
    first = [dl(client, _variant(i)) for i in range(2)]
    assert all(r.status_code == 200 for r in first)
    r = client.post("/downloads", json=_variant(2))
    assert r.status_code == 503 and r.json()["detail"]["code"] == "downloads_capacity"
    left = {d.name for d in results_root().iterdir() if d.is_dir() and not d.name.startswith(".")}
    assert left == {x.id for x in first}


def test_prefix_terms_in_sequences_query(seq_client):
    client, _ = seq_client
    body = {"query": {"terms": [{"q": "uniq*", "exclude": ["uniq4"]}], "alias": "off"},
            "action": "cds"}
    est = estimate(client, body).json()
    assert est["n_genes"] == 4
    t = est["terms"][0]
    assert t["kind"] == "prefix" and t["n_names"] == 4 and "uniq4" not in t["names"]
    # option order must not change the cached result
    same = {"query": {"alias": "off", "terms": [{"exclude": ["UNIQ4"], "q": "uniq*"}]},
            "action": "cds"}
    assert estimate(client, same).json()["id"] == est["id"]
    bad = dl(client, {"query": {"terms": [{"q": "tp53", "split": True}]}, "action": "cds"})
    assert bad.status_code == 422 and bad.json()["detail"]["code"] == "option_needs_prefix"
    short = dl(client, {"query": {"names": ["un*"]}, "action": "cds"})
    assert short.status_code == 422 and short.json()["detail"]["code"] == "prefix_too_short"


def test_transcripts_of_one_gene_by_locus_ref(seq_client):
    client, _ = seq_client
    aid = client.post("/transcripts", json=TP53).json()["genes"][0]["annotation_id"]
    g = client.get(f"/genes/{aid}:1/transcripts").json()
    assert [t["transcript_id"] for t in g["transcripts"]] == ["tx1a", "tx1b"]
    assert [t["selected"] for t in g["transcripts"]] == [True, False]          # longest CDS
    most = client.get(f"/genes/{aid}:1/transcripts", params={"transcript": "most_exons"}).json()
    assert [t["selected"] for t in most["transcripts"]] == [False, True]
    assert client.get(f"/genes/{aid}:999/transcripts").status_code == 404
    assert client.get("/genes/nope:1/transcripts").json()["detail"]["code"] == "annotation_not_found"
    bad = client.get(f"/genes/{aid}:1/transcripts", params={"action": "dna"})
    assert bad.status_code == 422
    assert client.get("/genes/no-colon/transcripts").json()["detail"]["code"] == "invalid_locus_ref"
