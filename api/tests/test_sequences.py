"""POST /sequences API + pipeline tests."""
from __future__ import annotations

import sqlite3
from pathlib import Path
from unittest.mock import patch

import httpx
import pysam
import pytest
from fastapi.testclient import TestClient

from helpers.gffread_runner import SequenceAction
from helpers.ncbi_eutils import NcbiEutilsClient, close_ncbi_client, reset_ncbi_client
from helpers.pagination import MAX_SEQUENCES_POST
from services.sequences import SequencesError, _merge_extra


def test_merge_extra_strips_reserved_keys():
    merged = _merge_extra(
        {"annotation_id": "a", "status": 999, "code": "x", "message": "y"},
        {"retry_after": "1", "seqid": "chr1"},
    )
    assert merged == {"annotation_id": "a", "seqid": "chr1"}
    # Must not TypeError when reserved keys were present in helper extras.
    err = SequencesError(
        "ncbi_error",
        "boom",
        status=502,
        retry_after="2",
        **_merge_extra({"status": 503, "retry_after": "9", "accession": "X"}),
    )
    assert err.status == 502
    assert err.retry_after == "2"
    assert err.extra == {"accession": "X"}


def _write_tabix_gff(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    plain = path.with_suffix(".gff")
    # tabix requires features sorted by seqid then start.
    headers: list[str] = []
    rows: list[tuple[str, int, str]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            headers.append(line)
            continue
        cols = line.split("\t")
        rows.append((cols[0], int(cols[3]), line))
    rows.sort(key=lambda r: (r[0], r[1]))
    plain.write_text(
        "\n".join([*headers, *[r[2] for r in rows]]) + "\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="gff", force=True)
    plain.unlink(missing_ok=True)
    return path


# Sorted by start for tabix (preset=gff).
TINY_GFF = """\
##gff-version 3
NC_000017.11	test	gene	100	500	.	+	.	ID=gene-TP53;Name=TP53
NC_000017.11	test	mRNA	100	400	.	+	.	ID=txA;Parent=gene-TP53;Name=txA
NC_000017.11	test	exon	100	200	.	+	.	ID=exA1;Parent=txA
NC_000017.11	test	CDS	120	200	.	+	0	ID=cdsA1;Parent=txA
NC_000017.11	test	mRNA	150	450	.	+	.	ID=txB;Parent=gene-TP53;Name=txB
NC_000017.11	test	exon	150	450	.	+	.	ID=exB1;Parent=txB
NC_000017.11	test	CDS	160	400	.	+	0	ID=cdsB1;Parent=txB
NC_000017.11	test	exon	300	400	.	+	.	ID=exA2;Parent=txA
NC_000017.11	test	CDS	300	380	.	+	2	ID=cdsA2;Parent=txA
"""


@pytest.fixture()
def sequences_env(corpus_dir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Corpus + annotations dir with a tiny tabix GFF for ann-human / gene 10."""
    ann_root = tmp_path / "annotations"
    gff = _write_tabix_gff(
        ann_root / "9606" / "GCF_000001405.40" / "RefSeq_test.gff.gz",
        TINY_GFF,
    )

    # Point the corpus annotation at the test GFF; shrink gene coords.
    db = corpus_dir / "gene_corpus.sqlite"
    conn = sqlite3.connect(str(db))
    conn.execute(
        """
        UPDATE annotation
        SET gff_path = ?, assembly_accession = 'GCF_000001405.40'
        WHERE annotation_id = 'ann-human'
        """,
        (str(gff),),
    )
    conn.commit()
    conn.close()

    shard = corpus_dir / "per_gff" / "ann-human" / "genes.sqlite"
    sconn = sqlite3.connect(str(shard))
    sconn.execute(
        """
        UPDATE gene
        SET seqid = 'NC_000017.11', start = 100, end = 500, strand = 1,
            source_gene_id = 'gene-TP53'
        WHERE local_id = 10
        """
    )
    sconn.commit()
    sconn.close()

    monkeypatch.setenv("ANNOTATIONS_DIR", str(ann_root))
    monkeypatch.setenv("SEQUENCES_CACHE_MB", "0")
    monkeypatch.setenv("NCBI_API_KEY", "")
    monkeypatch.setenv("NCBI_WORKERS", "1")
    import settings as settings_mod
    import helpers.session as session_mod

    settings_mod.settings.ANNOTATIONS_DIR = str(ann_root)
    settings_mod.settings.ANNOTATIONS_GFF_PATH_PREFIX = ""
    settings_mod.settings.SEQUENCES_CACHE_MB = 0
    settings_mod.settings.NCBI_API_KEY = ""
    settings_mod.settings.NCBI_WORKERS = 1
    settings_mod.settings.SEQUENCES_DEADLINE_S = 60.0
    session_mod.close_reader()
    close_ncbi_client()

    yield {"ann_root": ann_root, "gff": gff, "corpus_dir": corpus_dir}
    close_ncbi_client()
    session_mod.close_reader()


def _mock_ncbi_sequence(length: int = 10_000) -> NcbiEutilsClient:
    """Return ACGT-repeating sequence covering any requested window."""

    def handler(request: httpx.Request) -> httpx.Response:
        params = dict(request.url.params)
        start = int(params["seq_start"])
        stop = int(params["seq_stop"])
        n = stop - start + 1
        # Bound for safety in tests.
        n = min(n, 50_000)
        bases = ("ACGT" * ((n // 4) + 1))[:n]
        body = f">NC_000017.11\n{bases}\n"
        return httpx.Response(200, text=body)

    return NcbiEutilsClient(
        rate_per_sec=1000,
        cache_mb=0,
        transport=httpx.MockTransport(handler),
    )


def _fake_gffread(
    *,
    gff_path: Path,
    fasta_path: Path,
    out_path: Path,
    action: SequenceAction,
    buffer_bp: int = 0,
    include_exon_coords: bool = False,
    timeout_s: float | None = None,
) -> None:
    """Write one FASTA record per mRNA/transcript id in the input GFF."""
    ids: list[str] = []
    for line in gff_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        if cols[2] in {"mRNA", "transcript"}:
            attrs = dict(
                part.split("=", 1)
                for part in cols[8].split(";")
                if "=" in part
            )
            tid = attrs.get("ID")
            if tid:
                ids.append(tid)
    with open(out_path, "w", encoding="utf-8") as fh:
        for tid in ids:
            fh.write(f">{tid}\n")
            fh.write("ACGTACGT\n")


@pytest.fixture()
def client(sequences_env) -> TestClient:
    from main import create_app

    reset_ncbi_client(_mock_ncbi_sequence())
    app = create_app()
    with TestClient(app) as test_client:
        yield test_client


def test_too_many_items(client: TestClient):
    items = [
        {
            "annotation_id": "ann-human",
            "local_id": 10,
            "transcript_id": f"tx{i}",
        }
        for i in range(MAX_SEQUENCES_POST + 1)
    ]
    resp = client.post("/sequences", json={"action": "cds", "items": items})
    assert resp.status_code == 422
    assert resp.json()["detail"]["code"] == "too_many_items"


def test_items_required(client: TestClient):
    resp = client.post("/sequences", json={"action": "cds", "items": []})
    assert resp.status_code == 422
    code = resp.json()["detail"]["code"]
    assert code in {"items_required", "validation_error"}


def test_missing_transcript_id_validation(client: TestClient):
    resp = client.post(
        "/sequences",
        json={
            "action": "cds",
            "items": [{"annotation_id": "ann-human", "local_id": 10}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["code"] == "validation_error"


def test_buffer_not_supported_for_cds(client: TestClient):
    resp = client.post(
        "/sequences",
        json={
            "action": "cds",
            "buffer_bp": 100,
            "items": [
                {
                    "annotation_id": "ann-human",
                    "local_id": 10,
                    "transcript_id": "txA",
                }
            ],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["code"] == "buffer_not_supported"


def test_buffer_out_of_range(client: TestClient):
    resp = client.post(
        "/sequences",
        json={
            "action": "transcript",
            "buffer_bp": 100_000,
            "items": [
                {
                    "annotation_id": "ann-human",
                    "local_id": 10,
                    "transcript_id": "txA",
                }
            ],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["code"] == "validation_error"


def test_fail_first_unknown_transcript_no_ncbi(client: TestClient):
    calls = {"n": 0}
    real_fetch = NcbiEutilsClient.fetch_window

    def counting_fetch(self, window, **kwargs):
        calls["n"] += 1
        return real_fetch(self, window, **kwargs)

    with patch(
        "helpers.gffread_runner.run_gffread",
        side_effect=_fake_gffread,
    ), patch.object(NcbiEutilsClient, "fetch_window", counting_fetch):
        resp = client.post(
            "/sequences",
            json={
                "action": "cds",
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "does-not-exist",
                    }
                ],
            },
        )
    assert resp.status_code == 404
    detail = resp.json()["detail"]
    assert detail["code"] == "transcript_not_found"
    assert detail["transcript_id"] == "does-not-exist"
    assert calls["n"] == 0


def test_export_cds_success(client: TestClient):
    with patch("helpers.gffread_runner.run_gffread", side_effect=_fake_gffread):
        # Patch where the service imports it.
        with patch("services.sequences.run_gffread", side_effect=_fake_gffread):
            resp = client.post(
                "/sequences",
                json={
                    "action": "cds",
                    "items": [
                        {
                            "annotation_id": "ann-human",
                            "local_id": 10,
                            "transcript_id": "txA",
                        }
                    ],
                },
            )
    assert resp.status_code == 200, resp.text
    assert resp.headers.get("content-type", "").startswith("text/x-fasta")
    assert resp.headers.get("X-Sequence-Count") == "1"
    body = resp.text
    assert body.startswith(">txA|")
    assert "NC_000017.11:100-400(+)" in body.splitlines()[0]
    assert "buffer=0" in body.splitlines()[0]
    assert "ACGTACGT" in body


def test_export_two_isoforms_one_gene(client: TestClient):
    with patch("services.sequences.run_gffread", side_effect=_fake_gffread):
        resp = client.post(
            "/sequences",
            json={
                "action": "transcript",
                "buffer_bp": 10,
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txB",
                    },
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    },
                ],
            },
        )
    assert resp.status_code == 200, resp.text
    assert resp.headers.get("X-Sequence-Count") == "2"
    # Item order preserved: txB then txA
    headers = [line for line in resp.text.splitlines() if line.startswith(">")]
    assert headers[0].startswith(">txB|")
    assert headers[1].startswith(">txA|")
    assert "buffer=10" in headers[0]


def test_dedupe_items(client: TestClient):
    with patch("services.sequences.run_gffread", side_effect=_fake_gffread):
        resp = client.post(
            "/sequences",
            json={
                "action": "cds",
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    },
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    },
                ],
            },
        )
    assert resp.status_code == 200
    assert resp.headers.get("X-Sequence-Count") == "1"


def test_gene_not_found(client: TestClient):
    with patch("services.sequences.run_gffread", side_effect=_fake_gffread):
        resp = client.post(
            "/sequences",
            json={
                "action": "cds",
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 9999,
                        "transcript_id": "txA",
                    }
                ],
            },
        )
    assert resp.status_code == 404
    assert resp.json()["detail"]["code"] == "gene_not_found"


def test_build_request_unit():
    from services.sequences import SequencesError, build_request

    req = build_request(
        action="cds",
        items=[
            {
                "annotation_id": "a",
                "local_id": 1,
                "transcript_id": "t1",
            },
            {
                "annotation_id": "a",
                "local_id": 1,
                "transcript_id": "t1",
            },
        ],
    )
    assert len(req.items) == 1
    assert req.action is SequenceAction.cds

    with pytest.raises(SequencesError) as exc:
        build_request(
            action="protein",
            buffer_bp=5,
            items=[{"annotation_id": "a", "local_id": 1, "transcript_id": "t"}],
        )
    assert exc.value.code == "buffer_not_supported"


def test_sequences_busy(client: TestClient):
    from helpers import concurrency as concurrency_mod

    # Saturate the non-blocking gate.
    acquired = []
    while concurrency_mod._SEQUENCES_SEMAPHORE.acquire(blocking=False):
        acquired.append(True)
    try:
        resp = client.post(
            "/sequences",
            json={
                "action": "cds",
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    }
                ],
            },
        )
        assert resp.status_code == 503
        assert resp.json()["detail"]["code"] == "sequences_busy"
        assert resp.headers.get("Retry-After") == "1"
    finally:
        for _ in acquired:
            concurrency_mod._SEQUENCES_SEMAPHORE.release()


def _replace_gff(sequences_env: dict, text: str) -> None:
    gff: Path = sequences_env["gff"]
    _write_tabix_gff(gff, text)


def test_transcript_incomplete(client: TestClient, sequences_env: dict):
    _replace_gff(
        sequences_env,
        """\
##gff-version 3
NC_000017.11	test	gene	100	500	.	+	.	ID=gene-TP53;Name=TP53
NC_000017.11	test	mRNA	100	400	.	+	.	ID=txA;Parent=gene-TP53
""",
    )
    with patch("services.sequences.run_gffread", side_effect=_fake_gffread):
        resp = client.post(
            "/sequences",
            json={
                "action": "transcript",
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    }
                ],
            },
        )
    assert resp.status_code == 422
    assert resp.json()["detail"]["code"] == "transcript_incomplete"


def test_transcript_strand_unknown(client: TestClient, sequences_env: dict):
    _replace_gff(
        sequences_env,
        """\
##gff-version 3
NC_000017.11	test	gene	100	500	.	.	.	ID=gene-TP53;Name=TP53
NC_000017.11	test	mRNA	100	400	.	.	.	ID=txA;Parent=gene-TP53
NC_000017.11	test	exon	100	400	.	.	.	Parent=txA
""",
    )
    with patch("services.sequences.run_gffread", side_effect=_fake_gffread):
        resp = client.post(
            "/sequences",
            json={
                "action": "transcript",
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    }
                ],
            },
        )
    assert resp.status_code == 422
    assert resp.json()["detail"]["code"] == "transcript_strand_unknown"


def test_no_cds(client: TestClient, sequences_env: dict):
    _replace_gff(
        sequences_env,
        """\
##gff-version 3
NC_000017.11	test	gene	100	500	.	+	.	ID=gene-TP53;Name=TP53
NC_000017.11	test	mRNA	100	400	.	+	.	ID=txA;Parent=gene-TP53
NC_000017.11	test	exon	100	400	.	+	.	Parent=txA
""",
    )
    with patch("services.sequences.run_gffread", side_effect=_fake_gffread):
        resp = client.post(
            "/sequences",
            json={
                "action": "cds",
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    }
                ],
            },
        )
    assert resp.status_code == 422
    assert resp.json()["detail"]["code"] == "no_cds"


def test_window_too_large(
    client: TestClient,
    sequences_env: dict,
):
    import settings as settings_mod

    previous = settings_mod.settings.SEQUENCES_MAX_WINDOW_BP
    settings_mod.settings.SEQUENCES_MAX_WINDOW_BP = 50
    try:
        with patch("services.sequences.run_gffread", side_effect=_fake_gffread):
            resp = client.post(
                "/sequences",
                json={
                    "action": "cds",
                    "items": [
                        {
                            "annotation_id": "ann-human",
                            "local_id": 10,
                            "transcript_id": "txA",
                        }
                    ],
                },
            )
    finally:
        settings_mod.settings.SEQUENCES_MAX_WINDOW_BP = previous
    assert resp.status_code == 422
    assert resp.json()["detail"]["code"] == "window_too_large"


def test_synthetic_transcript_e2e(client: TestClient, sequences_env: dict):
    _replace_gff(
        sequences_env,
        """\
##gff-version 3
NC_000017.11	test	gene	100	200	.	+	.	ID=gene-TP53;Name=TP53
NC_000017.11	test	exon	100	140	.	+	.	Parent=gene-TP53
NC_000017.11	test	CDS	110	140	.	+	0	Parent=gene-TP53
NC_000017.11	test	exon	160	200	.	+	.	Parent=gene-TP53
NC_000017.11	test	CDS	160	190	.	+	0	Parent=gene-TP53
""",
    )
    # Shrink shard gene bounds to match synthetic feature span.
    shard = sequences_env["corpus_dir"] / "per_gff" / "ann-human" / "genes.sqlite"
    sconn = sqlite3.connect(str(shard))
    sconn.execute(
        "UPDATE gene SET start = 100, end = 200 WHERE local_id = 10"
    )
    sconn.commit()
    sconn.close()

    with patch("services.sequences.run_gffread", side_effect=_fake_gffread):
        resp = client.post(
            "/sequences",
            json={
                "action": "cds",
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "gene-TP53",
                    }
                ],
            },
        )
    assert resp.status_code == 200, resp.text
    assert resp.text.startswith(">gene-TP53|")


def test_coalesce_two_genes_one_fetch(
    client: TestClient,
    sequences_env: dict,
):
    # Add a second nearby gene on the same contig/accession.
    _replace_gff(
        sequences_env,
        """\
##gff-version 3
NC_000017.11	test	gene	100	200	.	+	.	ID=gene-TP53;Name=TP53
NC_000017.11	test	mRNA	100	200	.	+	.	ID=txA;Parent=gene-TP53
NC_000017.11	test	exon	100	200	.	+	.	Parent=txA
NC_000017.11	test	CDS	110	190	.	+	0	Parent=txA
NC_000017.11	test	gene	250	350	.	+	.	ID=gene-BRCA1;Name=BRCA1
NC_000017.11	test	mRNA	250	350	.	+	.	ID=txC;Parent=gene-BRCA1
NC_000017.11	test	exon	250	350	.	+	.	Parent=txC
NC_000017.11	test	CDS	260	340	.	+	0	Parent=txC
""",
    )
    shard = sequences_env["corpus_dir"] / "per_gff" / "ann-human" / "genes.sqlite"
    sconn = sqlite3.connect(str(shard))
    sconn.execute(
        "UPDATE gene SET start = 100, end = 200 WHERE local_id = 10"
    )
    sconn.execute(
        """
        UPDATE gene
        SET seqid = 'NC_000017.11', start = 250, end = 350, strand = 1,
            source_gene_id = 'gene-BRCA1'
        WHERE local_id = 20
        """
    )
    sconn.commit()
    sconn.close()

    calls = {"n": 0}
    real_fetch = NcbiEutilsClient.fetch_window

    def counting_fetch(self, window, **kwargs):
        calls["n"] += 1
        return real_fetch(self, window, **kwargs)

    with patch("services.sequences.run_gffread", side_effect=_fake_gffread), patch.object(
        NcbiEutilsClient, "fetch_window", counting_fetch
    ):
        resp = client.post(
            "/sequences",
            json={
                "action": "cds",
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    },
                    {
                        "annotation_id": "ann-human",
                        "local_id": 20,
                        "transcript_id": "txC",
                    },
                ],
            },
        )
    assert resp.status_code == 200, resp.text
    assert calls["n"] == 1
    assert resp.headers.get("X-Sequence-Count") == "2"


def test_truncation_buffer_ok_vs_body_fail(client: TestClient):
    from helpers.ncbi_eutils import FastaSlice

    # Window for txA is 100-400 (+ buffer). Returning a short sequence that
    # truncates inside the transcript body must 502; clipping only the 3'
    # buffer flank must succeed.
    def short_body_fetch(self, window, *, allow_truncation=False):
        # Cover only up to genomic 250 — truncates inside txA (ends 400).
        n = 250 - window.start + 1
        if n < 1:
            n = 1
        return FastaSlice(
            accession=window.accession,
            start=window.start,
            stop=window.start + n - 1,
            sequence="A" * n,
            expected_length=window.length,
        )

    with patch("services.sequences.run_gffread", side_effect=_fake_gffread), patch.object(
        NcbiEutilsClient, "fetch_window", short_body_fetch
    ):
        resp = client.post(
            "/sequences",
            json={
                "action": "transcript",
                "buffer_bp": 20,
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    }
                ],
            },
        )
    assert resp.status_code == 502
    assert resp.json()["detail"]["code"] == "ncbi_sequence_truncated"

    def flank_only_fetch(self, window, *, allow_truncation=False):
        # Cover through transcript end (400) but clip the 3' buffer.
        # Requested stop = 400 + 20 = 420; return through 410.
        stop = min(window.stop, 410)
        n = stop - window.start + 1
        return FastaSlice(
            accession=window.accession,
            start=window.start,
            stop=stop,
            sequence="A" * n,
            expected_length=window.length,
        )

    with patch("services.sequences.run_gffread", side_effect=_fake_gffread), patch.object(
        NcbiEutilsClient, "fetch_window", flank_only_fetch
    ):
        resp = client.post(
            "/sequences",
            json={
                "action": "transcript",
                "buffer_bp": 20,
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    }
                ],
            },
        )
    assert resp.status_code == 200, resp.text


def test_gffread_record_count_mismatch(client: TestClient):
    def bad_gffread(**kwargs):
        out = kwargs["out_path"]
        out.write_text(">only-one\nACGT\n", encoding="utf-8")

    with patch("services.sequences.run_gffread", side_effect=bad_gffread):
        resp = client.post(
            "/sequences",
            json={
                "action": "cds",
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    },
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txB",
                    },
                ],
            },
        )
    assert resp.status_code == 500
    detail = resp.json()["detail"]
    assert detail["code"] == "gffread_failed"
    assert "transcript_ids" in detail


def test_sequences_timeout(client: TestClient):
    import settings as settings_mod

    previous = settings_mod.settings.SEQUENCES_DEADLINE_S
    settings_mod.settings.SEQUENCES_DEADLINE_S = 0.0
    try:
        with patch("services.sequences.run_gffread", side_effect=_fake_gffread):
            resp = client.post(
                "/sequences",
                json={
                    "action": "cds",
                    "items": [
                        {
                            "annotation_id": "ann-human",
                            "local_id": 10,
                            "transcript_id": "txA",
                        }
                    ],
                },
            )
    finally:
        settings_mod.settings.SEQUENCES_DEADLINE_S = previous
    assert resp.status_code == 504
    assert resp.json()["detail"]["code"] == "sequences_timeout"


def test_rewrite_headers_exact_only(tmp_path: Path):
    from helpers.gffread_runner import GffreadError, rewrite_fasta_headers

    raw = tmp_path / "in.fa"
    out = tmp_path / "out.fa"
    raw.write_text(">txAextra\nACGT\n", encoding="utf-8")
    with pytest.raises(GffreadError) as exc:
        rewrite_fasta_headers(raw, out, header_for_id={"txA": "txA|gene|sp|ann"})
    assert exc.value.code == "gffread_failed"


def test_gid_unique_across_groups(
    client: TestClient,
    sequences_env: dict,
):
    _replace_gff(
        sequences_env,
        """\
##gff-version 3
NC_000017.11	test	gene	100	200	.	+	.	ID=gene-TP53;Name=TP53
NC_000017.11	test	mRNA	100	200	.	+	.	ID=txA;Parent=gene-TP53
NC_000017.11	test	exon	100	200	.	+	.	Parent=txA
NC_000017.11	test	CDS	110	190	.	+	0	Parent=txA
NC_000017.11	test	gene	250	350	.	+	.	ID=gene-BRCA1;Name=BRCA1
NC_000017.11	test	mRNA	250	350	.	+	.	ID=txC;Parent=gene-BRCA1
NC_000017.11	test	exon	250	350	.	+	.	Parent=txC
NC_000017.11	test	CDS	260	340	.	+	0	Parent=txC
""",
    )
    shard = sequences_env["corpus_dir"] / "per_gff" / "ann-human" / "genes.sqlite"
    sconn = sqlite3.connect(str(shard))
    sconn.execute(
        "UPDATE gene SET start = 100, end = 200 WHERE local_id = 10"
    )
    sconn.execute(
        """
        UPDATE gene
        SET seqid = 'NC_000017.11', start = 250, end = 350, strand = 1,
            source_gene_id = 'gene-BRCA1'
        WHERE local_id = 20
        """
    )
    sconn.commit()
    sconn.close()

    seen_gids: list[str] = []
    real = None

    def capturing_gffread(**kwargs):
        gff_path: Path = kwargs["gff_path"]
        seen_gids.append(gff_path.stem)
        return _fake_gffread(**kwargs)

    with patch("services.sequences.run_gffread", side_effect=capturing_gffread):
        resp = client.post(
            "/sequences",
            json={
                "action": "cds",
                "items": [
                    {
                        "annotation_id": "ann-human",
                        "local_id": 10,
                        "transcript_id": "txA",
                    },
                    {
                        "annotation_id": "ann-human",
                        "local_id": 20,
                        "transcript_id": "txC",
                    },
                ],
            },
        )
    assert resp.status_code == 200, resp.text
    assert seen_gids == ["g0", "g1"]


@pytest.mark.skipif(
    __import__("shutil").which("gffread") is None,
    reason="gffread binary not installed",
)
def test_real_gffread_e2e(client: TestClient):
    resp = client.post(
        "/sequences",
        json={
            "action": "cds",
            "items": [
                {
                    "annotation_id": "ann-human",
                    "local_id": 10,
                    "transcript_id": "txA",
                }
            ],
        },
    )
    assert resp.status_code == 200, resp.text
    assert resp.text.startswith(">txA|")
    assert len(resp.text.splitlines()) >= 2
