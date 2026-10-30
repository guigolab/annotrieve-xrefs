"""GFF transcript extraction / synthetic / shift tests."""
from __future__ import annotations

import pytest

from helpers.gff_slice import (
    load_features_from_lines,
    resolve_transcripts,
    shift_gff_lines,
)


GENE_GFF = """\
chr1	src	gene	100	500	.	+	.	ID=gene1;Name=GENE1
chr1	src	mRNA	100	400	.	+	.	ID=tx1;Parent=gene1
chr1	src	exon	100	200	.	+	.	ID=ex1;Parent=tx1
chr1	src	CDS	120	200	.	+	0	ID=cds1;Parent=tx1
chr1	src	exon	300	400	.	+	.	ID=ex2;Parent=tx1
chr1	src	CDS	300	380	.	+	0	ID=cds2;Parent=tx1
chr1	src	mRNA	150	450	.	+	.	ID=tx2;Parent=gene1
chr1	src	exon	150	450	.	+	.	ID=ex3;Parent=tx2
chr1	src	CDS	160	400	.	+	0	ID=cds3;Parent=tx2
"""

SYNTHETIC_GFF = """\
chr1	src	gene	10	100	.	-	.	ID=geneS;Name=SYN
chr1	src	exon	10	40	.	-	.	Parent=geneS
chr1	src	CDS	20	40	.	-	0	Parent=geneS
chr1	src	exon	60	100	.	-	.	Parent=geneS
chr1	src	CDS	60	90	.	-	1	Parent=geneS
"""


def test_resolve_normal_transcript():
    features = load_features_from_lines(GENE_GFF.splitlines())
    models = resolve_transcripts(
        features,
        source_gene_id="gene1",
        transcript_ids=["tx1"],
    )
    tx = models["tx1"]
    assert tx.start == 100
    assert tx.end == 400
    assert tx.strand == "+"
    assert tx.has_cds
    assert not tx.synthetic
    assert any("ID=tx1" in line for line in tx.lines)
    assert sum(1 for line in tx.lines if "\texon\t" in line) == 2


def test_resolve_multiple_isoforms():
    features = load_features_from_lines(GENE_GFF.splitlines())
    models = resolve_transcripts(
        features,
        source_gene_id="gene1",
        transcript_ids=["tx1", "tx2"],
    )
    assert set(models) == {"tx1", "tx2"}
    assert models["tx2"].start == 150


def test_missing_transcript():
    features = load_features_from_lines(GENE_GFF.splitlines())
    with pytest.raises(KeyError):
        resolve_transcripts(
            features,
            source_gene_id="gene1",
            transcript_ids=["missing"],
        )


def test_synthetic_transcript():
    features = load_features_from_lines(SYNTHETIC_GFF.splitlines())
    models = resolve_transcripts(
        features,
        source_gene_id="geneS",
        transcript_ids=["geneS"],
    )
    tx = models["geneS"]
    assert tx.synthetic
    assert tx.has_cds
    assert tx.start == 10
    assert tx.end == 100
    assert tx.strand == "-"
    assert any("\tmRNA\t" in line and "ID=geneS" in line for line in tx.lines)
    # Children re-parented to geneS
    assert all(
        "Parent=geneS" in line
        for line in tx.lines
        if "\texon\t" in line or "\tCDS\t" in line
    )


def test_shift_gff_lines():
    lines = [
        "chr1\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=tx1",
        "chr1\tsrc\texon\t100\t150\t.\t+\t.\tParent=tx1",
    ]
    shifted = shift_gff_lines(lines, new_seqid="g0", origin=100)
    assert shifted[0].startswith("g0\t")
    cols = shifted[0].split("\t")
    assert cols[3] == "1"
    assert cols[4] == "101"
    cols2 = shifted[1].split("\t")
    assert cols2[3] == "1"
    assert cols2[4] == "51"


def test_shift_gff_lines_fails_closed_on_garbage():
    from helpers.gff_slice import shift_gff_lines

    with pytest.raises(ValueError, match="unparseable"):
        shift_gff_lines(
            ["chr1\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=tx1", "not-a-gff-line"],
            new_seqid="g0",
            origin=100,
        )


def test_dedupe_gff_lines_same_id_same_locus_collapses():
    from helpers.gff_slice import dedupe_gff_lines

    lines = [
        "chr1\tsrc\texon\t100\t200\t.\t+\t.\tID=ex1;Parent=txA",
        "chr1\tsrc\texon\t100\t200\t.\t+\t.\tID=ex1;Parent=txA,txB",
        "chr1\tsrc\texon\t300\t400\t.\t+\t.\tID=ex2;Parent=txA",
    ]
    out = dedupe_gff_lines(lines)
    assert len(out) == 2
    assert "ID=ex1;Parent=txA" in out[0]
    assert "ID=ex2" in out[1]


def test_dedupe_keeps_shared_id_cds_segments():
    """NCBI/RefSeq reuses one CDS ID across exons — all segments must survive."""
    from helpers.gff_slice import dedupe_gff_lines

    lines = [
        "chr1\tsrc\tCDS\t120\t200\t.\t+\t0\tID=cds-NP1;Parent=tx1",
        "chr1\tsrc\tCDS\t300\t380\t.\t+\t2\tID=cds-NP1;Parent=tx1",
        "chr1\tsrc\tCDS\t400\t450\t.\t+\t1\tID=cds-NP1;Parent=tx1",
    ]
    out = dedupe_gff_lines(lines)
    assert len(out) == 3
    assert "120\t200" in out[0]
    assert "300\t380" in out[1]
    assert "400\t450" in out[2]


def test_dedupe_keeps_anonymous_exons_at_different_coords():
    from helpers.gff_slice import dedupe_gff_lines

    lines = [
        "chr1\tsrc\texon\t100\t200\t.\t+\t.\tParent=txA",
        "chr1\tsrc\texon\t300\t400\t.\t+\t.\tParent=txA",
    ]
    out = dedupe_gff_lines(lines)
    assert len(out) == 2


def test_dedupe_collapses_identical_anonymous_lines():
    from helpers.gff_slice import dedupe_gff_lines

    line = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tParent=txA"
    out = dedupe_gff_lines([line, line])
    assert len(out) == 1
