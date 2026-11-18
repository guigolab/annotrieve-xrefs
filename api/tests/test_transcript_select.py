"""Transcript selection modes: all / longest / canonical_or_longest / canonical_only."""
from __future__ import annotations

import pytest

from helpers.gff_slice import load_features_from_lines
from helpers.transcript_select import (
    TranscriptMode,
    feature_tags,
    list_candidates,
    select_transcripts,
)

# txA: spliced 200, CDS 150.  txB: spliced 301 (single exon), CDS 100, tag Ensembl_canonical.
# txC: non-coding, spliced 500.  txD: same spliced length as txA (tie), CDS 150.
GENE = """\
chr1	s	gene	1	1000	.	+	.	ID=g1
chr1	s	mRNA	100	400	.	+	.	ID=txA;Parent=g1
chr1	s	exon	100	199	.	+	.	Parent=txA
chr1	s	exon	300	399	.	+	.	Parent=txA
chr1	s	CDS	120	199	.	+	0	Parent=txA
chr1	s	CDS	300	369	.	+	0	Parent=txA
chr1	s	mRNA	100	400	.	+	.	ID=txB;Parent=g1;tag=basic,Ensembl_canonical
chr1	s	exon	100	400	.	+	.	Parent=txB
chr1	s	CDS	150	249	.	+	0	Parent=txB
chr1	s	lnc_RNA	500	999	.	+	.	ID=txC;Parent=g1
chr1	s	exon	500	999	.	+	.	Parent=txC
chr1	s	mRNA	100	400	.	+	.	ID=txD;Parent=g1
chr1	s	exon	100	199	.	+	.	Parent=txD
chr1	s	exon	300	399	.	+	.	Parent=txD
chr1	s	CDS	120	199	.	+	0	Parent=txD
chr1	s	CDS	300	369	.	+	0	Parent=txD
"""


def _cands(text=GENE, gene="g1"):
    return list_candidates(load_features_from_lines(text.splitlines()), gene)


def _ids(sel):
    return [c.transcript_id for c in sel.chosen]


def test_candidates_lengths_and_tags():
    by_id = {c.transcript_id: c for c in _cands()}
    assert set(by_id) == {"txA", "txB", "txC", "txD"}
    assert by_id["txA"].spliced_len == 200 and by_id["txA"].cds_len == 150
    assert by_id["txB"].spliced_len == 301 and by_id["txB"].canonical
    assert by_id["txC"].cds_len == 0 and not by_id["txC"].canonical


def test_longest_depends_on_action_and_breaks_ties_by_id():
    c = _cands()
    assert _ids(select_transcripts(c, "longest", "transcript")) == ["txC"]   # 500 bp
    assert _ids(select_transcripts(c, "longest", "cds")) == ["txA"]          # 150 tie -> lowest id
    assert _ids(select_transcripts(c, "longest", "protein")) == ["txA"]


def test_all_skips_transcripts_that_cannot_give_the_action():
    c = _cands()
    assert _ids(select_transcripts(c, "all", "transcript")) == ["txA", "txB", "txC", "txD"]
    assert _ids(select_transcripts(c, "all", "cds")) == ["txA", "txB", "txD"]


def test_canonical_or_longest_prefers_flag_then_falls_back():
    c = _cands()
    sel = select_transcripts(c, "canonical_or_longest", "cds")
    assert _ids(sel) == ["txB"] and sel.rule == "canonical"
    no_flag = _cands(GENE.replace(";tag=basic,Ensembl_canonical", ""))
    sel = select_transcripts(no_flag, "canonical_or_longest", "cds")
    assert _ids(sel) == ["txA"] and sel.rule == "longest"


def test_canonical_only_drops_genes_without_a_flag():
    sel = select_transcripts(_cands(), TranscriptMode.canonical_only, "transcript")
    assert _ids(sel) == ["txB"]
    no_flag = _cands(GENE.replace(";tag=basic,Ensembl_canonical", ""))
    sel = select_transcripts(no_flag, "canonical_only", "transcript")
    assert sel.chosen == [] and sel.skipped_reason == "no_canonical"


def test_canonical_flag_without_cds_is_not_eligible_for_cds():
    text = GENE.replace(";tag=basic,Ensembl_canonical", "").replace(
        "ID=txC;Parent=g1", "ID=txC;Parent=g1;tag=MANE_Select")
    sel = select_transcripts(_cands(text), "canonical_only", "cds")
    assert sel.chosen == [] and sel.skipped_reason == "no_canonical"
    sel = select_transcripts(_cands(text), "canonical_only", "transcript")
    assert _ids(sel) == ["txC"]


def test_several_canonical_picks_the_longest_flagged():
    text = GENE.replace("ID=txA;Parent=g1", "ID=txA;Parent=g1;tag=MANE Select")
    sel = select_transcripts(_cands(text), "canonical_only", "transcript")
    assert _ids(sel) == ["txB"]  # 301 bp beats txA 200 bp


def test_tag_normalisation_covers_ncbi_and_ensembl_spellings():
    feats = load_features_from_lines([
        "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=a;tag=RefSeq Select,MANE Select",
        "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=b;Tag=Ensembl_canonical",
    ])
    assert feature_tags(feats[0]) == ("refseq_select", "mane_select")
    assert feature_tags(feats[1]) == ("ensembl_canonical",)


SYNTHETIC = """\
chr1	s	gene	10	100	.	-	.	ID=geneS
chr1	s	exon	10	40	.	-	.	Parent=geneS
chr1	s	CDS	20	40	.	-	0	Parent=geneS
chr1	s	exon	60	100	.	-	.	Parent=geneS
"""


def test_gene_with_direct_structure_gives_one_synthetic_candidate():
    c = _cands(SYNTHETIC, "geneS")
    assert len(c) == 1 and c[0].synthetic and c[0].transcript_id == "geneS"
    assert _ids(select_transcripts(c, "longest", "transcript")) == ["geneS"]


def test_empty_and_noncoding_reasons():
    assert select_transcripts([], "longest", "cds").skipped_reason == "no_transcript"
    nc = [c for c in _cands() if c.transcript_id == "txC"]
    assert select_transcripts(nc, "longest", "cds").skipped_reason == "no_cds"
    with pytest.raises(ValueError):
        select_transcripts(nc, "bogus", "cds")
