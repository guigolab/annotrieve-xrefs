"""Transcript rules: presets, rank_by, require, keep, parsing."""
from __future__ import annotations

import pytest

from helpers.gff_slice import load_features_from_lines
from helpers.transcript_select import RuleError, list_candidates, parse_rule, select_transcripts

# txA: 2 exons, spliced 200, CDS 150.   txB: 1 exon, spliced 301, CDS 100, tag basic.
# txC: non-coding lncRNA, spliced 500.  txD: same shape as txA (tie with txA).
GENE = """\
chr1	s	gene	1	1000	.	+	.	ID=g1
chr1	s	mRNA	100	400	.	+	.	ID=txA;Parent=g1;biotype=protein_coding
chr1	s	exon	100	199	.	+	.	Parent=txA
chr1	s	exon	300	399	.	+	.	Parent=txA
chr1	s	CDS	120	199	.	+	0	Parent=txA
chr1	s	CDS	300	369	.	+	0	Parent=txA
chr1	s	mRNA	100	400	.	+	.	ID=txB;Parent=g1;tag=basic,Ensembl_canonical;biotype=protein_coding
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


def cands(text=GENE, gene="g1"):
    return list_candidates(load_features_from_lines(text.splitlines()), gene)


def pick(spec, action, text=GENE, gene="g1"):
    sel = select_transcripts(cands(text, gene), parse_rule(spec, action), action)
    return [c.transcript_id for c in sel.chosen], sel


def test_candidate_metrics():
    by_id = {c.transcript_id: c for c in cands()}
    assert set(by_id) == {"txA", "txB", "txC", "txD"}
    a = by_id["txA"]
    assert (a.spliced_len, a.cds_len, a.n_exons, a.span, a.utr_len) == (200, 150, 2, 301, 50)
    assert a.biotype == "protein_coding" and by_id["txC"].biotype == "lnc_rna"
    assert not by_id["txC"].coding and by_id["txB"].tags == ("basic", "ensembl_canonical")


def test_default_longest_is_action_aware():
    assert pick(None, "cds")[0] == ["txA"]            # CDS 150 tie with txD -> lowest id
    assert pick(None, "protein")[0] == ["txA"]
    assert pick(None, "transcript")[0] == ["txB"]     # coding first, then longest spliced (301)
    assert pick("longest_spliced", "transcript")[0] == ["txC"]   # literal longest, non-coding


def test_cds_and_exon_presets():
    assert pick("longest_cds", "cds")[0] == ["txA"]
    assert pick("shortest_cds", "cds")[0] == ["txB"]
    assert pick("most_exons", "transcript")[0] == ["txA"]    # 2 exons, tie with txD
    assert pick("fewest_exons", "transcript")[0] == ["txB"]  # 1 exon, tie with txC -> lowest id


def test_all_returns_every_eligible_transcript_by_id():
    assert pick("all", "transcript")[0] == ["txA", "txB", "txC", "txD"]
    assert pick("all", "cds")[0] == ["txA", "txB", "txD"]


def test_custom_rule_rank_require_keep():
    spec = {"rank_by": ["n_exons:desc", "span"], "require": ["coding"], "keep": 2}
    ids, sel = pick(spec, "transcript")
    assert ids == ["txA", "txD"] and sel.rule == "custom" and sel.n_candidates == 4
    assert [c.transcript_id for c in sel.ranked] == ["txA", "txD", "txB"]


def test_require_tag_and_biotype():
    assert pick({"require": ["tag:Basic"]}, "transcript")[0] == ["txB"]
    assert pick({"require": ["biotype:protein_coding"], "keep": "all"}, "transcript")[0] == ["txA", "txB"]
    ids, sel = pick({"require": ["tag:nothing"]}, "transcript")
    assert ids == [] and sel.skipped_reason == "filtered_out"


def test_skip_reasons():
    nc = [c for c in cands() if c.transcript_id == "txC"]
    assert select_transcripts([], parse_rule(None, "cds"), "cds").skipped_reason == "no_transcript"
    assert select_transcripts(nc, parse_rule(None, "cds"), "cds").skipped_reason == "no_cds"


SYNTHETIC = """\
chr1	s	gene	10	100	.	-	.	ID=geneS
chr1	s	exon	10	40	.	-	.	Parent=geneS
chr1	s	CDS	20	40	.	-	0	Parent=geneS
chr1	s	exon	60	100	.	-	.	Parent=geneS
"""


def test_gene_with_direct_structure_gives_one_synthetic_candidate():
    c = cands(SYNTHETIC, "geneS")
    assert len(c) == 1 and c[0].synthetic and c[0].transcript_id == "geneS" and c[0].n_exons == 2
    ids, _ = pick("longest", "transcript", SYNTHETIC, "geneS")
    assert ids == ["geneS"]


@pytest.mark.parametrize("spec,code", [
    ("canonical_or_longest", "unknown_transcript_preset"),
    ({"rank_by": ["canonical"]}, "unknown_rank_key"),
    ({"rank_by": ["cds_len:up"]}, "invalid_transcript_rule"),
    ({"require": ["canonical"]}, "unknown_require"),
    ({"require": ["tag:"]}, "unknown_require"),
    ({"keep": 0}, "invalid_transcript_rule"),
    ({"keep": True}, "invalid_transcript_rule"),
    ({"bogus": 1}, "invalid_transcript_rule"),
    (42, "invalid_transcript_rule"),
])
def test_invalid_rules_give_coded_errors(spec, code):
    with pytest.raises(RuleError) as exc:
        parse_rule(spec, "cds")
    assert exc.value.code == code
