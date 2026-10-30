"""Annotation GFF path + seqid→accession resolver tests."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from helpers.annotation_files import (
    AnnotationFilesError,
    clear_annotation_file_caches,
    is_insdc_accession,
    resolve_annotation_files,
    resolve_seqid_accession,
)


@pytest.fixture(autouse=True)
def _clear_caches():
    clear_annotation_file_caches()
    yield
    clear_annotation_file_caches()


def _write_gff_pair(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"##gff-version 3\n")
    (Path(str(path) + ".csi")).write_bytes(b"csi")
    return path


def test_is_insdc_accession():
    assert is_insdc_accession("NC_000017.11")
    assert is_insdc_accession("CM085869.1")
    assert not is_insdc_accession("1")
    assert not is_insdc_accession("chr1")
    assert not is_insdc_accession("NC_000017")  # needs version


def test_resolve_by_convention(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = tmp_path / "annotations"
    gff = _write_gff_pair(
        root / "9606" / "GCF_000001405.40" / "RefSeq_abc.gff.gz"
    )
    monkeypatch.setenv("ANNOTATIONS_DIR", str(root))
    import settings as settings_mod

    settings_mod.settings.ANNOTATIONS_DIR = str(root)
    settings_mod.settings.ANNOTATIONS_GFF_PATH_PREFIX = ""

    files = resolve_annotation_files(
        taxid=9606,
        assembly_accession="GCF_000001405.40",
        gff_path="/harvest/root/9606/GCF_000001405.40/RefSeq_abc.gff.gz",
        source_database="RefSeq",
    )
    assert files.gff_path == gff
    assert files.csi_path == Path(str(gff) + ".csi")


def test_resolve_by_rebase_prefix(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = tmp_path / "annotations"
    harvest = tmp_path / "harvest"
    gff = _write_gff_pair(
        root / "9606" / "GCA_1" / "Ensembl_x.gff.gz"
    )
    # Same relative layout under harvest (not used as the live file).
    _ = gff
    monkeypatch.setenv("ANNOTATIONS_DIR", str(root))
    import settings as settings_mod

    settings_mod.settings.ANNOTATIONS_DIR = str(root)
    settings_mod.settings.ANNOTATIONS_GFF_PATH_PREFIX = str(harvest)

    files = resolve_annotation_files(
        taxid=9606,
        assembly_accession="GCA_1",
        gff_path=str(harvest / "9606" / "GCA_1" / "Ensembl_x.gff.gz"),
        source_database="Ensembl",
    )
    assert files.gff_path.name == "Ensembl_x.gff.gz"
    assert files.gff_path.is_file()


def test_missing_annotations_dir(monkeypatch: pytest.MonkeyPatch):
    import settings as settings_mod

    settings_mod.settings.ANNOTATIONS_DIR = ""
    with pytest.raises(AnnotationFilesError) as exc:
        resolve_annotation_files(
            taxid=1,
            assembly_accession="GCA_1",
            gff_path="/x.gff.gz",
        )
    assert exc.value.code == "annotations_unavailable"


def test_missing_gff(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = tmp_path / "annotations"
    root.mkdir()
    import settings as settings_mod

    settings_mod.settings.ANNOTATIONS_DIR = str(root)
    settings_mod.settings.ANNOTATIONS_GFF_PATH_PREFIX = ""
    with pytest.raises(AnnotationFilesError) as exc:
        resolve_annotation_files(
            taxid=9606,
            assembly_accession="GCA_1",
            gff_path="/9606/GCA_1/missing.gff.gz",
        )
    assert exc.value.code == "annotation_file_not_found"


def test_seqid_already_insdc(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = tmp_path / "annotations"
    gff = _write_gff_pair(root / "9606" / "GCF_1" / "x.gff.gz")
    import settings as settings_mod

    settings_mod.settings.ANNOTATIONS_DIR = str(root)
    files = resolve_annotation_files(
        taxid=9606,
        assembly_accession="GCF_1",
        gff_path=str(gff),
        source_database="RefSeq",
    )
    assert resolve_seqid_accession(files, "NC_000017.11") == "NC_000017.11"


def test_seqid_via_chromosomes_json(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = tmp_path / "annotations"
    ass = root / "9606" / "GCA_041900085.1"
    gff = _write_gff_pair(ass / "Ensembl_x.gff.gz")
    (ass / "chromosomes.json").write_text(
        json.dumps(
            [
                {
                    "chr_name": "1",
                    "sequence_name": "scaffoldA",
                    "ucsc_style_name": "chr1",
                    "genbank_accession": "CM085869.1",
                    "refseq_accession": "",
                    "canonical_id": "CM085869.1",
                }
            ]
        ),
        encoding="utf-8",
    )
    import settings as settings_mod

    settings_mod.settings.ANNOTATIONS_DIR = str(root)
    files = resolve_annotation_files(
        taxid=9606,
        assembly_accession="GCA_041900085.1",
        gff_path=str(gff),
        source_database="Ensembl",
    )
    assert resolve_seqid_accession(files, "1") == "CM085869.1"
    assert resolve_seqid_accession(files, "chr1") == "CM085869.1"


def test_seqid_prefer_refseq(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = tmp_path / "annotations"
    ass = root / "9606" / "GCF_000001405.40"
    gff = _write_gff_pair(ass / "RefSeq_x.gff.gz")
    (ass / "chromosomes.json").write_text(
        json.dumps(
            [
                {
                    "chr_name": "17",
                    "genbank_accession": "CM000679.2",
                    "refseq_accession": "NC_000017.11",
                    "canonical_id": "NC_000017.11",
                }
            ]
        ),
        encoding="utf-8",
    )
    import settings as settings_mod

    settings_mod.settings.ANNOTATIONS_DIR = str(root)
    files = resolve_annotation_files(
        taxid=9606,
        assembly_accession="GCF_000001405.40",
        gff_path=str(gff),
        source_database="RefSeq",
    )
    assert resolve_seqid_accession(files, "17") == "NC_000017.11"


def test_seqid_unresolved(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = tmp_path / "annotations"
    ass = root / "9606" / "GCA_1"
    gff = _write_gff_pair(ass / "x.gff.gz")
    import settings as settings_mod

    settings_mod.settings.ANNOTATIONS_DIR = str(root)
    files = resolve_annotation_files(
        taxid=9606,
        assembly_accession="GCA_1",
        gff_path=str(gff),
        source_database="Ensembl",
    )
    with pytest.raises(AnnotationFilesError) as exc:
        resolve_seqid_accession(files, "chrZ")
    assert exc.value.code == "seqid_not_resolved"
