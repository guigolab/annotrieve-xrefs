"""Search index copy in RAM."""
from __future__ import annotations

import os
import shutil

import pytest

import helpers.presence as presence_mod
import settings as settings_mod

from tests.test_presence_api import presence_client  # noqa: F401  (reuse the fixture)


@pytest.fixture
def ram(tmp_path, monkeypatch):
    d = tmp_path / "ram"
    monkeypatch.setattr(settings_mod.settings, "PRESENCE_RAM_DIR", str(d))
    presence_mod._RAM_FAILED.clear()
    presence_mod.close_presence()
    yield d
    presence_mod.close_presence()
    presence_mod._RAM_FAILED.clear()


def test_copy_is_used_and_loci_path_stays(presence_client, ram):  # noqa: F811
    client, _ = presence_client
    reader = presence_mod.presence_reader()
    assert reader.db_path == ram / "v1" / "presence.sqlite"
    assert reader.db_path.is_file()
    assert reader.path.parent.name == "v1" and ram not in reader.path.parents
    body = client.get("/hits", params={"names": "TP53"}).json()
    assert body["terms"] and (ram / "v1" / "presence.sqlite").is_file()


def test_second_reader_reuses_copy(presence_client, ram):  # noqa: F811
    first = presence_mod.presence_reader().db_path
    mtime = os.stat(first).st_mtime_ns
    presence_mod.close_presence()
    assert presence_mod.presence_reader().db_path == first
    assert os.stat(first).st_mtime_ns == mtime


def test_too_large_falls_back_to_disk(presence_client, ram, monkeypatch):  # noqa: F811
    monkeypatch.setattr(settings_mod.settings, "PRESENCE_RAM_MAX_MB", 0)
    reader = presence_mod.presence_reader()
    assert reader.db_path == reader.path
    assert not ram.exists() or not any(ram.rglob("presence.sqlite"))


def test_not_enough_space_falls_back(presence_client, ram, monkeypatch):  # noqa: F811
    monkeypatch.setattr(shutil, "disk_usage",
                        lambda p: shutil._ntuple_diskusage(10, 10, 0))
    reader = presence_mod.presence_reader()
    assert reader.db_path == reader.path


def test_recheck_interval_skips_filesystem(presence_client, ram, monkeypatch):  # noqa: F811
    monkeypatch.setattr(settings_mod.settings, "PRESENCE_RECHECK_S", 60.0)
    presence_mod.close_presence()
    r1 = presence_mod.presence_reader()
    monkeypatch.setattr(presence_mod, "presence_db_path",
                        lambda: (_ for _ in ()).throw(AssertionError("stat called")))
    assert presence_mod.presence_reader() is r1
