"""Thread-local gene_corpus reader tests."""
from __future__ import annotations

import threading
from pathlib import Path

import pytest


def test_open_reader_is_thread_local(corpus_dir: Path) -> None:
    import helpers.session as session_mod

    conns: list[object] = []
    errors: list[BaseException] = []

    def worker() -> None:
        try:
            conn = session_mod.open_reader()
            row = conn.execute("SELECT 1 AS n").fetchone()
            assert int(row["n"]) == 1
            conns.append(conn)
        except BaseException as exc:  # noqa: BLE001
            errors.append(exc)

    threads = [threading.Thread(target=worker) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert not errors
    assert len(conns) == 2
    assert conns[0] is not conns[1]


def test_close_reader_clears_handles(corpus_dir: Path) -> None:
    import helpers.session as session_mod

    first = session_mod.open_reader()
    session_mod.close_reader()
    second = session_mod.open_reader()
    assert first is not second
    second.execute("SELECT 1").fetchone()
