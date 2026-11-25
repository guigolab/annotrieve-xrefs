"""Slots and streaming bodies are given back whatever the client does."""
from __future__ import annotations

import gc

from helpers.concurrency import GuardedBody


def _body(log):
    def gen():
        try:
            yield "a"
            yield "b"
        finally:
            log.append("generator")
    return gen()


def test_cleanup_runs_once_after_a_full_read():
    log = []
    body = GuardedBody(_body(log), lambda: log.append("cleanup"))
    assert list(body) == ["a", "b"]
    body.close()
    del body
    gc.collect()
    assert log.count("cleanup") == 1


def test_cleanup_runs_when_the_body_was_never_started():
    """A generator that never ran does not execute its ``finally`` when it is collected."""
    log = []
    body = GuardedBody(_body(log), lambda: log.append("cleanup"))
    del body
    gc.collect()
    assert log == ["cleanup"]


def test_cleanup_runs_when_the_client_leaves_half_way():
    log = []
    body = GuardedBody(_body(log), lambda: log.append("cleanup"))
    assert next(body) == "a"
    del body
    gc.collect()
    assert log.count("cleanup") == 1 and "generator" in log


def test_cleanup_runs_when_the_body_raises():
    log = []

    def failing():
        yield "a"
        raise RuntimeError("boom")

    body = GuardedBody(failing(), lambda: log.append("cleanup"))
    next(body)
    try:
        next(body)
    except RuntimeError:
        pass
    assert log == ["cleanup"]
