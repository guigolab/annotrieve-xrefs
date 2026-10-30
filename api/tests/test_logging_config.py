"""Ensure app INFO logs actually reach stderr (Docker / gunicorn)."""
from __future__ import annotations

import logging

from helpers.logging_config import configure_logging


def test_configure_logging_emits_info_to_stderr(capsys):
    configure_logging()
    logging.getLogger("routers.sequences").info(
        "POST /sequences start action=%s items=%d",
        "cds",
        2,
    )
    logging.getLogger("services.sequences").info(
        "sequences phase=resolve groups=%d wall_s=%.2f",
        1,
        0.01,
    )
    err = capsys.readouterr().err
    assert "POST /sequences start action=cds items=2" in err
    assert "sequences phase=resolve groups=1" in err
    assert "[routers.sequences]" in err
    assert "[services.sequences]" in err


def test_configure_logging_idempotent_no_duplicate_handlers(capsys):
    configure_logging()
    before = len(logging.getLogger().handlers)
    configure_logging()
    configure_logging()
    assert len(logging.getLogger().handlers) == before
    logging.getLogger("helpers.logging_config").info("idempotent-smoke")
    assert "idempotent-smoke" in capsys.readouterr().err
