"""Process-wide logging so app ``logger.info`` reaches gunicorn/Docker stdout."""
from __future__ import annotations

import logging
import os
import sys

_BANNER_EMITTED = False


def _stdio_handler(root: logging.Logger) -> logging.Handler | None:
    for handler in root.handlers:
        if not isinstance(handler, logging.StreamHandler):
            continue
        if getattr(handler, "stream", None) in (sys.stderr, sys.stdout):
            return handler
    return None


def configure_logging() -> None:
    """
    Ensure INFO+ app logs go to stderr (what ``docker logs`` shows).

    Idempotent: safe at import, in ``create_app``, and in lifespan. If a stdio
    handler is missing (e.g. after gunicorn reconfigures logging), one is
    added; if one already exists, its level/format are updated in place.
    Override level with ``LOG_LEVEL`` (DEBUG, INFO, WARNING, …).
    """
    global _BANNER_EMITTED

    level_name = (os.getenv("LOG_LEVEL") or "INFO").strip().upper()
    level = getattr(logging, level_name, logging.INFO)

    root = logging.getLogger()
    root.setLevel(level)

    formatter = logging.Formatter(
        fmt="%(asctime)s %(levelname)s [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    handler = _stdio_handler(root)
    if handler is None:
        handler = logging.StreamHandler(sys.stderr)
        root.addHandler(handler)
    handler.setLevel(level)
    handler.setFormatter(formatter)

    if level > logging.DEBUG:
        logging.getLogger("httpx").setLevel(logging.WARNING)
        logging.getLogger("httpcore").setLevel(logging.WARNING)

    if not _BANNER_EMITTED:
        _BANNER_EMITTED = True
        logging.getLogger(__name__).info(
            "logging configured level=%s handlers=%d",
            logging.getLevelName(root.level),
            len(root.handlers),
        )
