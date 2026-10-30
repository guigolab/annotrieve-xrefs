"""Lightweight process stats for request timing / memory logs."""
from __future__ import annotations

import resource


def rss_mb() -> float:
    """Current process RSS in MiB.

    Prefers ``VmRSS`` from ``/proc/self/status`` (Linux). Falls back to peak
    ``ru_maxrss`` (KiB on Linux). Returns ``-1.0`` if unavailable.
    """
    try:
        with open("/proc/self/status", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    parts = line.split()
                    if len(parts) >= 2:
                        return int(parts[1]) / 1024.0
    except (OSError, ValueError):
        pass
    try:
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    except (OSError, ValueError, AttributeError):
        return -1.0
