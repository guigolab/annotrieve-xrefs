"""Run gffread to extract CDS / transcript / protein FASTA."""
from __future__ import annotations

import os
import shutil
import subprocess
from enum import Enum
from pathlib import Path

from settings import settings


class SequenceAction(str, Enum):
    cds = "cds"
    transcript = "transcript"
    protein = "protein"


class GffreadError(Exception):
    def __init__(self, code: str, message: str, **extra) -> None:
        super().__init__(message)
        self.code = code
        self.extra = extra


_ACTION_FLAG = {
    SequenceAction.cds: "-x",
    SequenceAction.transcript: "-w",
    SequenceAction.protein: "-y",
}


def gffread_bin() -> str:
    configured = (settings.GFFREAD_BIN or "gffread").strip() or "gffread"
    if "/" in configured or configured.startswith("."):
        if Path(configured).is_file() and os.access(configured, os.X_OK):
            return configured
        raise GffreadError(
            "gffread_unavailable",
            f"gffread binary not found: {configured}",
        )
    found = shutil.which(configured)
    if not found:
        raise GffreadError(
            "gffread_unavailable",
            f"gffread binary not found: {configured}",
        )
    return found


def run_gffread(
    *,
    gff_path: Path,
    fasta_path: Path,
    out_path: Path,
    action: SequenceAction,
    buffer_bp: int = 0,
    include_exon_coords: bool = False,
    timeout_s: float | None = None,
) -> None:
    """
    Invoke gffread with an argument list (no shell).

    ``buffer_bp`` maps to ``--w-add`` and is only valid for ``transcript``.
    """
    binary = gffread_bin()
    flag = _ACTION_FLAG[action]
    cmd = [
        binary,
        str(gff_path),
        "-g",
        str(fasta_path),
        flag,
        str(out_path),
    ]
    if include_exon_coords:
        cmd.append("-W")
    if buffer_bp:
        if action != SequenceAction.transcript:
            raise GffreadError(
                "buffer_not_supported",
                "buffer_bp is only supported with action=transcript",
                buffer_bp=buffer_bp,
                action=action.value,
            )
        cmd.extend(["--w-add", str(int(buffer_bp))])

    timeout = (
        settings.GFFREAD_TIMEOUT_S if timeout_s is None else float(timeout_s)
    )
    try:
        proc = subprocess.run(
            cmd,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise GffreadError(
            "gffread_failed",
            f"gffread timed out after {timeout}s",
            action=action.value,
        ) from exc
    except FileNotFoundError as exc:
        raise GffreadError(
            "gffread_unavailable",
            f"gffread binary not found: {binary}",
        ) from exc

    if proc.returncode != 0:
        stderr = (proc.stderr or proc.stdout or "").strip()
        raise GffreadError(
            "gffread_failed",
            f"gffread exited with {proc.returncode}",
            action=action.value,
            stderr=stderr[-2000:],
        )

    if not out_path.is_file() or out_path.stat().st_size == 0:
        stderr = (proc.stderr or "").strip()
        raise GffreadError(
            "gffread_failed",
            "gffread produced empty output",
            action=action.value,
            stderr=stderr[-2000:] if stderr else None,
        )


def count_fasta_records(path: Path) -> int:
    count = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def rewrite_fasta_headers(
    in_path: Path,
    out_path: Path,
    *,
    header_for_id: dict[str, str],
) -> list[str]:
    """
    Stream ``in_path`` to ``out_path``, replacing deflines.

    ``header_for_id`` maps the first whitespace-delimited token of the
    gffread defline (transcript id) to the full new header (without ``>``).
    Returns the ordered list of transcript ids written.
    """
    written: list[str] = []
    current_id: str | None = None
    with open(in_path, encoding="utf-8", errors="replace") as src, open(
        out_path, "w", encoding="utf-8"
    ) as dst:
        for line in src:
            if line.startswith(">"):
                token = line[1:].strip().split()[0] if line[1:].strip() else ""
                # gffread may emit ``id`` or ``id seqid:coords``. Exact match only.
                tid = token.split("|")[0]
                header = header_for_id.get(tid)
                if header is None and token != tid:
                    header = header_for_id.get(token)
                    if header is not None:
                        tid = token
                if header is None:
                    raise GffreadError(
                        "gffread_failed",
                        f"unexpected gffread record id {token!r}",
                        transcript_id=token,
                    )
                dst.write(f">{header}\n")
                written.append(tid)
                current_id = tid
            else:
                if current_id is None:
                    continue
                dst.write(line)
    return written


def split_fasta_to_record_files(
    in_path: Path,
    out_dir: Path,
    *,
    header_for_id: dict[str, str],
) -> dict[str, Path]:
    """
    Stream ``in_path``, rewrite headers, and write one file per record.

    Returns ``{transcript_id: path}``. Peak memory is one record's sequence
    buffer only (deflines are rewritten immediately).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, Path] = {}
    current_id: str | None = None
    current_fh = None
    try:
        with open(in_path, encoding="utf-8", errors="replace") as src:
            for line in src:
                if line.startswith(">"):
                    if current_fh is not None:
                        current_fh.close()
                        current_fh = None
                    token = (
                        line[1:].strip().split()[0] if line[1:].strip() else ""
                    )
                    tid = token.split("|")[0]
                    header = header_for_id.get(tid)
                    if header is None and token != tid:
                        header = header_for_id.get(token)
                        if header is not None:
                            tid = token
                    if header is None:
                        raise GffreadError(
                            "gffread_failed",
                            f"unexpected gffread record id {token!r}",
                            transcript_id=token,
                        )
                    if tid in result:
                        raise GffreadError(
                            "gffread_failed",
                            f"duplicate gffread record id {tid!r}",
                            transcript_id=tid,
                        )
                    path = out_dir / f"{len(result):04d}.fa"
                    current_fh = open(path, "w", encoding="utf-8")
                    current_fh.write(f">{header}\n")
                    result[tid] = path
                    current_id = tid
                else:
                    if current_fh is None or current_id is None:
                        continue
                    current_fh.write(line)
    finally:
        if current_fh is not None:
            current_fh.close()
    return result
