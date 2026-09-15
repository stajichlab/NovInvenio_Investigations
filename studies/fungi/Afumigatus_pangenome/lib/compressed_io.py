"""Transparent gzip/zstd input support for this study's `bin/` scripts --
per the general storage-compression convention (large text intermediates
should default to compressed, and scripts that read them should accept
compressed input directly rather than requiring a manual decompression
step first).

`.gz` via the stdlib `gzip` module; `.zst` by shelling out to the `zstd`
CLI (no extra Python dependency needed -- `zstd` is already a pixi
dependency for this workspace). Plain, uncompressed files still work
unchanged.
"""
from __future__ import annotations

import gzip
import subprocess
from pathlib import Path
from typing import IO


def open_maybe_compressed(path: str | Path) -> IO[str]:
    """Open `path` for text reading, transparently decompressing based on
    its extension (.gz, .zst) or returning a plain text handle otherwise.

    The .zst case shells out to `zstd -dc` and returns its stdout as a
    text stream -- there is no line-iteration difference from a normal
    file handle, but the subprocess is not explicitly waited on here;
    for a short-lived script that reads the whole stream and exits, the
    OS reaps it on process exit, which is an acceptable simplification
    for how this helper is actually used (one big streaming read, not a
    long-lived server)."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith(".zst"):
        proc = subprocess.Popen(
            ["zstd", "-dc", path], stdout=subprocess.PIPE, text=True,
        )
        return proc.stdout
    return open(path)
