"""Locate the sibling NovInvenio checkout's ``lib/`` directory.

Several scripts in this study import NovInvenio's ``config_parser`` (the
``config.csv`` reader). NovInvenio is a *separate* repository that normally
sits next to NovInvenio_Investigations under a shared ``projects/`` directory::

    projects/
      NovInvenio/                    <- lib/config_parser.py lives here
      NovInvenio_Investigations/     <- this repo (or a git worktree of it)

A git worktree of this repo is NOT necessarily named
``NovInvenio_Investigations`` and is NOT necessarily inside ``projects/`` at the
same depth, so a fixed ``Path(__file__).resolve().parents[N]`` offset is wrong
(it was, before this module existed -- it pointed inside the worktree itself and
every import failed with a bare ``ModuleNotFoundError``).

Resolution order:

1. ``$NOVINVENIO_LIB`` -- a directory containing ``config_parser.py``.
2. ``$NOVINVENIO_ROOT`` -- a NovInvenio checkout root; its ``lib/`` is used.
3. Walk up from this file and, for each ancestor, try ``<ancestor>/NovInvenio/lib``
   (this finds the sibling checkout from an arbitrarily-named worktree).
4. Raise a clear, actionable ``RuntimeError``.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

#: Walking further up than this is pointless -- it leaves the projects directory.
_MAX_ANCESTORS = 8

#: A directory only counts as NovInvenio's lib/ if it actually has this module.
_SENTINEL = "config_parser.py"


def _is_novinvenio_lib(path: Path) -> bool:
    return (path / _SENTINEL).is_file()


def find_novinvenio_lib() -> Path:
    """Return the path to NovInvenio's ``lib/`` directory.

    Raises:
        RuntimeError: with instructions, if no checkout can be located.
    """
    tried: list[str] = []

    env_lib = os.environ.get("NOVINVENIO_LIB")
    if env_lib:
        candidate = Path(env_lib).expanduser().resolve()
        if _is_novinvenio_lib(candidate):
            return candidate
        tried.append(f"$NOVINVENIO_LIB={candidate}")

    env_root = os.environ.get("NOVINVENIO_ROOT")
    if env_root:
        candidate = (Path(env_root).expanduser() / "lib").resolve()
        if _is_novinvenio_lib(candidate):
            return candidate
        tried.append(f"$NOVINVENIO_ROOT/lib={candidate}")

    here = Path(__file__).resolve()
    for ancestor in list(here.parents)[:_MAX_ANCESTORS]:
        candidate = ancestor / "NovInvenio" / "lib"
        if _is_novinvenio_lib(candidate):
            return candidate
        tried.append(str(candidate))

    raise RuntimeError(
        "Could not locate the NovInvenio checkout's lib/ directory (needed for "
        "config_parser.py, the config.csv reader).\n"
        "Set NOVINVENIO_LIB to that lib/ directory, or NOVINVENIO_ROOT to the "
        "NovInvenio checkout root, e.g.:\n"
        "  export NOVINVENIO_ROOT=/path/to/projects/NovInvenio\n"
        "Locations tried:\n  " + "\n  ".join(tried)
    )


def add_novinvenio_lib_to_path() -> Path:
    """Prepend NovInvenio's ``lib/`` to ``sys.path`` (idempotent) and return it."""
    lib = find_novinvenio_lib()
    if str(lib) not in sys.path:
        sys.path.insert(0, str(lib))
    return lib
