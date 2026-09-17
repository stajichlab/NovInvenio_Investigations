#!/usr/bin/env python3
"""Downloads starbase's Zenodo-archived reference SQLite dump (Starship
annotations across 19,863 fungal genomes, incl. a curated set for
*A. fumigatus* Af293) into this study's `data/` (class-1 ephemeral,
gitignored, regenerable -- see studies/*/*/data/ in the repo root
.gitignore, added 2026-09-17). Used by bin/starbase_crossvalidation.py.

Source: starbase (starbase.serve.scilifelab.se / github.com/FungAGE/starbase),
Zenodo record DOI 10.5281/zenodo.17533381 (Forsythe, Gluck-Thaler & Vogan,
"starbase: A Database and Toolkit for Exploring Giant Cargo-Carrying Mobile
Elements in Fungi", v0.1.0-pre, published 2025-11-05, CC-BY-4.0).

Usage:
  fetch_starbase_reference.py --output_dir data/starbase
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path

ZENODO_RECORD = "17533381"
FILES = {
    "starbase_v0.1.0-pre.sqlite": {
        "url": f"https://zenodo.org/api/records/{ZENODO_RECORD}/files/starbase_v0.1.0-pre.sqlite/content",
        "md5": "fc946b8e5ac8c7487f04a85db4a7b1c1",
        "size": 382242816,
    },
    "starship-scan.py": {
        "url": f"https://zenodo.org/api/records/{ZENODO_RECORD}/files/starship-scan.py/content",
        "md5": "6f476c360073c5753405d6efde773a73",
        "size": 20212,
    },
}


def _md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_starbase_files(output_dir: Path, force: bool = False) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    fetched = []
    for name, meta in FILES.items():
        dest = output_dir / name
        if dest.exists() and not force and _md5sum(dest) == meta["md5"]:
            print(f"fetch_starbase_reference: {name} already present with correct checksum, skipping",
                  file=sys.stderr)
            fetched.append(dest)
            continue
        print(f"fetch_starbase_reference: downloading {name} ({meta['size']:,} bytes)...", file=sys.stderr)
        urllib.request.urlretrieve(meta["url"], dest)
        actual_md5 = _md5sum(dest)
        if actual_md5 != meta["md5"]:
            raise RuntimeError(
                f"{name}: checksum mismatch (expected {meta['md5']}, got {actual_md5}) -- "
                "download is corrupt or the Zenodo record changed; do not trust this file."
            )
        print(f"fetch_starbase_reference: {name} downloaded and checksum-verified", file=sys.stderr)
        fetched.append(dest)
    return fetched


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output_dir", default="data/starbase")
    ap.add_argument("--force", action="store_true", help="re-download even if a correct-checksum copy exists")
    args = ap.parse_args()
    fetch_starbase_files(Path(args.output_dir), force=args.force)


if __name__ == "__main__":
    main()
