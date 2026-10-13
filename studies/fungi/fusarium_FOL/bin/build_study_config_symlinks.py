#!/usr/bin/env python3
"""Run NII bin/build_study_config.py unchanged, but link instead of copy.

build_study_config.py copies every local_faa/local_genome/local_gff3 file
into data_dir/ (shutil.copyfile). For ~330 BFD strains that is ~30 GB of
duplicate data on /bigdata. This wrapper imports the tool and replaces its
module-level shutil.copyfile with an os.symlink to the source path, so
data_dir/pep, dna, gff3 hold symlinks. Everything else (stem naming,
config.csv, DATA_MANIFEST.yaml records) is the tool's own code. The
manifest checksum is sha256 of the link target, i.e. the source file; the
record's derived_by text still says "copied" -- it is a symlink.

Usage: build_study_config_symlinks.py --study-dir <study> --skip-fetch
"""
from __future__ import annotations

import os
import sys
import types
from pathlib import Path

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII / "bin"))
import build_study_config as bsc  # noqa: E402


def _link(src, dst, *a, **k):
    src, dst = Path(src).resolve(), Path(dst)
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    os.symlink(src, dst)
    return str(dst)


bsc.shutil = types.SimpleNamespace(**{n: getattr(bsc.shutil, n) for n in dir(bsc.shutil) if not n.startswith("__")})
bsc.shutil.copyfile = _link

if __name__ == "__main__":
    raise SystemExit(bsc.main())
