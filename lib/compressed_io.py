"""Transparent-gzip text opener, shared by scripts that read/write the per-species
annotation TSVs (bin/extract_dat_annotations.py output, consumed by
bin/domain_enrichment.py, bin/go_enrichment.py, bin/merge_uniprot_annotations.py).

Same opener pattern extract_dat_annotations.py already used for UniProt's own
.dat.gz input (`gzip.open if path.suffix == ".gz" else open`) -- centralized here
so writers and readers agree on it, instead of every consumer re-implementing the
suffix check.
"""
import gzip
from pathlib import Path


def open_text(path: Path, mode: str = "rt", **kwargs):
    """Open `path` as text, transparently gzip-decompressing/-compressing when its
    suffix is `.gz`. `newline`/`encoding`/etc. pass through as with plain `open()`.
    """
    opener = gzip.open if Path(path).suffix == ".gz" else open
    return opener(path, mode, **kwargs)
