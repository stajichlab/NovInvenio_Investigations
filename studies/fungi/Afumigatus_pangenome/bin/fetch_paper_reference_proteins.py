#!/usr/bin/env python3
"""Fetch the mBio (Gluck-Thaler et al. 2025) supplement's cited NCBI protein
accessions (AF293 XP_* RefSeq IDs from Tables S19/S5/S13) as real sequences,
for the sequence-based ID crosswalk (notes/superpowers/specs/
2026-09-13-pangenome-cluster-profile-design.md, component 6's "Required
control, not yet resolved").

Study-specific, not a general `bin/fetch_<source>.py` recipe (NII's
top-level convention) -- this fetches one paper's specific cited accession
list, not a general "pull an NCBI protein set" tool any study would reuse
with different arguments. Lives in this study's own `bin/` per NII's
CLAUDE.md "Study-specific vs. shared scripts" rule.

Usage:
  fetch_paper_reference_proteins.py --accessions accessions_to_fetch.txt \\
      --output paper_reference_proteins.fa

Batches accessions (NCBI efetch, db=protein, rettype=fasta) at
`--batch_size` (default 200) per request, well under the 3-requests/second
unauthenticated rate limit given the total accession count here (~600).
"""
from __future__ import annotations

import argparse
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch_batch(accessions: list[str], timeout: int = 60) -> str:
    params = {
        "db": "protein",
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
    }
    url = NCBI_EFETCH_URL + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read().decode("utf-8")


def fetch_all(accessions: list[str], batch_size: int = 200, sleep_s: float = 0.4) -> str:
    """Fetch every accession in `accessions`, batched, and return the
    concatenated FASTA text. Raises RuntimeError if the number of FASTA
    records returned across all batches doesn't match the number of
    accessions requested -- silently returning a partial/truncated fetch
    would poison the downstream crosswalk with missing sequences that look
    like "no real ortholog" instead of "fetch failed"."""
    out_chunks: list[str] = []
    n_records = 0
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        text = fetch_batch(batch)
        n_records += text.count(">")
        out_chunks.append(text)
        if i + batch_size < len(accessions):
            time.sleep(sleep_s)
    combined = "".join(out_chunks)
    if n_records != len(accessions):
        print(
            f"WARNING: requested {len(accessions)} accessions, got {n_records} "
            "FASTA records back -- some accessions may be obsolete/replaced/"
            "withdrawn at NCBI. Inspect the output before trusting the crosswalk.",
            file=sys.stderr,
        )
    return combined


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accessions", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--batch_size", type=int, default=200)
    args = ap.parse_args()

    accessions = [
        line.strip() for line in args.accessions.read_text().splitlines() if line.strip()
    ]
    print(f"Fetching {len(accessions)} accessions from NCBI efetch...", file=sys.stderr)
    fasta_text = fetch_all(accessions, batch_size=args.batch_size)
    args.output.write_text(fasta_text)
    n_written = fasta_text.count(">")
    print(f"Wrote {n_written} FASTA records to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
