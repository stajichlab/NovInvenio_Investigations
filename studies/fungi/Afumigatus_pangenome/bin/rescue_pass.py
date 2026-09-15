#!/usr/bin/env python3
"""Genome-level rescue pass: tblastn every tier-1 family representative
against each strain's own genome to catch coverage-failed protein-model
absences (fragmented/split gene models, draft-assembly artifacts) --
notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md,
component 1b, Fable review finding 2.

Usage:
  rescue_pass.py --matrix presence_matrix.tsv --tblastn_tsv all_vs_all.tblastn.tsv \
      --output presence_matrix.rescued.tsv

`--tblastn_tsv` accepts a plain, `.gz`, or `.zst` file, and may be repeated
to read several chunked tblastn outputs (e.g. from
run_rescue_pass_tblastn_chunked.sh's per-array-task .tsv.zst files)
without concatenating them first.

Upstream isoform collapse: before building the protein FASTA fed into
Task 4's clustering, run NovInvenio's existing bin/collapse_isoforms.py
per strain (see NovInvenio's CLAUDE.md "collapse_isoforms.py" entry) so
per-strain copy counts reflect genes, not alternative transcripts:

  NovInvenio/bin/collapse_isoforms.py \
      --protein-fasta data_dir/pep/<Short>.pep.fa \
      --feature-table <Short>_feature_table.txt.gz \
      --output data_dir/pep_collapsed/<Short>.pep.fa
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import ABSENT, GENOME_ONLY, PresenceMatrix  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402


def parse_tblastn_hits(
    lines: list[str], min_pident: float = 90.0, min_qcov: float = 80.0
) -> set[tuple[str, str]]:
    """Parse tblastn outfmt6 lines (with a trailing qcovs column) into the
    set of (family_rep, strain) pairs with a qualifying genomic hit.
    Subject IDs are expected as '<Short>|<contig>' (via renamed genome FASTA
    headers before makeblastdb, so the strain is recoverable from the hit)."""
    hits: set[tuple[str, str]] = set()
    for line in lines:
        line = line.rstrip("\n")
        # Skip blank lines and comment lines
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 13:  # Minimum columns for outfmt6 with qcovs
            continue
        try:
            family, subject, pident, qcovs = parts[0], parts[1], float(parts[2]), float(parts[-1])
        except (ValueError, IndexError):
            continue
        if pident >= min_pident and qcovs >= min_qcov:
            strain = subject.split("|", 1)[0]
            hits.add((family, strain))
    return hits


def apply_rescue(matrix: PresenceMatrix, rescue_hits: set[tuple[str, str]]) -> tuple[int, int]:
    """Upgrade ABSENT calls to GENOME_ONLY wherever a qualifying rescue
    hit exists. Never downgrades an existing PRESENT call.

    Returns: (num_applied, num_skipped) where skipped = hits with
    unrecognized strain or family."""
    applied = 0
    skipped = 0
    for family, strain in rescue_hits:
        if strain not in matrix.strains or family not in matrix.families:
            skipped += 1
            continue
        if matrix.call(family, strain) == ABSENT:
            matrix.set_call(family, strain, GENOME_ONLY)
            applied += 1
    return applied, skipped


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument(
        "--tblastn_tsv", required=True, action="append",
        help="plain, .gz, or .zst tblastn outfmt6+qcovs file; repeat to read "
        "several chunked outputs (e.g. run_rescue_pass_tblastn_chunked.sh's "
        "per-array-task files) without concatenating them first",
    )
    ap.add_argument("--min_pident", type=float, default=90.0)
    ap.add_argument("--min_qcov", type=float, default=80.0)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)
    hits: set[tuple[str, str]] = set()
    for tblastn_path in args.tblastn_tsv:
        with open_maybe_compressed(tblastn_path) as fh:
            hits |= parse_tblastn_hits(fh.readlines(), args.min_pident, args.min_qcov)

    applied, skipped = apply_rescue(matrix, hits)
    print(f"Rescue: {applied} ABSENT→GENOME_ONLY, {skipped} skipped (unrecognized strain/family)", file=sys.stderr)

    if hits and skipped == len(hits):
        print(f"ERROR: All {len(hits)} parsed tblastn hits were skipped (likely wrong genome-DB naming)", file=sys.stderr)
        sys.exit(1)

    matrix.to_tsv(args.output)


if __name__ == "__main__":
    main()
