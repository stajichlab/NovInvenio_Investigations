#!/usr/bin/env python3
"""Sequence-based crosswalk between the paper's gene/protein IDs (AFUB_*,
Afu*g*, per-strain g#, XP_* accessions) and this study's own protein IDs
-- notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-
design.md, component 6's required ID-crosswalk control. String matching
does not work (the two ID schemes are unrelated); this runs the paper's
cited protein sequences (fetched separately, e.g. via NCBI efetch for the
accessions in Table S19/S3-S16, into a FASTA passed as --paper_fasta)
through diamond blastp against this study's own proteomes and keeps the
best hit per paper ID.

Usage:
  diamond makedb --in all_ingroup.fa -d all_ingroup
  diamond blastp -q paper_reference_proteins.fa -d all_ingroup \\
      -o paper_vs_study.tsv --outfmt 6
  id_crosswalk.py --diamond_tsv paper_vs_study.tsv --output crosswalk.tsv
"""
from __future__ import annotations

import argparse


def parse_diamond_blastp_besthits(lines: list[str]) -> dict[str, str]:
    """Parse diamond blastp outfmt6 lines and keep, for each query
    (paper gene/protein ID), the subject (this study's protein ID) with
    the highest bitscore."""
    best_bitscore: dict[str, float] = {}
    best_subject: dict[str, str] = {}
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        query, subject, bitscore = parts[0], parts[1], float(parts[11])
        if bitscore > best_bitscore.get(query, -1.0):
            best_bitscore[query] = bitscore
            best_subject[query] = subject
    return best_subject


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--diamond_tsv", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open(args.diamond_tsv) as fh:
        crosswalk = parse_diamond_blastp_besthits(fh.readlines())
    with open(args.output, "w") as fh:
        fh.write("paper_id\tstudy_protein_id\n")
        for paper_id, study_id in sorted(crosswalk.items()):
            fh.write(f"{paper_id}\t{study_id}\n")


if __name__ == "__main__":
    main()
