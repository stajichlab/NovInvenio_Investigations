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


def parse_diamond_blastp_hits_per_strain(
    lines: list[str], id_sep: str = "|"
) -> dict[tuple[str, str], dict]:
    """Parse diamond blastp outfmt6 lines against a Short-prefixed multi-
    strain subject database, keeping the best hit PER (query, strain) pair
    rather than collapsing to one global best hit (`besthits` above).

    Needed for a single-copy reference gene's presence/absence CALL across
    every strain (e.g. hacA/hrmA vs all 293 strains) -- collapsing to one
    global best hit would silently discard every strain but the single
    closest match, the wrong question here (is it present in strain X, not
    which strain's copy is the single best match overall).

    Returns {(query, strain): {"subject", "pident", "align_len", "bitscore",
    "evalue"}},
    keeping the highest-bitscore hit per (query, strain). A subject ID
    without `id_sep` (not Short-prefixed) is skipped and counted; if EVERY
    subject is unrecognized, raises ValueError (same "your IDs are wrong"
    guard other scripts in this study use).
    """
    best: dict[tuple[str, str], dict] = {}
    n_total, n_unrecognized = 0, 0
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12:
            continue
        n_total += 1
        query, subject = parts[0], parts[1]
        if id_sep not in subject:
            n_unrecognized += 1
            continue
        strain = subject.split(id_sep, 1)[0]
        bitscore = float(parts[11])
        key = (query, strain)
        if bitscore > best.get(key, {}).get("bitscore", -1.0):
            best[key] = {
                "subject": subject,
                "pident": float(parts[2]),
                "align_len": int(parts[3]),
                "bitscore": bitscore,
                "evalue": float(parts[10]),
            }
    if n_total and n_unrecognized == n_total:
        raise ValueError(
            f"parse_diamond_blastp_hits_per_strain: none of {n_total} hits had "
            f"a {id_sep!r}-separated Short prefix in the subject ID -- the "
            "subject database was not built from a Short-prefixed FASTA"
        )
    return best


def parse_hmmsearch_tblout_hits(
    lines: list[str], id_sep: str = "|"
) -> dict[str, list[dict]]:
    """Parse `hmmsearch --tblout` output (Short-prefixed target sequence
    IDs) into {strain: [{"target", "evalue", "bitscore"}]} for every strain
    with at least one hit -- used for a gene FAMILY (e.g. PF11001), where a
    strain can carry zero, one, or several matching proteins, unlike the
    single-locus diamond crosswalk above."""
    hits: dict[str, list[dict]] = {}
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        target = parts[0]
        if id_sep not in target:
            continue
        strain = target.split(id_sep, 1)[0]
        hits.setdefault(strain, []).append({
            "target": target,
            "evalue": float(parts[4]),
            "bitscore": float(parts[5]),
        })
    return hits


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
