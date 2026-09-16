#!/usr/bin/env python3
"""Two-tier, swappable (mmseqs2 | diamond) clustering backend for the
pangenome cluster-profile analysis (notes/superpowers/specs/
2026-09-13-pangenome-cluster-profile-design.md, component 1). Tier 1 is
the allele/ortholog unit used by every downstream step (presence,
frequency, co-occurrence, synteny); tier 2 re-clusters tier-1
representatives loosely, purely as a superfamily annotation label -- it
is never used for presence/frequency calls.

Usage:
  cluster_backend.py mmseqs-tier1 --fasta all_ingroup.fa --out_prefix tier1
  cluster_backend.py mmseqs-tier2 --fasta tier1_rep_seq.fasta --out_prefix tier2
  cluster_backend.py diamond-tier1 --fasta all_ingroup.fa --out_prefix tier1
  cluster_backend.py diamond-tier2 --fasta tier1_rep_seq.fasta --out_prefix tier2

Input FASTA requirement: `all_ingroup.fa` is every strain's (isoform-collapsed)
proteome concatenated, with each header Short-PREFIXED as
`><Short>|<original_protein_id>`. The cluster TSV carries nothing but sequence
IDs, so that prefix is the only way build_presence_matrix.py can recover which
strain a family member came from -- see that script's docstring for the exact
convention and a worked prefixing one-liner.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import build_families, read_cluster_tsv  # noqa: E402


def run_mmseqs_cluster(
    fasta: Path, out_prefix: str, min_seq_id: float, cov: float, cluster_reassign: bool
) -> Path:
    cmd = [
        "mmseqs", "easy-cluster", str(fasta), out_prefix, "tmp_mmseqs",
        "--min-seq-id", str(min_seq_id), "-c", str(cov), "--cov-mode", "0",
    ]
    if cluster_reassign:
        cmd.append("--cluster-reassign")
    subprocess.run(cmd, check=True)
    return Path(f"{out_prefix}_cluster.tsv")


def run_diamond_cluster(
    fasta: Path, out_prefix: str, approx_id: float, member_cover: float, threads: int | None = None
) -> Path:
    cmd = [
        "diamond", "cluster", "-d", str(fasta), "-o", f"{out_prefix}_cluster.tsv",
        "--approx-id", str(approx_id), "--member-cover", str(member_cover),
    ]
    if threads:
        cmd += ["--threads", str(threads)]
    subprocess.run(cmd, check=True)
    return Path(f"{out_prefix}_cluster.tsv")


def two_tier_families(tier1_cluster_tsv: str | Path, tier2_cluster_tsv: str | Path) -> dict:
    """Combine tier-1 family membership with the tier-2 superfamily label
    for each tier-1 representative. Returns
    {tier1_rep: {"members": [...], "superfamily": tier2_rep_of_tier1_rep}}."""
    tier1_families = build_families(read_cluster_tsv(tier1_cluster_tsv))
    tier1_rep_to_tier2_rep = read_cluster_tsv(tier2_cluster_tsv)
    return {
        rep: {
            "members": members,
            "superfamily": tier1_rep_to_tier2_rep.get(rep, rep),
        }
        for rep, members in tier1_families.items()
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)
    for name in ("mmseqs-tier1", "mmseqs-tier2", "diamond-tier1", "diamond-tier2"):
        p = sub.add_parser(name)
        p.add_argument("--fasta", required=True, type=Path)
        p.add_argument("--out_prefix", required=True)
        if name.startswith("diamond"):
            p.add_argument("--threads", type=int, default=None)
    args = ap.parse_args()

    if args.command == "mmseqs-tier1":
        run_mmseqs_cluster(args.fasta, args.out_prefix, min_seq_id=0.9, cov=0.8, cluster_reassign=True)
    elif args.command == "mmseqs-tier2":
        run_mmseqs_cluster(args.fasta, args.out_prefix, min_seq_id=0.4, cov=0.8, cluster_reassign=False)
    elif args.command == "diamond-tier1":
        run_diamond_cluster(args.fasta, args.out_prefix, approx_id=90, member_cover=80, threads=args.threads)
    elif args.command == "diamond-tier2":
        run_diamond_cluster(args.fasta, args.out_prefix, approx_id=40, member_cover=80, threads=args.threads)


if __name__ == "__main__":
    main()
