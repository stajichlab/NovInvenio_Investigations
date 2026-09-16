#!/usr/bin/env python3
"""Tests whether a set of island member genes are real homologs of a named
reference gene cluster, rather than trusting a single SwissProt best-hit
name at face value.

Motivating case (2026-09-16): the `Asfu_HMR_AF_270` 8-gene island (187
strains) had a SwissProt best-hit named "Phomopsin biosynthesis cluster
protein D" (UniProt A0A8J9R8G6, PHOC1_DIALO, from the real Phomopsin
biosynthetic gene cluster in *Diaporthe leptostromiformis* described in
Ding et al. 2016 PNAS 113:3527 -- a ribosomally-synthesized, prenylated
cyclic peptide toxin pathway). A single best-hit name can look like a
specific, meaningful lead, but a biosynthetic gene cluster is defined by
co-inheritance of MULTIPLE genes together -- a real horizontally-shared or
orthologous cluster should show several of the island's genes hitting
several different cluster genes, not one coincidental domain-level match.
This script makes that distinction explicit and reusable for any future
named-cluster lead, instead of re-deriving the logic ad hoc each time.

Usage:
  cluster_homology_test.py --query_fasta island_genes.fasta \\
      --reference_fasta phomopsin_cluster_reference.fasta \\
      --evalue 1e-5 --qcov_threshold 0.5 \\
      --output cluster_homology_report.tsv
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

OUTFMT6_COLUMNS = (
    "qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
    "evalue bitscore qlen slen"
).split()


def parse_blast_outfmt6(lines: list[str]) -> list[dict]:
    """Parses `blastp -outfmt '6 qseqid sseqid pident length mismatch gapopen
    qstart qend sstart send evalue bitscore qlen slen'` output, adding a
    derived `qcov` (aligned length / query length) -- the coverage check
    that distinguishes a real full-length homolog from a single shared
    promiscuous domain."""
    hits = []
    for line in lines:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        row = dict(zip(OUTFMT6_COLUMNS, parts))
        row["pident"] = float(row["pident"])
        row["length"] = int(row["length"])
        row["evalue"] = float(row["evalue"])
        row["qlen"] = int(row["qlen"])
        row["slen"] = int(row["slen"])
        row["qcov"] = row["length"] / row["qlen"] if row["qlen"] else 0.0
        hits.append(row)
    return hits


def best_hit_per_query(hits: list[dict]) -> dict[str, dict]:
    """{qseqid: best_hit_row}, best = lowest e-value (blastp output is
    already sorted by bitscore descending per query, but this is explicit
    and doesn't depend on that ordering being preserved)."""
    best: dict[str, dict] = {}
    for hit in hits:
        q = hit["qseqid"]
        if q not in best or hit["evalue"] < best[q]["evalue"]:
            best[q] = hit
    return best


def summarize_cluster_homology(
    query_ids: list[str],
    best_hits: dict[str, dict],
    evalue_threshold: float,
    qcov_threshold: float,
) -> dict:
    """A query "qualifies" only if its best hit passes BOTH the e-value and
    query-coverage bar -- coverage is what rules out a single shared
    promiscuous domain (e.g. a Cupin or AAA-ATPase fold present in
    thousands of unrelated proteins) masquerading as a full-length hit.
    Verdict is "supported" only with >=2 qualifying queries hitting >=2
    distinct reference targets -- a single qualifying hit, however good,
    is exactly the "one coincidental match" case this script exists to
    catch, not evidence of a shared multi-gene cluster."""
    qualifying = {
        q: hit for q, hit in best_hits.items()
        if q in query_ids and hit["evalue"] <= evalue_threshold and hit["qcov"] >= qcov_threshold
    }
    n_distinct_targets = len({hit["sseqid"] for hit in qualifying.values()})
    verdict = "supported" if len(qualifying) >= 2 and n_distinct_targets >= 2 else "not_supported"
    return {
        "n_total": len(query_ids),
        "n_with_qualifying_hit": len(qualifying),
        "n_distinct_targets": n_distinct_targets,
        "verdict": verdict,
        "qualifying_hits": qualifying,
    }


def run_blastp(query_fasta: str, reference_fasta: str, evalue: float) -> list[str]:
    with tempfile.TemporaryDirectory() as tmp:
        db_prefix = str(Path(tmp) / "ref_db")
        subprocess.run(
            ["makeblastdb", "-in", reference_fasta, "-dbtype", "prot", "-out", db_prefix],
            check=True, capture_output=True,
        )
        result = subprocess.run(
            ["blastp", "-query", query_fasta, "-db", db_prefix,
             "-outfmt", "6 " + " ".join(OUTFMT6_COLUMNS),
             "-evalue", str(evalue), "-max_target_seqs", "5"],
            check=True, capture_output=True, text=True,
        )
        return result.stdout.splitlines()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query_fasta", required=True, help="island member representative sequences")
    ap.add_argument("--reference_fasta", required=True, help="named reference cluster's protein sequences")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--qcov_threshold", type=float, default=0.5)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    query_ids = []
    with open(args.query_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                query_ids.append(line[1:].strip().split()[0])
    print(f"cluster_homology_test: {len(query_ids)} query genes, running blastp...", file=sys.stderr)

    lines = run_blastp(args.query_fasta, args.reference_fasta, args.evalue)
    hits = parse_blast_outfmt6(lines)
    best_hits = best_hit_per_query(hits)
    summary = summarize_cluster_homology(query_ids, best_hits, args.evalue, args.qcov_threshold)

    print(f"cluster_homology_test: {summary['n_with_qualifying_hit']}/{summary['n_total']} "
          f"query genes have a qualifying hit (E<={args.evalue}, qcov>={args.qcov_threshold}), "
          f"{summary['n_distinct_targets']} distinct reference targets hit -- "
          f"verdict: {summary['verdict'].upper()}", file=sys.stderr)

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["query", "best_hit", "pident", "qcov", "evalue", "qualifies"])
        for q in query_ids:
            hit = best_hits.get(q)
            if hit is None:
                writer.writerow([q, "", "", "", "", "N"])
                continue
            qualifies = hit["evalue"] <= args.evalue and hit["qcov"] >= args.qcov_threshold
            writer.writerow([q, hit["sseqid"], f"{hit['pident']:.1f}", f"{hit['qcov']:.2f}",
                              f"{hit['evalue']:.2e}", "Y" if qualifies else "N"])
        writer.writerow([])
        writer.writerow(["VERDICT", summary["verdict"],
                          f"{summary['n_with_qualifying_hit']}/{summary['n_total']} qualify",
                          f"{summary['n_distinct_targets']} distinct targets"])


if __name__ == "__main__":
    main()
