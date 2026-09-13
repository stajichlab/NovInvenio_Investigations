#!/usr/bin/env python3
"""Score each clustering backend (mmseqs2 vs. diamond) against known
Starship positive controls (Tables S6/S21/S12/S13/S7/S14-16) plus
negative controls (conserved non-mobile SM clusters, random family
pairs) -- notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-
design.md, component 6.

Usage:
  benchmark_scorecard.py --crosswalk crosswalk.tsv \\
      --ground_truth_starships starship_ground_truth_starships.tsv \\
      --ground_truth_cargo starship_ground_truth_cargo.tsv \\
      --matrix_mmseqs presence_matrix.mmseqs.rescued.tsv \\
      --matrix_diamond presence_matrix.diamond.rescued.tsv \\
      --output benchmark_scorecard.tsv
"""
from __future__ import annotations

import argparse


def score_presence_recovery(predicted: dict[str, bool], truth: dict[str, bool]) -> dict:
    tp = fp = tn = fn = 0
    for key, truth_val in truth.items():
        pred_val = predicted.get(key, False)
        if truth_val and pred_val:
            tp += 1
        elif truth_val and not pred_val:
            fn += 1
        elif not truth_val and pred_val:
            fp += 1
        else:
            tn += 1
    total = tp + fp + tn + fn
    accuracy = (tp + tn) / total if total else 0.0
    union = tp + fp + fn
    jaccard = tp / union if union else 0.0
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn, "accuracy": accuracy, "jaccard": jaccard}


def score_cargo_grouping(
    predicted_family_of_gene: dict[str, str], truth_cargo_sets: dict[str, set[str]]
) -> dict[str, dict]:
    """For each true cargo set, find the predicted family containing the
    largest overlap with it, and report purity (fraction of that
    family's members that are true cargo members) and completeness
    (fraction of true cargo members recovered in that family)."""
    result = {}
    for starship_id, true_genes in truth_cargo_sets.items():
        family_hits: dict[str, set[str]] = {}
        for gene in true_genes:
            fam = predicted_family_of_gene.get(gene)
            if fam is not None:
                family_hits.setdefault(fam, set()).add(gene)
        if not family_hits:
            result[starship_id] = {"purity": 0.0, "completeness": 0.0}
            continue
        best_family = max(family_hits, key=lambda f: len(family_hits[f]))
        family_members = {
            g for g, f in predicted_family_of_gene.items() if f == best_family
        }
        overlap = family_hits[best_family]
        purity = len(overlap) / len(family_members) if family_members else 0.0
        completeness = len(overlap) / len(true_genes) if true_genes else 0.0
        result[starship_id] = {"purity": purity, "completeness": completeness}
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--crosswalk", required=True)
    ap.add_argument("--ground_truth_starships", required=True)
    ap.add_argument("--ground_truth_cargo", required=True)
    ap.add_argument("--matrix_mmseqs", required=True)
    ap.add_argument("--matrix_diamond", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    # The real run wires PresenceMatrix.from_tsv() for each backend, the
    # crosswalk TSV, and Task 9's ground-truth TSVs into calls to
    # score_presence_recovery/score_cargo_grouping per control Starship;
    # see Task 10 Step 5 in the plan for the exact sequence, since it
    # depends on the crosswalk actually being populated first (Task 9).
    print("See plan Task 10 Step 5 for the real-data scorecard run.", )


if __name__ == "__main__":
    main()
