#!/usr/bin/env python3
"""Score each clustering backend (mmseqs2 vs. diamond) against known
Starship positive controls (Tables S6/S21/S13, via
bin/build_ground_truth_tables.py) plus a negative control (Table S19's
conserved, non-Starship secondary-metabolite/virulence gene catalog) --
notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md,
component 6.

Usage:
  benchmark_scorecard.py --crosswalk crosswalk.tsv \\
      --ground_truth_starships ground_truth_starships_by_short.tsv \\
      --ground_truth_cargo ground_truth_cargo_by_nameid.tsv \\
      --negative_control_genes ground_truth_negative_control_conserved_genes.tsv \\
      --tier1_cluster_mmseqs tier1_cluster.tsv \\
      --matrix_mmseqs presence_matrix.rescued.tsv \\
      [--tier1_cluster_diamond tier1_diamond_cluster.tsv] \\
      [--matrix_diamond presence_matrix.diamond.rescued.tsv] \\
      --output benchmark_scorecard.tsv

The two `--tier1_cluster_*`/`--matrix_diamond` diamond arguments are
optional: if this study has not (yet) run a diamond-backend clustering at
full scale, pass only the mmseqs arguments and the scorecard is written for
mmseqs alone, with a note in its header rather than a fabricated diamond
comparison.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix, read_cluster_tsv  # noqa: E402
from frequency_bins import assign_bin  # noqa: E402


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


def load_crosswalk(path: str | Path) -> dict[str, str]:
    """{paper_id: study_protein_id} from id_crosswalk.py's output TSV."""
    crosswalk: dict[str, str] = {}
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                crosswalk[parts[0]] = parts[1]
    return crosswalk


def load_cargo_by_name(path: str | Path) -> dict[str, set[str]]:
    """{nameID: {geneID, ...}} from build_ground_truth_tables.py's
    ground_truth_cargo_by_nameid.tsv."""
    cargo: dict[str, set[str]] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cargo.setdefault(row["nameID"], set()).add(row["geneID"])
    return cargo


def load_starships_by_short(path: str | Path) -> dict[str, dict]:
    """{nameID: {"population_freq": float, "presence": {Short: bool}}} from
    build_ground_truth_tables.py's ground_truth_starships_by_short.tsv."""
    starships: dict[str, dict] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            starships[row["nameID"]] = {
                "population_freq": float(row["population_freq"]),
                "presence": json.loads(row["presence_by_short_json"]),
            }
    return starships


def build_protein_to_family(cluster_tsv_path: str | Path) -> dict[str, str]:
    """member protein_id -> tier-1 family (cluster representative) id --
    read_cluster_tsv already returns exactly this shape (rep\\tmember lines
    parsed into {member: rep})."""
    return read_cluster_tsv(cluster_tsv_path)


def score_backend(
    backend_name: str,
    crosswalk: dict[str, str],
    protein_to_family: dict[str, str],
    ground_truth_cargo: dict[str, set[str]],
    ground_truth_starships: dict[str, dict],
    matrix: PresenceMatrix,
) -> list[dict]:
    predicted_family_of_gene = {
        paper_gene: protein_to_family[study_protein]
        for paper_gene, study_protein in crosswalk.items()
        if study_protein in protein_to_family
    }
    matrix_strains = set(matrix.strains)
    rows: list[dict] = []

    cargo_scores = score_cargo_grouping(predicted_family_of_gene, ground_truth_cargo)
    for name_id, s in cargo_scores.items():
        n_crosswalked = sum(1 for g in ground_truth_cargo.get(name_id, ()) if g in predicted_family_of_gene)
        rows.append({
            "backend": backend_name,
            "control_type": "positive_cargo_grouping",
            "entity_id": name_id,
            "n_true_genes": len(ground_truth_cargo.get(name_id, ())),
            "n_crosswalked_genes": n_crosswalked,
            "purity": round(s["purity"], 4),
            "completeness": round(s["completeness"], 4),
        })

    for name_id, info in ground_truth_starships.items():
        cargo_genes = ground_truth_cargo.get(name_id, set())
        families = [predicted_family_of_gene[g] for g in cargo_genes if g in predicted_family_of_gene]
        if not families:
            rows.append({
                "backend": backend_name,
                "control_type": "positive_presence_recovery",
                "entity_id": name_id,
                "notes": "no crosswalked cargo gene -> cannot pick a diagnostic family",
            })
            continue
        diagnostic_family, n_votes = Counter(families).most_common(1)[0]
        truth = {s: v for s, v in info["presence"].items() if s in matrix_strains}
        predicted = {s: matrix.is_present(diagnostic_family, s) for s in truth}
        result = score_presence_recovery(predicted, truth)
        rows.append({
            "backend": backend_name,
            "control_type": "positive_presence_recovery",
            "entity_id": name_id,
            "population_freq_paper": info["population_freq"],
            "diagnostic_family": diagnostic_family,
            "diagnostic_family_votes": n_votes,
            "n_strains_scored": len(truth),
            **result,
        })

    return rows


def score_negative_control(
    backend_name: str,
    crosswalk: dict[str, str],
    protein_to_family: dict[str, str],
    negative_control_path: str | Path,
    matrix: PresenceMatrix,
) -> list[dict]:
    """Table S19's conserved, non-Starship gene catalog: expect HIGH
    recovered frequency (core/soft_core), not a Starship-style
    presence/absence split -- see build_ground_truth_tables.py's module
    docstring point 3 for why this is a defensible negative control."""
    rows: list[dict] = []
    with open(negative_control_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            accession = row["accession"]
            study_protein = crosswalk.get(accession)
            family = protein_to_family.get(study_protein) if study_protein else None
            if family is None:
                rows.append({
                    "backend": backend_name,
                    "control_type": "negative_control_conserved",
                    "entity_id": accession,
                    "notes": "not crosswalked" if study_protein is None else "crosswalked protein not in this backend's clusters",
                })
                continue
            freq = matrix.frequency(family)
            rows.append({
                "backend": backend_name,
                "control_type": "negative_control_conserved",
                "entity_id": accession,
                "gene": row["gene"],
                "cluster": row["cluster"],
                "diagnostic_family": family,
                "frequency": round(freq, 4),
                "freq_bin": assign_bin(freq, matrix.strain_count(family)),
            })
    return rows


FIELDNAMES = [
    "backend", "control_type", "entity_id", "gene", "cluster",
    "n_true_genes", "n_crosswalked_genes", "purity", "completeness",
    "population_freq_paper", "diagnostic_family", "diagnostic_family_votes",
    "n_strains_scored", "tp", "fp", "tn", "fn", "accuracy", "jaccard",
    "frequency", "freq_bin", "notes",
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--crosswalk", required=True)
    ap.add_argument("--ground_truth_starships", required=True)
    ap.add_argument("--ground_truth_cargo", required=True)
    ap.add_argument("--negative_control_genes", required=True)
    ap.add_argument("--tier1_cluster_mmseqs", required=True)
    ap.add_argument("--matrix_mmseqs", required=True)
    ap.add_argument("--tier1_cluster_diamond", default=None)
    ap.add_argument("--matrix_diamond", default=None)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    crosswalk = load_crosswalk(args.crosswalk)
    ground_truth_cargo = load_cargo_by_name(args.ground_truth_cargo)
    ground_truth_starships = load_starships_by_short(args.ground_truth_starships)

    all_rows: list[dict] = []

    print(f"Loading mmseqs tier-1 cluster membership from {args.tier1_cluster_mmseqs} ...", file=sys.stderr)
    mmseqs_protein_to_family = build_protein_to_family(args.tier1_cluster_mmseqs)
    print(f"Loading mmseqs presence matrix from {args.matrix_mmseqs} ...", file=sys.stderr)
    mmseqs_matrix = PresenceMatrix.from_tsv(args.matrix_mmseqs, read_copy_number=False)
    all_rows += score_backend(
        "mmseqs", crosswalk, mmseqs_protein_to_family, ground_truth_cargo, ground_truth_starships, mmseqs_matrix
    )
    all_rows += score_negative_control(
        "mmseqs", crosswalk, mmseqs_protein_to_family, args.negative_control_genes, mmseqs_matrix
    )

    ran_diamond = bool(args.tier1_cluster_diamond and args.matrix_diamond)
    if ran_diamond:
        print(f"Loading diamond tier-1 cluster membership from {args.tier1_cluster_diamond} ...", file=sys.stderr)
        diamond_protein_to_family = build_protein_to_family(args.tier1_cluster_diamond)
        print(f"Loading diamond presence matrix from {args.matrix_diamond} ...", file=sys.stderr)
        diamond_matrix = PresenceMatrix.from_tsv(args.matrix_diamond, read_copy_number=False)
        all_rows += score_backend(
            "diamond", crosswalk, diamond_protein_to_family, ground_truth_cargo, ground_truth_starships, diamond_matrix
        )
        all_rows += score_negative_control(
            "diamond", crosswalk, diamond_protein_to_family, args.negative_control_genes, diamond_matrix
        )
    else:
        print(
            "No --tier1_cluster_diamond/--matrix_diamond given -- scorecard covers "
            "mmseqs only, not a backend comparison.",
            file=sys.stderr,
        )

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    print(f"Wrote {len(all_rows)} scorecard rows to {args.output} "
          f"({'mmseqs + diamond' if ran_diamond else 'mmseqs only'})", file=sys.stderr)


if __name__ == "__main__":
    main()
