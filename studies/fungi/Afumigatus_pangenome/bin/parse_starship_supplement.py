#!/usr/bin/env python3
"""Extract structured ground truth from mbio.01092-25-s0002.xlsx (Gluck-
Thaler et al. 2025, doi:10.1128/mbio.01092-25) for the benchmark
scorecard -- notes/superpowers/specs/2026-09-13-pangenome-cluster-
profile-design.md, component 6.

Usage:
  parse_starship_supplement.py --xlsx mbio.01092-25-s0002.xlsx \\
      --output_prefix starship_ground_truth
"""
from __future__ import annotations

import argparse

import pandas as pd


def high_confidence_starships(table_s6: pd.DataFrame, table_s21: pd.DataFrame) -> dict:
    """Combine Table S6 (nameID, freq) with Table S21 (isolateID, nameID,
    presence/absence) into {name_id: {population_freq, presence: {isolate: bool}}}."""
    result = {}
    freq_by_name = dict(zip(table_s6["nameID"], table_s6["freq"]))
    for name_id, freq in freq_by_name.items():
        subset = table_s21[table_s21["nameID"] == name_id]
        # Explicitly compare to 1 to avoid NaN→True coercion; NaN becomes False naturally
        presence = dict(zip(subset["isolateID"], subset["presence/absence"] == 1))
        result[name_id] = {"population_freq": freq, "presence": presence}
    return result


def cargo_gene_sets(table_s12_or_s13: pd.DataFrame) -> dict[str, set[str]]:
    """Group Table S12/S13's (starshipID, geneID) rows into
    {starship_id: {gene_id, ...}}."""
    result: dict[str, set[str]] = {}
    for starship_id, gene_id in zip(table_s12_or_s13["starshipID"], table_s12_or_s13["geneID"]):
        result.setdefault(starship_id, set()).add(gene_id)
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--output_prefix", required=True)
    args = ap.parse_args()

    xl = pd.ExcelFile(args.xlsx)
    table_s6 = xl.parse("Table S6", header=1)
    table_s21 = xl.parse("Table S21", header=1)
    table_s13 = xl.parse("Table S13", header=1)

    starships = high_confidence_starships(table_s6, table_s21)
    cargo = cargo_gene_sets(table_s13)

    with open(f"{args.output_prefix}_starships.tsv", "w") as fh:
        fh.write("nameID\tpopulation_freq\tpresence_json\n")
        import json
        for name_id, info in starships.items():
            fh.write(f"{name_id}\t{info['population_freq']}\t{json.dumps(info['presence'])}\n")

    with open(f"{args.output_prefix}_cargo.tsv", "w") as fh:
        fh.write("starshipID\tgeneID\n")
        for starship_id, genes in cargo.items():
            for gene in sorted(genes):
                fh.write(f"{starship_id}\t{gene}\n")


if __name__ == "__main__":
    main()
