#!/usr/bin/env python3
"""Assemble a single Markdown summary-tables report from this study's
already-computed per-step outputs -- pangenome composition (frequency
bins), assembly completeness (BUSCO), the HAC/hacA targeted screen, and
co-occurring-pair classification. Pure aggregation: does not recompute
anything, just reads existing TSVs and tabulates them for a human to read
in one place instead of five separate files.

Each input is optional except --frequency_table -- pass only the ones you
have (e.g. run this against pre-rescue outputs today, rerun unchanged
against post-rescue outputs once available). --label stamps the report
with which run it summarizes, since this study has had both a known-wrong
pre-fix rescue result and a corrected one on the same day (2026-09-15) --
see PANGENOME_CLUSTER_PROFILE_NOTES.md.

Usage:
  build_summary_report.py --frequency_table frequency_table.rescued.tsv \\
      --busco_summary busco_completeness_summary.tsv \\
      --hac_screen hac_crosswalk/hac_reference_screen.tsv \\
      --pair_classification pair_classification.rescued.tsv \\
      --label "post-rescue, corrected (2026-09-15)" \\
      --output SUMMARY.md
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path


def read_tsv(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def band_counts(frequency_table: list[dict]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in frequency_table:
        counts[row["bin"]] = counts.get(row["bin"], 0) + 1
    return counts


def band_table_markdown(counts: dict[str, int]) -> str:
    order = ["core", "soft_core", "shell", "cloud", "singleton"]
    total = sum(counts.values())
    lines = ["| Bin | Families | % of total |", "|---|---:|---:|"]
    for bin_name in order:
        n = counts.get(bin_name, 0)
        pct = 100 * n / total if total else 0.0
        lines.append(f"| {bin_name} | {n:,} | {pct:.1f}% |")
    for bin_name, n in counts.items():
        if bin_name not in order:
            pct = 100 * n / total if total else 0.0
            lines.append(f"| {bin_name} | {n:,} | {pct:.1f}% |")
    lines.append(f"| **Total** | **{total:,}** | **100.0%** |")
    return "\n".join(lines)


def busco_table_markdown(busco_rows: list[dict]) -> str:
    complete = [float(r["Complete_pct"]) for r in busco_rows if r.get("Complete_pct")]
    contigs = [int(r["Contigs"]) for r in busco_rows if r.get("Contigs")]
    n_draft = sum(1 for c in contigs if c > 100)
    n_near_chromosome = sum(1 for c in contigs if c <= 20)
    lines = [
        "| Metric | Value |", "|---|---:|",
        f"| Strains | {len(busco_rows):,} |",
        f"| Complete BUSCO % (min / median / max) | "
        f"{min(complete):.1f} / {statistics.median(complete):.1f} / {max(complete):.1f} |",
        f"| Strains with >100 contigs (draft-assembly-consistent) | {n_draft:,} |",
        f"| Strains with <=20 contigs (near-chromosome-consistent) | {n_near_chromosome:,} |",
    ]
    return "\n".join(lines)


def hac_table_markdown(hac_rows: list[dict]) -> str:
    n = len(hac_rows)
    hac_a_present = sum(1 for r in hac_rows if r.get("hacA_present") == "Y")
    hrm_a_present = sum(1 for r in hac_rows if r.get("hrmA_present") == "Y")
    lines = [
        "| Locus | Present in | % of strains |", "|---|---:|---:|",
        f"| hacA (Afu3g04070, UPR transcription factor) | {hac_a_present}/{n} | "
        f"{100 * hac_a_present / n:.1f}% |",
        f"| hrmA ortholog (subtelomeric HAC, Starship-mobile) | {hrm_a_present}/{n} | "
        f"{100 * hrm_a_present / n:.1f}% |",
    ]
    return "\n".join(lines)


def pair_classification_table_markdown(pc_rows: list[dict]) -> str:
    counts: dict[str, int] = {}
    for row in pc_rows:
        label = row["classification"]
        counts[label] = counts.get(label, 0) + 1
    total = sum(counts.values())
    lines = ["| Classification | Pairs | % of total |", "|---|---:|---:|"]
    for label, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        pct = 100 * n / total if total else 0.0
        lines.append(f"| {label} | {n:,} | {pct:.1f}% |")
    lines.append(f"| **Total FDR-significant pairs** | **{total:,}** | **100.0%** |")
    return "\n".join(lines)


def build_report(
    frequency_table: list[dict],
    busco_rows: list[dict] | None,
    hac_rows: list[dict] | None,
    pair_classification_rows: list[dict] | None,
    label: str,
) -> str:
    sections = [f"# Pangenome summary -- {label}\n"]

    sections.append("## Pangenome composition\n")
    sections.append(band_table_markdown(band_counts(frequency_table)))

    if busco_rows is not None:
        sections.append("\n## Assembly completeness (BUSCO, genome mode)\n")
        sections.append(busco_table_markdown(busco_rows))

    if hac_rows is not None:
        sections.append("\n## HAC / hacA targeted screen\n")
        sections.append(hac_table_markdown(hac_rows))

    if pair_classification_rows is not None:
        sections.append("\n## Co-occurring pair classification\n")
        sections.append(pair_classification_table_markdown(pair_classification_rows))

    return "\n".join(sections) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--busco_summary")
    ap.add_argument("--hac_screen")
    ap.add_argument("--pair_classification")
    ap.add_argument("--label", default="")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    frequency_table = read_tsv(args.frequency_table)
    busco_rows = read_tsv(args.busco_summary) if args.busco_summary else None
    hac_rows = read_tsv(args.hac_screen) if args.hac_screen else None
    pair_classification_rows = read_tsv(args.pair_classification) if args.pair_classification else None

    report = build_report(frequency_table, busco_rows, hac_rows, pair_classification_rows, args.label)
    Path(args.output).write_text(report)
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
