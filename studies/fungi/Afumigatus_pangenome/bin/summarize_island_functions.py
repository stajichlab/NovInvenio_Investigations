#!/usr/bin/env python3
"""Functional (Pfam-domain) summary of the significant accessory islands
found by find_accessory_islands.py, plus a domain-enrichment test of
island-member families against the CORRECT background: all families
actually eligible for co-occurrence testing (shell+cloud bins,
cooccurrence.py's own selection), never the whole genome -- core genes
were never eligible for testing in the first place, so a whole-genome
background would spuriously enrich for "accessory-genome-typical"
domains regardless of which specific island is being tested (the same
methodological point raised before this analysis was built).

Two outputs:
  1. --output_islands: significant_islands.tsv (or its sorted variant)
     with an added `pfam_domains` column -- the union of Pfam domains
     found among that island's member families, for browsing.
  2. --output_enrichment: one row per Pfam domain that appears in ANY
     island-member family, with a one-sided Fisher's exact test
     (enrichment among island-member families vs. the eligible
     background) and BH-FDR correction across all domains tested.

Usage:
  summarize_island_functions.py \\
      --significant_islands significant_islands.tsv \\
      --domtblout island_family_reps_vs_pfam.domtblout \\
      --domtblout extra_background_vs_pfam.domtblout \\
      --frequency_table frequency_table.rescued.tsv \\
      --output_islands significant_islands.with_domains.tsv \\
      --output_enrichment domain_enrichment.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact, false_discovery_control


def parse_domtblout(paths: list[str], max_ievalue: float = 1e-3) -> dict[str, set[str]]:
    """{family_id: {pfam_domain_name, ...}} from one or more hmmscan
    --domtblout files (query=protein/family rep, target=Pfam profile).
    Keeps a domain hit only if its independent E-value clears
    `max_ievalue` (redundant with hmmscan's own -E cutoff at the
    full-sequence level, but domtblout's domain-level i-Evalue can be
    looser than the sequence-level one, so this re-checks at the
    domain level explicitly)."""
    hits: dict[str, set[str]] = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 13:
                    continue
                target_name, query_name = parts[0], parts[3]
                try:
                    i_evalue = float(parts[12])
                except ValueError:
                    continue
                if i_evalue > max_ievalue:
                    continue
                hits.setdefault(query_name, set()).add(target_name)
    return hits


def load_eligible_background(frequency_table_path: str) -> set[str]:
    """All families in the shell or cloud bin -- cooccurrence.py's own
    eligibility selection (`row["bin"] in ("shell", "cloud")`), NOT
    singleton and NOT core. This is the correct enrichment background:
    the population of families that could possibly have been tested for
    co-occurrence in the first place."""
    background = set()
    with open(frequency_table_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["bin"] in ("shell", "cloud"):
                background.add(row["family"])
    return background


def domain_enrichment(
    island_member_families: set[str],
    background_families: set[str],
    family_domains: dict[str, set[str]],
) -> list[dict]:
    """One-sided Fisher's exact test per Pfam domain: is this domain more
    common among island-member families than among the full eligible
    background? `background_families` must be a SUPERSET of
    `island_member_families` (the treatment group is drawn from the same
    population as the background, not a disjoint comparison group) --
    the 2x2 table is (in-island vs. not-in-island) x (has-domain vs.
    not), both restricted to the eligible background population.
    BH-FDR corrected across every domain that appears in at least one
    island-member family (the actual number of hypotheses tested).

    `island_member_families` is intersected with `background_families`
    before any counting: `accessory_islands()` merges runs of ANY
    non-core gene (shell, cloud, OR singleton), but singleton families
    were never eligible for co-occurrence testing in the first place
    (cooccurrence.py only tests shell/cloud) -- a singleton island member
    can't sensibly count as either "in" or "out" of a population it was
    never part of, so it's excluded from the enrichment test entirely
    (not an error; a real, expected, and reported occurrence, not a data
    bug)."""
    excluded = island_member_families - background_families
    island_member_families = island_member_families & background_families
    if excluded:
        print(
            f"domain_enrichment: {len(excluded)} island-member families are outside "
            "the eligible (shell+cloud) background -- most likely singleton families "
            "an island happened to include; excluded from the enrichment test "
            f"({len(island_member_families)} remain)", file=sys.stderr,
        )

    all_domains: set[str] = set()
    for family in island_member_families:
        all_domains |= family_domains.get(family, set())

    n_background = len(background_families)
    n_island = len(island_member_families)
    rows = []
    pvalues = []
    for domain in sorted(all_domains):
        families_with_domain = {f for f in background_families if domain in family_domains.get(f, ())}
        both = len(families_with_domain & island_member_families)
        island_only = n_island - both
        domain_only = len(families_with_domain) - both
        neither = n_background - n_island - domain_only
        _, p = fisher_exact([[both, island_only], [domain_only, neither]], alternative="greater")
        rows.append({
            "domain": domain,
            "n_with_domain_in_islands": both,
            "n_with_domain_in_background": len(families_with_domain),
            "n_island_families": n_island,
            "n_background_families": n_background,
            "fisher_p": p,
        })
        pvalues.append(p)

    if pvalues:
        qvalues = false_discovery_control(np.asarray(pvalues), method="bh")
        for row, q in zip(rows, qvalues):
            row["fdr_q"] = float(q)
    rows.sort(key=lambda r: r["fisher_p"])
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--significant_islands", required=True)
    ap.add_argument("--domtblout", required=True, action="append")
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--output_islands", required=True)
    ap.add_argument("--output_enrichment", required=True)
    args = ap.parse_args()

    family_domains = parse_domtblout(args.domtblout)
    print(f"summarize_island_functions: {len(family_domains)} families with >=1 Pfam domain hit",
          file=sys.stderr)

    background = load_eligible_background(args.frequency_table)
    print(f"summarize_island_functions: {len(background)} eligible (shell+cloud) background families",
          file=sys.stderr)

    with open(args.significant_islands, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    island_member_families: set[str] = set()
    for row in rows:
        island_member_families.update(row["member_families"].split(","))

    with open(args.output_islands, "w", newline="") as out:
        fieldnames = list(rows[0].keys()) + ["pfam_domains"] if rows else []
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in sorted(rows, key=lambda r: -int(r["island_size"])):
            members = row["member_families"].split(",")
            domains: set[str] = set()
            for m in members:
                domains |= family_domains.get(m, set())
            row["pfam_domains"] = ",".join(sorted(domains)) if domains else "-"
            writer.writerow(row)

    enrichment = domain_enrichment(island_member_families, background, family_domains)
    with open(args.output_enrichment, "w") as out:
        out.write(
            "domain\tn_with_domain_in_islands\tn_with_domain_in_background\t"
            "n_island_families\tn_background_families\tfisher_p\tfdr_q\n"
        )
        for r in enrichment:
            out.write(
                f"{r['domain']}\t{r['n_with_domain_in_islands']}\t{r['n_with_domain_in_background']}\t"
                f"{r['n_island_families']}\t{r['n_background_families']}\t"
                f"{r['fisher_p']:.3e}\t{r['fdr_q']:.3e}\n"
            )

    print(f"summarize_island_functions: {len(enrichment)} domains tested, "
          f"{sum(1 for r in enrichment if r['fdr_q'] < 0.05)} significant at FDR<0.05",
          file=sys.stderr)


if __name__ == "__main__":
    main()
