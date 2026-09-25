#!/usr/bin/env python3
"""Functional (Pfam-domain) summary of the Leiden trans-network modules
found by cooccurrence.py's community-detection step
(results/full_293run/modules_preview/family_modules_r*.tsv), plus a
domain-enrichment test of non-singleton-module member families against
the same eligible background used for significant_islands' enrichment
test (shell+cloud bins, cooccurrence.py's own selection).

Direct mirror of summarize_island_functions.py, generalized from "island
member families" to "trans-module member families" -- imports that
script's parse_domtblout/load_eligible_background/domain_enrichment
functions UNCHANGED rather than reimplementing the same math, since the
enrichment test itself doesn't care whether a "group of co-varying
families" came from physical adjacency or trans co-occurrence.

Why this exists: a family's Pfam/SwissProt annotation coverage was
already complete for every family that ever appears in a significant
cooccurring pair (verified directly: 8,940 distinct families appear in
pair_classification.rescued.tsv, all 8,940 already covered by the
existing island_family_reps.fa + extra_background_reps.fa Pfam scans) --
so no new hmmscan/diamond run is needed. What was missing was simply an
analysis that USES that existing coverage at the trans-module level the
way summarize_island_functions.py already does at the island level. A
hand-checked worked example motivated this: the family module rooted at
Asfu_08190230|KAK9559653.1 (22 families, jaccard=1.0 across 53/295
strains) is itself a real, physically contiguous 53-gene genomic block
(confirmed against family_positions.rescued.tsv) enriched for Ankyrin/
NACHT/NPHP3_N domains (fungal NLR/heterokaryon-incompatibility genes) --
this script makes that kind of lookup available for every module, not
just ones found by hand.

Two outputs, matching summarize_island_functions.py's shape exactly so
the SAME downstream join scripts (annotate_islands_with_enrichment.py,
add_swissprot_to_islands.py -- both generic over any table with a
`pfam_domains` / `member_families` column) work unchanged on either:
  1. --output_modules: one row per non-singleton module (module_id,
     module_size, member_families, pfam_domains) -- module_size >=
     --min_module_size (default 2; a singleton module has no internal
     co-occurrence structure to test).
  2. --output_enrichment: one row per Pfam domain found in >=1
     non-singleton-module member family, one-sided Fisher's exact test
     (enrichment among ALL non-singleton-module members, pooled, vs. the
     eligible background), BH-FDR corrected.

Usage:
  summarize_trans_module_functions.py \\
      --family_modules results/full_293run/modules_preview/family_modules_r5.0.tsv \\
      --domtblout results/accessory_islands/island_family_reps_vs_pfam.domtblout \\
      --domtblout results/accessory_islands/extra_background_vs_pfam.domtblout \\
      --frequency_table results/full_293run/frequency_table.rescued.tsv \\
      --output_modules results/trans_modules/family_modules_r5.0.with_domains.tsv \\
      --output_enrichment results/trans_modules/trans_module_domain_enrichment.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from summarize_island_functions import (  # noqa: E402
    domain_enrichment,
    load_eligible_background,
    parse_domtblout,
)


def load_modules(path: str, min_module_size: int) -> dict[str, list[str]]:
    """{module_id: [family, ...]} for every module with >= min_module_size
    members (singleton modules, module_size == 1, are excluded by
    default -- there is no internal co-occurrence structure inside a
    module of one family)."""
    modules: dict[str, list[str]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if int(row["module_size"]) < min_module_size:
                continue
            modules.setdefault(row["module_id"], []).append(row["family"])
    return modules


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family_modules", required=True)
    ap.add_argument("--domtblout", required=True, action="append")
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--min_module_size", type=int, default=2)
    ap.add_argument("--output_modules", required=True)
    ap.add_argument("--output_enrichment", required=True)
    args = ap.parse_args()

    family_domains = parse_domtblout(args.domtblout)
    print(f"summarize_trans_module_functions: {len(family_domains)} families with >=1 Pfam domain hit",
          file=sys.stderr)

    background = load_eligible_background(args.frequency_table)
    print(f"summarize_trans_module_functions: {len(background)} eligible (shell+cloud) background families",
          file=sys.stderr)

    modules = load_modules(args.family_modules, args.min_module_size)
    print(f"summarize_trans_module_functions: {len(modules)} modules with >={args.min_module_size} members",
          file=sys.stderr)

    module_member_families: set[str] = set()
    for members in modules.values():
        module_member_families.update(members)

    with open(args.output_modules, "w", newline="") as out:
        fieldnames = ["module_id", "module_size", "member_families", "pfam_domains"]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for module_id, members in sorted(modules.items(), key=lambda kv: -len(kv[1])):
            domains: set[str] = set()
            for m in members:
                domains |= family_domains.get(m, set())
            writer.writerow({
                "module_id": module_id,
                "module_size": len(members),
                "member_families": ",".join(members),
                "pfam_domains": ",".join(sorted(domains)) if domains else "-",
            })

    enrichment = domain_enrichment(module_member_families, background, family_domains)
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

    print(f"summarize_trans_module_functions: {len(enrichment)} domains tested, "
          f"{sum(1 for r in enrichment if r['fdr_q'] < 0.05)} significant at FDR<0.05", file=sys.stderr)


if __name__ == "__main__":
    main()
