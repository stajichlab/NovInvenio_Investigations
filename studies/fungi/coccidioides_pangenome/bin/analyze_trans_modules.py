#!/usr/bin/env python3
"""Cross-reference Leiden trans-module family assignments with their Pfam
domain hits: per-module domain summary, largest modules, and a targeted
search for NACHT/HET (heterokaryon incompatibility) and secondary
metabolite biosynthesis (PKS/NRPS) domains among chained families.

Usage:
  analyze_trans_modules.py --family_modules family_modules_r2.0.tsv \\
      --domtblout pfam.domtblout --output_prefix trans_modules/r2.0
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict

NACHT_HET_KEYWORDS = ("NACHT", "HET", "Het-", "Ankyrin", "Ank_")
SM_KEYWORDS = (
    "PKS", "Ketoacyl-synt", "AMP-binding", "Condensation", "PP-binding",
    "Thioesterase", "Epimerase", "NRPS", "Polyketide", "TE_",
)


def parse_domtblout(path: str) -> dict[str, set[tuple[str, str]]]:
    """{query_protein: {(domain_name, accession), ...}}"""
    hits: dict[str, set[tuple[str, str]]] = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            domain_name, accession, query = parts[0], parts[1], parts[3]
            hits[query].add((domain_name, accession))
    return hits


def load_family_modules(path: str) -> dict[str, tuple[str, int]]:
    """{family: (module_id, module_size)}"""
    out = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            family, module_id, module_size = line.rstrip("\n").split("\t")
            out[family] = (module_id, int(module_size))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family_modules", required=True)
    ap.add_argument("--domtblout", required=True)
    ap.add_argument("--output_prefix", required=True)
    ap.add_argument("--min_module_size", type=int, default=2)
    args = ap.parse_args()

    family_module = load_family_modules(args.family_modules)
    domain_hits = parse_domtblout(args.domtblout)
    print(f"{len(family_module)} families in modules, {len(domain_hits)} proteins with >=1 Pfam hit",
          file=sys.stderr)

    module_families: dict[str, list[str]] = defaultdict(list)
    for family, (module_id, size) in family_module.items():
        if size >= args.min_module_size:
            module_families[module_id].append(family)

    # Per-module domain summary
    with open(f"{args.output_prefix}.module_domains.tsv", "w") as out:
        out.write("module_id\tmodule_size\tn_families_with_domain\tdomain_names\n")
        rows = []
        for module_id, families in module_families.items():
            domains: Counter = Counter()
            n_with = 0
            for fam in families:
                fam_domains = {d for d, _ in domain_hits.get(fam, set())}
                if fam_domains:
                    n_with += 1
                domains.update(fam_domains)
            top_domains = ",".join(f"{d}({c})" for d, c in domains.most_common(15))
            rows.append((module_id, len(families), n_with, top_domains))
        for module_id, size, n_with, top_domains in sorted(rows, key=lambda r: -r[1]):
            out.write(f"{module_id}\t{size}\t{n_with}\t{top_domains or '-'}\n")

    # Targeted keyword search: NACHT/HET and secondary-metabolite domains
    def keyword_search(keywords: tuple[str, ...]) -> list[tuple[str, str, str, str]]:
        rows = []
        for family, (module_id, size) in family_module.items():
            for domain_name, accession in domain_hits.get(family, set()):
                if any(k.lower() in domain_name.lower() for k in keywords):
                    rows.append((module_id, size, family, f"{domain_name}({accession})"))
        return sorted(rows, key=lambda r: -r[1])

    with open(f"{args.output_prefix}.nacht_het_hits.tsv", "w") as out:
        out.write("module_id\tmodule_size\tfamily\tdomain\n")
        for row in keyword_search(NACHT_HET_KEYWORDS):
            out.write("\t".join(map(str, row)) + "\n")

    with open(f"{args.output_prefix}.secondary_metabolite_hits.tsv", "w") as out:
        out.write("module_id\tmodule_size\tfamily\tdomain\n")
        for row in keyword_search(SM_KEYWORDS):
            out.write("\t".join(map(str, row)) + "\n")

    n_nacht_het = len(keyword_search(NACHT_HET_KEYWORDS))
    n_sm = len(keyword_search(SM_KEYWORDS))
    print(f"NACHT/HET-like hits: {n_nacht_het} family rows; secondary-metabolite-like hits: {n_sm} family rows",
          file=sys.stderr)


if __name__ == "__main__":
    main()
