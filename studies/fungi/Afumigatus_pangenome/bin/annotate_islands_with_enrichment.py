#!/usr/bin/env python3
"""Join significant_islands.with_domains.tsv (per-island Pfam domain
lists, from summarize_island_functions.py) with domain_enrichment.tsv
(per-domain FDR q-values, from the same script's enrichment test) into
one file that answers directly "which islands carry a significantly
enriched domain, and at what FDR" -- neither source file has both
pieces together (the per-island file has domain NAMES only, no q-values;
the per-domain file has q-values but no island membership).

Usage:
  annotate_islands_with_enrichment.py \\
      --islands significant_islands.with_domains.tsv \\
      --domain_enrichment domain_enrichment.tsv \\
      --fdr_alpha 0.05 \\
      --output significant_islands.with_enrichment.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys


def load_domain_qvalues(path: str) -> dict[str, float]:
    qvalues: dict[str, float] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            qvalues[row["domain"]] = float(row["fdr_q"])
    return qvalues


def annotate_islands(
    island_rows: list[dict], domain_qvalues: dict[str, float], fdr_alpha: float,
) -> list[dict]:
    """Adds `n_significant_domains`, `significant_domains` (name:q pairs,
    sorted by q ascending), and `min_fdr_q` (the single most significant
    domain's q-value, or empty if none) to each island row. Rows are
    returned in the SAME order they were given (callers control sorting;
    this function only annotates, matching the source file's existing
    island_size-descending order by default)."""
    annotated = []
    for row in island_rows:
        domains = [d for d in row.get("pfam_domains", "").split(",") if d and d != "-"]
        hits = sorted(
            ((d, domain_qvalues[d]) for d in domains if d in domain_qvalues and domain_qvalues[d] < fdr_alpha),
            key=lambda dq: dq[1],
        )
        new_row = dict(row)
        new_row["n_significant_domains"] = str(len(hits))
        new_row["significant_domains"] = ",".join(f"{d}:{q:.2e}" for d, q in hits)
        new_row["min_fdr_q"] = f"{hits[0][1]:.2e}" if hits else ""
        annotated.append(new_row)
    return annotated


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--islands", required=True)
    ap.add_argument("--domain_enrichment", required=True)
    ap.add_argument("--fdr_alpha", type=float, default=0.05)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    domain_qvalues = load_domain_qvalues(args.domain_enrichment)
    with open(args.islands, newline="") as fh:
        island_rows = list(csv.DictReader(fh, delimiter="\t"))

    annotated = annotate_islands(island_rows, domain_qvalues, args.fdr_alpha)
    n_with_hit = sum(1 for r in annotated if int(r["n_significant_domains"]) > 0)
    print(
        f"annotate_islands_with_enrichment: {n_with_hit}/{len(annotated)} islands carry "
        f">=1 domain significant at FDR<{args.fdr_alpha}", file=sys.stderr,
    )

    fieldnames = list(annotated[0].keys()) if annotated else []
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(annotated)


if __name__ == "__main__":
    main()
