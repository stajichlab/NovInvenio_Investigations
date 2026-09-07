#!/usr/bin/env python3
"""Pfam/InterPro over-representation analysis (ORA) for a candidate protein set.

ORA only, hypergeometric test -- no GSEA (DESIGN.md Sec 7/10: candidates are a discrete
set-membership call, not a continuous ranking). This is the custom equivalent of what
goatools does for GO (bin/go_enrichment.py); Pfam/InterPro have no GO-DAG to propagate
through, so this is a flat per-domain test, not per-vocabulary-hierarchy.

Addresses two points from the Fable-model design review (DESIGN.md Sec 10 addendum):
  - The ORA universe (background) is restricted to proteins carrying >=1 annotation in
    the tested vocabulary -- candidates/background are disproportionately unannotated
    by construction, and including unannotated proteins in the denominator biases
    every test toward "depleted."
  - Counting happens at the level given in --candidates/--background as-is: pass
    one row per gene FAMILY (not per individual candidate protein) if family-level
    counting is wanted -- this script doesn't cluster on its own, it tests exactly the
    membership you give it. See DESIGN.md Sec 10 on why per-protein counting inflates
    significance for expanded families.

Input: one or more --annotations TSVs from bin/extract_dat_annotations.py (accession,
taxon_id, gene_name, go_ids, pfam_ids, interpro_ids), a --candidates file (one
accession per line), and a --background file (one accession per line; every candidate
must also appear in background, since candidates are the "successes" drawn from it).
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

from scipy.stats import hypergeom

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from stats_util import bh_fdr  # noqa: E402


def load_annotations(paths: list[Path], field: str) -> dict[str, set[str]]:
    """accession -> set of domain IDs (Pfam or InterPro) from one or more TSVs."""
    out: dict[str, set[str]] = {}
    for p in paths:
        with open(p, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                ids = {x for x in row[field].split("|") if x}
                if ids:
                    out[row["accession"]] = ids
    return out


def load_ids(path: Path) -> set[str]:
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotations", nargs="+", required=True, type=Path, help="One or more extract_dat_annotations.py TSVs")
    ap.add_argument("--candidates", required=True, type=Path, help="File of candidate accessions, one per line")
    ap.add_argument("--background", required=True, type=Path, help="File of background accessions, one per line")
    ap.add_argument("--domain-field", choices=["pfam_ids", "interpro_ids"], required=True)
    ap.add_argument("--min-term-size", type=int, default=2, help="Skip domains present in fewer than this many background proteins (default: 2)")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    annot = load_annotations(args.annotations, args.domain_field)
    candidates = load_ids(args.candidates)
    background = load_ids(args.background)

    missing = candidates - background
    if missing:
        sys.exit(f"ERROR: {len(missing)} candidate accession(s) not in background (candidates must be a subset), e.g. {sorted(missing)[:3]}")

    # restrict to the annotatable universe (Fable review point: don't let
    # unannotated proteins bias the denominator toward "depleted")
    background_annot = background & annot.keys()
    candidates_annot = candidates & annot.keys()
    n_total = len(background_annot)
    n_candidates = len(candidates_annot)

    print(
        f"Universe: {len(background)} background ({n_total} annotatable, "
        f"{n_total/len(background):.0%}), {len(candidates)} candidates "
        f"({n_candidates} annotatable, {n_candidates/len(candidates) if candidates else 0:.0%})",
        file=sys.stderr,
    )

    domain_bg_count: dict[str, int] = defaultdict(int)
    domain_cand_count: dict[str, int] = defaultdict(int)
    for acc in background_annot:
        for d in annot[acc]:
            domain_bg_count[d] += 1
    for acc in candidates_annot:
        for d in annot[acc]:
            domain_cand_count[d] += 1

    rows = []
    for domain, n_bg in domain_bg_count.items():
        if n_bg < args.min_term_size:
            continue
        k = domain_cand_count.get(domain, 0)
        # one-sided over-representation: P(X >= k) given M=n_total, n=n_bg, N=n_candidates
        pval = hypergeom.sf(k - 1, n_total, n_bg, n_candidates) if k > 0 else 1.0
        expected = n_bg * n_candidates / n_total if n_total else 0
        rows.append({
            "domain": domain,
            "n_candidates_with_domain": k,
            "n_candidates_annotatable": n_candidates,
            "n_background_with_domain": n_bg,
            "n_background_annotatable": n_total,
            "expected": round(expected, 3),
            "fold_enrichment": round(k / expected, 3) if expected > 0 else float("nan"),
            "pvalue": pval,
        })

    fdrs = bh_fdr([r["pvalue"] for r in rows])
    for r, fdr in zip(rows, fdrs):
        r["fdr_bh"] = fdr
    rows.sort(key=lambda r: r["pvalue"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "domain", "n_candidates_with_domain", "n_candidates_annotatable",
            "n_background_with_domain", "n_background_annotatable",
            "expected", "fold_enrichment", "pvalue", "fdr_bh",
        ], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_sig = sum(1 for r in rows if r["fdr_bh"] < 0.05)
    print(f"Wrote {args.output} ({len(rows)} domains tested, {n_sig} at FDR<0.05)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
