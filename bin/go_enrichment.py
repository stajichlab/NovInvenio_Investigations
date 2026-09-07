#!/usr/bin/env python3
"""GO-term over-representation analysis (ORA) for a candidate protein set, via goatools.

ORA only (one-sided, over-represented terms), hypergeometric + BH-FDR, no GSEA --
same rationale as bin/domain_enrichment.py (DESIGN.md Sec 7/10). goatools supplies
GO-DAG true-path propagation (--propagate-counts, on by default) given go-basic.obo
(bin/fetch_go_obo.py), which the Pfam/InterPro test has no equivalent of.

Evidence-code choice (Fable-model review, DESIGN.md Sec 10 addendum): TrEMBL `DR GO`
annotations are mostly IEA (electronic, InterPro2GO/UniRule-derived) -- GO is largely
*derived from* Pfam/InterPro for exactly the less-studied species this pipeline cares
about. Excluding IEA would leave species like Ztri/Cimm with near-zero GO evidence.
Default here is to INCLUDE IEA; --exclude-iea is available but changes results a lot
for exactly the species that need GO evidence most -- state which was used in any
report built from this output.

Same annotatable-universe restriction as domain_enrichment.py: background/candidates
are restricted to proteins carrying >=1 GO term before testing.
"""
import argparse
import csv
import sys
from pathlib import Path

from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from stats_util import bh_fdr  # noqa: E402
from uniprot_ids import bare_accession  # noqa: E402


def load_go_assoc(paths: list[Path], exclude_iea: bool) -> dict[str, set[str]]:
    assoc: dict[str, set[str]] = {}
    for p in paths:
        with open(p, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                terms = set()
                for entry in row["go_ids"].split("|"):
                    if not entry:
                        continue
                    # entry format is "GO:0016020:IEA" -- rpartition on the last colon,
                    # since GO IDs themselves contain a colon ("GO:0016020")
                    go_id, _, evidence = entry.rpartition(":")
                    if not go_id:
                        go_id, evidence = entry, ""
                    if exclude_iea and evidence == "IEA":
                        continue
                    terms.add(go_id)
                if terms:
                    assoc[row["accession"]] = terms
    return assoc


def load_ids(path: Path) -> set[str]:
    with open(path) as fh:
        return {bare_accession(line.strip()) for line in fh if line.strip()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotations", nargs="+", required=True, type=Path, help="One or more extract_dat_annotations.py TSVs")
    ap.add_argument("--obo", required=True, type=Path, help="go-basic.obo (bin/fetch_go_obo.py)")
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--background", required=True, type=Path)
    ap.add_argument("--exclude-iea", action="store_true", help="Exclude IEA-evidence GO terms (see module docstring -- changes results a lot for less-studied species)")
    ap.add_argument("--no-propagate", action="store_true", help="Disable GO-DAG true-path propagation (default: propagate)")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    print(f"Loading GO DAG from {args.obo}...", file=sys.stderr)
    obodag = GODag(str(args.obo))

    assoc = load_go_assoc(args.annotations, args.exclude_iea)
    candidates = load_ids(args.candidates)
    background = load_ids(args.background)

    missing = candidates - background
    if missing:
        sys.exit(f"ERROR: {len(missing)} candidate accession(s) not in background, e.g. {sorted(missing)[:3]}")

    background_annot = background & assoc.keys()
    candidates_annot = candidates & assoc.keys()
    print(
        f"Universe: {len(background)} background ({len(background_annot)} with >=1 GO term), "
        f"{len(candidates)} candidates ({len(candidates_annot)} with >=1 GO term). "
        f"IEA evidence {'EXCLUDED' if args.exclude_iea else 'included'}.",
        file=sys.stderr,
    )

    # methods=[] -- goatools' own "fdr_bh" method requires statsmodels'
    # sandbox.stats.multicomp.multipletests, removed in current statsmodels releases.
    # Apply BH-FDR ourselves instead (lib/stats_util.bh_fdr), the same implementation
    # domain_enrichment.py uses, so GO/Pfam/InterPro results are corrected identically.
    study = GOEnrichmentStudy(
        list(background_annot),
        assoc,
        obodag,
        propagate_counts=not args.no_propagate,
        alpha=0.05,
        methods=[],
    )
    results = study.run_study(list(candidates_annot))

    # ORA-only: keep over-represented terms only (enrichment == 'e'), matching the
    # one-sided decision in DESIGN.md Sec 7/10
    kept = [r for r in results if r.enrichment == "e"]
    fdrs = bh_fdr([r.p_uncorrected for r in kept])
    rows = []
    for r, fdr in zip(kept, fdrs):
        rows.append({
            "GO_id": r.GO,
            "namespace": r.NS,
            "name": r.name,
            "n_candidates_with_term": r.study_count,
            "n_candidates_annotatable": r.study_n,
            "n_background_with_term": r.pop_count,
            "n_background_annotatable": r.pop_n,
            "pvalue": r.p_uncorrected,
            "fdr_bh": fdr,
        })
    rows.sort(key=lambda r: r["pvalue"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "GO_id", "namespace", "name", "n_candidates_with_term",
            "n_candidates_annotatable", "n_background_with_term",
            "n_background_annotatable", "pvalue", "fdr_bh",
        ], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_sig = sum(1 for r in rows if r["fdr_bh"] < 0.05)
    print(f"Wrote {args.output} ({len(rows)} over-represented terms tested, {n_sig} at FDR<0.05)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
