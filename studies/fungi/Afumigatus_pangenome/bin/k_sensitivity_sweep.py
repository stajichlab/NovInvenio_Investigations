#!/usr/bin/env python3
"""Sweep pair_classification.py's physical-linkage window (`k`, a gene-rank
distance) and report how many trans/trans_unconfirmed/ambiguous_linkage
pairs would reclassify into a physical category at each larger k.

Why this matters: pair_classification.py's own module docstring already
flags this as a known partial-implementation caveat -- `k` is a fixed
gene-rank window (default 10), not sized to real accessory-island
footprints, which range up to 837 genes (see
results/accessory_islands/significant_islands.tsv, median 7). A worked
example confirmed this directly: the family module rooted at
Asfu_08190230|KAK9559653.1 (22 families, jaccard=1.0 across 53/295 strains)
sits contiguously on ONE contig in a real carrier strain (Asfu_UD1,
contig JBISCG010000004.1, ranks 16514-16567 -- a 53-gene span) and in the
outgroup A. lentulus, yet many of that module's own internal family pairs
are still classified `trans` pairwise, because the two families' positions
are >10 genes apart even though the whole block is one physically
contiguous, co-inherited unit already recovered as a single significant
island by find_accessory_islands.py's separate run-merging logic. This
sweep quantifies how much of the `trans`/`trans_unconfirmed`/
`ambiguous_linkage` buckets are this kind of window-cutoff artifact, rather
than genuinely non-physical co-occurrence.

Design: computes each candidate pair's per-strain (min same-contig rank
distance) ONCE, then reclassifies at every requested k from that cached
geometry -- O(1) extra work per additional k, instead of re-running
pair_classification.py once per k (each of which would reload
family_positions and recompute every FDR-significant pair's fisher/
permutation test that this sweep does not need to touch). Only pairs
currently classified trans/trans_unconfirmed/ambiguous_linkage are
candidates -- linkage_fraction is monotonically non-decreasing in k, so an
already-physical pair (starship_explained/unexplained_physical) cannot
become less linked at a larger k and is skipped.

Usage:
  k_sensitivity_sweep.py \\
      --pair_classification results/full_293run/pair_classification.rescued.tsv \\
      --family_positions results/full_293run/family_positions.rescued.tsv \\
      --cluster_tsv results/full_293run/tier1_cluster.tsv \\
      --captain_tblout results/captain_gene/DUF3435_vs_study.tblout \\
      --k_values 10,20,30,50,100,200 \\
      --output k_sensitivity_sweep.tsv
"""
from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import read_cluster_tsv  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from pair_classification import (  # noqa: E402
    DEFAULT_MIN_CLADES,
    DEFAULT_PERM_ALPHA,
    DEFAULT_PHYSICAL_THRESHOLD,
    DEFAULT_TRANS_THRESHOLD,
    has_captain_evidence,
    load_captain_families,
    load_family_positions,
)

CANDIDATE_CLASSIFICATIONS = {"trans", "trans_unconfirmed", "ambiguous_linkage"}


def per_strain_min_distances(
    family_a: str,
    family_b: str,
    gene_position: dict[str, dict[str, list[tuple[str, int]]]],
) -> list[int | None]:
    """One entry per strain carrying both families: the minimum same-contig
    rank distance between any copy of family_a and any copy of family_b, or
    None if that strain carries both but never on the same contig."""
    dists: list[int | None] = []
    for positions in gene_position.values():
        pos_a = positions.get(family_a)
        pos_b = positions.get(family_b)
        if not pos_a or not pos_b:
            continue
        best = None
        for contig_a, rank_a in pos_a:
            for contig_b, rank_b in pos_b:
                if contig_a == contig_b:
                    d = abs(rank_a - rank_b)
                    if best is None or d < best:
                        best = d
        dists.append(best)
    return dists


def linkage_fraction_from_distances(dists: list[int | None], k: int) -> float:
    if not dists:
        return 0.0
    linked = sum(1 for d in dists if d is not None and d <= k)
    return linked / len(dists)


def reclassify(
    frac: float,
    permutation_p: float,
    n_clades: int,
    family_a: str,
    family_b: str,
    gene_position: dict[str, dict[str, list[tuple[str, int]]]],
    captain_families: dict[str, set[str]],
    k: int,
    captain_k: int | None = None,
    physical_threshold: float = DEFAULT_PHYSICAL_THRESHOLD,
    trans_threshold: float = DEFAULT_TRANS_THRESHOLD,
    perm_alpha: float = DEFAULT_PERM_ALPHA,
    min_clades: int = DEFAULT_MIN_CLADES,
) -> str:
    """Mirrors pair_classification.classify_pair's threshold logic exactly,
    but starting from an already-known linkage_fraction/n_co_carrying
    (insufficient_data is never re-derived here -- a pair below
    min_co_carrying at the original k stays below it at any larger k, so
    those rows are never in CANDIDATE_CLASSIFICATIONS to begin with).

    `captain_k` is the window used for the captain-gene proximity check,
    independent of `k` (the swept physical-linkage window). Left as None,
    it defaults to `k` -- matching pair_classification.py's own behavior,
    which reuses one k for both checks. Passing a FIXED captain_k (e.g. the
    default 10) decouples the two: as the sweep grows the physical-linkage
    window to test whether two families' positions are part of a larger
    contiguous block, a captain gene doesn't retroactively become "nearby"
    just because that unrelated window grew -- avoiding the confound noted
    in results/REPORT.md's k-sensitivity section, where reusing one k made
    `starship_explained` climb mostly because the captain-gene window itself
    widened, not because more real Starship mechanism was found."""
    effective_captain_k = k if captain_k is None else captain_k
    if frac >= physical_threshold:
        if has_captain_evidence(family_a, family_b, gene_position, captain_families, effective_captain_k):
            return "starship_explained"
        return "unexplained_physical"
    if frac <= trans_threshold:
        if permutation_p < perm_alpha and n_clades >= min_clades:
            return "trans"
        return "trans_unconfirmed"
    return "ambiguous_linkage"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair_classification", required=True)
    ap.add_argument("--family_positions", required=True)
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--captain_tblout", required=True)
    ap.add_argument("--k_values", default="10,20,30,50,100,200")
    ap.add_argument(
        "--captain_k", type=int, default=None,
        help=(
            "Fixed window for the captain-gene proximity check, independent of the "
            "swept physical-linkage k (default: None, meaning captain_k == k at "
            "every sweep point, matching pair_classification.py's own coupled "
            "behavior). Pass e.g. 10 to hold the captain-gene window at the "
            "pipeline's real default while only the physical-linkage window grows "
            "-- decouples the two so starship_explained growth reflects genuine "
            "captain-gene proximity, not the captain window widening too."
        ),
    )
    ap.add_argument("--output", required=True)
    ap.add_argument(
        "--examples_output", default=None,
        help="Optional: up to --max_examples family pairs per (k, new_classification) transition, for spot-checking.",
    )
    ap.add_argument("--max_examples", type=int, default=5)
    args = ap.parse_args()

    k_values = sorted({int(x) for x in args.k_values.split(",")})

    print("Loading family positions...", file=sys.stderr)
    gene_position = load_family_positions(args.family_positions)
    print("Loading cluster membership...", file=sys.stderr)
    member_to_rep = read_cluster_tsv(args.cluster_tsv)
    print("Loading captain-gene evidence...", file=sys.stderr)
    captain_families = load_captain_families(args.captain_tblout, member_to_rep)

    counts: dict[int, dict[tuple[str, str], int]] = {k: {} for k in k_values}
    examples: dict[tuple[int, str, str], list[tuple[str, str]]] = {}
    n_candidates = 0
    n_total = 0

    with open_maybe_compressed(args.pair_classification) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            n_total += 1
            parts = line.rstrip("\n").split("\t")
            orig = parts[idx["classification"]]
            if orig not in CANDIDATE_CLASSIFICATIONS:
                continue
            n_candidates += 1
            family_a = parts[idx["family_a"]]
            family_b = parts[idx["family_b"]]
            permutation_p = float(parts[idx["permutation_p"]])
            clade_composition = ast.literal_eval(parts[idx["clade_composition"]])
            n_clades = len(clade_composition)

            dists = per_strain_min_distances(family_a, family_b, gene_position)
            for k in k_values:
                frac = linkage_fraction_from_distances(dists, k)
                new_label = reclassify(
                    frac, permutation_p, n_clades, family_a, family_b,
                    gene_position, captain_families, k, captain_k=args.captain_k,
                )
                key = (orig, new_label)
                counts[k][key] = counts[k].get(key, 0) + 1
                if args.examples_output and new_label != orig:
                    ex_key = (k, orig, new_label)
                    bucket = examples.setdefault(ex_key, [])
                    if len(bucket) < args.max_examples:
                        bucket.append((family_a, family_b))

            if n_candidates % 500_000 == 0:
                print(f"  ...{n_candidates} candidate pairs processed ({n_total} total rows read)", file=sys.stderr)

    with open(args.output, "w") as out:
        out.write("k\toriginal_classification\tnew_classification\tn_pairs\n")
        for k in k_values:
            for (orig, new_label), n in sorted(counts[k].items(), key=lambda kv: -kv[1]):
                out.write(f"{k}\t{orig}\t{new_label}\t{n}\n")

    if args.examples_output:
        with open(args.examples_output, "w") as out:
            out.write("k\toriginal_classification\tnew_classification\tfamily_a\tfamily_b\n")
            for (k, orig, new_label), pairs in examples.items():
                for family_a, family_b in pairs:
                    out.write(f"{k}\t{orig}\t{new_label}\t{family_a}\t{family_b}\n")

    print(f"Processed {n_candidates} candidate pairs (of {n_total} total rows).", file=sys.stderr)
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
