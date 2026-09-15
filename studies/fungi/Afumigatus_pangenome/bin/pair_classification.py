#!/usr/bin/env python3
"""Classify each FDR-significant co-occurring family pair as physically
linked (Starship-explained or not) vs. trans (candidate non-physical
interaction/co-evolution) -- design spec component 8, addendum
2026-09-13, revised per the Opus review (must-fix M1-M5).

Labels (a pair gets exactly one):
  - insufficient_data:    fewer than `min_co_carrying` strains carry both
                           families with a resolvable genomic position.
  - starship_explained:   linkage_fraction >= `physical_threshold` AND at
                           least one co-carrying strain has a DUF3435
                           (Starship captain gene) family within `k` genes
                           of either family.
  - unexplained_physical: linkage_fraction >= `physical_threshold`, no
                           captain-gene evidence found -- a real candidate
                           for an unannotated/novel mobile element, not
                           discarded as "expected Starship."
  - ambiguous_linkage:    `trans_threshold` < linkage_fraction <
                           `physical_threshold` -- reported, not discarded.
  - trans:                linkage_fraction <= `trans_threshold` AND the
                           pair ALSO clears `permutation_p < perm_alpha`
                           AND spans >= `min_clades` distinct clade labels
                           (component 4's own "independent-contrasts proxy",
                           now actually enforced here, not just prose) --
                           must-fix M5.
  - trans_unconfirmed:    linkage_fraction <= `trans_threshold` but fails
                           the permutation/clade gate above -- physically
                           non-linked but not statistically robust enough
                           to call a confident interaction candidate. Not
                           one of the addendum's original six labels; added
                           here because a pair that is confidently NON-
                           physical but statistically shaky has nowhere
                           honest to go in that taxonomy otherwise (folding
                           it into "trans" would smuggle back exactly the
                           unguarded case M5 was written to close).

Known partial-implementation caveats (recorded, not hidden):
  - `linkage_fraction` here is rank-window only (component 5's original
    "within k genes" statistic) -- the fuller must-fix M3 fix (an
    additional same-accessory-island check and a real bp-distance window
    sized to an actual Starship footprint, not a gene-rank count) is NOT
    yet implemented. A physical block much larger than `k` genes that this
    function would still recognize via accessory-island membership is not
    yet caught by this script.
  - Starship evidence is DUF3435-hit-adjacency only; Table S5/S21
    paper-coordinate cross-referencing (available for only 254/293 strains,
    and not cross-walkable to this study's own assembly versions per the
    Asfu_A1163 case in the study notes) is not joined in here.
  - `taxon_scheme` restriction (must-fix M5's fix for TaxonGroup mixing
    three incompatible clustering schemes) is NOT applied -- `min_clades`
    counts distinct entries in `clade_composition` as given, which may mix
    Barber/DAPC/Mash-derived labels. Treat `trans` calls as provisional
    until that restriction is added.

Usage:
  pair_classification.py --cooccurring_pairs cooccurring_pairs.tsv \\
      --family_positions family_positions.tsv \\
      --cluster_tsv tier1_cluster.tsv \\
      --captain_tblout DUF3435_vs_study.tblout \\
      --output pair_classification.tsv
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
from synteny_windows import linkage_fraction  # noqa: E402

DEFAULT_K = 10
DEFAULT_PHYSICAL_THRESHOLD = 0.5
DEFAULT_TRANS_THRESHOLD = 0.05
DEFAULT_MIN_CO_CARRYING = 5
DEFAULT_PERM_ALPHA = 0.05
DEFAULT_MIN_CLADES = 2


def load_family_positions(path: str) -> dict[str, dict[str, list[tuple[str, int]]]]:
    gene_position: dict[str, dict[str, list[tuple[str, int]]]] = {}
    with open_maybe_compressed(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            short, family, contig, rank = parts[0], parts[1], parts[2], int(parts[3])
            gene_position.setdefault(short, {}).setdefault(family, []).append((contig, rank))
    return gene_position


def load_captain_families(
    tblout_path: str, member_to_rep: dict[str, str], id_sep: str = "|"
) -> dict[str, set[str]]:
    """{Short: {family_id, ...}} for every strain with at least one DUF3435
    (Starship captain gene) hit, resolved to that protein's OWN tier-1
    family -- every captain-hit protein was part of the same all-strains
    clustering input, so it always has a family assignment already."""
    captain_families: dict[str, set[str]] = {}
    with open_maybe_compressed(tblout_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            target = line.split()[0]
            if id_sep not in target:
                continue
            short, protein_id = target.split(id_sep, 1)
            member_id = f"{short}{id_sep}{protein_id}"
            family = member_to_rep.get(member_id)
            if family is None:
                continue
            captain_families.setdefault(short, set()).add(family)
    return captain_families


def has_captain_evidence(
    family_a: str,
    family_b: str,
    gene_position: dict[str, dict[str, list[tuple[str, int]]]],
    captain_families: dict[str, set[str]],
    k: int,
) -> bool:
    """True if, in ANY strain, a captain-gene family sits within k genes of
    family_a or family_b (reuses linkage_fraction's own proximity test by
    treating each captain family as one more family to check linkage
    against)."""
    for short, captains in captain_families.items():
        positions = gene_position.get(short, {})
        if family_a not in positions and family_b not in positions:
            continue
        for captain_family in captains:
            if captain_family not in positions:
                continue
            single_strain_positions = {short: positions}
            for fam in (family_a, family_b):
                if fam not in positions:
                    continue
                if linkage_fraction(fam, captain_family, single_strain_positions, k=k) > 0:
                    return True
    return False


def classify_pair(
    family_a: str,
    family_b: str,
    permutation_p: float,
    clade_composition: dict[str, int],
    gene_position: dict[str, dict[str, list[tuple[str, int]]]],
    captain_families: dict[str, set[str]],
    k: int = DEFAULT_K,
    physical_threshold: float = DEFAULT_PHYSICAL_THRESHOLD,
    trans_threshold: float = DEFAULT_TRANS_THRESHOLD,
    min_co_carrying: int = DEFAULT_MIN_CO_CARRYING,
    perm_alpha: float = DEFAULT_PERM_ALPHA,
    min_clades: int = DEFAULT_MIN_CLADES,
) -> tuple[str, float]:
    """Returns (classification, linkage_fraction)."""
    n_co_carrying = sum(
        1 for positions in gene_position.values()
        if positions.get(family_a) and positions.get(family_b)
    )
    if n_co_carrying < min_co_carrying:
        return "insufficient_data", 0.0

    frac = linkage_fraction(family_a, family_b, gene_position, k=k)

    if frac >= physical_threshold:
        if has_captain_evidence(family_a, family_b, gene_position, captain_families, k):
            return "starship_explained", frac
        return "unexplained_physical", frac

    if frac <= trans_threshold:
        n_clades = len(clade_composition)
        if permutation_p < perm_alpha and n_clades >= min_clades:
            return "trans", frac
        return "trans_unconfirmed", frac

    return "ambiguous_linkage", frac


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cooccurring_pairs", required=True)
    ap.add_argument("--family_positions", required=True)
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--captain_tblout", required=True)
    ap.add_argument("--k", type=int, default=DEFAULT_K)
    ap.add_argument("--physical_threshold", type=float, default=DEFAULT_PHYSICAL_THRESHOLD)
    ap.add_argument("--trans_threshold", type=float, default=DEFAULT_TRANS_THRESHOLD)
    ap.add_argument("--min_co_carrying", type=int, default=DEFAULT_MIN_CO_CARRYING)
    ap.add_argument("--perm_alpha", type=float, default=DEFAULT_PERM_ALPHA)
    ap.add_argument("--min_clades", type=int, default=DEFAULT_MIN_CLADES)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    print("Loading family positions...", file=sys.stderr)
    gene_position = load_family_positions(args.family_positions)
    print("Loading cluster membership...", file=sys.stderr)
    member_to_rep = read_cluster_tsv(args.cluster_tsv)
    print("Loading captain-gene evidence...", file=sys.stderr)
    captain_families = load_captain_families(args.captain_tblout, member_to_rep)
    print(
        f"{sum(len(v) for v in captain_families.values())} captain-gene family "
        f"assignments across {len(captain_families)} strains", file=sys.stderr,
    )

    with open_maybe_compressed(args.cooccurring_pairs) as fh, open(args.output, "w") as out:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        out.write(
            "family_a\tfamily_b\tclassification\tlinkage_fraction\tjaccard\t"
            "fisher_p\tfdr_q\tpermutation_p\tdirection_a\tclade_composition\n"
        )
        n = 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            family_a, family_b = parts[idx["family_a"]], parts[idx["family_b"]]
            permutation_p = float(parts[idx["permutation_p"]])
            clade_composition = ast.literal_eval(parts[idx["clade_composition"]])
            classification, frac = classify_pair(
                family_a, family_b, permutation_p, clade_composition,
                gene_position, captain_families,
                k=args.k, physical_threshold=args.physical_threshold,
                trans_threshold=args.trans_threshold,
                min_co_carrying=args.min_co_carrying,
                perm_alpha=args.perm_alpha, min_clades=args.min_clades,
            )
            out.write(
                f"{family_a}\t{family_b}\t{classification}\t{frac:.4f}\t"
                f"{parts[idx['jaccard']]}\t{parts[idx['fisher_p']]}\t"
                f"{parts[idx['fdr_q']]}\t{parts[idx['permutation_p']]}\t"
                f"{parts[idx['direction_a']]}\t{parts[idx['clade_composition']]}\n"
            )
            n += 1
            if n % 100_000 == 0:
                print(f"pair_classification: classified {n} pairs", file=sys.stderr)

    print(f"pair_classification: {n} pairs classified", file=sys.stderr)


if __name__ == "__main__":
    main()
