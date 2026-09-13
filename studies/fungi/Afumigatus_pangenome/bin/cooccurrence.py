#!/usr/bin/env python3
"""Co-occurrence (co-loss/co-gain) statistics between shell+cloud gene
families -- notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-
design.md, component 4 (revised per the Fable bioinformatics review:
frequency floor + Fisher/BH-FDR instead of a raw Jaccard threshold,
outgroup-based gain/loss polarization, within-clade permutation null as a
cheap stand-in for a full phylogenetic correction).

Usage:
  cooccurrence.py --matrix presence_matrix.rescued.tsv \\
      --frequency_table frequency_table.tsv --config config.csv \\
      --output cooccurring_pairs.tsv
"""
from __future__ import annotations

import argparse
import itertools
import random
import sys
from pathlib import Path

from scipy.stats import fisher_exact, false_discovery_control

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix  # noqa: E402


def jaccard(a: list[bool], b: list[bool]) -> float:
    intersection = sum(1 for x, y in zip(a, b) if x and y)
    union = sum(1 for x, y in zip(a, b) if x or y)
    return intersection / union if union else 0.0


def fisher_pvalue(a: list[bool], b: list[bool]) -> float:
    both = sum(1 for x, y in zip(a, b) if x and y)
    a_only = sum(1 for x, y in zip(a, b) if x and not y)
    b_only = sum(1 for x, y in zip(a, b) if y and not x)
    neither = sum(1 for x, y in zip(a, b) if not x and not y)
    _, p = fisher_exact([[both, a_only], [b_only, neither]], alternative="greater")
    return p


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    if not pvalues:
        return []
    return list(false_discovery_control(pvalues, method="bh"))


def polarize_direction(present_in_outgroup: bool, freq_in_ingroup: float) -> str:
    """A family present in the outgroup that most ingroup strains also
    carry implies the strains lacking it lost it; a family absent from
    the outgroup implies the strains carrying it gained it."""
    if present_in_outgroup:
        return "loss"
    return "gain"


def clade_composition(strains_present: list[str], clade_of_strain: dict[str, str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for s in strains_present:
        clade = clade_of_strain.get(s, "unknown")
        counts[clade] = counts.get(clade, 0) + 1
    return counts


def permutation_null_pvalue(
    a: list[bool], b: list[bool], clades: list[str], n_perms: int, rng: random.Random
) -> float:
    """Shuffle `b` WITHIN each clade group (not globally) n_perms times,
    recompute Fisher's p-value each time, and return the fraction of
    permutations at least as extreme as the observed statistic -- a
    pair significant only because both families mark the same clade
    will look unremarkable under this null."""
    observed = fisher_pvalue(a, b)
    indices_by_clade: dict[str, list[int]] = {}
    for i, clade in enumerate(clades):
        indices_by_clade.setdefault(clade, []).append(i)

    at_least_as_extreme = 0
    for _ in range(n_perms):
        b_perm = list(b)
        for indices in indices_by_clade.values():
            values = [b_perm[i] for i in indices]
            rng.shuffle(values)
            for i, v in zip(indices, values):
                b_perm[i] = v
        if fisher_pvalue(a, b_perm) <= observed:
            at_least_as_extreme += 1
    return at_least_as_extreme / n_perms


def find_cooccurring_pairs(
    matrix: PresenceMatrix,
    frequency_table: list[dict],
    clade_of_strain: dict[str, str],
    outgroup_presence: dict[str, bool],
    min_strain_count: int = 5,
    fdr_alpha: float = 0.05,
    n_perms: int = 1000,
    seed: int = 0,
) -> list[dict]:
    eligible = [
        row["family"] for row in frequency_table
        if row["bin"] in ("shell", "cloud") and matrix.strain_count(row["family"]) >= min_strain_count
    ]
    clades = [clade_of_strain.get(s, "unknown") for s in matrix.strains]

    candidates = []
    for fam_a, fam_b in itertools.combinations(eligible, 2):
        vec_a = matrix.presence_vector(fam_a)
        vec_b = matrix.presence_vector(fam_b)
        p = fisher_pvalue(vec_a, vec_b)
        candidates.append((fam_a, fam_b, vec_a, vec_b, p))

    if not candidates:
        return []
    qvalues = benjamini_hochberg([c[4] for c in candidates])

    rng = random.Random(seed)
    results = []
    for (fam_a, fam_b, vec_a, vec_b, p), q in zip(candidates, qvalues):
        if q >= fdr_alpha:
            continue
        strains_present_a = [s for s, present in zip(matrix.strains, vec_a) if present]
        results.append({
            "family_a": fam_a,
            "family_b": fam_b,
            "jaccard": jaccard(vec_a, vec_b),
            "fisher_p": p,
            "fdr_q": q,
            "permutation_p": permutation_null_pvalue(vec_a, vec_b, clades, n_perms, rng),
            "direction_a": polarize_direction(outgroup_presence.get(fam_a, False), sum(vec_a) / len(vec_a)),
            "clade_composition": clade_composition(strains_present_a, clade_of_strain),
        })
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--min_strain_count", type=int, default=5)
    ap.add_argument("--fdr_alpha", type=float, default=0.05)
    ap.add_argument("--n_perms", type=int, default=1000)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    NOVINVENIO_LIB = Path(__file__).resolve().parents[4] / "NovInvenio" / "lib"
    sys.path.insert(0, str(NOVINVENIO_LIB))
    from config_parser import parse_config  # noqa: E402

    samples = parse_config(args.config)
    clade_of_strain = {s.short: s.taxon_group for s in samples}
    outgroup_shorts = {s.short for s in samples if s.group == "OUT"}

    matrix = PresenceMatrix.from_tsv(args.matrix)
    with open(args.frequency_table) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        frequency_table = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]

    outgroup_presence = {
        fam: any(matrix.is_present(fam, s) for s in outgroup_shorts) for fam in matrix.families
    }

    pairs = find_cooccurring_pairs(
        matrix, frequency_table, clade_of_strain, outgroup_presence,
        args.min_strain_count, args.fdr_alpha, args.n_perms,
    )
    with open(args.output, "w") as fh:
        fh.write("family_a\tfamily_b\tjaccard\tfisher_p\tfdr_q\tpermutation_p\tdirection_a\tclade_composition\n")
        for row in pairs:
            fh.write(
                f"{row['family_a']}\t{row['family_b']}\t{row['jaccard']:.4f}\t"
                f"{row['fisher_p']:.2e}\t{row['fdr_q']:.2e}\t{row['permutation_p']:.4f}\t"
                f"{row['direction_a']}\t{row['clade_composition']}\n"
            )


if __name__ == "__main__":
    main()
