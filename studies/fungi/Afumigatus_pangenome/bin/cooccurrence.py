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
from novinvenio_path import add_novinvenio_lib_to_path  # noqa: E402
from pangenome_matrix import PresenceMatrix  # noqa: E402


def presence_vector_within(matrix: PresenceMatrix, family: str, strains: list[str]) -> list[bool]:
    """Like PresenceMatrix.presence_vector, but restricted to a strain
    subset (e.g. ingroup-only) instead of always using matrix.strains."""
    return [matrix.is_present(family, s) for s in strains]


def strain_count_within(matrix: PresenceMatrix, family: str, strains: list[str]) -> int:
    return sum(presence_vector_within(matrix, family, strains))


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


def polarize_direction(outgroup_present_count: int, outgroup_total: int) -> str:
    """A family present in EVERY outgroup strain implies the ingroup
    strains lacking it lost it ("loss"); a family absent from every
    outgroup strain implies the ingroup strains carrying it gained it
    ("gain"). With no outgroup strains to check, or with outgroup
    strains disagreeing (present in some but not all), the polarity
    can't be called -- "ambiguous"."""
    if outgroup_total == 0 or 0 < outgroup_present_count < outgroup_total:
        return "ambiguous"
    if outgroup_present_count == outgroup_total:
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
    # Phipson & Smyth (2010): never report a Monte-Carlo p-value of exactly
    # 0.0 -- the true p-value is bounded below by 1/(n_perms+1), not 0.
    return (at_least_as_extreme + 1) / (n_perms + 1)


def find_cooccurring_pairs(
    matrix: PresenceMatrix,
    frequency_table: list[dict],
    clade_of_strain: dict[str, str],
    outgroup_presence: dict[str, tuple[int, int]],
    min_strain_count: int = 5,
    fdr_alpha: float = 0.05,
    n_perms: int = 1000,
    seed: int = 0,
    strains: list[str] | None = None,
) -> list[dict]:
    """`strains` restricts every ingroup statistic (the strain-count floor,
    Fisher presence vectors, clade composition) to that subset -- callers
    doing a real ingroup-vs-outgroup study MUST pass the ingroup-only strain
    list here, or outgroup strains silently inflate strain counts and
    contaminate clade composition. Defaults to matrix.strains for backward
    compatibility with callers that have no outgroup (or already pre-filtered
    the matrix). `outgroup_presence` maps family -> (outgroup_present_count,
    outgroup_total), always computed over the FULL matrix (it needs to see
    the outgroup strains) -- see polarize_direction for how the counts are
    used."""
    strains = matrix.strains if strains is None else strains

    eligible = [
        row["family"] for row in frequency_table
        if row["bin"] in ("shell", "cloud")
        and strain_count_within(matrix, row["family"], strains) >= min_strain_count
    ]
    clades = [clade_of_strain.get(s, "unknown") for s in strains]

    candidates = []
    for fam_a, fam_b in itertools.combinations(eligible, 2):
        vec_a = presence_vector_within(matrix, fam_a, strains)
        vec_b = presence_vector_within(matrix, fam_b, strains)
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
        strains_present_a = [s for s, present in zip(strains, vec_a) if present]
        out_count_a, out_total_a = outgroup_presence.get(fam_a, (0, 0))
        results.append({
            "family_a": fam_a,
            "family_b": fam_b,
            "jaccard": jaccard(vec_a, vec_b),
            "fisher_p": p,
            "fdr_q": q,
            "permutation_p": permutation_null_pvalue(vec_a, vec_b, clades, n_perms, rng),
            "direction_a": polarize_direction(out_count_a, out_total_a),
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

    add_novinvenio_lib_to_path()
    from config_parser import parse_config  # noqa: E402

    samples = parse_config(args.config)
    clade_of_strain = {s.short: s.taxon_group for s in samples}
    outgroup_shorts = [s.short for s in samples if s.group == "OUT"]
    ingroup_shorts = [s.short for s in samples if s.group == "IN"]

    matrix = PresenceMatrix.from_tsv(args.matrix)
    with open(args.frequency_table) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        frequency_table = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]

    outgroup_total = len(outgroup_shorts)
    outgroup_presence = {
        fam: (
            sum(1 for s in outgroup_shorts if matrix.is_present(fam, s)),
            outgroup_total,
        )
        for fam in matrix.families
    }

    pairs = find_cooccurring_pairs(
        matrix, frequency_table, clade_of_strain, outgroup_presence,
        args.min_strain_count, args.fdr_alpha, args.n_perms,
        strains=ingroup_shorts,
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
