#!/usr/bin/env python3
"""Real-data run of `synteny_windows.py`'s `accessory_islands()` (built and
tested during the design phase, 2026-09-13, never yet run against real
data -- see PANGENOME_CLUSTER_PROFILE_NOTES.md) against the corrected
293-strain pangenome, then cross-references the resulting islands against
the statistically-significant physically-linked co-occurring pairs
(starship_explained / unexplained_physical / ambiguous_linkage in
pair_classification.rescued.tsv) to answer: which of those pairs are part
of a REAL multi-gene genomic island (not just adjacent by k-window
coincidence), and of those islands, which have a Starship captain gene
(DUF3435) nearby vs. a secondary-metabolite backbone gene (PKS
ketosynthase / NRPS condensation domain) nearby vs. neither.

Reuses already-computed per-strain gene order directly from
family_positions.rescued.tsv's `rank` column (a strain's full,
contig-then-start-sorted gene order, already built by
build_family_positions.py from each strain's own GFF3) instead of
re-parsing 295 GFF3 files from scratch -- `accessory_islands()` only
consumes (gene_id, contig) pairs in traversal order, so a rank-sorted
family_positions.tsv row list is a drop-in equivalent to a freshly-parsed
GFF3 gene_order list for this purpose. Family ID stands in for gene ID
throughout (this study's family/protein-ID convention:
PANGENOME_CLUSTER_PROFILE_NOTES.md's pipeline-conventions item 1-2).

Usage:
  find_accessory_islands.py \\
      --family_positions results/full_293run/family_positions.rescued.tsv \\
      --frequency_table results/full_293run/frequency_table.rescued.tsv \\
      --pair_classification results/full_293run/pair_classification.rescued.tsv \\
      --captain_tblout results/captain_gene/DUF3435_vs_study.tblout \\
      --sm_backbone_tblout results/sm_backbone/SM_backbone_vs_study.tblout \\
      --cluster_tsv results/full_293run/tier1_cluster.tsv \\
      --output results/accessory_islands/significant_islands.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import read_cluster_tsv  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from synteny_windows import accessory_islands  # noqa: E402


PHYSICAL_CLASSIFICATIONS = frozenset(
    {"starship_explained", "unexplained_physical", "ambiguous_linkage"}
)


def load_is_core(frequency_table_path: str) -> dict[str, bool]:
    """{family: True if core/soft_core, False if shell/cloud/singleton}."""
    is_core: dict[str, bool] = {}
    with open(frequency_table_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            is_core[row["family"]] = row["bin"] in ("core", "soft_core")
    return is_core


def load_strain_gene_orders(family_positions_path: str) -> dict[str, list[tuple]]:
    """{Short: [(family, contig, rank, rank), ...]}, sorted by rank -- rank
    is already a per-strain, contig-then-start-sorted global index (see
    build_family_positions.py), so sorting by it recovers genomic order
    without re-reading any GFF3. The tuple shape (family, contig, rank,
    rank) matches accessory_islands()'s expected 4-tuple even though only
    the first two fields are ever read."""
    by_strain: dict[str, list[tuple]] = {}
    with open(family_positions_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            short, family, contig, rank = row["Short"], row["family"], row["contig"], int(row["rank"])
            by_strain.setdefault(short, []).append((family, contig, rank, rank))
    for rows in by_strain.values():
        rows.sort(key=lambda r: r[2])
    return by_strain


def load_significant_physical_pairs(pair_classification_path: str) -> dict[frozenset, str]:
    """{frozenset({family_a, family_b}): classification} for every pair
    classified as physically linked (any of the 3 PHYSICAL_CLASSIFICATIONS)."""
    pairs: dict[frozenset, str] = {}
    with open(pair_classification_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["classification"] in PHYSICAL_CLASSIFICATIONS:
                pairs[frozenset({row["family_a"], row["family_b"]})] = row["classification"]
    return pairs


def load_hit_families(tblout_path: str, member_to_rep: dict[str, str], id_sep: str = "|") -> set[str]:
    """Families with at least one qualifying hmmsearch hit anywhere in the
    cohort -- same tblout parsing convention as
    pair_classification.py's load_captain_families, generalized to any
    hmmsearch --tblout (works identically for DUF3435 or the PKS/NRPS
    backbone search)."""
    families: set[str] = set()
    with open_maybe_compressed(tblout_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            target = line.split()[0]
            if id_sep not in target:
                continue
            short, protein_id = target.split(id_sep, 1)
            family = member_to_rep.get(f"{short}{id_sep}{protein_id}")
            if family is not None:
                families.add(family)
    return families


def build_pair_index(significant_pairs: dict[frozenset, str]) -> dict[str, list[tuple[frozenset, str]]]:
    """{family: [(pair, classification), ...]} for every significant pair
    touching that family -- lets find_significant_islands() check only the
    pairs that could possibly be relevant to a given island (those
    touching one of its member families) instead of scanning all
    significant pairs for every island. At real scale (thousands of
    islands x tens of thousands of pairs), the naive full-scan version was
    O(islands x total_pairs); this makes it O(islands x island_size x
    avg_pairs_per_family), which is what actually finished in reasonable
    time -- verified by killing a run of the naive version after 12+
    minutes with no output (2026-09-15)."""
    index: dict[str, list[tuple[frozenset, str]]] = {}
    for pair, classification in significant_pairs.items():
        for family in pair:
            index.setdefault(family, []).append((pair, classification))
    return index


def find_significant_islands(
    strain_gene_orders: dict[str, list[tuple]],
    is_core: dict[str, bool],
    significant_pairs: dict[frozenset, str],
) -> list[dict]:
    """For every strain, build accessory islands, then keep only islands
    that contain at least one significant physically-linked pair (both
    members of the pair present in the SAME island in THAT strain) --
    the real, adjacency-confirmed evidence this whole script exists to
    surface, as opposed to accessory_islands() alone (which would also
    report every non-core run regardless of any statistical support)."""
    pairs_by_family = build_pair_index(significant_pairs)
    results: list[dict] = []
    for short, gene_order in strain_gene_orders.items():
        islands = accessory_islands(gene_order, is_core)
        for island in islands:
            if len(island) < 2:
                continue
            members = [g[0] for g in island]
            member_set = set(members)
            candidates: dict[frozenset, str] = {}
            for family in member_set:
                for pair, classification in pairs_by_family.get(family, ()):
                    if pair <= member_set:
                        candidates[pair] = classification
            if not candidates:
                continue
            results.append({
                "strain": short,
                "contig": island[0][1],
                "size": len(island),
                "members": members,
                "supporting_pairs": list(candidates.items()),
            })
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family_positions", required=True)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--pair_classification", required=True)
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--captain_tblout")
    ap.add_argument("--sm_backbone_tblout")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    is_core = load_is_core(args.frequency_table)
    strain_gene_orders = load_strain_gene_orders(args.family_positions)
    significant_pairs = load_significant_physical_pairs(args.pair_classification)
    print(f"find_accessory_islands: {len(significant_pairs)} significant physically-linked "
          f"pairs, {len(strain_gene_orders)} strains", file=sys.stderr)

    islands = find_significant_islands(strain_gene_orders, is_core, significant_pairs)
    print(f"find_accessory_islands: {len(islands)} significant islands found across all strains",
          file=sys.stderr)

    captain_families: set[str] = set()
    sm_families: set[str] = set()
    if args.captain_tblout or args.sm_backbone_tblout:
        member_to_rep = read_cluster_tsv(args.cluster_tsv)
        if args.captain_tblout:
            captain_families = load_hit_families(args.captain_tblout, member_to_rep)
        if args.sm_backbone_tblout:
            sm_families = load_hit_families(args.sm_backbone_tblout, member_to_rep)

    with open(args.output, "w") as out:
        out.write(
            "n_strains\texample_strain\tisland_size\tmember_families\t"
            "n_supporting_pairs\tclassifications\thas_captain_gene\thas_sm_backbone_gene\n"
        )
        # Dedupe by (frozenset(members), classifications) so the same
        # recurring island in many strains collapses to one row with a
        # strain-support count, rather than one row per strain.
        by_key: dict[tuple, dict] = {}
        for island in islands:
            member_set = frozenset(island["members"])
            classifications = frozenset(c for _, c in island["supporting_pairs"])
            key = (member_set, classifications)
            entry = by_key.setdefault(key, {
                "strains": [], "size": island["size"], "members": island["members"],
                "classifications": classifications,
                "n_supporting_pairs": len(island["supporting_pairs"]),
            })
            entry["strains"].append(island["strain"])

        for (member_set, classifications), entry in sorted(
            by_key.items(), key=lambda kv: -len(kv[1]["strains"])
        ):
            has_captain = bool(member_set & captain_families)
            has_sm = bool(member_set & sm_families)
            out.write(
                f"{len(entry['strains'])}\t{entry['strains'][0]}\t{entry['size']}\t"
                f"{','.join(entry['members'])}\t{entry['n_supporting_pairs']}\t"
                f"{','.join(sorted(classifications))}\t"
                f"{'Y' if has_captain else 'N'}\t{'Y' if has_sm else 'N'}\n"
            )

    print(f"find_accessory_islands: {len(by_key)} distinct significant islands "
          f"(deduped across strains) written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
