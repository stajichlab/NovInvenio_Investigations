# Coccidioides pangenome species tree (for Dollo polarization)

Built 2026-09-25. These trees are the input for `bin/pangenome_cooccurrence.py --species_tree` (issue #140 / PR #141).
They are not wired into any run yet.

## Files

| File | Use |
|---|---|
| `coccidioides_genus_with_Uree.full.nwk` | 530 tips. For `rescue_structural_genus_vs_ureesii` (and `mmseqs_genus_vs_ureesii`). |
| `coccidioides_genus_with_Uree.collapsed_ufboot70.nwk` | The same tree with nodes of UFboot < 70 collapsed. |
| `coccidioides_coccidioides_only.full.nwk` | 529 tips, no U. reesii. For the two `rescue_freqpol_*` runs. |
| `coccidioides_coccidioides_only.collapsed_ufboot70.nwk` | The same tree, collapsed. |
| `tip_name_map.tsv` | Every CDS-tree tip mapped to its pangenome `Short`, or marked UNMAPPED. |
| `grafts.tsv` | Every strain that is missing from the CDS tree and was grafted in by mash distance. |
| `build_species_tree.py`, `build_log.txt` | Build script and its log. |
| `validate_species_tree.py`, `validation.txt` | Validation script and its results. |

## Sources

- **Topology and support.** The existing 50-BUSCO-locus CDS tree: 82,437 nt, IQ-TREE, partitioned, UFboot.
  Path: `/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Phylogeny/results/msa_filter_cds_ascomycota-buildtree/Cocci_cds.488taxa_ascomycota.fa.part.aicc.contree`.
  "ascomycota" in the name refers to the BUSCO lineage set. **The tree has no non-Coccidioides tips.**
- **Grafts.** The pipeline's own mash distance matrix (k=21, 1000 hashes). It is not published.
  Path: `.nf_launch/genus_vs_ureesii/work/9b/04c1b520f258337221235ee5b85e9d/all_strains.mash_dist.tsv`. Work dir deleted 2026-09-26; a byte-identical gzipped copy is at `results/preserved_from_work/genus_vs_ureesii/all_strains.mash_dist.tsv.gz` (SHA256SUMS alongside).
- **Names and species.** `config_genus_vs_ureesii.csv`, which carries the species corrections in `CORRECTIONS.md`.

## Name mapping

- 478 of the 488 CDS tips map to pangenome strains.
- The mapping strips the `Coccidioides_<species>_` prefix, then matches on `Strain`/`Short`.
- **Spelling variant.** `Coccidioides_posadasii_Coahuilla_2` was mapped to `Coahuila_2`. Its mash nearest neighbour is not in its CDS parent clade. This is normal: for the 478 mapped strains, the mash nearest neighbour lies in the ancestor-1/2/3 clade only 33% / 44% / 50% of the time. The mash data therefore give no evidence against this mapping, and none for it.
- **Species-prefix conflicts, kept.** Two tips have a species prefix that disagrees with the config:
  - `Coccidioides_immitis_485B-1_L_OLD_CPA0023` maps to `485B-1_L_OLD_CPA0023`, which the config calls Cp.
  - `Coccidioides_posadasii_B3476` maps to `B3476`, which the config calls Ci.

  These are the two corrections in `CORRECTIONS.md`. In the CDS tree, each tip sits inside a clade of the corrected species. So the tree agrees with the correction, and the prefix is only a stale file name.
- **Unmapped, pruned (10).**
  - `CimmitisRS`, `CimmitisWA211`, `CposadasiiSilveira2022`. These are reference-assembly duplicates of `RS`, `WA_211` and `UCSF_Cp_Silveira`, which map from their own tips.
  - `GT-153`, `NM_459`, `SOIL_582-1_S_NEW_CPA0065`, `VFC054_7SA`, `M158` (Cp); `SD7`, `SJV_3` (Ci). These are not in the pangenome. `SOIL_582-1_S_NEW_CPA0065` is not `582-1_L_NEW_CPA0064`: it is a different isolate code.

## Rooting

- The CDS tree has no Ascomycota outgroup tips, so it cannot be rooted on one.
- It was rooted on the Ci | Cp bipartition. This bipartition has UFboot 100 on both sides.
- U. reesii (`Uree`) was then added as the sister of all Coccidioides at a new root, so the root has two children: (Coccidioides, Uree).
- This is the position an outgroup rooting would give if Coccidioides is monophyletic relative to Uncinocarpus. That is expected biology, but this tree does not test it.
- The Uree branch length and the Coccidioides stem length are arbitrary (0.05). Dollo does not use branch lengths.
- Bio.Phylo left a trifurcating root, with Cp split in two. The two pure-Cp children were wrapped back into one Cp node (support 100). The unrooted topology is unchanged.

## Grafts (51)

- The pangenome has 529 Coccidioides strains. 478 are in the CDS tree and 51 are not: 28 Cp and 23 Ci.
- Each missing strain was attached with zero-length branches as the sister of its mash-nearest strain among the original CDS-tree strains, within the same species in every case.
- Several missing strains attached to the same sister strain. They are nested, so the full tree stays binary.
- Mash distance to the sister: min 0.00029, median 0.0031, max 0.0046. This is similar to typical within-species pair distances (median 0.0032 Ci, 0.0047 Cp).
- For 17 of the 51 grafts, the nearest strain overall is another missing strain.
- **Graft positions below the species level are therefore weak.**

## Collapse

- Nodes with UFboot < 70 were collapsed. The root and its two children are protected. Graft nodes carry no support and are not collapsed.
- 111 nodes were collapsed: internal nodes went from 529 to 418.
- The result has 30 polytomies. The largest has 26 children.

## Validation (`validation.txt`)

- **Tip sets.** Each tree was checked against both `presence_matrix.tsv` and `presence_matrix.rescued.tsv` of its run(s). In every case: 0 missing tips, 0 extra tips, 0 duplicates.
- **Parsing.** All trees parse with `Bio.Phylo.read(..., 'newick')`, the parser `bin/pangenome_cooccurrence.py:533` uses.
- **Root.** The root has 2 children, as `lib/ancestral_states.py:87` requires.
- **Monophyly.** Ci (169) and Cp (360) are each monophyletic in all four trees.
- **Dollo smoke test.** `dollo_polarize` was run on 1,000 random families from `rescue_structural_genus_vs_ureesii/presence_matrix.rescued.tsv`, with presence = `present` or `genome_only`, as in `lib/pangenome_matrix.py:113-114`.
  - 0 errors.
  - About 0.09 s per family per tree.
  - Collapsed vs full tree: 425 families have more loss events, 0 have fewer, 575 are equal. The total is 39,511 vs 31,497 (+25%).

## Caveat: polytomies inflate Dollo loss counts

`lib/ancestral_states.py:111-119` walks every child of a node. Each all-absent child counts as one loss. At a polytomy, k absent children therefore give k losses, although one resolution of the polytomy could explain them with one loss. The collapsed tree therefore gives upper-bound loss counts, not minimum ones. The full tree gives counts for one arbitrary resolution of poorly supported nodes. Neither count is the Dollo minimum over resolutions.
