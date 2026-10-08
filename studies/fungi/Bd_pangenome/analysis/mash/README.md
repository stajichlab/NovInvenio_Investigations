# Bd lineage groups from Mash

`run_mash.sh`: `mash sketch -k 21 -s 10000` of the 344 genomes, `mash triangle`.
`cluster_mash.py`: average-linkage clustering of the distances; cut at 0.004. Four singletons
(40.OZ, AP15_IT, CJB4, JEL427-P39, TST75 is the fifth) join the nearest group by mean distance (all
above 0.004, so these are weak placements). Output: `bd_taxongroups.tsv`, `bd_mash_upgma.nwk`.

| Group | Strains | Anchor strains inside |
|---|---|---|
| Mash_A | 283 | JEL423 |
| Mash_B | 26 | CLFT024-02 |
| Mash_C | 24 | none of the anchors checked |
| Mash_D | 7 | KBO_317, KRBOOR_317 |
| Mash_E | 4 | none |

Pairwise distances span 0.00006 to 0.0117. At this resolution BdBRAZIL strains (CJB5-2, CJB7) sit inside
Mash_A. The groups are NOT named lineages; the maintainer must confirm which, if any, match BdGPL,
BdBRAZIL, BdASIA.
