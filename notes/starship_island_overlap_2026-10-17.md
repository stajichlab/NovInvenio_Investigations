# Starship / captain / SM overlap with pangenome islands (2026-10-17)

Abbreviations: `AF` = `studies/fungi/Afumigatus_pangenome`, `CP` = `studies/fungi/coccidioides_pangenome`,
`SIO` = `analysis/starship_island_overlap`. Outputs are not committed.

## Annotation sources

| Run | Starship elements (coordinates) | Captain (DUF3435, PF11917.14, E <= 1e-5) | SM annotation |
|---|---|---|---|
| AF full_v072 | Table S7, Gluck-Thaler 2025: 8 Af293 elements, 7 high confidence (`AF/SIO/inputs/af293_starships_table_s7.elements.tsv`). Af293 only. | `AF/results/captain_gene/DUF3435_vs_study.tblout`, all 295 strains, 2,233 hits | No antiSMASH for the study proteins. Proxy: `AF/results/sm_backbone/SM_backbone_vs_study.tblout` (PKS KS PF00109 + NRPS C PF00668), all strains |
| CP immitis_in_posadasii_out_v5 | not found (`CP/analysis/locus_summary/README.md:14`) | new: `CP/SIO/inputs/DUF3435_vs_study.tblout.zst`, 529 strains, 6,292 hits | antiSMASH 8.0.1 regions, 479 of 529 strains, own funannotate annotation (not the 20260112 freeze) |
| CP posadasii_in_immitis_out_v2 | not found | same tblout (2,323 hits in strains not in this run) | same, 335 of 362 strains |

## Island counts (report.md)

AF 11,676 (`AF/results/full_v072/config/pangenome/report/report.md:63`); CP v5 11,968 (`.../report.md:87`);
CP posadasii-in v2 16,014 (`.../report.md:65`). The 8,226 figure is the N=121 sweep, not full_v072.

## Headline (AF full_v072, `AF/SIO/full_v072/`)

- Starship families (>= 1 Af293 gene midpoint inside an element): 167; 163 in >= 1 island (summary.tsv).
  Non-core universe: 163/163 vs 11,683/40,122 (29.1%); hypergeometric p = 2.04e-88 (enrichment.tsv).
- Islands with >= 1 Starship family: 984 of 11,676 (8.4%). Island families that are Starship families: 163 of 11,683 (1.4%).
- Coordinate level in Af293: 9 islands have an Af293 instance overlapping an element; all 7 high-confidence elements
  overlapped; each element's genes lie 94-100% inside Af293 island instances (per_starship.tsv). All 9 are single-strain islands.
- Gene-level metric of the sweep: reproduced on full_v070 (21.28% vs 5.23%); full_v072 21.35% vs 5.23% (reproduce_sweep_metric/).
- Captain families: 217 non-core, 122 in islands (56.2%), p = 6.78e-17.
- Caveat: 961 of 9,647 Af293 members (Swiss-Prot 'sp|' entries) have no row in gene_positions; 941 recovered
  (with_unplaced_af293_coords/). The UniProt table used for recovery has no provenance record.

## CP

- v5: captain families 354 non-core, 238 in islands (67.2%) vs 23.3%, p = 1.93e-69. v2: 270, 211 (78.1%) vs 38.0%, p = 1.37e-41.
- SM: v5 12,528 regions, 113 contain a captain gene; captain genes in SM regions 155 vs 329.21 expected (depletion p = 3.71e-28).
  C. immitis T1PKS: 89 of 861 regions contain a captain gene (89 strains, one region each; 88 on RagTag scaffold GG704915).
  v2: 8,739 regions, 22 with a captain gene; 26 vs 212.62 expected (depletion p = 2.69e-61).

## SM in AF (with_unplaced_af293_coords)

32 Af293 PKS/NRPS backbone genes; 2 inside Starships (Q4WFD7 in AF293_s00032 Osiris h4; FGNA Q4WKW9 in AF293_e00029 Lamia h4);
expected 0.59, p = 0.117. Families: 90 of 308 non-core backbone families in islands, expected 89.69, p = 0.506.
