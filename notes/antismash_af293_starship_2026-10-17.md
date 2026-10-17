# antiSMASH SM clusters vs Starships in A. fumigatus Af293 (2026-10-17)

Folder: `studies/fungi/Afumigatus_pangenome/analysis/starship_island_overlap/full_v072_antismash/` (`O/` below).
Outputs are not committed. Follows `notes/starship_island_overlap_2026-10-17.md` (proxy: PKS KS + NRPS C genes).

## Run

- antiSMASH 8.0.4 (HPCC module `antismash/8.0.4`; the `8.0.4` package dir is a symlink to `8.0.1/`). Pfam 35.0, knownclusterblast 4.0.
  `--check-prereqs`: all satisfied. The 2025_All_Cocci runs used 8.0.1 (`index.html`).
- Input: study Af293 genome `data_dir/dna/UP000002530_330879.dna.fa` + NCBI GFF3 `data_dir/gff3/UP000002530_330879.gff3` (gene models, no gene finding).
- Options: `--taxon fungi --fullhmmer --clusterhmmer --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go -c 8`.
- SLURM 30035698, partition short (stajichlab: Priority; epyc: AssocGrpMemLimit), r07, 8 cpus. Runtime 4,145 s (01:10:10).
- ID mapping: 9,630 antiSMASH CDS (GFF3 locus tags). 9,623 map to pangenome members by exact gene coordinates,
  4 by `GN=` in the protein FASTA, 3 unmapped. 9,627 have a tier-1 family (`O/summary.tsv`).

## Results (`O/summary.tsv`, `O/sm_by_type.tsv`, `O/sm_starship_overlap.tsv`, `O/sm_permutation.tsv`)

- 44 regions, 873 genes, 66 core genes. Types: terpene 12, T1PKS 8, NRPS 6, NRPS-like 4, 2 each T1PKS+terpene,
  isocyanide, terpene-precursor, 1 each of 8 other types or hybrids.
- Overlap with the 8 Table S7 elements: >= 1 gene inside: 3 of 44; >= 50% of genes inside: 1 of 44.
  - CM000169.1.region001 T1PKS, 13/13 genes in AF293_e00029 (Lamia h4).
  - CM000171.1.region002 T1PKS+NRPS-like, 12/30 genes in AF293_s00032 (Osiris h4); top knownclusterblast hit BGC0002267.2 (azasperpyranone).
  - CM000173.1.region009 terpene, 3/7 genes in AF293_s00037 (Nebuchadnezzar h1).
- Clusters with a DUF3435 captain gene: 3 (the same 3 clusters; AFUA_1G00890, AFUA_3G02550, AFUA_5G15040).
- Elements with >= 1 SM cluster: 3 of 8. The other 5 have none.
- Element-shift permutation (1,000, seed 20261017): clusters >= 1 gene: observed 3, null mean 1.88, p(>=) = 0.283.
  Clusters >= 50%: observed 1, null mean 1.21, p(>=) = 0.715. Cluster genes inside: 28, null mean 24.02, p(>=) = 0.383.
- Hypergeometric, all Af293 CDS: cluster genes 28 of 873 inside (expected 16.05, p = 0.00247);
  core genes 2 of 66 inside (expected 1.21, p = 0.343). The gene test treats adjacent genes as independent; the permutation does not.
- Families (non-core universe 40,122; island families 11,683): core-gene families 10 of 15 in islands (expected 4.37, p = 0.00286);
  all cluster-gene families 212 of 310 (expected 90.27, p = 1.39e-46). 50 of 65 core-gene families are core or soft_core.

## Proxy comparison

Proxy: 2 of 32 backbone genes inside Starships, expected 0.59, p = 0.117. antiSMASH core genes: 2 of 66, expected 1.21, p = 0.343.
The 2 genes are the same (AFUA_1G01010 = Q4WKW9 FGNA; AFUA_3G02570 = Q4WFD7). 28 of 32 proxy genes lie in an antiSMASH region;
27 are antiSMASH core genes. Both analyses find no significant SM excess inside Af293 Starships.

## Not verified

- antiSMASH regions are predictions; boundaries were not checked by hand.
- Starship coordinates exist for Af293 only; no other strain was tested.
- The 21 MB `Af293.json.gz` size; whether it should be committed or kept as a release asset is not decided.
