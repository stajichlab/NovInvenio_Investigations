# HrmA-domain genes in the Coccidioides pangenome: presence/absence and genomic neighbourhoods

Date: 2026-09-24. Read-only analysis. Nothing in the pipeline or `results/` was changed. Nothing was committed.

Run analysed: `results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome/` (529 strains: 169 *C. immitis* (Ci) and 360 *C. posadasii* (Cp)). The other-direction run was not used. Family IDs from the two runs were never joined.

All intermediate tables and scripts are in `analysis/hrmA_2026-09-24/` next to this file.

## Plain-language summary

- The Pfam "HrmA" model (PF28515, the small "subtelomeric hrmA-associated cluster protein") hits 2,142 proteins in 42 tier-1 families. Every strain has at least one. Ci strains have a median of 6 per genome (range 4-8). Cp strains have a median of 3 (range 1-7). The difference is significant (Mann-Whitney U, p = 1.8e-69). The count does not correlate with assembly contig number in either species (Spearman rho = -0.06 and 0.06, p > 0.2).
- The two domains of the true hrmA protein (HrmA_N, PF28647, and HrmA_middle, PF28646) are rare. They occur in 26 proteins in 23 strains (7 Ci, 16 Cp). No protein carries both domains. In A. fumigatus Af293, hrmA (Afu5g14900) carries both on one protein.
- The 42 families group into 24 loci by shared conserved flanking genes. 15 loci could be scored per strain. Most large loci are fixed in one species and empty at the same site in the other. Ci-only loci: L04, L06, L07, L11. Cp-only or Cp-biased loci: L01, L02, L08, L13. L03 is in both species, but the clustering splits it into species-specific families (F05 in Ci, F11 in Cp).
- 1,839 of 2,168 HrmA-hit genes (84.8%) have a gene with an APH (PF01636, phosphotransferase) domain as their immediate neighbour. The genome-wide rate for other genes is 1.0% (Fisher OR = 535, p < 1e-300). The same APH-HrmA gene pair exists in A. fumigatus: Afu5g14880 (APH) next to Afu5g14890 (HrmA), and Afu1g00750 (APH) next to Afu1g00760 (HrmA). This pair is the clearest shared feature with the A. fumigatus HAC.
- The full A. fumigatus HAC arrangement (Pkinase, APH, HrmA, hrmA, cgnA, PF11001 protein, HAC6 protein) was not found intact in Coccidioides. The closest match is family F01 (HrmA_middle). In its neighbourhood there is a protein kinase, a Zn_clus protein, a PF11001 protein, BTB proteins and a DUF3723 protein. No cgnA (collagen-like) homolog was found.
- HrmA genes are not strongly enriched at contig ends. 41.9% lie within 20 kb of a contig end. For all genes the value is 40.1% (OR 1.08, p = 0.080). The assemblies are too fragmented to test subtelomeric position: only 647 of 737,258 contig ends carry a telomeric repeat array.
- Several apparent "species-specific families" and many "genome_only" cells are artefacts. Orthologs at one locus are split into separate families (for example F05 and F11), and rescue calls fall on top of annotated paralogs (for example F20 and F12). Read the family-level presence/absence table only together with the locus-level table.

## Methods

Tools: HMMER 3.4 and DIAMOND from the NII pixi environment. Python scripts ran with the NII pixi python (pandas 3.0.3, scipy 1.18.1, statsmodels 0.15.0). All runs were on node c05, inside the existing interactive SLURM allocation. Temporary FASTA files were written to `$SCRATCH` and deleted afterwards.

1. HMMs. `hmmfetch -f Pfam-A.hmm` (Pfam-A 38.2) for HrmA (PF28515.1, GA 27.0), HrmA_middle (PF28646.1, GA 27.0), HrmA_N (PF28647.1, GA 27.0) and, for context, AFUB_07903_YDR124W_hel (PF11001.15, GA 31.3). File: `hrmA_models.hmm`.
2. Census against all tier-1 representatives: `hmmsearch --cut_ga` against `cluster/tier1_rep_seq.fasta` (47,745 reps). Output: `reps.tblout`, `reps.domtblout`.
3. Census against all proteins: all 529 proteomes were joined with a `Short|` prefix (4,504,703 proteins). Then `hmmsearch --cut_ga --cpu 6` was run (`run_allprot.sh`; output `allprot.domtblout`). Every hit protein was mapped to its tier-1 family through `tier1_cluster.tsv`. All 2,168 hit proteins were in the cluster map.
4. Presence/absence: rows of `presence_matrix.rescued.tsv` for the 42 families. Two-sided Fisher exact test per family, Ci vs Cp. Test (a): carried = present or genome_only. Test (b): carried = present only. Benjamini-Hochberg correction was applied across the 42 families for each test.
5. Location: gene coordinates came from `gene_positions.tsv.zst`. Contig lengths came from the `data_dir/dna` FASTA files (`contig_lengths.tsv.gz`). Distance to the nearest contig end = min(start-1, contig_len-end). Telomere check: a contig end is "telomeric" if its terminal 500 bp contain >= 4 tandem TTAGGG or CCCTAA copies (`telomere_ends.py`, `contig_telomere_ends.tsv.gz`). Fisher tests compare HrmA genes against all other genes.
6. Neighbourhoods: genes on each contig were ordered by start coordinate. For every HrmA-hit gene the +/-10 genes were recorded with their family, frequency bin and Pfam domains (`hrmA_neighbourhoods.tsv.gz`). Pfam domains for the 1,403 neighbour families came from `hmmscan --cut_ga` of their representatives against full Pfam-A (`neighbour_reps_pfam.domtblout`). Conservation per family: (i) the modal ordered +/-2 family signature, made orientation-independent, as a fraction of anchors with a complete +/-2 window; (ii) the mean pairwise Jaccard index of +/-10 neighbour-family sets (at most 5,000 random pairs). Exemplars: the two best-assembled carriers (lowest `n_contigs`) per species for each family with >= 10 members (`hrmA_exemplar_neighbourhoods.tsv`). GFF3 `product=` values were collected. All of them are "hypothetical protein", so they give no gene names.
7. Loci (`hrmA_step3_loci.py`): a family's markers are neighbour families in the core or soft_core bin that flank >= 50% of its anchors. Families that share >= 2 markers (or >= 1 marker when one family has only one) were joined into one locus. For each locus and strain the up-to-6 top markers were located. The strain was then scored as follows. `occupied` = an HrmA-hit gene lies within 12 genes of a marker gene. `empty_site` = a marker gene has >= 12 genes on both sides and no HrmA gene is in that window. `unresolved` = markers are present but too close to a contig end. `no_marker` = no marker gene is present. `+rescue` = a TBLASTN rescue position of any HrmA family falls in the marker window. Fisher exact test (occupied vs empty_site, Ci vs Cp), with BH correction across loci.
8. Islands and pairs: HrmA families were looked up in the `member_families` column of `report_tables/islands_with_domains.tsv` (the same 12,197 islands as `significant_islands.tsv`), and in `family_a`/`family_b` of `pair_classification.tsv.zst`.
9. Leiden trans-modules: `family_modules.tsv`, `module_summary.tsv` and `module_domains.tsv` in this run contain only a header. This run has no modules, so module membership cannot be reported.
10. A. fumigatus comparison: Af293 proteins Afu5g14850-Afu5g14950 and Afu1g00700-Afu1g00800 were taken from `FungiDB-68_AfumigatusAf293_AnnotatedProteins.fasta` (`af293_hac_region.fa`). They were then (a) scanned with `hmmscan --cut_ga` against Pfam-A and (b) searched with `diamond blastp --sensitive -e 1e-5` against all tier-1 reps. In addition, `hmmsearch --cut_ga` with APH, HAC6_N and HAC6_helical was run against all 529 proteomes (`run_aph_hac6.sh`, `hrmA_step4_aph_hac6.py`).

Note on frequency bins: `frequency_table.tsv` bins count only the 133 dereplicated representative Ci (IN) strains (max `strain_count` = 133). A Cp-only family therefore shows bin "singleton" and strain_count 0, even when it is in 98% of Cp strains. Use the per-species columns, not the bin, for Cp.

## 1. Census

Rep-level search: 38 reps hit HrmA, 3 hit HrmA_N and 1 hits HrmA_middle. Together that is 42 families. For context, 30 reps hit PF11001. The all-protein search found 2,142 HrmA, 13 HrmA_N and 13 HrmA_middle hits, i.e. 2,168 proteins. The all-protein search added one family whose rep has no GA hit: F01 (`574-0_S_OLD_CPA0039|C27CBD2_004020-T1`, 667 aa). In F01, 11 of 23 members hit HrmA_middle. In every other family, all members hit the same model as the rep. No protein hits more than one HrmA model. No HrmA-hit protein also carries PF11001.

Below-GA homologs (found by DIAMOND with Af293 queries, not counted in the census):
- `B15145|CD39A0E_006880-T1` (69 aa): 73.5% identity over 68 aa to Afu5g14890 (HrmA). Its HrmA score is 21.6 (GA 27). It is present in 71 Ci and 0 Cp strains.
- `VFC085_CPA0068|CB235D1_004510-T1` (344 aa): 38.8% identity to hrmA (Afu5g14900). Its HrmA_middle score is 21.0. It is present in 19 Ci and 91 Cp strains.
- `CA26|C27950D_008239-T1` did not score against any HrmA model. It co-occurs with F42 (see Section 4d).

The GA cutoffs therefore miss at least two homologous families. The census below undercounts HrmA-related genes.

Per-family table. Bin = C. immitis-based bin (see the note above). Present = annotated protein. genome_only = TBLASTN rescue only. q-values are BH-adjusted two-sided Fisher exact tests, Ci vs Cp. Full columns, including copy number, are in `hrmA_family_census.tsv`.

| Fam | Representative | Domain | Rep aa | Members | Bin (Ci) | Ci present / genome_only (of 169) | Cp present / genome_only (of 360) | BH q (any) | BH q (annotated only) | Locus |
|---|---|---|---:|---:|---|---|---|---:|---:|---|
| F01 | `574-0_S_OLD_CPA0039\|C27CBD2_004020-T1` | (rep none; members HrmA_middle:11) | 667 | 23 | cloud | 7 / 0 | 16 / 0 | 1.000 | 1.000 | L12 |
| F02 | `578-1_L_NEW_CPA0049\|C79DB71_000776-T1` | HrmA | 124 | 354 | singleton | 0 / 0 | 354 / 0 | 1.6e-131 | 1.6e-131 | L01 |
| F03 | `578-1_L_NEW_CPA0049\|C79DB71_008063-T1` | HrmA | 122 | 343 | singleton | 0 / 0 | 343 / 0 | 1.4e-118 | 1.4e-118 | L02 |
| F04 | `GT-126\|C8C79ED_001931-T1` | HrmA | 125 | 231 | core | 162 / 0 | 0 / 0 | 5.6e-128 | 5.6e-128 | L04 |
| F05 | `5SD\|CCCC1E1_001005-T1` | HrmA | 117 | 191 | core | 168 / 0 | 22 / 0 | 2.2e-111 | 2.6e-111 | L03 |
| F06 | `5SD\|CCCC1E1_003562-T1` | HrmA | 170 | 167 | core | 167 / 0 | 0 / 0 | 8.7e-137 | 8.7e-137 | L05 |
| F07 | `5SD\|CCCC1E1_000191-T1` | HrmA | 129 | 165 | core | 165 / 0 | 0 / 0 | 4.8e-133 | 4.8e-133 | L06 |
| F08 | `GT-146\|CDEB87D_007880-T1` | HrmA | 117 | 135 | shell | 132 / 0 | 3 / 0 | 8.2e-85 | 1.1e-84 | L07 |
| F09 | `5SD\|CCCC1E1_000093-T1` | HrmA | 124 | 114 | cloud | 11 / 0 | 103 / 0 | 2.2e-09 | 2.9e-09 | L08 |
| F10 | `574-0_S_OLD_CPA0039\|C27CBD2_008386-T1` | HrmA | 119 | 109 | singleton | 0 / 0 | 109 / 0 | 5.3e-21 | 7.7e-21 | L09 |
| F11 | `TX11\|CEC36E2_008626-T1` | HrmA | 87 | 102 | singleton | 0 / 0 | 102 / 0 | 1.8e-19 | 2.5e-19 | L03 |
| F12 | `GT-126\|C8C79ED_008251-T1` | HrmA | 178 | 57 | shell | 57 / 88 | 0 / 0 | 1.4e-104 | 1.5e-31 | L10 |
| F13 | `CA2\|CD8B963_008506-T1` | HrmA | 117 | 37 | shell | 37 / 0 | 0 / 0 | 7.1e-20 | 9.9e-20 | L11 |
| F14 | `730332_Guatemala\|C45F23D_004938-T1` | HrmA | 124 | 36 | singleton | 1 / 0 | 32 / 0 | 1.3e-04 | 1.6e-04 | L05 |
| F15 | `B3314\|CBC6125_000776-T1` | HrmA | 146 | 18 | singleton | 0 / 0 | 18 / 0 | 0.003 | 0.003 | L13 |
| F16 | `582-1_L_NEW_CPA0064\|CF9574D_004546-T1` | HrmA | 120 | 17 | cloud | 3 / 0 | 14 / 0 | 0.407 | 0.463 | L14 |
| F17 | `B3353\|C210CF2_008619-T1` | HrmA | 108 | 14 | singleton | 0 / 0 | 14 / 0 | 0.013 | 0.017 | L15 |
| F18 | `B15142\|CE5AC85_006064-T1` | HrmA | 62 | 13 | cloud | 13 / 0 | 0 / 0 | 6.4e-07 | 8.4e-07 | L07 |
| F19 | `B15146\|CF83403_004486-T1` | HrmA | 103 | 12 | cloud | 12 / 0 | 0 / 0 | 2.0e-06 | 2.6e-06 | L07 |
| F20 | `B14298\|C71D630_004309-T1` | HrmA | 225 | 3 | soft_core | 3 / 151 | 0 / 0 | 5.2e-116 | 0.075 | L06 |
| F21 | `B3417\|C574D2D_005626-T1` | HrmA | 145 | 3 | singleton | 0 / 0 | 3 / 0 | 0.685 | 0.728 | L05 |
| F22 | `4545-MICE\|CB4077B_003725-T1` | HrmA | 180 | 2 | singleton | 0 / 0 | 2 / 156 | 7.8e-33 | 1.000 | L13 |
| F23 | `B3489\|CCD713F_001533-T1` | HrmA | 134 | 2 | singleton | 0 / 0 | 2 / 0 | 1.000 | 1.000 | L02 |
| F24 | `GT-111\|C5A2EA4_008657-T1` | HrmA | 62 | 2 | cloud | 2 / 0 | 0 / 0 | 0.171 | 0.225 | L04 |
| F25 | `UTAH_23742X195\|C3F3A6B_006523-T1` | HrmA | 98 | 2 | singleton | 0 / 0 | 2 / 0 | 1.000 | 1.000 | L01 |
| F26 | `4M3\|CB16A6A_002786-T1` | HrmA | 151 | 1 | shell | 1 / 86 | 0 / 0 | 7.8e-52 | 0.463 | L20 |
| F27 | `B11057\|CCF9097_002456-T1` | HrmA | 165 | 1 | singleton | 1 / 0 | 0 / 136 | 1.8e-25 | 0.463 | L03 |
| F28 | `B15146\|CF83403_001699-T1` | HrmA | 103 | 1 | singleton | 1 / 0 | 0 / 0 | 0.407 | 0.463 | L06 |
| F29 | `B3348\|C8D1221_008493-T1` | HrmA | 140 | 1 | singleton | 0 / 0 | 1 / 0 | 1.000 | 1.000 | L21 |
| F30 | `CA29\|CCF2728_004905-T1` | HrmA | 100 | 1 | singleton | 1 / 0 | 0 / 0 | 0.407 | 0.463 | L07 |
| F31 | `COCSP-8945\|C2B0656_003994-T1` | HrmA | 90 | 1 | singleton | 1 / 0 | 0 / 0 | 0.407 | 0.463 | L22 |
| F32 | `GT-118\|C574A51_008475-T1` | HrmA | 154 | 1 | singleton | 1 / 0 | 0 / 0 | 0.407 | 0.463 | L05 |
| F33 | `RS\|C3CED3A_008636-T1` | HrmA | 202 | 1 | cloud | 1 / 1 | 0 / 0 | 0.171 | 0.463 | L23 |
| F34 | `SDVA14\|C979C11_006396-T1` | HrmA | 71 | 1 | singleton | 0 / 0 | 1 / 0 | 1.000 | 1.000 | L05 |
| F35 | `SDVA14\|C979C11_007553-T1` | HrmA | 126 | 1 | singleton | 0 / 0 | 1 / 0 | 1.000 | 1.000 | L02 |
| F36 | `SJV_2\|CEBEDBD_008473-T1` | HrmA | 62 | 1 | singleton | 1 / 0 | 0 / 0 | 0.407 | 0.463 | L07 |
| F37 | `S_675\|CD889B5_008191-T1` | HrmA | 53 | 1 | singleton | 0 / 0 | 1 / 0 | 1.000 | 1.000 | L03 |
| F38 | `UTAH_23742X170\|C9AA314_008204-T1` | HrmA | 90 | 1 | singleton | 1 / 0 | 0 / 0 | 0.407 | 0.463 | L24 |
| F39 | `M199\|C82F880_004216-T1` | HrmA_N | 290 | 7 | cloud | 3 / 0 | 4 / 5 | 0.912 | 0.872 | L16 |
| F40 | `730332_Guatemala\|C45F23D_000045-T1` | HrmA_N | 186 | 3 | singleton | 0 / 0 | 3 / 1 | 0.407 | 0.728 | L17 |
| F41 | `B3313\|CDF0487_002293-T1` | HrmA_N | 366 | 3 | singleton | 0 / 0 | 3 / 10 | 0.022 | 0.728 | L18 |
| F42 | `730333_Guatemala\|C537722_008849-T1` | HrmA_middle | 427 | 2 | shell | 0 / 23 | 2 / 102 | 2.5e-04 | 1.000 | L19 |

Copy number: 39 of 42 families have at most one member per strain. The exceptions are F04 (69 Ci strains with 2 copies; mean 1.43 per carrier), F05 (1 strain with 2 copies) and F14 (3 Cp strains with 2 copies; mean 1.09). Summed over families, HrmA copy number per strain differs by species: median 6 in Ci, median 3 in Cp.

Families that differ between species (BH q < 0.05, both tests), with the dominant species:
- Cp only or Cp-biased: F02 (354/360 Cp, 0 Ci), F03 (343 Cp, 0 Ci), F09 (103 Cp vs 11 Ci), F10 (109 Cp, 0 Ci), F11 (102 Cp, 0 Ci), F14 (32 Cp vs 1 Ci), F15 (18 Cp), F17 (14 Cp). F16 (3 Ci, 14 Cp) is not significant (q = 0.41).
- Ci only or Ci-biased: F04 (162 Ci, 0 Cp), F05 (168 Ci vs 22 Cp), F06 (167 Ci), F07 (165 Ci), F08 (132 Ci vs 3 Cp), F12 (57 Ci annotated + 88 genome_only), F13 (37 Ci), F18 (13 Ci), F19 (12 Ci).
- F20, F22, F26, F27, F41 and F42 differ only when genome_only cells are counted. Most of their calls are rescue-only (Section 5).

## 2. Loci: presence/absence at the genomic site

Family-level calls mix two things: a real gene gain or loss, and ortholog splitting by the clustering. The locus table resolves this. It asks whether the conserved flanking genes are present and whether an HrmA gene sits between them.

Pairwise identity between family reps (`hrmA_family_rep_pairwise.tsv`) supports these locus groups. Examples: F05-F27 99.1%; F06-F14 100% over 124 aa; F07-F20 98.2%; F08-F19 88.0%; F04-F12 84.8% (F12 is a separate locus); F39-F41 99.6%; F01-F42 97.1% over 339 aa.

| Locus | Families | Markers | Ci occupied | Ci empty site | Ci unresolved | Ci rescue-only | Cp occupied | Cp empty site | Cp unresolved | Cp rescue-only | Ci occupied / resolved | Cp occupied / resolved | BH q |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| L01 | F02,F25 | 6 | 1 | 88 | 79 | 0 | 351 | 1 | 7 | 0 | 1/89 | 351/352 | 8.9e-91 |
| L02 | F03,F23,F35 | 6 | 1 | 112 | 55 | 0 | 344 | 8 | 6 | 0 | 1/113 | 344/352 | 3.3e-96 |
| L03 | F05,F11,F27,F37 | 6 | 169 | 0 | 0 | 0 | 124 | 3 | 120 | 112 | 169/169 | 124/127 | 0.11 |
| L04 | F04,F24 | 6 | 146 | 13 | 10 | 0 | 0 | 294 | 64 | 0 | 146/159 | 0/294 | 1.5e-103 |
| L05 | F06,F14,F21,F32,F34 | 6 | 169 | 0 | 0 | 0 | 33 | 19 | 305 | 2 | 169/169 | 33/52 | 1.2e-13 |
| L06 | F07,F20,F28 | 6 | 169 | 0 | 0 | 0 | 0 | 269 | 89 | 0 | 169/169 | 0/269 | 5.1e-125 |
| L07 | F08,F18,F19,F30,F36 | 6 | 159 | 3 | 7 | 0 | 3 | 196 | 158 | 0 | 159/162 | 3/199 | 1.1e-94 |
| L08 | F09 | 6 | 7 | 118 | 43 | 0 | 63 | 59 | 237 | 0 | 7/125 | 63/122 | 2.0e-16 |
| L11 | F13 | 1 | 33 | 7 | 129 | 0 | 0 | 27 | 136 | 0 | 33/40 | 0/27 | 2.8e-12 |
| L13 | F15,F22 | 6 | 0 | 113 | 56 | 0 | 21 | 123 | 89 | 124 | 0/113 | 21/144 | 4.6e-06 |
| L14 | F16 | 1 | 3 | 35 | 131 | 0 | 15 | 89 | 238 | 0 | 3/38 | 15/104 | 0.51 |
| L16 | F39 | 1 | 2 | 0 | 3 | 0 | 9 | 3 | 38 | 8 | 2/2 | 9/12 | 1.00 |
| L18 | F41 | 1 | 0 | 50 | 117 | 0 | 2 | 109 | 246 | 1 | 0/50 | 2/111 | 1.00 |
| L20 | F26 | 1 | 1 | 0 | 165 | 0 | 0 | 0 | 29 | 0 | 1/1 | 0/0 |  |
| L22 | F31 | 6 | 1 | 67 | 101 | 0 | 1 | 100 | 257 | 0 | 1/68 | 1/101 | 1.00 |

Only loci with >= 1 marker are shown. L09 (F10), L10 (F12), L12 (F01), L15 (F17), L17 (F40), L19 (F42), L21 (F29), L23 (F33) and L24 (F38) have no marker: their anchors sit at contig ends, or the flanking genes are not core. Their locus-level state cannot be scored. Full table: `hrmA_loci.tsv`. Per-strain states: `hrmA_locus_occupancy.tsv.gz`. "Rescue-only" means that no annotated HrmA gene is present, but a TBLASTN rescue position of an HrmA family falls in the marker window. "Resolved" = occupied + empty_site.

Reading the table:
- L01 (F02) and L02 (F03) are Cp loci. In Ci, 88 and 112 strains have the flanking core genes on a contig long enough to show that the site is empty. Only 1 Ci strain has an HrmA gene there.
- L04 (F04), L06 (F07) and L07 (F08 and its satellite families) are Ci loci. In Cp, 294, 269 and 196 strains show an empty site.
- L03 is present in all 169 Ci strains and in 124 of 127 resolved Cp strains (q = 0.11, not significant). The family split is an artefact of the clustering: the Ci allele is in F05 (168 Ci + 22 Cp strains), and most Cp alleles are in F11 (102 Cp). F11 reps are shorter (87 aa vs 117 aa). Another 112 Cp strains carry only a rescue call (F27) at this locus. So "F05 is Ci-specific" and "F11 is Cp-specific" are both wrong at the locus level.
- L05 (F06 in Ci, F14 in Cp) is in all 169 Ci strains and in 33 of 52 resolved Cp strains. 305 Cp strains are unresolved, so the Cp frequency is poorly known.
- L08 (F09) is Cp-biased (63/122 resolved vs 7/125).
- L13 (F15/F22): 21 Cp strains have an annotated gene. 124 Cp strains have only a rescue call. Ci has 113 empty sites.

## 3. Location on contigs

Per family (families with >= 10 members). "APH gene adjacent" = an APH (PF01636) domain gene at offset -1 or +1. Pfam for neighbours comes from the family rep.

| Fam | Genes | Median bp to contig end | <=20 kb | <=50 kb | Median contig length | Median genes to contig end | Anchors with full +/-10 window | Modal +/-2 order fraction | Mean Jaccard +/-10 | APH gene adjacent |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| F01 | 11 | 13,468 | 1.00 | 1.00 | 39,170 | 6 | 0 | 0.67 (n=9) | 0.341 | 0/11 |
| F02 | 354 | 29,569 | 0.44 | 0.73 | 120,508 | 10 | 181 | 0.38 (n=309) | 0.472 | 335/354 |
| F03 | 343 | 42,379 | 0.31 | 0.57 | 165,293 | 11 | 188 | 0.37 (n=229) | 0.573 | 342/343 |
| F04 | 231 | 20,762 | 0.48 | 0.67 | 111,059 | 3 | 92 | 0.49 (n=127) | 0.32 | 179/231 |
| F05 | 191 | 36,181 | 0.34 | 0.86 | 125,904 | 7 | 18 | 0.83 (n=186) | 0.598 | 188/191 |
| F06 | 167 | 50,492 | 0.24 | 0.49 | 126,585 | 10 | 86 | 0.99 (n=138) | 0.475 | 167/167 |
| F07 | 165 | 47,590 | 0.33 | 0.52 | 167,011 | 17 | 106 | 0.57 (n=160) | 0.634 | 165/165 |
| F08 | 135 | 17,931 | 0.55 | 0.76 | 97,005 | 4 | 38 | 0.63 (n=130) | 0.398 | 0/135 |
| F09 | 114 | 50,132 | 0.32 | 0.47 | 155,069 | 10 | 65 | 0.30 (n=86) | 0.483 | 113/114 |
| F10 | 109 | 13,913 | 0.63 | 0.98 | 52,335 | 0 | 0 |  | 0.526 | 107/109 |
| F11 | 102 | 34,368 | 0.32 | 0.72 | 98,124 | 1 | 2 | 0.68 (n=44) | 0.609 | 100/102 |
| F12 | 57 | 12,358 | 0.82 | 0.91 | 32,149 | 0 | 5 | 0.33 (n=6) | 0.707 | 56/57 |
| F13 | 37 | 5,787 | 1.00 | 1.00 | 18,417 | 0 | 0 |  | 0.688 | 33/37 |
| F14 | 36 | 15,073 | 0.53 | 0.94 | 71,216 | 0 | 0 | 0.33 (n=6) | 0.218 | 35/36 |
| F15 | 18 | 35,160 | 0.11 | 0.56 | 249,391 | 17 | 14 | 0.33 (n=15) | 0.584 | 0/18 |
| F16 | 17 | 30,809 | 0.29 | 0.71 | 214,888 | 0 | 3 | 0.50 (n=4) | 0.329 | 17/17 |
| F17 | 14 | 9,070 | 1.00 | 1.00 | 26,013 | 0 | 0 |  | 0.052 | 0/14 |
| F18 | 13 | 19,553 | 0.54 | 0.92 | 100,982 | 5 | 3 | 0.75 (n=12) | 0.554 | 0/13 |
| F19 | 12 | 18,408 | 0.50 | 0.75 | 86,510 | 6.5 | 5 | 0.91 (n=11) | 0.467 | 0/12 |

All HrmA-hit genes compared with the rest of the genome (`hrmA_location_tests.tsv`; two-sided Fisher exact test):

| Test | HrmA genes yes / no | Other genes yes / no | OR | p |
|---|---|---|---:|---:|
| <= 20 kb from a contig end | 909 / 1,259 (41.9%) | 1,804,189 / 2,698,346 (40.1%) | 1.08 | 0.080 |
| <= 50 kb from a contig end | 1,499 / 669 (69.1%) | 2,928,472 / 1,574,063 (65.0%) | 1.20 | 5.5e-05 |
| Nearest contig end has a telomeric repeat | 1 / 2,167 | 596 / 4,501,939 | 3.49 | 0.25 |

Background by frequency bin (`hrmA_location_vs_background.tsv`): <= 20 kb is 38.1% for core genes, 45.7% for shell, 50.5% for cloud and 42.1% for singletons. The median contig length that carries any gene is 153 kb.

Interpretation:
- The assemblies are highly fragmented. Contig N50 values are mostly 100-400 kb, and 1M0 has 1,943 contigs. Even core genes are within 20 kb of a contig end 38% of the time. "Near a contig end" therefore does not mean "subtelomeric".
- Across all 2,168 genes, HrmA genes are only slightly closer to contig ends than other genes (OR 1.20 at 50 kb).
- Some families sit at contig ends almost always: F10 (98% <= 50 kb; median 0 genes to the end; no anchor has a complete +/-2 window), F12, F13, F17 and F42 (on 2.7 kb contigs). F01 is 100% <= 20 kb. For these families, the assembly breaks next to the HrmA gene. This pattern fits a repeat-rich or hard-to-assemble region. It also means that their neighbourhood and presence calls are least reliable.
- Telomere repeats are almost never assembled. Only 647 of 737,258 contig ends (0.09%) have >= 4 TTAGGG/CCCTAA copies. 281 of 529 strains have none. I cannot test from these assemblies whether HrmA genes are subtelomeric. I also did not verify that TTAGGG is the Coccidioides telomere repeat. The test would need chromosome-level assemblies.

## 4. What clusters and neighbourhoods they are in

### 4a. Tier-1 families

See the census table. Sizes range from 1 to 354 members. The biggest families are F02 (354), F03 (343), F04 (231), F05 (191), F06 (167), F07 (165), F08 (135), F09 (114), F10 (109) and F11 (102). 20 families have <= 3 members. Most of these small families are fragments or divergent alleles of a larger family at the same locus (Section 2).

### 4b. Neighbourhoods

The most frequent neighbour families per HrmA family are in `hrmA_top_neighbours.tsv` (8 per family, with bin and Pfam). The full +/-10 windows for the best-assembled strains are in `hrmA_exemplar_neighbourhoods.tsv`.

The main pattern is an **APH-HrmA gene pair**. 1,839 of 2,168 HrmA-hit genes (84.8%) have an APH-domain gene as the immediate neighbour. For non-APH, non-HrmA genes the rate is 46,276 of 4,475,098 (1.0%). Fisher OR = 535, p < 1e-300 (`hrmA_aph_hac6_summary.tsv`). The genomes have 27,437 APH genes (domain-GA; median 52 per strain), so only 1,840 of them (6.7%) sit next to an HrmA gene. The APH partner belongs to a different tier-1 family at each locus: for example `Beeville|C6949CF_002794-T1` at F02, `B14135|CD69CCA_003448-T1` at F03, `578-1_L_NEW_CPA0049|C79DB71_007395-T1` at F05/F11 and `4SD|C8AA998_007437-T1` at F06/F14. The unit therefore looks duplicated across loci. It is not one mobile family. Families without the APH neighbour: F01 and F39-F42 (the HrmA_N/HrmA_middle families), F08/F18/F19 (L07), F15, F17, F22, F24 and five single-member families (F29, F30, F33, F36, F38).

Neighbourhood conservation:
- High (same flanking families in the same order in most strains): F06 (modal +/-2 order in 99% of 138 anchors), F19 (91%), F05 (83%), F18 (75%), F11 (68%), F08 (63%). The flanks are mostly core genes. These are stable chromosomal positions, not mobile elements.
- Lower: F02 (38%), F03 (37%), F09 (30%). Here the +/-2 window includes singleton or cloud neighbours that vary between strains. For example, the F02 window has a DUF3435 singleton at -2/-3 and an APH shell gene at -1. The outer core flanks (DNA_methylase, WD40, Mtf2, COQ9 on one side; AA_permease, FmdA_AmdA on the other) are the same in both exemplar strains.
- F09 (L08) sits in a block of cloud-frequency genes (Glyco_tranf_2, Pkinase, BTB, ABC_tran) next to core genes (Sld7_C, Vps55, Proteasome, Med6). This is the most "accessory-cluster-like" HrmA neighbourhood.
- F01 (HrmA_middle, the most hrmA-like family) sits in cloud/shell genes: Pkinase_fungal/Pkinase at -1, Zn_clus at +3, a PF11001 protein at +4, BTB, DUF3723 and AAA. In 3 of 4 exemplar strains a PF11001 protein is 4 genes away. F01 contigs are short (median 39 kb). No anchor has a complete +/-10 window.

### 4c. Accessory islands and trans-modules

14 of 42 families occur in at least one island in `islands_with_domains.tsv` (302 island-family rows; `hrmA_island_membership.tsv`). The main ones:

| Fam | Islands | Max strains sharing an island | Island size (genes), median (range) |
|---|---:|---:|---|
| F42 | 122 | 2 | 25 (2-57) |
| F09 | 93 | 7 | 20 (2-34) |
| F01 | 23 | 1 | 19 (5-63) |
| F02 | 20 | 84 | 6.5 (4-8) |
| F12 | 16 | 50 | 3 (2-7) |
| F26 | 10 | 42 | 3 (2-7) |
| F41 | 6 | 1 | 53.5 (22-74) |
| F39 | 5 | 1 | 30 (22-38) |

F03, F08, F13, F14, F22 and F31 are in 1-2 islands each. All island classifications are `unexplained_physical` and/or `ambiguous_linkage`. The HrmA_N/HrmA_middle families (F01, F39, F41, F42) sit in large, strain-private islands (up to 74 genes). The APH-HrmA core loci (F04-F07) are not in islands, because they are core in Ci.

Leiden trans-modules: this run has none (all module files contain only a header).

### 4d. Co-occurring partners (pair_classification)

Only 9 of 42 families are in `pair_classification.tsv.zst`: F01, F08, F09, F12, F13, F18, F19, F26 and F42. All 75,849 rows list only *C. immitis* in `clade_composition`. I infer that the pair tests were run on the Ci (IN) strains only, so Cp-only families cannot appear. I did not check the pipeline code to confirm this. Counts for the 302 pairs that include an HrmA family (all fdr_q <= 0.049):

| Fam | trans_unconfirmed | unexplained_physical | ambiguous_linkage | insufficient_data |
|---|---:|---:|---:|---:|
| F01 | 25 | 7 | 2 | 21 |
| F08 | 11 | 0 | 1 | 0 |
| F09 | 24 | 11 | 2 | 0 |
| F12 | 3 | 1 | 0 | 0 |
| F13 | 75 | 0 | 0 | 0 |
| F18 | 3 | 1 | 1 | 0 |
| F19 | 1 | 1 | 1 | 0 |
| F26 | 1 | 2 | 2 | 0 |
| F42 | 77 | 20 | 9 | 0 |

The physical partners are the neighbours from 4b:
- F09 pairs with its APH neighbour `5SD|CCCC1E1_000094-T1` (linkage 1.0, Jaccard 1.0, q = 2.1e-11) and with the cloud genes in its block.
- F12 pairs with its APH neighbour `5SD|CCCC1E1_006841-T1` (linkage 1.0, Jaccard 0.95).
- F01 pairs with the Zn_clus neighbour `B1249_Guatemala|C40D72F_008025-T1` (linkage 0.94) and with the Pkinase neighbour `CA17|C6A2BBC_006283-T1`.
- F42 pairs with the F01 neighbourhood genes (`578-1_L_NEW_CPA0049|C79DB71_000424-T1`, `..._000433-T1` Pkinase, `GT-139|CB29A3D_006564-T1`) and with the below-GA hrmA homologs `VFC085_CPA0068|CB235D1_004510-T1` and `CA26|C27950D_008239-T1`. F42's own annotated genes are on 2.7 kb contigs, so `partner_frac_anchors_within10` is 0 for these pairs. The linkage comes from the rescue positions.

So F01, F42 and the two below-GA families form one hrmA-like (HrmA_middle) neighbourhood in about 20-23 Ci strains. This neighbourhood contains Pkinase, Zn_clus, PF11001, BTB and DUF3723 genes.

The trans partners of F13 (75) and F42 (77) have median Jaccard 0.35 and 0.38. I did not interpret them further.

### 4e. Comparison with the A. fumigatus HAC

Pfam domains of the Af293 genes at the HAC (`af293_hac_region_pfam.domtblout`):

| Af293 gene | FungiDB product | Pfam (GA) |
|---|---|---|
| Afu5g14870 | putative Ser/Thr kinase | Pkinase |
| Afu5g14880 | mitochondrion localization | APH |
| Afu5g14890 | ortholog of Afu1g00760 | HrmA (PF28515) |
| Afu5g14900 | hrmA | HrmA_N + HrmA_middle |
| Afu5g14910 | cgnA, collagen-like | Collagen |
| Afu5g14920 | unknown | PF11001 |
| Afu5g14930 | unknown | HAC6_N + HAC6_helical |
| Afu1g00750 / 00760 / 00780 / 00790 | paralogous chr1 copy | APH / HrmA / PF11001 / DUF3723 |

Homologs in Coccidioides (DIAMOND vs tier-1 reps, e <= 1e-5; `af293_hac_region_hits_annotated.tsv`):
- Afu5g14880 / Afu1g00750 (APH): many hits (the capped limit of 200 reps). The APH neighbour of HrmA is shared (Section 4b).
- Afu5g14890 / Afu1g00760 (HrmA): best hit `B15145|CD39A0E_006880-T1`, 73.5% over 68 aa (below GA; Ci only). Other hits are the PF28515 families F05, F27, F03, F15, F22, F23, F35, F11 and F37 (31.8-60% identity, partial length).
- Afu5g14900 (hrmA): F42 (40.9% over 215 aa), F01 and `VFC085_CPA0068|CB235D1_004510-T1` (38.8% over 85 aa), F39/F41 (50% over 78 aa). No Coccidioides protein matches hrmA over its full length with both domains.
- Afu5g14910 (cgnA): no hit.
- Afu5g14920 / Afu1g00780 (PF11001): 30 rep families. The two PF11001 families next to F01 in the exemplars (`TX14|CEF4951_008063-T1`, `CA17|C6A2BBC_006277-T1`) are among them (36.2% over 298 aa and 35.0% over 280 aa).
- Afu5g14930 (HAC6): hits `5SD|CCCC1E1_008641-T1` (40.8% over 201 aa) and others. `hmmsearch` finds 1,180 HAC6-domain genes in 15 families (833 in Cp, 347 in Ci). Only 4 lie within 10 genes of an HrmA gene. 82 genes of family `B3315|C362512_005075-T1` lie near F09, at a median of 15 genes. The median distance for all 92 HAC6 genes on a contig with an HrmA gene is 37.7 kb (`hac6_to_hrmA_distance.tsv`).

Conclusion: the APH-HrmA(PF28515) pair is conserved and repeated at many loci. The hrmA-like HrmA_middle genes (F01/F42) keep a partial HAC-like context (Pkinase, PF11001, DUF3723 nearby). cgnA was not found. HAC6 is not adjacent. The A. fumigatus study screened `hacA` (Afu3g04070). FungiDB names this gene "Transcriptional activator hacA, putative", the UPR bZIP. It is not a HAC component, so it was not used here.

## 5. Artefacts and caveats

1. **Ortholog splitting by clustering.** Alleles at one locus fall into different families by species or length: L03 = F05 (Ci) + F11 (Cp) + F27 + F37; L05 = F06 (Ci) + F14 (Cp) + others; L07 = F08 + F18 + F19 + F30 + F36. Family-level "species-specific" results for F05, F11 and F14 are artefacts of this split. Use the locus table.
2. **Rescue-only (genome_only) calls on top of paralogs** (`hrmA_rescue_overlap.tsv`; a rescue start within +/-500 bp of an annotated gene):
   - F12: 87 of 88 genome_only cells overlap annotated F04 genes.
   - F20: 151 of 151 overlap F07 genes.
   - F26: 86 of 86 overlap F04 genes.
   These cells are the same locus counted twice. They make F20 look "soft_core" (strain_count 120) with only 3 annotated members.
   - F42: 123 of 125 genome_only cells overlap an annotated gene. For 108 of them this is the below-GA hrmA-like gene `VFC085_CPA0068|CB235D1_004510-T1`, and for 12 it is an F01 gene.
   - F22 (156 Cp) and F27 (136 Cp): mostly overlap other genes or fall in unannotated sequence at L13/L03. They may be real unannotated or divergent alleles. I could not decide from these data.
   This matches the run-level diagnostic: 77.5% of rescuable cells overlap a gene of another family.
3. **Assembly fragmentation.** Median contig length at HrmA genes is 118 kb. F10, F12, F13, F17, F42 and F01 are nearly always at a contig end. The locus test excludes strains where the site is not resolved. Resolved strains are 52-352 per locus in the main loci, so the species calls rest on subsets.
4. **Split hrmA gene models.** HrmA_N and HrmA_middle never occur on one protein. 3 strains (730332_Guatemala, 730333_Guatemala, B1249_Guatemala) have both domains, on different contigs. In 2 of them the HrmA_middle gene is on a 2.7 kb contig. The data cannot tell a fragmented assembly or split gene model apart from a real two-gene arrangement.
5. **GA cutoffs miss homologs** (Section 1). The census is a lower bound.
6. **Bins are C. immitis-only.** Bins and strain_count ignore Cp.
7. **Pair and island tables** appear to be Ci-driven (4d). Cp-specific loci cannot show pairs.
8. **Gene orientation** was not analysed. `gene_positions.tsv.zst` has no strand column. I did not determine whether the APH-HrmA pair is divergent, convergent or tandem.
9. **Telomere test** is not informative with these assemblies (Section 3).
10. **GFF3 product names** are all "hypothetical protein". They give no gene names.

## Files (in `analysis/hrmA_2026-09-24/`)

- Scripts: `run_allprot.sh`, `run_aph_hac6.sh`, `telomere_ends.py`, `hrmA_step1_families.py`, `hrmA_step2_context.py`, `hrmA_step3_loci.py`, `hrmA_step4_aph_hac6.py`. Run them with `pixi run --manifest-path <NII>/pixi.toml python <script>` in the order 1, 2, 3, 4. Step 2 needs `neighbour_reps_pfam.domtblout`, `af293_hac_region_vs_reps.tsv` and `contig_telomere_ends.tsv.gz` first.
- HMM outputs: `hrmA_models.hmm`, `reps.domtblout`, `allprot.domtblout`, `aph_hac6_allprot.domtblout.gz`, `neighbour_reps_pfam.domtblout`, `af293.domtblout`, `af293_hac_region_pfam.domtblout`.
- Tables: `hrmA_family_census.tsv`, `hrmA_presence_matrix.tsv.gz`, `hrmA_gene_locations.tsv`, `hrmA_family_locations.tsv`, `hrmA_location_tests.tsv`, `hrmA_location_vs_background.tsv`, `hrmA_location_telomere_background.tsv`, `hrmA_neighbourhoods.tsv.gz`, `hrmA_top_neighbours.tsv`, `hrmA_neighbourhood_conservation.tsv`, `hrmA_exemplar_neighbourhoods.tsv`, `hrmA_loci.tsv`, `hrmA_locus_occupancy.tsv.gz`, `hrmA_island_membership.tsv`, `hrmA_pairs.tsv`, `hrmA_in_pair_table.tsv`, `hrmA_rescue_overlap.tsv`, `hrmA_APH_adjacency.tsv`, `hrmA_aph_hac6_summary.tsv`, `hac6_to_hrmA_distance.tsv`, `hrmA_genes_per_strain.tsv`, `hrmA_family_rep_pairwise.tsv`, `af293_hac_region_hits_annotated.tsv`, `contig_lengths.tsv.gz`, `contig_telomere_ends.tsv.gz`.
