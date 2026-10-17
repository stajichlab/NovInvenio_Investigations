# Genomic context of the domains shared by A. fumigatus and Coccidioides (2026-10-17)

Question: the Pfam domains that are over-represented in A. fumigatus (full_v072, "Af") and in both Coccidioides runs (immitis v5, posadasii v2) — is their genomic context the same in the two genera?

Directory: `studies/fungi/coccidioides_pangenome/analysis/shared_domains_context/` (`SDC/` below). Nothing was committed, pushed or published. `docs/` was not touched.

## 0. Domains (Pfam 38.2, `Pfam-A.hmm.dat`)

| label | ID | AC | DE | GA |
|---|---|---|---|---|
| tier 1 | APH | PF01636.30 | Phosphotransferase enzyme family | 23.2 |
| tier 1 | DUF3435 | PF11917.14 | Protein of unknown function (DUF3435) | 25.8 |
| tier 1 | DUF3723 | PF12520.15 | Protein of unknown function (DUF3723) central domain | 25.0 |
| tier 1 | DUF3723_N | PF29418.1 | Protein of unknown function (DUF3723) N-terminal domain | 25.0 |
| tier 2 | AFUB_07903_YDR124W_hel ("PF11001") | PF11001.15 | YDR124W-like, helical bundle domain | 31.3 |
| HrmA (one label) | HrmA / HrmA_middle / HrmA_N | PF28515.1 / PF28646.1 / PF28647.1 | Subtelomeric hrmA-associated cluster protein / HrmA middle domain / HrmA N-terminal domain | 27 |
| secondary | NPHP3_N, Condensation, PP-binding, AMP-binding, SPOK2_N | PF24883.3, PF00668.26, PF00550.32, PF00501.35, PF27653.1 | Nephrocystin 3, N-terminal; Condensation domain; Phosphopantetheine attachment site; AMP-binding enzyme; Spore killer 2 N-terminal domain | |

## 1. Methods (short)

- **Gene level.** `hmmsearch --cut_ga` of the 13 models (`SDC/hmm/shared_domains.hmm`, hmmfetch) over the PF11001 manifest: 295 Af-study, 529 Coccidioides, 144 survey proteomes. Each model is scored as a separate query, so this equals 13 single-profile runs. SLURM 30042262 (05:26). Check: PF11001 gives 1,133 Af-study and 913 Coccidioides proteins, the same as `notes/pf11001_2026-10-17.md` Sec 4.
- Per-gene tables: `SDC/tables/genes.{afum,cocci,survey}.tsv.gz`. One row per (protein, model): domain accession, envelope and HMM coordinates, scores, protein architecture, contig/start/end/rank (`gene_positions.tsv.zst`; Af from full_v072 plus 941 unplaced Af293 Swiss-Prot coordinates; Coccidioides from v5), family, family bin, strain count, number of islands, per run.
- **Family/island level.** Family labels come from `analysis/gain_loss_enrichment/<run>/family_pfam.tsv.gz` (representative only, `--cut_ga`). Islands = rows of `significant_islands.tsv`.
- **Permutations.** 1,000, seed 20261017. p = (1 + n_perm >= obs) / 1001, so the smallest possible p is 0.001. BH within run (islands) or within reference x window (windows).

## 2. Task 1: genes per genome (`SDC/tables/genes_per_genome_summary.tsv`)

Median (min-max) proteins per proteome.

| label | A. fumigatus (293) | C. immitis (169) | C. posadasii (360) |
|---|---|---|---|
| APH | 12 (8-20) | 54 (48-64) | 50 (40-64) |
| DUF3435 | 7 (2-23) | 15 (9-24) | 11 (6-22) |
| DUF3723 | 3 (0-11); 275 with >= 1 | 2 (0-7) | 2 (0-11) |
| DUF3723_N | 2 (0-9); 224 with >= 1 | 1 (0-4); 107 with >= 1 | 1 (0-6); 246 with >= 1 |
| PF11001 | 3 (3-12) | 1 (0-4) | 2 (0-6) |
| HrmA | 1 (0-6); 171 with >= 1 | 6 (4-8) | 3 (1-7) |
| NPHP3_N | 16 (11-25) | 6 (5-7) | 6 (4-8) |
| Condensation | 17 (13-22) | 7 (6-12) | 7 (6-9) |
| PP-binding | 40 (33-45) | 20 (16-25) | 20 (17-22) |
| AMP-binding | 52 (46-69) | 40 (37-43) | 40 (37-43) |
| SPOK2_N | 3 (0-6); 262 with >= 1 | 10 (8-18) | 9 (5-17) |

Family classes of the hit proteins (`SDC/tables/gene_family_class_split.tsv`; distinct families, core / soft_core / shell / cloud / singleton / nonrep_only / outgroup_only):

| label | Af | v5 | v2 |
|---|---|---|---|
| APH | 7/1/19/44/23/11/23 (128) | 40/3/38/58/43/7/197 (386) | 15/10/67/140/53/7/15 (307) |
| DUF3435 | 1/0/20/84/86/30/5 (226) | 3/1/33/78/54/9/202 (380) | 1/2/19/154/83/16/11 (286) |
| DUF3723 | 0/0/14/26/41/38/5 (124) | 0/0/3/33/15/5/43 (99) | 1/0/2/46/11/4/2 (66) |
| DUF3723_N | 0/0/9/20/25/9/2 (65) | 0/0/1/17/7/4/18 (47) | 0/0/1/24/5/1/1 (32) |
| PF11001 | 3/0/1/22/18/7/4 (55) | 1/0/1/7/4/2/14 (29) | 1/0/1/13/8/1/0 (24) |
| HrmA | 0/0/3/16/9/2/0 (30) | 4/1/5/8/8/1/16 (43) | 1/1/5/13/4/1/6 (31) |

- APH: 1,470 of 3,626 Af proteins are in core families. In v5, 11,678 of 27,437 are core.
- DUF3435, DUF3723 and DUF3723_N have no or almost no core family in any run.
- HrmA: Af has 0 core HrmA proteins. v5 has 754 of 2,168, and v2 has 354 of 1,163.

## 3. Task 2: co-occurrence

### (i) Same protein (`SDC/tables/same_protein_combinations.tsv`)

| combination | A. fumigatus (293): proteins (genomes) | Coccidioides (529) |
|---|---|---|
| DUF3723 + DUF3723_N | 519 (221) | 519 (339) |
| Condensation + PP-binding + AMP-binding | 4,556 (293) | 3,654 (529) |
| PP-binding + AMP-binding | 2,077 (293) | 2,040 (529) |
| APH + HrmA | 13 (13) | 0 |
| APH + SPOK2_N | 0 | 3 (3) |

- DUF3723_N occurs with DUF3723 on 519 of 555 Af proteins (93.5%) and on 519 of 546 Coccidioides proteins (95.1%).
- No protein carries DUF3435 together with another label.
- No protein carries PF11001 together with another label.
- APH occurs with any other tier 1 label on 0 proteins.

### (ii) Same island (`SDC/tables/island_pair_cooccurrence.tsv`, `SDC/figures/fig1_*`)

Universe: island-resident families (union of member_families; Af 11,683, v5 9,394, v2 10,934). Null: the label set of each family is shuffled among the island-resident families. Two counts:

- **islands** (`distinct_families`): islands that have label a and label b on different families.
- **family pairs** (`family_pairs`): (a-family, b-family) pairs, f != g, that share at least 1 island.

Results for the tier 1 labels, PF11001 and HrmA. Each cell is observed / permutation mean, then p.

| pair | Af islands | v5 islands | v2 islands | Af family pairs | v5 family pairs | v2 family pairs |
|---|---|---|---|---|---|---|
| APH-DUF3435 | 199/60, p 0.010 | 297/299, p 0.48 | 352/261, p 0.20 | 103/28, 0.002 | 282/103, 0.001 | 217/72, 0.001 |
| DUF3435-DUF3723 | 242/57, 0.002 | 140/195, 0.77 | 103/83, 0.26 | 129/27, 0.001 | 166/66, 0.001 | 99/22, 0.001 |
| DUF3435-DUF3723_N | 201/38, 0.001 | 54/95, 0.77 | 48/49, 0.35 | 82/17, 0.001 | 50/31, 0.050 | 45/12, 0.001 |
| DUF3435-PF11001 | 110/29, 0.024 | 51/49, 0.37 | 48/33, 0.20 | 67/13, 0.001 | 51/16, 0.001 | 61/8, 0.001 |
| APH-HrmA | 90/15, 0.018 | 133/17, 0.002 | 211/19, 0.009 | 51/6, 0.001 | 17/5, 0.009 | 22/4, 0.001 |
| PF11001-HrmA | 96/7, 0.005 | 111/3, 0.001 | 90/2, 0.003 | 41/3, 0.001 | 10/1, 0.001 | 8/0.5, 0.001 |
| DUF3723-PF11001 | 80/19, 0.033 | 94/21, 0.022 | 77/9, 0.017 | 47/7, 0.001 | 38/5, 0.001 | 32/2, 0.001 |
| DUF3723-HrmA | 42/13, 0.060 | 68/12, 0.016 | 50/6, 0.016 | 53/5, 0.001 | 20/3, 0.001 | 18/1, 0.001 |

- **Family-pair count.** Every pair of the 6 labels is above the null in all three runs (BH q < 0.05). There are two exceptions: DUF3435-DUF3723_N in v5 (50 vs 31, p 0.050, q 0.11) and HrmA-HrmA in v2 (1 vs 0.15, q 0.25).
- **Island count.** The null is wide (for example, Af APH-DUF3435 has a 95% range of 15-161), because some families sit in many overlapping islands. At this level:
  - Af: all DUF3435 pairs are above the null.
  - v5 and v2: DUF3435 is at or below the null with APH, DUF3723 and DUF3723_N (obs/exp 0.57-1.35).
  - In all three runs, PF11001-HrmA and APH-HrmA are above the null.
  - BH q < 0.05 (island count, all 11 labels): Af: APH-DUF3435, DUF3435 with DUF3435/DUF3723/DUF3723_N, PF11001-HrmA, HrmA-HrmA, and 8 pairs with SPOK2_N or NPHP3_N (14 in total). v5: APH-APH, APH-HrmA, PF11001-HrmA. v2: none.

### (iii) Same window in the reference genomes (`SDC/tables/reference_window_pairs.tsv`)

The test counts gene pairs (two different genes) with labels a and b on the same contig, within 5, 10 or 20 genes, or within 50 kb. Null: each label's gene ranks are rotated by a random offset within each contig.

Reference assemblies (from `data_dir/dna`, esl-seqstat):
- Af293: 8 chromosomes.
- RS: 463 sequences, N50 141,890 bp. 349 of them carry genes.
- Silveira: 346 sequences, N50 250,466 bp. 202 carry genes.

On these fragmented Coccidioides assemblies, pairs on different contigs cannot be counted. Many Coccidioides pairs therefore have an expected value of 0.

| pair | Af293 obs/exp, p (BH q) | RS | Silveira |
|---|---|---|---|
| APH-HrmA <= 5 genes | 3/0.04, 0.001 (0.028) | 6/1.8, 0.003 (0.083) | 3/0.86, 0.013 (0.36) |
| APH-PF11001 <= 10 genes | 3/0.12, (0.011) | 0/0.4, 1.0 | 0/0 |
| PF11001-HrmA <= 5 genes | 3/0.03, 0.001 (0.028) | 0/0 | 0/0 |
| DUF3723-PF11001 <= 10 genes | 2/0.06 (0.011) | 0 (exp 0.55 at <= 5) | 0/0 |
| DUF3435-DUF3723 <= 20 genes | 6/0.31, 0.001 (0.014) | 2/1.5, 0.56 | 0/0 |
| DUF3435-PF11001 <= 50 kb | 2/0.15, 0.007 (0.043) | 0/0.67, 1.0 | 0/0 |
| DUF3435-HrmA <= 20 genes | 3/0.12, 0.001 (0.014) | 0/0.5, 1.0 | 1/0.59, 0.59 |
| APH-DUF3435 <= 5 genes | 2/0.15, 0.005 (0.055) | 9/3.3, 0.002 (0.083) | 3/0.34, 0.007 (0.36) |

- Af293: 28 pair x window rows have q < 0.05. The pairs are APH, DUF3435, DUF3723, DUF3723_N, PF11001 and HrmA, plus DUF3723-NPHP3_N.
- RS and Silveira: no row has q < 0.05. The nominal p < 0.05 rows are APH-DUF3435, APH-HrmA and APH-SPOK2_N (both references), plus three Silveira rows with AMP-binding/PP-binding/SPOK2_N.

**All strains, descriptive** (`SDC/tables/allstrain_neighbour_rates.tsv`). This is the fraction of label-a genes with a label-b gene within 5 genes. No null was computed.

| a -> b | A. fumigatus | Coccidioides |
|---|---|---|
| PF11001 -> HrmA | 104/1,125 (9.2%) | 8/913 (0.9%) |
| PF11001 -> APH | 5.1% | 3.6% |
| PF11001 -> DUF3723 | 5.1% | 0/913 |
| PF11001 -> DUF3435 | 7.9% | 1.8% |
| DUF3723 -> DUF3435 | 18.4% | 5.4% |
| DUF3435 -> APH | 3.7% | 37.3% |
| HrmA -> APH | 50.7% | 85.2% |
| HrmA -> DUF3723 | 6.3% | 0/2,168 |

**Answer to 2(iii).** In Af293, PF11001 lies near APH, HrmA, DUF3723 and DUF3435 more often than the rotation null (q 0.011-0.043). This pattern is not detected in RS or Silveira:
- 0 PF11001-HrmA and 0 PF11001-DUF3723 pairs.
- Over all 529 Coccidioides proteomes, 0 PF11001 genes have a DUF3723 gene within 5 genes.
- What the two genera share at gene level is APH next to HrmA and APH next to DUF3435. These pairs are significant in Af293 (q <= 0.055 for <= 5 genes) and nominal only in RS/Silveira.

## 4. Task 3: island content (`SDC/tables/island_content.<af|v5|v2>.tsv.gz`, `island_tier1_count_summary.tsv`, `island_composition_summary.tsv`, `island_largest_multidomain.tsv`, `island_other_pfam_top.tsv`, `SDC/figures/fig2_island_tier1_counts.png`)

| | Af | v5 | v2 |
|---|---|---|---|
| islands (all) | 11,676 | 11,968 | 16,014 |
| with >= 1 tier 1 label | 1,642 | 2,614 | 3,144 |
| with 1 / 2 / 3 / 4 tier 1 labels | 1,016 / 390 / 176 / 60 | 1,922 / 498 / 167 / 27 | 2,520 / 491 / 110 / 23 |
| single-strain | 1,336 | 2,112 | 2,510 |
| strains: median (max) | 1 (67) | 1 (122) | 1 (247) |
| member families: median (max) | 12 (160) | 10 (79) | 8 (77) |
| with captain family (SIO `n_captain_families` > 0) | 901 | 828 | 929 |
| with DUF3435 representative | 903 | 839 | 953 |
| fraction with APH / DUF3435 / DUF3723 / DUF3723_N | 0.42 / 0.55 / 0.34 / 0.25 | 0.65 / 0.32 / 0.24 / 0.13 | 0.73 / 0.30 / 0.13 / 0.08 |
| fraction with PF11001 / HrmA / NPHP3_N / SPOK2_N | 0.08 / 0.09 / 0.16 / 0.25 | 0.05 / 0.08 / 0.0004 / 0.12 | 0.04 / 0.08 / 0 / 0.11 |
| family slots shell / cloud / outgroup_only | 0.66 / 0.30 / 0.002 | 0.31 / 0.44 / 0.17 | 0.54 / 0.43 / 0.003 |

- Af only: 590 of the 1,642 tier 1 islands contain a family with an Af293 Table S7 element member (family level). 9 overlap an element by Af293 coordinates (`with_unplaced_af293_coords/per_island.tsv.gz`).
- **Largest tier 1 multi-label islands** (`island_largest_multidomain.tsv`):
  - The 10 largest in each run are single-strain islands: Af 103-160 families, v5 68-79, v2 63-77.
  - In islands with >= 2 strains, the largest Af islands have 38-43 families and 2-9 strains. All 10 carry DUF3723+DUF3723_N, and all 10 have a family in AF293_s00042 (family level). None overlaps by coordinate.
  - The largest multi-strain v5 islands have 6-13 families. One is in 85 strains (island 170: APH, DUF3435, 6 families).
  - The largest multi-strain v2 islands have 12-29 families and 2-6 strains. 8 of 10 are APH+DUF3435.
- **Other Pfams in tier 1 islands** (families, top shared): Pkinase is first in all three runs (Af 26%, v5 25%, v2 22% of islands).
  - Af next: Ank/Ank_2/Ank_4 (17-24%), FAD_binding_8, NAD_binding_6, Kdo, HSF_DNA-bind.
  - v5/v2 next: PK_Tyr_Ser-Thr (17%), DUF6540 (16% / 12%), DUF6914, Pkinase_fungal, BTB, ABC_tran.
  - Ankyrin repeats are in 5% (v5) and 4% (v2) of the Coccidioides tier 1 islands. Transposase-like Pfams (DDE_1, HTH_Tnp_Tc5, DDE_3) are in the Af top 25 but not in the v5/v2 top 25.
- **Comparison.** Coccidioides tier 1 islands are dominated by APH (65-73%). Af islands are dominated by DUF3435 (55%) and carry DUF3723/DUF3723_N more often (34/25% vs 13-24/8-13%). Af islands are larger (median 12 vs 8-10 families). The v5 family slots include 17% outgroup_only families (C. posadasii-only families in a C. immitis-IN run).

## 5. Task 4: cross-genus best-hit families (`SDC/tables/crossgenus_pairs_by_domain.tsv`, `crossgenus_island_families_by_domain.tsv`)

Source: `cross_species_overrep/family_pairs.tsv.gz` (mmseqs, both coverages >= 0.5). The table counts Af families with the label and a Coccidioides partner that also carries the label.

| label | Af families | best hit, partner has label (v5 / v2) | RBH, partner has label (v5 / v2) | RBH identity range (median), v5 | both in islands (best hit, v5 / v2) |
|---|---|---|---|---|---|
| APH | 121 | 81 of 109 / 83 of 109 | 24 of 27 / 25 of 28 | 0.34-0.93 (0.59) | 14 / 27 |
| DUF3435 | 224 | 203 of 207 / 200 of 207 | 20 of 21 / 22 of 23 | 0.36-0.76 (0.55) | 61 / 77 |
| DUF3723 | 122 | 90 of 92 / 86 of 88 | 6 of 6 / 6 of 6 | 0.30-0.56 (0.48) | 45 / 40 |
| DUF3723_N | 62 | 49 of 50 / 48 of 50 | 3 of 3 / 5 of 5 | 0.30-0.56 (0.48) | 33 / 33 |
| PF11001 | 55 | 37 of 37 / 36 of 36 | 4 of 4 / 3 of 3 | 0.41-0.56 (0.53) | 7 / 5 |

- The best-hit partners collapse onto fewer Coccidioides families: DUF3435 has 203 best-hit pairs onto 81 v5 families, and APH has 81 pairs onto 41 families.
- Coverage ranges, best-hit pairs, v5: DUF3435 qcov 0.50-1.00, tcov 0.50-1.00; DUF3723 qcov 0.51-1.00, tcov 0.50-0.99; PF11001 qcov 0.50-0.94, tcov 0.51-0.88.
- Partner classes (best hit, v5): DUF3435 partners are outgroup_only 124, cloud 36, singleton 31, shell 12. APH partners include 30 core families.
- Cross-genus island families (111 rows): DUF3435 4, APH 3, DUF3723 3, DUF3723_N 3 (the same 3 families), PF11001 0, HrmA 0. In every case the v5 and v2 partners also carry the label.
- Of the 4 DUF3435 families, 1 is in AF293_e00039. Of the 3 DUF3723 families, 2 are in AF293_s00031/s00037 (family level).

## 6. Task 5: DUF3435 domain tree (`SDC/tree/`, `SDC/figures/fig3_duf3435_tree.png`)

- **Input.** 89 DUF3435 envelopes: Af293 10 (8 proteins), RS 18, Silveira 8, and 53 from 12 other Dikarya proteomes (A. megalospora, A. nidulans, C. globosum, C. rosea, C. higginsianum, C. ligniaria, C. parasitica, M. anisopliae, V. dahliae, Z. tritici, T. caries, V. palustris).
- **Method.** hmmalign --trim → FastTree 2.2.0 `-lg -gamma` (log-likelihood -39,125.66) → midpoint root. FastTree reports SH-like local support.
- **Separate groups.** Af293 copies form 8 maximal Af293-only clades (sizes 3, 1×7). Coccidioides copies form 12 maximal Coccidioides-only clades (sizes 5, 4, 3, 3, 2, 2, 2, 1×5).
- **Mixed clades.** There are 5 smallest clades that hold both Af293 and Coccidioides tips. Sizes and support: 5 tips, 1.0; 4, 0.946; 2, 0.981; 3, 0.211; and 8, 0.507 (with 4 other-species tips).
- **Nearest tip.** 2 Af293 domains have a Coccidioides domain as the nearest tip (A4DA40 → RS 005222; Q4WMN5 → RS 002663, distance 0.59).
- **Summary.** The Af293 and Coccidioides copies do not fall in one Af293 clade and one Coccidioides clade. They are spread over several clades, and 3 clades with support >= 0.94 contain both. This is a descriptive result; it says nothing about orthology.
- Some envelopes are short or low-scoring (domain score 28-42). Branch lengths reach 3 substitutions/site, so deep nodes are poorly resolved.

## 7. Summary answer

- **Shared.** In both genera:
  - The DUF3723+DUF3723_N fusion is the same, on 94-95% of DUF3723_N proteins.
  - The labelled families share islands with each other more than random island families do (family-pair count; all six labels, all runs).
  - APH lies next to HrmA, and APH next to DUF3435, at gene level.
  - Best-hit partners almost always carry the same domain.
- **Different.**
  - In A. fumigatus, DUF3435 islands also carry DUF3723/DUF3723_N/PF11001 (island count above null; Af293 windows q <= 0.043). In Coccidioides, DUF3435 co-occurs with these at or below the island-count null.
  - Coccidioides has 0 PF11001-DUF3723 and 0 HrmA-DUF3723 neighbours within 5 genes over 529 proteomes.
  - The Af293 PF11001-HrmA-APH-DUF3723 neighbourhood is not seen in RS/Silveira.
  - Coccidioides tier 1 islands are mostly APH islands. Af tier 1 islands are mostly DUF3435 islands.
  - Copy numbers differ: APH ~50 vs 12, HrmA 3-6 vs 1, NPHP3_N 6 vs 16.

## 8. Caveats

- Island co-occurrence uses family-representative labels only (`family_pfam.tsv.gz`). A family "has" a domain when its representative has it.
- v5 and v2 share genomes. They are not independent replicates.
- Starship coordinates exist only for Af293 (Table S7). Coccidioides has captain (DUF3435) annotation only.
- Island calls and family clustering were not rechecked. Islands overlap and share families, so island counts are not independent. This is why the island-count null is wide and why two count types are reported. Island co-membership is not physical adjacency: v5 PF11001-HrmA share 111 islands, but only 8 of 913 PF11001 genes have an HrmA gene within 5 genes.
- Label shuffles among island-resident families do not control for family frequency class, protein length or copy number.
- RS and Silveira assemblies are fragmented (N50 142 kb and 250 kb). Window tests on them have low power, and between-contig pairs cannot be counted.
- The rotation null moves each label independently. It splits fused domains (DUF3723/DUF3723_N), and the same-label diagonal is not tested.
- The Coccidioides gene positions are from v5 for all 529 proteomes. The HrmA label merges 3 Pfam models.
- The hmmsearch job ran on partition short (node c09, not milan) after AssocGrpMemLimit on epyc.

## 9. Not verified

- No check of whether the Af293 / Coccidioides DUF3435 clades are orthologous. No outgroup rooting (midpoint only). No bootstrap: FastTree SH-like support only.
- The tree uses one proteome per species for the survey species. Taxon sampling is small.
- I did not recompute the tier lists. They are taken from `notes/cross_species_overrep_2026-10-17.md`.
- The Pfam names and descriptions were not checked against the literature.
- Starship overlap for islands is taken from the `starship_island_overlap` tables as they are.

## 10. Files (all new, uncommitted; `SDC/`, 30 MB total)

- Scripts (`scripts/`):
  - `run_hmmsearch.sbatch`, `run_genes.sbatch`
  - `sdc_genes.py`, `sdc_islands.py`, `sdc_gene_analysis.py`, `sdc_duf3435_tree.py`, `sdc_figures.py`
- `PROVENANCE.yaml`
- `hmm/shared_domains.hmm` (1.5 MB) and `domain_names.txt`.
- `hmmsearch/`:
  - `manifest.tsv` (copy of the pf11001 manifest)
  - `all_proteomes.shared_domains.domtbl.idx.zst` (4.3 MB)
  - `domtbl_header.txt`, `seq_counts.pf11001_run.tsv`
- Large per-gene tables (`tables/`): `genes.afum.tsv.gz` (2.5 MB), `genes.cocci.tsv.gz` (4.7 MB), `genes.survey.tsv.gz` (0.36 MB), `reference_gene_order.tsv.gz` (0.37 MB).
- Small tables (`tables/`, < 100 kB each): the other 18 TSV files named in Sec 2-5.
- `tree/` (0.5 MB): env FASTA, Stockholm/afa alignments, Newick trees, tip/clade/summary tables, FastTree log, DUF3435.hmm.
- `figures/`: fig1 (2 versions), fig2, fig3, as PNG and PDF (1.1 MB). All PNGs were viewed.
- `logs/`.

**Suggested commit split:**
- PR 1: `scripts/` and `PROVENANCE.yaml` (about 60 kB).
- PR 2: small tables, `tree/`, figures and this note (about 2 MB).
- Decide separately on the large regenerable files: `genes.*.tsv.gz` (7.6 MB), `all_proteomes.shared_domains.domtbl.idx.zst` (4.3 MB), `hmm/shared_domains.hmm` (1.5 MB, hmmfetch output). My suggestion is to keep them untracked or as release data, because the scripts regenerate them in about 7 minutes of SLURM time.
