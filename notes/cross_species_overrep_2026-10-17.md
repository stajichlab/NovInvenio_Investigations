# Domains and families over-represented in both A. fumigatus and Coccidioides (2026-10-17)

Question: are any Pfam domains or gene families over-represented in C. immitis AND C. posadasii AND ALSO in A. fumigatus?

Runs: A. fumigatus `full_v072` (Af), C. immitis `immitis_in_posadasii_out_v5` (v5), C. posadasii `posadasii_in_immitis_out_v2` (v2). v5 and v2 share genomes, so they are not independent.

Directory: `studies/fungi/coccidioides_pangenome/analysis/cross_species_overrep/` (abbreviated `XS/` below). Nothing was committed, pushed or published.

## 1. Definitions and tests

- **A, island.** The pipeline test, recomputed. Background = families with bin shell or cloud. A family has a domain when one row of `pfam.domtblout` (hmmscan `--domE 1e-3`, background representatives only) has i-Evalue <= 1e-3. Island set = union of `significant_islands.tsv` member_families, intersected with the background. Test = one-sided Fisher. BH is applied over the domains that occur in at least one island family. My counts match `island_pfam_enrichment.tsv` for every domain in all three runs: 1,977 / 928 / 1,325 domains, 0 count mismatches. The largest relative difference in p is 4.8e-4, which comes from 4-digit rounding in the published file (`XS/logs/check_island_recompute.tsv`). The clan level is new: a family has a clan when one of its counted domains maps to that clan (Pfam-A.hmm.dat `CL`). The test is the same.
- **B, gain/loss.** Read from `analysis/gain_loss_enrichment/<run>/{shell_soft_core,all_accessory}/enrichment_domains.tsv.gz` (`bin/pangenome_gain_loss_enrichment.py`). Test = one-sided hypergeometric. BH is applied within each test set and vocabulary. A term is tested when K >= 5. Annotation = `hmmsearch --cut_ga` on all family representatives. Test sets: gain and loss. Vocabularies: Pfam and Pfam clan.
- **C, accessory vs core.** New, in `XS/scripts/xsp_domains.py`. Annotation = `family_pfam.tsv.gz` (`--cut_ga`, all representatives). The table exists for all three runs. Test = one-sided hypergeometric. BH is applied over the tested terms. A term is tested when K >= 5.
  - C1: test set = soft_core+shell+cloud+singleton. Universe = core+soft_core+shell+cloud+singleton families with >= 1 Pfam.
  - C2: test set = shell+cloud. Universe = core+shell+cloud families with >= 1 Pfam.
- **Tiers.** Only terms tested in >= 2 runs are included. The tiers are exclusive.
  - Tier 1: q < 0.05 in Af, v5 and v2.
  - Tier 2: q < 0.05 in Af and in exactly one of v5/v2.
  - Tier 3: p < 0.05 in all three runs, and the term is not in tier 1 or tier 2.
- **Overlap test.** Universe = terms tested in Af and in the Coccidioides run(s). Hypergeometric P(overlap >= observed). Expected overlap = n_Af_sig × n_Cocci_sig / universe.
- **Rank correlation.** Spearman of -log10 p over the terms tested in both runs.

## 2. Pfam results (TASK 1)

Counts are from `XS/summary_domains.tsv` lines 2-13. The tier lists are in `XS/tiers.tsv`.

| definition | terms q<0.05 Af / v5 / v2 | tested >=2 runs | tier 1 | tier 2 | tier 3 | overlap Af×v5: obs / exp / p | Af×v2: obs / exp / p | Af×(v5 and v2): obs / exp / p | Spearman Af~v5 rho (n, p) | Af~v2 rho (n, p) |
|---|---|---|---|---|---|---|---|---|---|---|
| A Pfam | 21 / 12 / 17 | 1,063 | 4 | 1 | 12 | 4 / 0.26 / 6.0e-5 | 5 / 0.36 / 1.1e-5 | 4 / 0.22 / 2.2e-5 | 0.236 (549, 2.3e-8) | 0.222 (715, 2.0e-9) |
| A clan | 9 / 6 / 8 | 250 | 3 | 1 | 2 | 3 / 0.27 / 1.1e-3 | 4 / 0.35 / 1.2e-4 | 3 / 0.26 / 8.6e-4 | 0.292 (179, 7.5e-5) | 0.247 (204, 3.6e-4) |
| B shell+soft-core gain, Pfam | 8 / 2 / 0 | 103 | 0 | 1 | 0 | 1 / 0.09 / 0.091 | 0 / 0 / - | 0 | 0.215 (77, 0.061) | 0.263 (83, 0.016) |
| B shell+soft-core gain, clan | 3 / 2 / 1 | 94 | 1 | 1 | 0 | 2 / 0.08 / 1.0e-3 | 1 / 0.04 / 0.038 | 1 / 0.05 / 0.047 | 0.182 (77, 0.11) | 0.346 (78, 1.9e-3) |
| B shell+soft-core loss, Pfam / clan | 0 / 0 / 0 | 103 / 94 | 0 | 0 | 0 | - | - | - | 0.318 / 0.082 | 0.235 / 0.201 |
| B all accessory gain, Pfam | 10 / 5 / 7 | 555 | 3 | 0 | 4 | 3 / 0.13 / 1.3e-4 | 3 / 0.15 / 2.1e-4 | 3 / 0.16 / 2.4e-4 | 0.109 (383, 0.033) | 0.225 (406, 4e-6) |
| B all accessory gain, clan | 4 / 4 / 4 | 231 | 3 | 0 | 1 | 3 / 0.08 / 1.3e-5 | 3 / 0.08 / 1.2e-5 | 3 / 0.09 / 1.8e-5 | 0.153 (196, 0.032) | 0.349 (201, 3.7e-7) |
| B all accessory loss, Pfam / clan | 0 / 0 / 0 (clan: Af 1) | 555 / 231 | 0 | 0 | 0 | - | - | - | 0.090 / 0.054 | 0.147 / 0.203 |
| C1 accessory vs all | 38 / 21 / 27 | 959 | 9 | 5 | 15 | 10 / 0.73 / 8.9e-11 | 13 / 0.92 / 3.4e-14 | 9 / 0.69 / 1.1e-9 | 0.351 (765, 1.2e-23) | 0.429 (778, 3.1e-36) |
| C2 shell+cloud vs core | 52 / 28 / 30 | 713 | 9 | 4 | 16 | 9 / 1.50 / 3.4e-6 | 13 / 1.62 / 1.2e-10 | 9 / 1.21 / 2.9e-7 | 0.381 (566, 5.2e-21) | 0.430 (608, 9.0e-29) |

So the overlap of significant terms between A. fumigatus and each Coccidioides run is larger than chance under every definition that has significant terms in both genera. The overlaps are small in absolute terms (1-13 terms). The rank correlations are positive but weak to moderate (rho 0.11-0.43). No term is significant for LOSS in more than one run.

### Tier 1 and tier 2 terms

The format is k/K (q). K = families with the term in the test universe. The number in [] is the number of families with the term among all families of that run (`--cut_ga`, all bins), which gives copy-number context. Rows are `XS/tiers.tsv` lines 2-80, as listed.

**A, island (Pfam)**, tiers.tsv lines 2-6:

| tier | domain | DE | Af | v5 | v2 |
|---|---|---|---|---|---|
| 1 | APH PF01636 | Phosphotransferase enzyme family | 67/72 (1.3e-10) [121] | 86/114 (5.8e-7) [365] | 180/244 (2.3e-7) [292] |
| 1 | DUF3435 PF11917 | Protein of unknown function (DUF3435) | 91/106 (3.8e-10) [224] | 98/113 (4.7e-15) [374] | 147/181 (6.6e-11) [281] |
| 1 | DUF3723 PF12520 | DUF3723 central domain | 37/39 (3.2e-6) [122] | 36/36 (1.5e-9) [99] | 40/49 (9.7e-3) [66] |
| 1 | DUF3723_N PF29418 | DUF3723 N-terminal domain | 27/28 (8.8e-5) [62] | 18/18 (2.8e-4) [46] | 23/24 (1.6e-3) [31] |
| 2 | AFUB_07903_YDR124W_hel PF11001 | YDR124W-like, helical bundle domain | 21/24 (0.036) [55] | 8/8 (0.138) [29] | 14/14 (0.018) [24] |

**A, island (clan)**, lines 19-22:
- Tier 1: PKinase CL0016, 226/308 / 281/401 / 460/631, q 1.3e-11 / 1.7e-17 / 2.1e-19.
- Tier 1: DNA-mend CL0382 (contains DUF3435), 91/108 / 98/113 / 147/181.
- Tier 1: ParB_Srx CL0248 (contains DUF3723_N), 27/28 / 18/18 / 23/24.
- Tier 2: POZ CL0033 (contains BTB, HAC6_N), Af 46/51 q 1e-6, v2 28/32 q 3.5e-3, v5 12/19 q 0.95.

**B, gained families**, lines 25-37:
- Shell+soft-core scope:
  - Tier 2: DUF3435, Af 20/20 q 4.2e-3, v5 30/33 q 2e-7, v2 13/21 q 0.95.
  - Tier 1: clan PKinase, 62/74 / 79/146 / 101/178, q 4.6e-3 / 0.026 / 0.032.
  - Tier 2: clan DNA-mend.
- All-accessory scope:
  - Tier 1: AMP-binding PF00501, 209/217 / 168/177 / 310/322, q 2.4e-5 / 7.9e-9 / 5.7e-13.
  - Tier 1: DUF3435, 181/188 / 151/164 / 231/253.
  - Tier 1: Condensation PF00668, 123/128 / 147/147 / 270/274.
  - Tier 1 clans: ANL CL0378, DNA-mend CL0382, CoA-acyltrans CL0149.
  - The fold values are small (1.10-1.30), because most all-accessory families are "gain" (rule_freq): Af 26,267/29,137 (0.90), v5 15,534/19,201 (0.81), v2 21,474/24,579 (0.87) (`analysis/gain_loss/<run>/family_gain_loss.tsv.gz`, column 10).

**C1, accessory vs all**, lines 39-52:
- Tier 1 (9):
  - DUF3435: 188/189, 164/167, 253/254.
  - DUF3723: 80/80, 51/51, 59/60.
  - NPHP3_N PF24883 (Nephrocystin 3, N-terminal): 128/136, 22/24, 28/30.
  - Condensation.
  - DUF3723_N.
  - PP-binding PF00550.
  - APH.
  - AMP-binding.
  - SPOK2_N PF27653 (Spore killer 2 N-terminal domain): 34/34, 48/55, 85/89.
- Tier 2 (5): Fung_rhodopsin PF20684, zf-CCHC PF00098, HNH_2 PF13391, PF11001 (Af 41/44 q 0.032; v2 22/23 q 0.021; v5 12/13 q 0.25), Chitin_bind_1 PF00187.

**C2, shell+cloud vs core**, lines 68-80:
- Tier 1 (9): DUF3435, NPHP3_N, Fung_rhodopsin, DUF3723, APH, DUF3723_N, SPOK2_N, zf-CCHC, DUF7104 PF23397.
- Tier 2 (4): HNH_2, AA_permease_2 PF13520, Chitin_bind_1, PF11001 (Af 23/26 q 0.038; v2 14/15 q 0.024; v5 8/9 q 0.097).

Tier 3 (nominal p only) adds, among others: Pkinase, GloB_C, Glyco_hydro_18 (A); PP-binding, ABC_tran, SPOK2_N (B); Ank, Ank_2, Ank_3, Ank_5, HAC6_N/HAC6_helical, HrmA, Pkinase_fungal, NB-ARC (C). See `XS/tiers.tsv` (tier 3 rows).

**Terms found under all three definitions.** DUF3435 is tier 1 in A, B (all accessory) and C1/C2. APH, DUF3723 and DUF3723_N are tier 1 in A and in C1/C2. These are the only domains that are tier 1 in more than one definition, apart from AMP-binding/Condensation (B and C1).

**PF11001.** It is tier 2 in A, C1 and C2: q < 0.05 in Af and v2, not in v5 (v5 q 0.138 / 0.25 / 0.097). It is not significant in B (Af all-accessory gain q 0.094; v5 q 1.0; v2 q 1.0; not tested in the shell+soft-core scope because K < 5). This matches `notes/pf11001_2026-10-17.md` Sec 3-4.

**Generic names and copy number.** PKinase (clan) has 648 / 1,244 / 979 families with the clan in Af / v5 / v2. AMP-binding has 379 / 643 / 433 families, and Condensation has 222 / 552 / 354. DUF3435 has 224 / 374 / 281. Ankyrin-repeat Pfams (tier 3 in C) are large multi-family groups. A significant result for such a term can reflect many paralogous families. It does not show that one gene is shared. The best-hit components in Sec 3 show the same: the largest component in v5 has 134 Af and 528 Coccidioides families with AMP-binding/Condensation/PP-binding (`XS/family_besthit_components_top15.tsv` line 2).

## 3. Families (TASK 2)

- **Search.** mmseqs 18.8cc5c `easy-search -e 1e-5 -c 0.5 --cov-mode 0` (query AND target coverage >= 50%), default sensitivity. Searches: Af reps vs v5 reps and v2 reps, both directions. SLURM 30041523 (stajichlab, milan), COMPLETED in 00:36:53 (`sacct`).
- **Best hit and RBH.** Best hit = the highest bits. RBH = mutual best hit. v5 and v2 are handled separately. I call these pairs "best-hit homologous families". I make no orthology claim.
- **Partners** (`XS/family_partner_coverage.tsv` lines 9, 20; `XS/logs/families.log`):
  - v5: 22,981 of 47,954 Af families (47.9%) have a best hit, and 7,518 (15.7%) have an RBH partner.
  - v2: 22,331 (46.6%) have a best hit, and 7,267 (15.2%) have an RBH partner.
  - RBH coverage by Af class (v5): core 55.9%, soft_core 31.8%, shell 17.0%, cloud 8.1%, singleton 5.2%.
  - Af shell+cloud island families: 10.7% have an RBH partner (v5), 10.4% (v2). Shell+cloud non-island families: 10.8% (v5), 10.2% (v2).
  - So RBH is restrictive for accessory families: most have no partner.
- **(i) Island concordance** (`XS/family_tests.tsv`; Fisher; OR with 95% CI):

| set | run | Af island: partner in island | Af non-island: partner in island | OR (95% CI) | p (two-sided) | line |
|---|---|---|---|---|---|---|
| all RBH pairs | v5 | 188/1,057 (0.178) | 471/6,461 (0.073) | 2.75 (2.28-3.32) | 2.7e-24 | 2 |
| all RBH pairs | v2 | 225/1,011 (0.223) | 617/6,256 (0.099) | 2.62 (2.20-3.11) | 9.5e-26 | 10 |
| both shell or cloud | v5 | 120/185 (0.649) | 63/172 (0.366) | 3.18 (2.03-5.04) | 1.0e-7 | 5 |
| both shell or cloud | v2 | 154/242 (0.636) | 100/211 (0.474) | 1.94 (1.31-2.88) | 6.2e-4 | 13 |
| Af best hit (non-reciprocal), both shell or cloud | v5 | 570/869 (0.656) | 306/687 (0.445) | 2.37 (1.92-2.93) | 9.4e-17 | 9 |
| same | v2 | 838/1,242 (0.675) | 508/916 (0.555) | 1.67 (1.39-2.00) | 1.4e-8 | 17 |
| RBH in both runs; partners island in both | v5+v2 | 111/967 (0.115) | 193/6,086 (0.032) | 3.96 (3.07-5.08) | 9.8e-25 | 18 |
| same, all three shell or cloud | v5+v2 | 61/102 (0.598) | 14/63 (0.222) | 5.15 (2.43-11.5) | 2.5e-6 | 19 |

  The "all pairs" rows mix frequency classes: core families are rarely in islands in either genus. The shell/cloud rows remove most of this effect, and the association remains. Within shell+cloud, carrier frequency was not controlled.
  - Gain direction (both shell or cloud): partner gain 140/242 vs 28/115, OR 4.25, p 2.4e-9 (v5, line 6). For v2: 207/309 vs 46/144, OR 4.31, p 2.8e-12 (line 14).
  - Loss direction: v5 OR 1.73, p 0.11 (line 7). v2 OR 6.50, p 4.0e-13 (line 15).
  - The p values printed as 0 in the all-pairs gain/loss rows are below float precision.
- **(ii) Cross-genus island families** (`XS/cross_genus_island_families.tsv`): 111 Af island families whose RBH partners are island families in both v5 and v2.
  - 61 of these are shell or cloud in all three runs.
  - Af classes: cloud 52, shell 34, singleton 13, outgroup_only 7, nonrep_only 5.
  - 43 of 111 Af representatives have no `--cut_ga` Pfam.
  - Most frequent Pfams (number of families): DUF3435 4, NAD_binding_6 4, FAD_binding_8 4, Pkinase 4, DUF3723 3, DUF3723_N 3, PK_Tyr_Ser-Thr 3, PGAP1 3, DUF676 3, APH 3.
  - Starship: 8 Af families lie in an Af293 Table S7 element (AF293_e00029, e00039, s00031, s00037, s00042).
  - Captain: 4 have a DUF3435 captain hit (`captain` = yes). 100 are not listed in `starship_island_overlap/full_v072/annotated_families.tsv.gz`. That table holds only element members and captain hits (380 rows).
  - PF11001 has no family in this list.
- **(iii) Paralog check** (`XS/logs/families.log`):
  - RBH pairs are 1:1. They involve 7,518 distinct Af and 7,518 distinct v5 families, and 7,267 and 7,267 for v2.
  - Every best-hit graph component holds at most one RBH pair (maximum 1 per component). This follows from one best hit per node.
  - The largest components have 662 nodes (v5: 134 Af / 528 Coccidioides, AMP-binding/Condensation/PP-binding) and 477 nodes (v2). Each contributes one RBH pair.
  - The largest Af->Coccidioides best-hit in-degree is 68 Af families to one Coccidioides family.
  - So big paralog groups are under-represented by RBH, not over-represented. The non-reciprocal best-hit rows above count each Af family once and can reuse one Coccidioides family many times. They give the same direction of result.

## 4. Figures (all viewed)

- `XS/figures/fig1_tiers_A_B.png`: -log10 q heat plot of tier 1 and tier 2 terms, definitions A and B.
- `XS/figures/fig1_tiers_C.png`: the same plot for definition C.
- `XS/figures/fig2_upset_A.png`: UpSet-style counts of q < 0.05 Pfam domains for definition A.
- `XS/figures/fig3_family_island.png`: fraction of Af island vs non-island families whose RBH partner is an island family.
- PDF copies exist for all figures.

## 5. Caveats

- v5 and v2 share genomes. They are not two independent replicates of Coccidioides. The tier 1 requirement ("both runs") is therefore weaker than it looks.
- Annotation is by family representative only. A family "has" a domain when its representative has it.
- The pipeline background (A) excludes core, soft_core, singleton, nonrep_only and outgroup_only families. Definition C uses a different annotation (`--cut_ga`) and includes core.
- **Selection effect.** Per run, 103-1,063 terms were compared across runs, in 12 definition/vocabulary combinations. BH is applied within each run and definition only, not across definitions or runs. The tier lists are not corrected for the number of definitions.
- The overlap hypergeometric test assumes that terms are independent. They are not: DUF3723 and DUF3723_N sit on the same proteins, and Pfam and clan levels are nested.
- Islands include families with bin outgroup_only and nonrep_only (seen in Sec 3 (ii)).
- RBH requires 50% coverage of both sequences. This drops multi-domain fusions and fragments.

## 6. Not verified

- The island calls and the gain/loss calls themselves were not re-examined.
- I made no phylogeny, so "best-hit homologous" is not orthology, and shared ancestry vs independent acquisition is unknown.
- No mapping between v5 and v2 families was done, so an Af family can pair with unrelated v5 and v2 families.
- Within shell+cloud, the island association was not controlled for carrier frequency or protein length.
- Starship status exists only for Af293 elements. Coccidioides partners have no Starship annotation.
- mmseqs used default sensitivity (-s 5.7) and default --max-seqs 300. I did not test sensitivity to these settings.

## 7. Files (all new, uncommitted)

- **Scripts:** `XS/scripts/run_mmseqs_search.sbatch`, `xsp_domains.py`, `xsp_families.py`, `xsp_figures.py`, `xsp_summary.py`, `write_provenance.py`.
- **Outputs:**
  - `XS/domains_<definition>.tsv` (12 files)
  - `XS/tiers.tsv`, `XS/summary_domains.tsv`, `XS/summary.tsv`
  - `XS/family_pairs.tsv.gz`, `XS/family_tests.tsv`, `XS/family_partner_coverage.tsv`, `XS/family_besthit_components_top15.tsv`
  - `XS/cross_genus_island_families.tsv`
  - `XS/PROVENANCE.yaml`
  - `XS/logs/`
  - `XS/mmseqs/*.m8.zst` (21 MB total)
  - `XS/figures/*`
- **Suggested commit split:**
  - PR 1 (code): `XS/scripts/`.
  - PR 2 (outputs and figures): the TSVs, `PROVENANCE.yaml`, figures and this note.
  - Decide whether `XS/mmseqs/*.m8.zst` (about 21 MB, regenerable by job script) should be committed or kept as untracked/release data.
