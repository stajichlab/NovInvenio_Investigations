# HET and fungal NLR genes in the Coccidioides pangenome: polymorphic families and loci

Date: 2026-09-25. Read-only analysis. Nothing in the pipeline or `results/` was changed. Nothing was committed.

Run analysed: `results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome/` (529 strains: 169 *C. immitis* (Ci), 360 *C. posadasii* (Cp); species labels from `samplesheet.with_clades.csv`). This is the same run as the HrmA report (`2026-09-24-hrmA-presence-and-neighbourhoods.md`). Family IDs from other runs were not used.

All scripts and tables are in `analysis/het_nlr_2026-09-25/`.

## Plain-language summary

- The Pfam panel found 6,462 HET/NLR-class proteins in 86 tier-1 families. Every strain has 10-16 of them (median 12 in both species). Coccidioides has few classical HET-domain genes: 807 proteins, 1-3 per strain.
- 39 of 86 families are polymorphic at the family level (30 within a species, 35 between species, 26 both). But the family level misleads here, as it did for HrmA.
- At the locus level, the 11 main HET/NLR sites are occupied in 99-100% of resolved strains in both species. Presence/absence of the site is not polymorphic. Only two small loci (L012, L013) show real presence/absence, at low frequency (1-9%).
- The polymorphism is in *which* sequence occupies a fixed site. Two kinds exist:
  1. **Divergent allele classes (the HET-like pattern).** At L005, three NLR allele classes share only 54-65% amino-acid identity. They sit between the same conserved core genes. All three occur in both species at 20-57%. Each class carries its own N-terminal domain (PGAP1-like, Abhydrolase_6/DUF676, or a C-terminal ANAPC4_WD40 hit). Each class is in complete linkage with its own divergent fructosamine-kinase neighbour gene (61-65% identity between the kinase classes). L006 shows the same pattern with two classes (77% identity), both species. L002 has four divergent Goodbye-domain classes (74-92% identity) next to an NPHP3_N gene. Its classes are species-specific.
  2. **Species-specific gene-model forms of one gene (95-100% identity).** At L004 (PNP_UDP-NPHP3_N-NACHT-WHD-Ank), the Cp gene is fused. In Ci, the Ank part is a separate, adjacent gene. At L008 (Patatin-NB-ARC-TPR), Cp has the fused gene. Most Ci have NB-ARC-TPR only, with a separate Patatin gene next to it, and 18 Ci carry the Cp form. L001a (short NB-ARC fragments) and L007 (HET-Het6_barrel) differ only in gene length. I cannot tell from these data whether these are real gene structures or annotation differences.
- One HET-domain family (H002, 835 aa) sits in an accessory block next to an APH gene in 272 of 275 copies. It is in 84% of Ci and 37% of Cp (BH q = 2.2e-24). It has no conserved core flanks, so its site cannot be scored as occupied or empty.
- The pipeline's own island Pfam enrichment shows no NLR/HET domain at FDR < 0.05 (best: Patatin, p = 0.013, q = 0.40). DUF676, which is on the L005-B and L006-B allele classes, is enriched (q = 0.029).
- HET/NLR genes are not near HrmA genes more often than other genes. They are depleted (OR 0.49, p = 3.8e-4). The one exception: the L001a NB-ARC site lies about 27 kb (10-11 genes) from HrmA locus L07 in Ci strains.
- HET/NLR copy number per strain does not track assembly contig count in a consistent way. The largest correlation is for HET-domain genes in Cp (rho = 0.32, p = 7e-10).

## Definitions

- **HET/NLR-class protein**: a protein with one of these domain combinations at Pfam GA (sequence and domain GA both passed, `--cut_ga`):
  - `NLR_NOD`: NACHT or NB-ARC.
  - `NLR_like_AAA`: AAA_16/AAA_22 plus an NLR-accessory or HET N-terminal domain, and no NACHT/NB-ARC. (No protein fell in this class.)
  - `NLR_accessory_only`: NPHP3_N, NPHP3_hel, TPR_NPHP3, NACHT_N, NACHT_sigma, Beta-prop_NWD2_C or an NLR-type WHD (WHD_GPIID, WHD_NACHT_C, WHD_NWD1, WHD_APAF1, WHD_CED4), with no NOD hit.
  - `HET_domain`: HET (PF06985) with no NOD hit.
  - `HeLo_SesA_Goodbye_HETs`: HeLo, HET-S, HET-s_218-289, SesA, Goodbye or Het-6_barrel without HET or NOD. (Only Goodbye was found.)
  - `Het-C`: Het-C (PF07217). This is the Neurospora het-c type. It is not an NLR. It is listed separately.
  - Effector domains (PNP_UDP_1, Patatin, CHAT, Peptidase_C14, TIR) and repeats (Ank, WD40, TPR) were recorded, but they do not make a protein HET/NLR-class alone. They are too common. For example, 3,172 of 3,567 Patatin proteins and 52,087 of 52,087 WD40 proteins carry no NOD or NLR-accessory domain.
- **Polymorphic**: a family, site or allele class carried by >= 5% and <= 95% of strains within at least one species ("within"), or significantly different between Ci and Cp (two-sided Fisher exact test, BH q < 0.05; "between"). Both kinds are reported separately. At the locus level a species frequency is used only when >= 10 strains of that species are resolved.
- **Resolved** strain at a site = occupied + empty_site.
- **Allele class**: families at one site whose representatives share >= 95% amino-acid identity are one class (length or gene-model variants). Families below that are separate allele classes. This 95% cut-off is my choice. It is not a standard.
- Frequency bins in `frequency_table.tsv` count only the 133 Ci ingroup representatives. I did not use the bins for species frequencies. All per-species counts below are my own.

## Methods

Tools: HMMER 3.4 and DIAMOND from the NII pixi environment. Python ran with the NII pixi python (pandas, scipy, statsmodels). The node had 6 CPUs; `--cpu 6` was used. Temporary FASTA went to `$SCRATCH/het_nlr` and was deleted at the end.

1. **Models.** 44 Pfam-A 38.2 models were fetched by name with `hmmfetch -f` (`models.txt`, `het_nlr_models.hmm`). Accessions and GA were read from the file. Examples: NACHT PF05729.19 (GA 26.7), NB-ARC PF00931.29 (23.5), HET PF06985.18 (22.5), HeLo PF14479.12 (23.2), HET-s_218-289 PF11558.14 (32.2), NPHP3_N PF24883.3 (24.8), WHD_GPIID PF22939.3 (27.0), Goodbye PF17109.11 (24.8), SesA PF17107.12 (26.3), PNP_UDP_1 PF01048.27 (25.1), Patatin PF01734.28 (27.7), CHAT PF12770.13 (24.0), Het-C PF07217.17 (27.0), Ank PF00023.37 (21.1), WD40 PF00400.39 (27.0), TPR_1 PF00515.35 (27.8). The full list with GA is in `models.txt` and the HMM file.
2. **Reuse check.** The run's `pfam.domtblout` covers only 6,117 representatives and was run without GA cutoffs. It was not used for the census. The run has no Leiden modules (module files hold only a header), as in the HrmA report.
3. **Census.** All 529 proteomes (4,504,703 proteins, `Short|` prefix) were searched with `hmmsearch --cut_ga` (`run_allprot.sh`). Architecture = domain labels ordered by envelope start (`het_step1_census.py`). All 102,474 hit proteins were in `tier1_cluster.tsv`.
4. **Families.** A family is HET/NLR if >= 1 member is HET/NLR-class. Family class = highest-priority class carried by >= 20% of classified members. Presence per species came from `presence_matrix.rescued.tsv`. Fisher tests Ci vs Cp: (a) present + genome_only, (b) present only; BH over 86 families.
5. **Locations and neighbourhoods.** Gene positions: `gene_positions.tsv.zst`. Contig lengths and telomere ends were reused from `hrmA_2026-09-24/contig_lengths.tsv.gz` and `contig_telomere_ends.tsv.gz`. +/-10-gene neighbourhoods were built for every HET/NLR gene.
6. **Loci** (`het_step2_loci.py`, same rules as `hrmA_step3_loci.py`): markers are core/soft_core neighbour families that flank >= 50% of a family's anchors. Families sharing >= 2 markers are joined (union-find). Up to 6 markers per locus. Per strain: `occupied` = a HET/NLR gene within 12 genes of a marker; `empty_site` = a marker gene with >= 12 genes on both sides and no HET/NLR gene; `unresolved`; `no_marker`; `+rescue` when a TBLASTN rescue position falls in the window. Step 2 was changed from the HrmA version in one way: it records all HET/NLR genes found in any marker window, not only the first window.
7. **Own-family re-scoring** (`het_step5_ownfamily_loci.py`): the same scoring, but only genes of the locus's own families count as occupants. Step 2 joined two physical sites into L001 through markers about 20 genes apart. Step 5 splits them into L001a (NB-ARC fragments) and L001b (DUF7104-NPHP3_N). The results tables below use step 5.
8. **Allele checks.** DIAMOND `--ultra-sensitive` all-vs-all of the 86 family representatives (`het_rep_pairwise.tsv`). Allele classes: `het_step6_allele_classes.py`. Full Pfam-A (`hmmscan --cut_ga`) on the 86 representatives and on 492 neighbour-family representatives not already scanned in the HrmA work (`run_small_scans.sh`, `reps_pfam.domtblout.gz`). A relaxed scan (no GA; E <= 1e-3) of NACHT, NB-ARC, AAA_16 and AAA_22 on all 6,462 HET/NLR proteins checked for a weak NTPase.
9. **Context** (`het_step3_context.py`, `het_step4_detail.py`): contig-end tests (Fisher, vs all other genes); copy number per strain (Mann-Whitney Ci vs Cp; Spearman vs `n_contigs` and assembly length within each species; all strains and dereplicated representatives only); distance to HrmA genes (`hrmA_gene_locations.tsv`); APH adjacency (APH hits from `hrmA_2026-09-24/aph_hac6_allprot.domtblout.gz`); islands (`report_tables/islands_with_domains.tsv`); pairs (`pair_classification.tsv.zst`); rescue overlap (`rescue_positions.tsv`, rescue start within +/-500 bp of an annotated gene).

## 1. Census

| Class | Proteins | Families | Ci proteins | Cp proteins | Main architectures (proteins) |
|---|---:|---:|---:|---:|---|
| NLR_NOD | 1,564 | 21 | 503 | 1,061 | NB-ARC (532); Patatin-NB-ARC-TPR (380); PNP_UDP-NPHP3_N-NACHT-WHD-Ank (309); NB-ARC-TPR (145); PNP_UDP-NPHP3_N-NACHT-WHD (144) |
| NLR_accessory_only | 2,639 | 41 | 839 | 1,800 | NPHP3_N-WHD (1,577); NPHP3_N (1,008) |
| HET_domain | 807 | 8 | 312 | 495 | HET-Het6_barrel (529); HET (278) |
| Goodbye (HeLo/SesA/Goodbye/HET-s class) | 277 | 9 | 95 | 182 | Goodbye (277) |
| Het-C | 1,175 | 7 | 341 | 834 | Het-C |
| NLR_like_AAA | 0 | 0 | | | |

- Models with no GA-passing domain in any protein: HeLo, HET-S, HET-s_218-289, SesA, NACHT_N, NACHT_sigma, NPHP3_hel, TPR_NPHP3, Beta-prop_NWD2_C, WHD_NACHT_C, WHD_NWD1, WHD_APAF1, WHD_CED4, TIR, TIR_2, WD40_like.
- No protein carries a WD40 repeat together with a NOD or NLR-accessory domain. The NLR repeat types found are TPR (L008) and Ank (L004).
- CHAT (1,032 proteins) and Peptidase_C14 (528, one per strain) never occur with a NOD or NLR-accessory domain.
- Relaxed NTPase check: 1,327 of 2,639 `NLR_accessory_only` proteins have a sub-GA NACHT or NB-ARC hit (E <= 1e-3). In the key families (H050, H052, H054, H047, H046) the weak NACHT hit overlaps the NPHP3_N envelope. This suggests that NPHP3_N covers part of a divergent NTPase region. That is an interpretation, not a test. No HET, Goodbye or Het-C protein has even a weak NTPase hit.
- Sequence-only hits: 1,177 proteins pass the NACHT sequence GA but no single domain passes the domain GA. 1,044 of them carry no other panel domain. They fall mostly in six core single-copy families (`seqlevel_only_nod_hits.tsv.gz`). They were excluded, following Pfam's two-threshold rule.
- Copy number: the largest HET/NLR families have at most one copy per strain.

## 2. Family-level polymorphism

Full table: `het_family_census.tsv` (86 families; labels H001-H086).

| Class | Families | Polymorphic within a species | Between species (BH q < 0.05) | Either |
|---|---:|---:|---:|---:|
| NLR_NOD | 21 | 9 | 11 | 11 |
| NLR_accessory_only | 41 | 12 | 12 | 16 |
| HET_domain | 8 | 4 | 4 | 4 |
| Goodbye | 9 | 4 | 7 | 7 |
| Het-C | 7 | 1 | 1 | 1 |
| Total | 86 | 30 | 35 | 39 |

(Annotated presence only. With genome_only cells counted, the totals change by at most 1 per class; see columns `poly_*_any`.)

Most of these family-level calls are not locus-level polymorphism (Section 3). Examples:
- H026 (Cp 96%, Ci 0.6%) and H027 (Ci 95%, Cp 1.9%) are the same gene at L004 (100% identity over 832 aa).
- H029-H035 are seven "species-specific" NB-ARC families at one site (L001a). Their representatives are 95-100% identical. They differ in length (150-398 aa).
- H001 (HET-Het6_barrel, 2,183 aa) is in 99% of Ci and 70% of Cp. The Cp "absences" are Cp strains with shorter models of the same gene (H003, H004, H005; 98.6-99.8% identity).

## 3. Loci: occupancy at the genomic site

Own-family scoring (step 5). `het_sites_ownfamily.tsv`; per-strain states: `het_sites_ownfamily_occupancy.tsv.gz`.

| Site | Families (main) | Architecture | Ci occupied / empty / unresolved | Cp occupied / empty / unresolved | Ci occ/resolved | Cp occ/resolved | BH q |
|---|---|---|---|---|---:|---:|---:|
| L001a | H029-H035, H039 | NB-ARC (150-398 aa) | 160 / 1 / 8 | 351 / 1 / 5 | 0.994 | 0.997 | 1.0 |
| L001b | H048, H053 | NPHP3_N + DUF7104 repeats | 155 / 1 / 13 | 321 / 10 / 27 | 0.994 | 0.970 | 1.0 |
| L002 | H009-H017 + H049, H055, H056, H058 | Goodbye gene + adjacent NPHP3_N gene | 167 / 0 / 2 | 356 / 2 / 2 | 1.000 | 0.994 | 1.0 |
| L003 | H018, H020 | Het-C | 167 / 0 / 1 | 356 / 0 / 1 | 1.000 | 1.000 | 1.0 |
| L004 | H026, H027 | PNP_UDP-NPHP3_N-NACHT-WHD(-Ank) | 168 / 0 / 1 | 360 / 0 / 0 | 1.000 | 1.000 | 1.0 |
| L005 | H050, H052, H054, H057, H059 | (PGAP1 / Abhydrolase_6+DUF676)-NPHP3_N-WHD_GPIID(-ANAPC4_WD40) | 169 / 0 / 0 | 358 / 0 / 1 | 1.000 | 1.000 | 1.0 |
| L006 | H047, H051 | (Abhydrolase_6+DUF676)-NPHP3_N-WHD_GPIID | 167 / 0 / 1 | 358 / 0 / 1 | 1.000 | 1.000 | 1.0 |
| L007 | H001, H003-H006 | DUF8348-HET-Het6_barrel | 166 / 0 / 3 | 358 / 0 / 1 | 1.000 | 1.000 | 1.0 |
| L008 | H025, H028 | Patatin-NB-ARC-TPR / NB-ARC-TPR | 158 / 0 / 7 | 256 / 1 / 103 | 1.000 | 0.996 | 1.0 |
| L009 | H019 | Het-C | 167 / 0 / 2 | 345 / 0 / 14 | 1.000 | 1.000 | 1.0 |
| L010 | H046 | PGAP1-NPHP3_N-WHD_GPIID | 167 / 0 / 1 | 355 / 0 / 2 | 1.000 | 1.000 | 1.0 |
| L012 | H036, H038 | Patatin-NB-ARC / NB-ARC | 5 / 105 / 59 | 17 / 169 / 81 | 0.045 | 0.091 | 1.0 |
| L013 | H040, H045 | NB-ARC | 1 / 96 / 71 | 4 / 169 / 187 | 0.010 | 0.023 | 1.0 |

- By the site definition, only L012 is polymorphic (Cp 9.1% of 186 resolved; within-species). L013 is below 5%. No site differs between species (all BH q = 1.0).
- L011 (H002, HET) and 71 small families have no core marker. Their site state cannot be scored.
- The step-2 scoring (any HET/NLR gene counts as occupant; `het_loci.tsv`) gives the same conclusion.

## 4. Allele classes at fixed sites

`het_allele_classes.tsv` (classes), `het_sites_ownfamily_alleles.tsv` (per family), `het_locus_allele_pairwise.tsv` (identities). Fractions are of resolved strains.

| Site | Allele class | Families | Ci | Cp | BH q | Within | Between | Both species >= 5% |
|---|---|---|---|---|---:|:---:|:---:|:---:|
| L005 | A: PGAP1-NPHP3_N-WHD | H050, H059, H069, H072, H077 | 52/169 (0.31) | 205/358 (0.57) | 2.9e-08 | Y | Y | Y |
| L005 | B: Abhydrolase_6/DUF676-NPHP3_N-WHD | H052, H057, H041 | 80/169 (0.47) | 82/358 (0.23) | 5.3e-08 | Y | Y | Y |
| L005 | C: NPHP3_N-WHD-ANAPC4_WD40 | H054 | 37/169 (0.22) | 71/358 (0.20) | 0.69 | Y | N | Y |
| L006 | A: NPHP3_N-WHD | H047, H065 | 100/167 (0.60) | 266/358 (0.74) | 1.3e-03 | Y | Y | Y |
| L006 | B: Abhydrolase_6/DUF676-NPHP3_N-WHD | H051, H075 | 67/167 (0.40) | 92/358 (0.26) | 1.3e-03 | Y | Y | Y |
| L002 | G1 Goodbye | H009, H015 | 1/167 | 136/358 (0.38) | 1.4e-25 | Y | Y | N |
| L002 | G2 Goodbye | H010 | 57/167 (0.34) | 0/358 | 6.3e-32 | Y | Y | N |
| L002 | G3 Goodbye | H011, H013 | 0/167 | 45/358 (0.13) | 2.9e-08 | Y | Y | N |
| L002 | G4 Goodbye | H012 | 20/167 (0.12) | 1/358 | 1.7e-09 | Y | Y | N |
| L008 | fused Patatin-NB-ARC-TPR | H025 | 18/158 (0.11) | 256/257 (1.00) | 1.7e-88 | Y | Y | Y |
| L008 | NB-ARC-TPR only | H028 | 140/158 (0.89) | 0/257 | 1.7e-90 | Y | Y | N |
| L004 | fused (+Ank) | H026, H062 | 1/168 | 353/360 (0.98) | 4.8e-127 | N | Y | N |
| L004 | no Ank | H027, H037 | 164/168 (0.98) | 14/360 (0.04) | 1.1e-112 | N | Y | N |
| L001a | one class (length variants) | H029-H044 | 160/161 | 351/352 | 0.61 | N | N | Y |
| L007 | one class (length variants) | H001, H003-H006 | 166/166 | 358/358 | 1.0 | N | N | Y |

Identity between allele classes (representatives, DIAMOND):
- L005: A-B 55.6% over 1,561 aa; A-C 53.7% over 1,576 aa; B-C 62.9% over 1,568 aa. Within a class: H050-H059 97.3%, H052-H057 98.8%, H052-H041 99.2%. For comparison, the L005 classes are 21-35% identical to the other NPHP3_N-WHD paralogs (H046, H047, H051, H048).
- L006: A-B 77.2% over 1,607 aa.
- L002 Goodbye: G1-G2 80.4%, G1-G3 86.5%, G1-G4 87.9%, G2-G3 79.7%, G2-G4 91.6%, G3-G4 74.2% (over 237-245 aa). The adjacent NPHP3_N genes (H049, H055, H056, H058) are 94.6-100% identical over 56-383 aa. They are short (105-690 aa) and have no NTPase hit, even relaxed. The Goodbye gene is present at L002 in 92/167 Ci (55%) and 182/358 Cp (51%) resolved strains (Fisher p = 0.40). So the Goodbye gene itself is polymorphic within both species at this site. The NPHP3_N gene is present in all.
- L008: H025-H028 90.1% over 624 aa.
- L004: H026-H027 100% over 832 aa.
- L001a: 94.9-100% among the NB-ARC families.

Same site, same flanks (exemplar neighbourhoods, `het_exemplar_loci.tsv`, `het_top_neighbours.tsv`):
- L005: every class sits between core genes 2OG-FeII_Oxy (`..._007783`) / Metallophos (`..._007501`) on one side and p450 (`..._008188`) / PAN (`..._006428`) on the other, in both species.
- L006: every class sits between core Aldo_ket_red (`..._004741`) / MFS (`..._002525`) and 1-cysPrx (`..._006974`).
- L002: between core LANC_like / CBF and Aa_trans / GCV_T genes.

Linked partner genes (`het_allele_partner_linkage.tsv`; presence from `presence_matrix.rescued.tsv`):

| Allele class | Partner family (adjacent) | Pfam of partner rep | Ci carriers with partner | Cp carriers with partner |
|---|---|---|---|---|
| L005-A (H050) | `B11518\|CF63038_008059-T1` | APH, Fructosamin_kin | 48/52 | 186/187 |
| L005-B (H052) | `730332_Guatemala\|C45F23D_005983-T1` | Fructosamin_kin | 45/45 | 82/82 |
| L005-B (H057) | same | | 33/33 | - |
| L005-C (H054) | `B11058\|C2FAE97_002404-T1` (+ `BFF340\|CF945FD_002404-T1`) | Fructosamin_kin | 37/37 | 71/71 |
| L006-A (H047) | `485B-1_L_OLD_CPA0023\|C4025D1_006853-T1` (+ `GT-134\|CE60774_006853-T1`) | Fructosamin_kin | 100/100 | 253/262 |
| L006-B (H051) | none | | 0/66 | 2/92 |

- The three L005 kinase partners are 60.9-64.6% identical to each other. They are not the partners of the other L005 classes. The pair table shows the same links: H047-006853 (linkage 1.0, Jaccard 1.0, q = 3.5e-30), H054-002404 (q = 3.6e-27), H050-008059 (q = 5.8e-25), H052-005983 (q = 9.3e-11).
- The L005/L006 kinase partners have no DIAMOND hit (e <= 1e-5) to the APH partners of HrmA (`kinase_partner_pairwise.tsv`).

Split versus fused forms (L004, L008):
- L004: in 158/161 Ci strains with H027, the adjacent Ank-repeat gene family `B0858-Guatemala|CBD0EE9_002373-T1` is present. It is present in 0/346 Cp strains with the fused H026 and in 7/7 Cp strains with H027. So the "no Ank" form always has the Ank part as a separate adjacent gene.
- L008: in 70/140 Ci strains with H028, a separate adjacent Patatin gene (`B11035|C5FC1D8_000795-T1`) is present. It is in 0/18 Ci strains with H025. 68 of 69 H025 genome_only (rescue) calls in Ci fall on this Patatin gene.
- I did not check the nucleotide sequence, so I cannot say whether the split is a stop codon/frameshift in Ci or an annotation difference. The pattern is consistent across strains within each species.

## 5. Copy number per strain

`het_genes_per_strain.tsv`, `het_copy_number_tests.tsv`. All 529 strains. Results with the 411 dereplicated representatives (133 Ci, 278 Cp) are in the same file and give the same direction.

| Count | Ci median (range) | Cp median (range) | Mann-Whitney p | Ci Spearman vs n_contigs (p) | Cp Spearman vs n_contigs (p) |
|---|---|---|---:|---|---|
| All HET/NLR | 12 (10-16) | 12 (10-15) | 3.9e-03 | 0.14 (0.071) | 0.07 (0.20) |
| NLR_NOD | 3 (2-5) | 3 (2-4) | 0.49 | 0.23 (0.0026) | 0.01 (0.83) |
| NLR_accessory_only | 5 (3-8) | 5 (4-7) | 0.32 | 0.05 (0.49) | -0.08 (0.12) |
| HET_domain | 2 (1-3) | 1 (0-2) | 2.3e-23 | 0.16 (0.038) | 0.32 (7.0e-10) |
| Goodbye | 1 (0-2) | 1 (0-1) | 0.25 | -0.07 (0.38) | -0.02 (0.73) |
| Het-C | 2 (1-3) | 2 (1-4) | 2.7e-13 | -0.09 (0.26) | -0.12 (0.025) |

- Means for all HET/NLR: Ci 12.37, Cp 12.14. The difference is small, though significant.
- The HET_domain difference comes mainly from H002 (Ci 84%, Cp 37%).
- NLR_NOD counts rise with contig number in Ci (rho 0.23). This fits gene models that split at fragmented sites (for example L008, where 26% of H025 and H028 genes are the last gene on their contig). HET_domain counts rise with contig number in Cp (rho 0.32). I cannot separate biology from assembly effects here.

## 6. Location, islands, HrmA/APH

Contig ends (`het_location_tests.tsv`; all other genes: 40.1% within 20 kb):

| Set | Genes | <= 20 kb of contig end | OR | p |
|---|---:|---:|---:|---:|
| All HET/NLR | 6,462 | 36.4% | 0.86 | 9.8e-10 |
| NLR_NOD | 1,564 | 28.4% | 0.59 | 6.5e-22 |
| NLR_accessory_only | 2,639 | 39.6% | 0.98 | 0.59 |
| HET_domain | 807 | 52.5% | 1.66 | 1.1e-12 |
| Goodbye | 277 | 11.6% | 0.20 | 1.0e-25 |
| Het-C | 1,175 | 34.6% | 0.79 | 1.1e-04 |

- H002 (HET) is the contig-end family: 68.7% within 20 kb, median contig 63 kb. H025, H028 (L008), H018 (Het-C) and H048 (L001b) are the last gene on their contig in 21-27% of copies. That is well above the other families (0-3%) and points to assembly breaks at these genes.
- No HET/NLR gene has a telomeric repeat at its nearest contig end (0 of 6,462). As in the HrmA report, only 647 of 737,258 contig ends carry a telomere array. Subtelomeric position cannot be tested.

Islands (`het_island_membership.tsv`): 37 HET/NLR families occur in 384 of 12,197 islands. All are `unexplained_physical` or `ambiguous_linkage`. H002 is in 200 islands (up to 14 strains per island; median 12 genes). H048 (45 islands), H050 (20), H053 (16), H036 (15) and H052 (10) follow. The islands with L005/L006 families have a median size of 3-5 genes.

Pipeline island Pfam enrichment (`island_pfam_enrichment.tsv`): 14 domains have FDR < 0.05. None is a NOD, HET, Ank, TPR or WD40 domain. NLR/HET-related rows: Patatin 12/15 families in islands, p = 0.013, q = 0.40; NPHP3_N 11/14, q = 0.54; CHAT q = 0.54; Goodbye q = 0.80; NACHT q = 0.80; NB-ARC q = 0.80; WHD_GPIID q = 0.80; HET q = 0.92; Ank q = 0.82. DUF676 is enriched (16/18, q = 0.029). DUF676 is on the L005-B and L006-B classes. APH (q = 7.3e-6) is also enriched. I did not test which island families drive these two domains.

Pairs (`het_pairs.tsv.gz`): 18 HET/NLR families occur in 133 significant pairs. The strongest physical pairs are the kinase partners in Section 4 and H002 with its APH neighbour (`578-1_L_NEW_CPA0049|C79DB71_008497-T1`, linkage 1.0, Jaccard 1.0, q = 9.0e-24). As the HrmA report noted, all pair rows list only Ci in `clade_composition`. Cp-only families cannot appear.

HrmA and APH (`het_context_summary.tsv`, `het_to_hrmA_distance.tsv.gz`):
- 21 of 6,462 HET/NLR genes (0.33%) lie within 10 genes of an HrmA gene. For other genes the rate is 0.66% (OR 0.49, p = 3.8e-4). HET/NLR genes are not in the HrmA loci.
- Exception: the L001a NB-ARC gene lies 10-11 genes (median about 27.5 kb) from HrmA locus L07 (F08/F18/F19) in 88 Ci strains, and in 2 Cp strains. L07 is a Ci-specific HrmA locus. L001a itself is fixed in both species.
- 438 of 6,462 HET/NLR genes (6.8%) have an APH-domain gene as immediate neighbour. Other genes: 1.06% (OR 6.8, p < 1e-300). This comes from four families: H002 (272/275), H028 (104/145), H050 (47/240; its partner family has both APH and Fructosamin_kin domains) and H025 (15/380). The H002 APH partner is only 24-29% identical to the HrmA APH partners.

## 7. Artefacts and caveats

1. **Ortholog splitting by length and gene model.** Many "species-specific" or "polymorphic" families are the same gene at one site with different lengths or split/fused models: L001a (H029-H044), L004 (H026/H027), L007 (H001/H003-H006), L003 (H020 is an 81-aa fragment next to H018), L002 NPHP3_N genes, L001b (H048/H053). Section 4 lists them. Use `het_allele_classes.tsv`, not family presence.
2. **Rescue on paralogs** (`het_rescue_overlap.tsv`): 361 of 368 rescue positions of HET/NLR families overlap an annotated gene; 258 overlap another HET/NLR gene. H029: 102/102 fall on H030 genes. H022: 121/121 fall on H018/H020 (Het-C). H025: 68 of 69 fall on the separate Ci Patatin gene. These genome_only cells are the same locus counted again.
3. **Contig ends.** L008 is unresolved in 103 Cp strains. H025/H028, H018 and H048 are often the last gene on a contig. H002 cannot be scored at the site level.
4. **Domain calls.** The census uses Pfam GA. It excludes sequence-only hits (1,177 NACHT, 3,042 NPHP3_N proteins) and weak NTPase domains. The relaxed check shows that many NPHP3_N-WHD proteins may be divergent NLRs. H046 (L010, core in both species, PGAP1-NPHP3_N-WHD_GPIID, 934 aa) may be a GPI inositol-deacylase (PGAP1/Bst1) rather than an NLR. WHD_GPIID is named after that enzyme. I did not test this. L005-A also carries a PGAP1 hit.
5. **Classification limits.** HeLo, HET-s, SesA and TIR were not found at GA. The screen may miss divergent N-terminal HET domains. No WD40-repeat NLR (HET-E/HET-D type) was found at GA. The C-terminal regions of the L005/L006 proteins (about 700-1,600 aa) have no GA repeat hit except one ANAPC4_WD40 hit in L005-C.
6. **Allele-class cut-off.** The 95% identity cut-off is my choice. The L002 fragment H014 is >= 98.9% identical to two classes that are 91.6% identical to each other. I left it unassigned.
7. **Split vs fused genes** (L004, L008) are consistent within species. I cannot say whether they are real gene structures or annotation differences without checking nucleotide sequence.
8. **Function is not shown.** Divergent, trans-species allele classes at a fixed site fit the known pattern of fungal het loci. This is an inference from the pattern. No incompatibility phenotype data exist here. I did not compare these loci with the A. fumigatus hetA-hetE loci by synteny.
9. **Pair and island tables are Ci-driven** (see the HrmA report). Cp-only families cannot show pairs.
10. **Strand.** `gene_positions.tsv.zst` has no strand column. I did not determine the orientation of the NLR-kinase pairs.

## 8. Comparison with the A. fumigatus classification

Using the categories of `Afumigatus_pangenome/results/GENE_CLUSTER.md`:

| A. fumigatus category | Found in Coccidioides (this run) | Where |
|---|---|---|
| NLR central domain (NACHT, NB-ARC) | Yes. 1,564 proteins, 21 families | L001a, L004, L008, L012, L013 |
| NLR C-terminal repeat, Ankyrin | Yes, on 311 NOD proteins and 43 NOD-less NPHP3_N/WHD proteins | L004 (Cp form) |
| NLR C-terminal repeat, WD40 | No NOD/NLR-accessory protein with WD40 at GA | L005-C has one ANAPC4_WD40 hit |
| TPR repeat (not seen in A. fumigatus) | Yes, 525 NB-ARC proteins | L008 |
| NPHP3_N, WHD_GPIID | Yes, the most common NLR signal (NPHP3_N on 2,625 NOD-less and 492 NOD proteins) | L002, L004, L005, L006, L010, L001b |
| Patatin effector (hetC-type) | Yes, on NB-ARC proteins | L008 (fused in Cp), L012 |
| PNP_UDP effector (hetA/hetD-type) | Yes, on NACHT proteins | L004 |
| CHAT (hetB-type) | CHAT exists (1,032 proteins) but never on an NLR | - |
| Island enrichment of these domains | Not significant in this run (all q >= 0.40) | Section 6 |

## 9. Ranked candidate loci for a locus view or manual inspection

Ranking criteria (my judgement): divergent allele classes at a fixed site > species-specific divergent classes > split/fused forms > presence/absence at low frequency. Numbers are occupied/resolved from step 5.

1. **L005: NPHP3_N-(weak NACHT)-WHD_GPIID NLR, three trans-species allele classes.** Families H050 (A), H052/H057 (B), H054 (C). Domains: PGAP1 (A) or Abhydrolase_6+DUF676 (B) N-terminal; NPHP3_N; WHD_GPIID; ANAPC4_WD40 (C). Site occupied Ci 169/169, Cp 358/358. Classes: A Ci 52 / Cp 205; B 80 / 82; C 37 / 71. Each class has its own linked fructosamine-kinase neighbour. Exemplars: A `UTAH_20380X16:scaffold_8:876856-881929` (Ci); B `B3411:scaffold_81:73238-78364` (Cp); C `UTAH_23742X123:scaffold_19:165851-170987` (Ci).
2. **L006: NPHP3_N-WHD_GPIID NLR, two trans-species allele classes (77% identity).** Families H047 (A), H051 (B). Domains: NPHP3_N, WHD_GPIID; B adds Abhydrolase_6+DUF676. Occupied Ci 167/167, Cp 358/358. A Ci 100 / Cp 266; B 67 / 92. A has a linked fructosamine-kinase neighbour; B has none. Exemplar: A `UTAH_20380X16:scaffold_37:228968-234107` (Ci); B `COCPO_Galgiani:scaffold_52:195565-200703` (Cp).
3. **L002: Goodbye-domain gene + adjacent NPHP3_N gene, four species-specific Goodbye classes (74-92% identity).** Families H009-H015 (Goodbye), H049/H055/H056/H058 (NPHP3_N). Site occupied Ci 167/167, Cp 356/358. Goodbye gene present in Ci 92/167, Cp 182/358. Exemplar: `M221:scaffold_52:47312-49811` (Cp; H009 + H049); Ci `UTAH_20380X10:scaffold_42:131502-132787` (H010 + H055).
4. **L008: Patatin-NB-ARC-TPR (hetC-type effector), fused in Cp, split in most Ci.** Families H025, H028 (90.1% identity). Occupied Ci 158/158, Cp 256/257 (103 Cp unresolved). Fused form Ci 18 / Cp 256. Exemplar: `COCPO_Galgiani:scaffold_4:87004-90882` (Cp, H025); Ci `UTAH_20380X10:scaffold_5:581687-584290` (H028).
5. **H002 (L011): HET-domain gene in an APH-containing accessory block.** Domain: HET (835 aa). Site not scorable (no core flanks). Family present Ci 142/169, Cp 133/360 (BH q = 2.2e-24). Exemplar: `UTAH_20380X16:scaffold_60:52132-54639` (Ci); `COCPO_Galgiani:scaffold_25:330503-333010` (Cp).
6. **L004: PNP_UDP-NPHP3_N-NACHT-WHD_GPIID(-Ank) NLR (hetA/hetD-type effector); Ank fused in Cp, separate gene in Ci.** Families H026, H027 (100% over 832 aa). Occupied Ci 168/168, Cp 360/360. Exemplar: `COCPO_Galgiani:scaffold_39:185645-189718` (Cp, H026); Ci `UTAH_20380X16:scaffold_8:396820-399336` (H027).
7. **L012: Patatin-NB-ARC, low-frequency presence/absence.** Families H036, H038. Occupied Ci 5/110, Cp 17/186. Exemplar: `UTAH_20380X1:scaffold_78:67309-68949` (Cp).
8. **L001a: short NB-ARC gene (150-398 aa) about 27 kb from HrmA locus L07 in Ci.** Families H029-H035 (95-100% identity). Occupied Ci 160/161, Cp 351/352. Exemplar: `UTAH_20380X16:scaffold_34:113059-114058` (Ci).

## Files (in `analysis/het_nlr_2026-09-25/`)

- Scripts, in run order: `run_allprot.sh`, `het_step1_census.py`, `het_step2_loci.py`, `het_step3_context.py`, `het_extract_seqs.py`, `run_small_scans.sh`, `het_step4_detail.py`, `het_step5_ownfamily_loci.py`, `het_step6_allele_classes.py`. Run the Python scripts with `pixi run --manifest-path <NII>/pixi.toml python <script>`. `het_extract_seqs.py` and `run_small_scans.sh` need `$SCRATCH`.
- Models and HMM outputs: `models.txt`, `het_nlr_models.hmm`, `allprot.domtblout.gz`, `allprot.tblout.gz`, `allprot.nseq`, `het_prots_nod_relaxed.domtblout.gz`, `reps_pfam.domtblout.gz`, `seqlevel_only_nod_hits.tsv.gz`.
- Census and families: `het_protein_domains.tsv.gz` (all 102,474 panel-hit proteins with architecture and class), `het_family_census.tsv`, `het_presence_matrix.tsv.gz`, `het_relaxed_nod_by_family.tsv`, `het_relaxed_nod_by_class.tsv`.
- Loci and alleles: `het_loci.tsv`, `het_locus_occupancy.tsv.gz`, `het_locus_alleles.tsv`, `het_locus_allele_diversity.tsv`, `het_locus_allele_spectrum_tests.tsv` (step 2); `het_sites_ownfamily.tsv`, `het_sites_ownfamily_alleles.tsv`, `het_sites_ownfamily_occupancy.tsv.gz` (step 5); `het_allele_classes.tsv` (step 6); `het_allele_partner_linkage.tsv`.
- Identity: `het_rep_pairwise.tsv`, `het_rep_pairwise_annotated.tsv`, `het_locus_allele_pairwise.tsv`, `kinase_partner_pairwise.tsv`.
- Context: `het_gene_locations.tsv.gz`, `het_family_locations.tsv`, `het_location_tests.tsv`, `het_neighbourhoods.tsv.gz`, `het_top_neighbours.tsv`, `het_exemplar_loci.tsv`, `het_exemplar_summary.tsv`, `het_genes_per_strain.tsv`, `het_copy_number_tests.tsv`, `het_to_hrmA_distance.tsv.gz`, `het_island_membership.tsv`, `het_pairs.tsv.gz`, `het_rescue_overlap.tsv`, `het_context_summary.tsv`, `neighbour_families.txt`.
- Reused from `analysis/hrmA_2026-09-24/`: `contig_lengths.tsv.gz`, `contig_telomere_ends.tsv.gz`, `hrmA_gene_locations.tsv`, `hrmA_family_census.tsv`, `aph_hac6_allprot.domtblout.gz`, `neighbour_reps_pfam.domtblout`.
