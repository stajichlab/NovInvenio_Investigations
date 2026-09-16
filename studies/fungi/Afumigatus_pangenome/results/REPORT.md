# Afumigatus_pangenome — full corrected results report

**295 genomes (293 *Aspergillus fumigatus* ingroup + 2 outgroup: *A. lentulus*, *A. fischeri*)**
**Status: corrected, 2026-09-15.** Two real bugs found and fixed the same day (tblastn
truncation, missing rescue-pass positions) — see "Correction history" below before citing
any number in this report against an earlier version of the pipeline's output.

Figures referenced below are in `figures/*.png` (raster, for viewing) and
`figures_pdf/*.pdf` (vector, for print/manuscript use — identical content). Full
per-metric tables in `SUMMARY.md`. Raw data in `full_293run/*.rescued.tsv`,
`id_crosswalk/*.tsv`, `full_293run/modules_preview/*.tsv`.

---

## 1. Correction history (read this first)

Two independent, verified bugs were found and fixed on 2026-09-15, both changing the
headline numbers substantially:

1. **tblastn rescue-pass truncation.** The genome-level rescue pass (recovering gene
   families missed by protein annotation) used `-max_target_seqs 5` against one combined
   295-strain BLAST database — a cap on hits *per query across the whole database*, not
   per strain. Verified: 89.4% of hit queries hit exactly the 5-strain cap. Fixed by
   redesigning to one BLAST database per strain (`run_rescue_pass_per_strain.sh`).
2. **Missing rescue-pass positions.** After fix #1, `pair_classification.py` couldn't
   resolve ~96% of pairs, because a rescue-pass hit has no annotated gene model, so it
   had no position in `family_positions.tsv`. Fixed by deriving positions directly from
   each hit's own genomic coordinates (`bin/extract_rescue_positions.py`).

| Metric | Pre-fix | Post-fix (this report) |
|---|---:|---:|
| Core families | 13.8% | **53.0%** |
| Singleton families | 57.6% | **22.6%** |
| Physically-linked co-occurring pairs | 10,687 | **83,415** |
| `insufficient_data` (pair classification) | 6.0% | **0.3%** |

Full technical detail: `PANGENOME_CLUSTER_PROFILE_NOTES.md`.

---

## 2. Pangenome composition

| Bin | Families | % of total |
|---|---:|---:|
| Core | 25,436 | 53.0% |
| Soft-core | 643 | 1.3% |
| Shell | 6,998 | 14.6% |
| Cloud | 4,054 | 8.4% |
| Singleton | 10,852 | 22.6% |
| **Total** | **47,983** | **100.0%** |

![Pangenome composition](figures/core_shell_cloud_pie.png)
![Frequency distribution](figures/frequency_distribution.png)

**The pangenome is open, not saturated.** The accumulation curve (20 random-order
permutations) shows neither the pangenome-size curve nor the core-size curve has
plateaued at n=295 — both are still moving in the same direction they started in.

![Accumulation curve](figures/accumulation_curve.png)

---

## 3. Assembly/annotation quality control

| Metric | Value |
|---|---:|
| Strains assessed | 295 |
| Complete BUSCO % (min / median / max) | 96.7 / 98.8 / 99.4 |
| Strains with >100 contigs (draft-consistent) | 279 |
| Strains with ≤20 contigs (near-chromosome) | 13 |

Completeness is uniformly high and narrow — not a confound for any "missing family"
call in this panel.

**One unresolved outlier**: `Asfu_H1106` (strain `H-1-10-6`, `GCA_020501995.1`) shows an
elevated shell/cloud/singleton family count (6,197 vs. panel median 3,591). Traced to
protein-level annotation (not the rescue pass). BUSCO (98.7%), contig count (678,
non-extreme), total protein count (10,961, normal range), and annotation source (NCBI,
same as 292/295 strains) were all checked and ruled out as explanations. **Kept in the
panel** — no positive evidence of a technical artifact. Flagged as an open QC item.

---

## 4. HAC / hacA targeted screen

| Locus | Present in | % of strains |
|---|---:|---:|
| hacA (Afu3g04070) — core UPR transcription factor | 293/293 | 100.0% |
| hrmA ortholog — subtelomeric, Starship-mobile HAC element | 8/293 | 2.7% |

A broad Pfam-family screen (PF11001) is **not usable** as a presence/absence marker —
it hits every strain (3-9 hits each), because it's an ancient, non-mobile paralogous
domain family, not specific to the Starship-mobile copy. The direct sequence-identity
screen resolves this: 8 strains carry the true hrmA ortholog, one (`Asfu_Af293`)
confirmed exactly against the published Starship `Nebuchadnezzar-h1` at the same base
coordinate.

---

## 5. Co-occurrence and pair classification — is loss/gain limited to Starships?

**Short answer: no.** There are 4,243,692 statistically robust (FDR<0.05, exact
clade-stratified test) co-occurring family pairs, and Starship-explained pairs are a
small minority even among the physically-linked ones.

| Classification | Pairs | % of total | What it means |
|---|---:|---:|---|
| trans | 2,613,303 | 61.6% | Not physically clustered; co-occurrence survives the clade-permutation null (not just population structure) |
| trans_unconfirmed | 1,533,331 | 36.1% | Not physically clustered; fails the permutation null (likely explained by shared clade membership) |
| **unexplained_physical** | **44,576** | **1.1%** | **Physically clustered — a real co-located, coordinately-varying gene block — with NO Starship captain gene nearby** |
| ambiguous_linkage | 29,949 | 0.7% | Partial physical clustering |
| insufficient_data | 13,643 | 0.3% | Too few resolvable co-carrying strains |
| **starship_explained** | **8,890** | **0.2%** | **Physically clustered AND a Starship captain gene sits nearby — mechanism identified** |

![Pair classification](figures/pair_classification_summary.png)

**There are 5x more `unexplained_physical` pairs than `starship_explained` ones.**
Among the 83,415 physically-linked pairs (starship_explained + unexplained_physical +
ambiguous_linkage), the *majority* have no identified mechanism. A concrete example —
one of the strongest `unexplained_physical` pairs (jaccard = 1.0, adjacent genes,
identical strain distribution, consistent "loss" direction across clades):

```
Neofi_ref|tr|A1DK62|A1DK62_ASPF1  ×  Neofi_ref|tr|A1DK67|A1DK67_ASPF1
  classification: unexplained_physical | linkage_fraction: 1.00 | jaccard: 1.00
  fisher_p: 1.58e-11 | fdr_q: 1.31e-09 | permutation_p: 0.0000 | direction: loss
```

This is real, tight, coordinated gene loss — not a statistical artifact (fdr_q ~1e-9,
permutation_p essentially 0) — with no Starship captain gene detected nearby. Whether
this reflects a different mobile-element family, simple deletion of a co-regulated
operon-like block, or some other mechanism is an **open question this pipeline
surfaces but does not resolve**. `results/full_293run/pair_classification.rescued.tsv`
has the full list, filterable on `classification == "unexplained_physical"`, sorted by
`jaccard`/`fdr_q` for the strongest, most confident candidates to follow up on.

**The `trans` category (61.6% of all significant pairs) is a separate phenomenon.**
These are NOT physically linked at all — two genes at different genomic locations
whose presence/absence correlates across strains. This is not "pattern of physical
loss/gain," but could reflect shared selective pressure, epistasis, or an unmodeled
confound. It is the largest category by far and not explored further by this pipeline
version.

### Are the physically-linked pairs part of real, larger multi-gene islands?

Yes — and most of them are unexplained by either mechanism this study checked. Ran
`bin/synteny_windows.py`'s `accessory_islands()` (tested during design, never
previously run on real data) against all 295 strains via a new script,
`bin/find_accessory_islands.py`, merging maximal runs of consecutive non-core genes
into islands, then keeping only islands that contain a statistically significant
physically-linked pair (not just any non-core run). Cross-referenced each island
against the existing Starship captain-gene (DUF3435) screen and a new secondary-
metabolite backbone screen (PKS ketosynthase `PF00109`, NRPS condensation domain
`PF00668`; `results/sm_backbone/`).

**Result: 12,861 distinct significant islands** (2–837 genes, median 7 — the largest
are very likely extended subtelomeric/repeat-rich regions, not compact gene clusters;
realistic biosynthetic-gene-cluster candidates sit in the 3–30 gene range).

| Captain gene (Starship) | SM backbone gene (PKS/NRPS) | Distinct islands |
|---|---|---:|
| No | No | 11,741 (91.3%) |
| Yes | No | 904 (7.0%) |
| No | Yes | 151 (1.2%) |
| Yes | Yes | 65 (0.5%) |

**91.3% of significant physical-linkage islands have neither mechanism identified.**
One concrete, checkable example in the realistic size range: a 30-family island in
strain `Asfu_G2141` (246 supporting pairs) contains the reviewed Swiss-Prot entry
`Asfu_Af293|sp|Q4WKX2|FGND_ASPFU`, carries a PKS/NRPS backbone hit, and has **no**
Starship captain gene — a real candidate secondary-metabolite-associated accessory
island independent of Starship mobilization. Full list:
`results/accessory_islands/significant_islands.tsv`.

*Caveats*: the captain-gene/SM-backbone flags are cohort-level (any strain, any copy
of that family) — a "Yes" doesn't confirm the specific copy in that specific island
carries the domain. The 12,861 count is "distinct exact member sets," not
non-overlapping loci — a larger island in one strain and a smaller nested subset in
another both count separately.

*Performance note*: the naive version of this cross-reference (checking all 83,415
significant pairs against every island) was killed after 12+ minutes with zero
output — an O(islands × pairs) blowup. Fixed with a family-indexed lookup
(Fable-reviewed, confirmed correct and empirically safe at this dataset's scale);
the real run then took 2 minutes.

### What functions are these islands enriched for?

Ran a real Pfam-A hmmscan (30,134 profiles) against all 11,052 families eligible for
co-occurrence testing (the 9,224 in significant islands plus the 2,620 additional
eligible families needed for a correctly-scoped background), then tested each Pfam
domain for enrichment among island-member families vs. that background (one-sided
Fisher's exact, BH-FDR corrected) — **901 domains tested, 71 significant at FDR<0.05**.

| Domain | Islands / background | Interpretation |
|---|---:|---|
| DUF3435 | 199/205 | The Starship captain domain itself — internal-consistency check |
| DUF3723 | 109/109 | Documented as part of the broader Starship "backbone" beyond the captain gene — some of the "unexplained" islands may still be Starship-associated via a signal outside this study's narrower DUF3435-only screen |
| Ank / Ank_2–5 (ankyrin repeat) | up to 215/231 | Documented Starship-cargo-associated domain family in other fungi |
| **DDE_1 (DDE transposase)** | 68/68 | **Independent evidence of a different mobile-element family** — direct support for the hypothesis that some unexplained_physical islands reflect a non-Starship transposon |
| NACHT | 70/73 | Fungal heterokaryon-incompatibility/innate-immune-like genes — classic accessory-genome category |
| Glyco_hydro_71, Patatin | ~55 each | Secreted cell-wall-modifying/lipase enzymes — another classic horizontally-variable category |

![Top enriched domains](figures/island_domain_enrichment.png)

**Full data file, with FDR included, sorted by island size**:
`results/accessory_islands/significant_islands.with_enrichment.tsv` — every distinct
significant island, sorted largest-first, with its full Pfam domain list, which of
those domains clear FDR<0.05, and each one's q-value
(`bin/annotate_islands_with_enrichment.py`, joining `significant_islands.with_domains
.tsv`'s per-island domain lists against `domain_enrichment.tsv`'s per-domain q-values
— neither source file alone has both). **5,712 of 12,861 islands (44.4%) carry at
least one FDR-significant domain.**

**12 concrete examples**, spanning the full range of strain support (1–187 strains)
and island size (3–30 genes):

| Island (strain, size) | Strains supporting | Classification(s) | Captain gene | Top domains (FDR q) |
|---|---:|---|:---:|---|
| `Asfu_HMR_AF_270`, 8 genes | 187 | unexplained_physical | N | NPHP3_N (2.2e-06), WHD_GPIID, NACHT, NB-ARC |
| `Asfu_CNMCM8686`, 6 genes | 170 | unexplained_physical | N | Patatin (1.9e-04) |
| `Asfu_niveus`, 3 genes | 158 | unexplained_physical | N | NPHP3_N (2.2e-06) |
| `Asfu_HMR_AF_270`, 4 genes | 132 | unexplained_physical | N | Ank_2, Ank_4, NPHP3_N, Ank_5, WHD_GPIID, PNP_UDP_1 (2.7e-09) |
| `Asfu_HMR_AF_270`, 5 genes | 116 | unexplained_physical | N | WD40_CDC20-Fz, WD40_Prp19 (2.7e-04) |
| `Asfu_NRZ2018236`, 30 genes | 1 | ambiguous_linkage, starship_explained, unexplained_physical | Y | DUF3435, DUF3723, Ank×3, NPHP3_N, NACHT +5 more (4.4e-14) |
| `Asfu_B170s1`, 30 genes | 2 | ambiguous_linkage, unexplained_physical | N | Ank_2, Ank, Ank_4, Ank_5, Ank_3, Myb_DNA-bind_6 (2.7e-09) |
| `Asfu_NRZ2017339`, 30 genes | 2 | ambiguous_linkage, unexplained_physical | N | Helicase_C, DEAD, GloB_C, Pkinase (1.8e-08) |
| `Asfu_E12510`, 30 genes | 3 | ambiguous_linkage, starship_explained, unexplained_physical | N | NPHP3_N, HSF_DNA-bind, WHD_GPIID, NACHT, APH (2.2e-06) |
| `Asfu_E166s1`, 30 genes | 5 | ambiguous_linkage, starship_explained, unexplained_physical | Y | DUF3435, GloB_C, APH (4.4e-14) |
| `Asfu_NRZ2018649`, 30 genes | 5 | ambiguous_linkage, unexplained_physical | N | Pkinase, Myb_DNA-bind_6 (1.1e-03) |
| `Asfu_E174s1`, 30 genes | 7 | ambiguous_linkage, unexplained_physical | N | HSF_DNA-bind (1.6e-05) |

The first five rows are the most striking: **highly recurring** (116–187 of 295
strains carry the exact same island), small-to-moderate (3–8 genes),
`unexplained_physical` (no Starship captain gene at all), yet strongly enriched for
specific, interpretable domains — NPHP3_N and NACHT (both linked to fungal
non-self-recognition/heterokaryon-incompatibility systems), ankyrin repeats, WD40
repeats, and Patatin (secreted lipase). These read as real, conserved,
non-Starship-driven accessory gene modules, not one-off noise. `Asfu_NRZ2018236`
(row 6) is the richest single example — captain gene, DUF3723, ankyrin repeats, and
NACHT all in the same 30-gene island in one strain, a genuine candidate for a hybrid
or compound accessory locus worth manual inspection.

*Caveat found while spot-checking real output, not a bug*: an island's own Pfam
domains and why one of its pairs was called `starship_explained` don't always
overlap — `pair_classification.py`'s captain-gene check uses a fixed window around
the *pair*, not the island's own boundary, so a captain gene can be "nearby" by the
pair-level definition while sitting just outside the island itself.

### Trans-network module structure (Leiden community detection)

Collapsing all 2,613,303 `trans`-classified pairs into communities (Leiden, resolution
sweep 0.5–10.0, `results/full_293run/modules_preview/`):

| Resolution | Modules | Largest module | Singleton modules |
|---|---:|---:|---:|
| 0.5 | 6 | 3,288 | 0 |
| 1.0 | 10 | 2,105 | 0 |
| 2.0 | 27 | 1,595 | 5 |
| 5.0 | 724 | 1,405 | 569 |
| 10.0 | 1,537 | 1,140 | 1,359 |

A persistent, densely-interconnected core structure across resolutions (no obviously
"right" resolution — module count and singleton fraction both grow steeply and
continuously). Whether this reflects real shared biology or a residual
population-structure confound not fully caught by the clade-stratified null is **not
resolved**.

---

## 6. Benchmark scorecard (validation against the reference paper)

Built a sequence-based ID crosswalk (614 NCBI-fetched reference sequences from
mbio.01092-25-s0002.xlsx Tables S6/S7/S12/S13/S19/S21) and scored both this study's
mmseqs clustering (rescued) and a newly-run diamond clustering (unrescued — no rescue
pass has been run for the diamond backend, so any mmseqs-vs-diamond gap below is
confounded by that, not a clean backend comparison) against it.

**Real coverage limit found**: 12 of 13 "reference-quality" strains in the paper's own
Tables S12/S13 use the paper's internal locus numbering, not public accessions — not
resolvable to a sequence from the supplement alone. The crosswalk is effectively
**AF293-anchored**: of 20 high-confidence Starships, only `Nebuchadnezzar-h1`/`-h2`
(the hrmA/HAC locus already covered in section 4) have a resolvable cargo gene.

| Control | Backend | Result |
|---|---|---|
| Cargo grouping (Nebuchadnezzar-h1/h2) | both | purity 1.0, completeness 0.05 — correctly does NOT over-merge the ~20 PF11001 paralogs into one family |
| Presence recovery (Nebuchadnezzar-h1) | both | accuracy 0.992, Jaccard 0.600 |
| Presence recovery (Nebuchadnezzar-h2) | both | accuracy 1.0, Jaccard 1.0 |
| Negative control (645 conserved genes, Table S19) | mmseqs (rescued) | 85.9% core (vs. 53.0% genome-wide) |
| Negative control | diamond (unrescued) | 79.5% core (vs. diamond's own baseline) |

Both backends correctly recover the one resolvable positive control and correctly
enrich the negative control toward core, relative to their own genome-wide baseline.
**Not yet extended** to the other 18 Starships — would need either external sequence
fetching beyond the paper's own supplement, or running `starfish`
(github.com/egluckthaler/starfish, the reference paper's own tool) directly on this
study's genomes. Considered, not done this session (see `PANGENOME_CLUSTER_PROFILE_NOTES.md`).

---

## 7. Dereplication threshold sensitivity

| Mash threshold | Representative strains (of 295) |
|---|---:|
| 0.0005 | 207 |
| 0.001 (default, used throughout this report) | 123 |
| 0.002 | 14 |
| 0.005 | 3 |

**Extreme sensitivity** — the default sits in the steep part of the curve, not a stable
plateau. Any claim that specifically leans on "123 representative strains" should carry
this caveat. Not resolved further (would need an independent criterion, e.g. a
published clonal/outbreak-cluster ground truth, to pick a threshold on grounds other
than "the value used so far").

---

## 8. Synthesis

Three findings reinforce each other. First, the pangenome is genuinely open despite
being core-dominated (53%) — a real, actively-generated accessory tail that keeps
expanding as more strains are sampled, not an artifact of incomplete sampling
saturating soon. Second, that accessory tail is *mostly not* explained by physical
clustering at all (98% of significant co-occurrence is `trans`, not physically
linked), and even within the physically-linked minority, **Starships explain a
minority of a minority** — 8,890 of 83,415 physically-linked pairs (10.7%). Third, the
HAC/hrmA screen shows why "mobile" and "common" aren't the same thing: the element
travels via Starships, but a specific Starship-riding gene is a minority event (2.7%
of strains).

Together this points to a species whose accessory genome is real, large, and driven by
at least two distinct, only-partly-understood mechanisms: a small, well-characterized
Starship-driven physical-linkage signal, a **larger unexplained physical-linkage
signal** worth its own targeted follow-up, and a much larger *trans* (non-physical)
co-occurrence signal whose biological meaning this pipeline surfaces but does not
explain.

---

## 9. Open items

**To finish the *A. fumigatus* profile:**
- Investigate the 44,576 `unexplained_physical` pairs directly — cluster them by
  genomic region, check for a non-Starship mobile-element signature or a shared
  regulatory context.
- Extend the benchmark scorecard past the single AF293-anchored control (starfish run,
  or external sequence fetching for the other 18 Starships).
- Resolve or bound the dereplication-threshold and Leiden-resolution open questions
  with an independent criterion.
- Follow up the `Asfu_H1106` outlier (allele-level divergence / clade placement check).

**To generalize into an NI Nextflow module** (parallel track, `NovInvenio` repo,
branch `pangenome-profiling-module`, not yet merged):
- Real DSL2 subworkflow built, reviewed (Opus), all findings fixed, **and now
  confirmed end-to-end on real data** (SLURM job 28428780: 57/57 processes, 0
  failures, real non-trivial output).
- Being tested live against a second species (Coccidioides) in parallel — found and is
  fixing a real GFF3-parsing gap (funannotate-style `Parent=` vs. NCBI-style
  `protein_id=` attributes).
- Diamond-backend clustering intentionally hard-disabled (unvalidated ID-matching
  path) until someone builds a safety net equivalent to the mmseqs path's
  `restore_mmseqs_cluster_ids.py`.

---

## 10. File manifest

| File | Contents |
|---|---|
| `figures/*.png`, `figures_pdf/*.pdf` | The 5 summary figures (raster + vector) |
| `SUMMARY.md` | Machine-generated aggregate tables (subset of this report) |
| `full_293run/presence_matrix.rescued.tsv` | Corrected presence/absence/genome-only calls, 47,983 × 295 |
| `full_293run/frequency_table.rescued.tsv` | Per-family bin + frequency |
| `full_293run/cooccurring_pairs.rescued.tsv` | 4,243,692 FDR-significant pairs (exact test) |
| `full_293run/pair_classification.rescued.tsv` | Same pairs + physical-linkage classification |
| `full_293run/modules_preview/` | Leiden trans-network modules, resolution sweep |
| `id_crosswalk/` | Paper ID crosswalk, ground-truth tables, benchmark scorecard |
| `busco_genome/busco_completeness_summary.tsv` | Per-strain BUSCO stats |
| `hac_crosswalk/hac_reference_screen.tsv` | Per-strain hacA/hrmA screen |
