# yeast_convergence: analysis plan

Status (2026-10-22): designed, not run. No data pulled.
Comparisons: `COMPARISONS.md` (generated from `batches/*.yaml`).
How the files fit together and how to run: `README.md`.

## 1. Question and hypotheses

Nagy et al. 2014 (*Nat Commun* 5:4471) concluded that yeast growth evolved
repeatedly from filamentous ancestors by re-deploying an old cell-separation
toolkit through regulatory change, with parallel Zn-cluster TF expansion and
convergent gene loss. They tested this once, on 59 genomes, with each yeast clade
compared against the whole tree. This plan re-tests it with one independent
pairwise run per yeast origin, against that origin's own filamentous sister, plus
reversal and control runs.

| | Hypothesis | Prediction in this design | Where it is read |
|---|---|---|---|
| H1 | Convergent loss is part of becoming a yeast | Loss terms recur across origins more than chance **and** are not called by the filamentous controls, deep ones especially | `term_recurrence.tsv`, class `yeast_convergent`, `n_pseudo_deep_sig = 0` |
| H1-null | The loss signal is generic genome streamlining | The same terms are called by the controls | class `streamlining_generic` |
| H2 | Yeast form reuses a latent toolkit (cell-separation genes present in filamentous ancestors) | Toolkit families are present on both sides of every origin run: not origin novelties, not lost in the filamentous outgroups | `focus_terms.tsv` (`retained` rows) and step 5b |
| H3 | Regulation, not new effectors, drives each origin | Few novelty terms recur across origins; any that do are TF or regulatory families | `term_recurrence.tsv`, novelty direction |
| H4 | Reversals to filamentous growth do not regain lost hyphal genes | No overlap between an origin's losses and its reversal's novelties | `reversal_overlap.tsv` |
| H5 | Zn-cluster TF families expanded in parallel in yeast clades | Higher Zn2Cys6 copy number in yeast ingroups than in their sisters, origin by origin | Not testable with presence/absence; step 5a |

## 2. Design

16 runs in three sets, one batch YAML per set:

| Set | YAML | Runs | Role in the test |
|---|---|---|---|
| Origins | `batches/yeastconv_origins.yaml` | 9 | Each yeast or yeast-dominant lineage vs its closest filamentous or hypha-dominant sister. Three test dominance (dimorphic outgroup), one tests capacity (Mucor) |
| Reversals | `batches/yeastconv_reversals.yaml` | 2 | Filamentous lineages nested in yeast clades (Eremothecium, Neolecta), each paired with an origin run (`compare_with`) |
| Controls (pseudo-origins) | `batches/yeastconv_pseudo_origins.yaml` | 5 | Filamentous lineages with large ecological gene loss and no yeast phase. Two are deep (subphylum, class) to match the deepest origins |

Rules the design follows, all enforced by `bin/yc_batches.py check`:

- Every origin has a yeast-dominant ingroup (or, for the capacity run, a
  yeast-capable one) and a filamentous or hypha-dominant outgroup.
  Controls are filamentous on both sides. Reversal ingroups are filamentous.
- Growth form uses three states: filamentous, dimorphic with hyphal dominance,
  and yeast-dominant (`growth_form` + `dominant_phase` per species).
- Every run records the rank of its IN/OUT split (`split_rank`) so controls
  can be compared at matching depth.
- Every run has at least two ingroup species, except `yc_rev_neolecta` (one
  genome exists). Single-species novelty lists are not used.
- The two deep controls share no species with any other run.

## 3. Steps

Each step lists its command, its output, and the check that must pass before
the next step.

### Step 0. Settle the open decisions (section 6)

### Step 1. Resolve and fetch data

```bash
F=studies/fungi/yeast_convergence
python3 $F/bin/yc_batches.py check --master-pool <NovInvenio>/config_support/master_pool.csv
for d in $F/yc_*/; do pixi run python bin/ni resolve --study-dir "$d"; done
```

- 31 of 106 species are prefilled from proteomes other studies already use.
  The rest go through `ni resolve` (UniProt reference proteome first, NCBI
  fallback). Where it finds nothing, use the lab's local copy
  (`config_support/1KFG_samples.csv`, `Chaetothyriales_samples.csv`) as
  `Protein_Source=local_faa`. If no genome exists, replace the species in the
  YAML and re-render that run only.
- **Check:** BUSCO completeness from the resolve output. Drop or replace any
  proteome below about 85% complete (fungi_odb10), since missing genes look
  like losses. Record the threshold actually used in each run's notes.
- **Check:** for lichen runs, mycobiont-only proteomes; for AMF, BUSCO before
  reading losses.
- `bin/build_study_config.py` writes each run's `DATA_MANIFEST.yaml`
  (provenance) when it fetches.

### Step 2. Pipeline runs

```bash
DRY_RUN=1 $F/bin/run_all.sh && $F/bin/run_all.sh
```

- One SLURM job per run, NovInvenio `main-d03f5a2`, diamond, the same settings as
  `yeast_filamentous`. `run_all.sh` skips any run with unresolved rows.
- Pfam annotation jobs (`NII_ANNOTATE_PFAM=1`) start for non-UniProt
  proteomes and must finish before step 3.
- **Check:** each run has `results/<id>`, and every species has an
  `annotations/` file.

### Step 3. Enrichment per run

```bash
sbatch $F/bin/run_enrichment_all.sh
```

- Novelty and loss, proteome universe, clustered background. A novelty is a
  family in at least 75% of ingroup species and no outgroup species; a loss
  is a family in at least 75% of outgroup species and no ingroup species.
- **Check:** `tests.tsv` status `ok` for each run and direction; annotated
  fraction of candidates above about 0.8 for Pfam.

### Step 4. Cross-run comparison

```bash
python3 $F/bin/convergence_compare.py            # -> $F/analysis/convergence/
```

- Primary call on Pfam and Pfam clan. GO/InterPro are a cross-check over
  all-UniProt runs only.
- Recurrence: Poisson-binomial null per term, BH across terms, at least 3 of
  the 9 origins.
- Sensitivity runs to report alongside the primary one:
  `--min-origins 2` and `4`; `--level protein`; `--include-mixed-annotation`.

### Step 5. Follow-up analyses (not built yet)

- **5a. Zn-cluster expansion (H5).** Per-family copy numbers for PF00172 and
  PF04082, from the same proteomes, for every run's species. Test origin by
  origin (yeast ingroup vs sister) and with a birth-death model (CAFE5) on a
  species tree of all 106 species.
- **5b. Toolkit age (H2).** Presence of the cell-separation orthologs Nagy et
  al. used (Chs1/2/3, Eng1, Agn1, Cts1, Cps1/Bgs1, Mok1, Ace2) across all 106
  species, mapped on the species tree. The early-diverging species now in the
  set (Mucorales, Umbelopsis, Mortierellomycotina, Glomeromycotina) are what can
  move the toolkit's origin.
- **5c. Regulatory signal (H3).** Ace2/Swi5-type motif scan in promoters of
  the 5b genes, yeast vs sister in each origin run.

## 4. How results will be read

| Outcome | Reading |
|---|---|
| Recurrent loss terms, absent from deep controls | Supports H1: convergent loss specific to yeast form |
| Recurrent loss terms, also in deep controls | Streamlining; H1 not supported for those terms |
| Recurrent loss terms, only in shallow controls | Inconclusive; check depth before reading |
| Few or no recurrent novelty terms | Consistent with H3; lineage-specific gains |
| Recurrent novelty terms that are TF families | Supports H3 (regulatory) and points to 5a |
| Toolkit families present on both sides of every origin | Supports H2 |
| Origin losses re-appear as reversal novelties | Against H4, or the family was never lost at that origin's base |
| Signal depends on the dominance runs only | Dominance, not capacity, drives the signal |

## 5. Risks

| Risk | Mitigation |
|---|---|
| Missing species proteomes or poor assemblies read as losses | BUSCO filter (step 1); a loss needs absence from every ingroup species, so one weak proteome is not enough |
| Annotation source differs between runs | Pfam primary; GO/InterPro only on all-UniProt runs; `annotation_coverage.tsv` |
| Runs not independent (shared outgroups) | Deep controls share nothing; remaining overlap in section 6, item A |
| Control lists shorter for depth reasons | Deep controls; `split_rank`; deep counts per term |
| Growth-form codings wrong for a few species | Codings listed per species in `COMPARISONS.md`; confirm the flagged ones (section 6, item B) |
| NCBI/JGI proteomes lack GO/InterPro | Pfam via HMMER for every proteome |

## 6. Open decisions

A. **Shared outgroups between `yc_saccharomycotina` and `yc_taphrinomycotina`**
   (*Neurospora crassa*, *Arthrobotrys megalospora*). Decision (2026-10-24): keep
   them, disclose, and do not run a leave-shared-out rerun.
   Reason: in the first comparison (`analysis/convergence/term_recurrence.tsv`) neither
   run is among the origins that make a recurrent call. F-box (novelty) is called in
   `yc_dothideales`, `yc_pucciniomycotina` and `yc_mucoromycota`. The RNase_H clan
   (loss) is called in `yc_tremellomycetes`, `yc_ustilaginomycotina`,
   `yc_pucciniomycotina` and `yc_mucoromycota`. A rerun could add calls in the two
   shared-outgroup runs; it could not remove these four. It would shift the null
   probabilities of those two runs slightly, and that was not measured.
   If the check is wanted later, `enrich_novelties_v2.py` would need an option such as
   `--exclude-species Amega,Ncra`. What I read in the script: it builds the target
   and other species lists from `config.csv` GROUP (lines 369-370) and recomputes
   candidates from the presence matrix with `rule_candidates` (line 431), so dropping
   species from those two lists would change `n_other` and the counts. Not verified:
   the pipeline's own `candidates.txt`, tblastn validation and the loss tables were
   made with the full species set, so they would also need to be rebuilt or
   re-filtered. The script is shared with `zoosporic_dikarya` and `yeast_filamentous`,
   so a new option needs a test and a default that leaves their output unchanged.

B. **Species and codings to confirm** (default: keep, replace any without a
   genome after step 1):
   - Genome availability: *Mucor racemosus*, *Hortaea thailandica*,
     *Geosiphon pyriformis*, *Lobosporangium transversale*, *Umbilicaria pustulata*.
     Already replaced from BFD: *Golovinomyces orontii* -> *Podosphaera xanthii* YZU573,
     *Mucor indicus* -> *Thamnidium elegans* (predicted dimorphic, untested), *Evernia prunastri*
     -> *Letharia lupina*, *Cladonia grayi* -> *C. rubrotincta*.
   - Growth form: *Eremothecium cymbalariae*, *E. sinecaudum* (traits.csv:
     yeast, from a genus-level call); *Exophiala oligosperma*, *E. mesophila*
     (yeast-dominant?); *Protomyces inouyei*, *Taphrina deformans* (coded
     hypha-dominant); *Yarrowia lipolytica* (traits.csv: yeast; here dimorphic).

C. **BUSCO threshold** for step 1. Default about 85% complete.

## 7. Status

- [x] Design, batch YAMLs, per-run `species.csv`, `run_params.txt`, `publish.yaml`
- [x] Checker, renderer, run and enrichment drivers, comparison script, tests
- [x] Nested-folder publishing in the shared scripts
- [ ] Decisions A–C
- [ ] Step 1: resolve, fetch, BUSCO
- [ ] Step 2: pipeline runs
- [ ] Step 3: enrichment
- [ ] Step 4: comparison
- [ ] Step 5: follow-ups (scripts to write)
