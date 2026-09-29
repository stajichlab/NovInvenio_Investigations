# Session handoff — 2026-09-28

Covers one long session (2026-09-26 to 2026-09-28). Earlier handoff:
`2026-09-23-session-handoff.md`. A parallel session worked on #212 / #214 and the
Afumigatus test sets at the same time; its state is in the memory note
`afumigatus_test_sets` and in NII commits b3444c7, 83c164a.

## 1. Code landed

NovInvenio (all merged to main, now at f87fd1e):

| PR | What |
|---|---|
| #201, #204 | Island locus view: loci, DNA presence check (blastn, flank-inclusive targets, rescue TBLASTN spans), five rankings (C default; A, B, D, E), 100 drawn loci; clinker synteny panel (clinker 0.0.32 from PyPI in the container). Released as v0.7.0 (a394ace); image built by CI. |
| #209 | Spec only: `docs/superpowers/specs/2026-09-28-outgroup-evidence-coverage-design.md` for #208. |
| #213 | Island steps (ISLAND_DNA_TARGETS, ISLAND_LOCI, ISLAND_SYNTENY) read the merged Pfam table (`pfam_domtblout_ch`) instead of the per-chunk FAMILY_PFAM_SCAN channel. Cause of both full-run failures (section 3). |

NII (all pushed to main):

| commit | What |
|---|---|
| NII PR #8 (merged earlier) | `sync_pangenome_report.py` stages `clinker/` for the top 25 loci (`--clinker_publish_top`). |
| 85c5911 | Preview run `cocci_locusview_preview` on the site (built by hand from `rescue_freqpol_immitis_in_posadasii_out`). |
| fbcfe63 | `pangenome_runs.yaml` for Afumigatus and Cocci at v0.7.0 (a394ace). |
| 3e5e64c, d5a66db | `studies/fungi/coccidioides_pangenome/analysis/sans/` (SANS test, tree comparison, prototype page) and `analysis/species_tree/` (BUSCO CDS tree + 51 mash grafts, first committed here). |
| 3219de7, ee563a0 | #132 priorities 2-4 notes in this folder. |

## 2. Results

**SANS and phylogeny (Cocci, `analysis/sans/README.md`).** SANS 2.8_2 (built at
`/bigdata/stajichlab/jstajich/software/sans`, maxN=600, static libstdc++) runs fast on
529 genomes (DNA k=31: 4 min, 15.6 GB; protein: 38 s). SANS has no gene-content mode;
`gene_content_splits.py` feeds family carrier sets to `SANS -s`. Gene presence/absence has
little tree structure inside a species: mean changes per family on the SNP tree vs a
label-shuffled tree are 10.88 vs 14.34 (Ci) and 17.65 vs 20.74 (Cp); pooled 19.06 vs 31.77.
Only 716 of 28,830 families (2.5%) fit one clade of the SNP IQ-TREE. The SNP tree explains
gene content better than the BUSCO tree and the SANS trees.

**Prototype page** (private): https://claude.ai/artifact/Wz53NNHTojwx4gCbfdAy8y — tree-ordered
heatmap (SNP or BUSCO tree, clade-grouped columns, homoplasy strip, clade filter) and split
networks (SANS splits on a NeighborNet ordering; DNA 663/1,028 splits drawn, 97% of weight;
gene content 104/566, 33%). Build: `analysis/sans/prototype/`. Not integrated into the
pipeline report; the user has not decided.

**Novelty case A0A090CJI0_PODAN (sordariales_shallow).** Called novel although TBLASTN hits
all 8 outgroups. Cause: filter 2 (paralog competition) removes all 33 outgroup protein hits
because the paralog B2B454_PODAN beats the query on the same targets; all TBLASTN hits cover
only ~20% of the query (Patatin PF01734). Same class as HEX-1. Design for recording this
evidence: #208 (not implemented).

**#132 priorities 2-4** (this folder, `2026-09-28-p2-*`, `-p3-*`, `-p4-*`; summary comment on
#132): core 0.95 agrees with a binomial-mixture core boundary (0.94-0.98); 0.90 and 0.15 are
conventions; cloud cannot exist below 14 strains. Pair-class trans call stable over k 5-20,
physical count doubles. Leiden resolution 1.0 not validated; stability-only selection picks
trivial partitions for pooled Cocci.

## 3. Runs

- **Failed 2026-09-28 at ISLAND_LOCI** (fixed by #213): Afumigatus `full_v070` (job 29171115,
  18 h) and Cocci `immitis_in_posadasii_out_v2` (job 29171116, 21 h). Both still pin a394ace.
  Plan (user, via the parallel session): iterate on the Afumigatus test sets first, then
  relaunch the full runs **under new run names** (`ni` refuses a commit change for an existing
  run name). Do not resume or resubmit without the user's OK.
- **Running at handoff:** `Afumigatus_test10` run `in10_v1` (job 29196843, f87fd1e), from the
  parallel session.
- **Published Cocci runs:** `cocci_locusview_preview`, `rescue_freqpol_immitis_in_posadasii_out`,
  `rescue_freqpol_posadasii_in_immitis_out`. Afumigatus: nothing published.

## 4. Open work

| item | state |
|---|---|
| #208 other-group evidence + coverage class (novelty tables/pages) | spec merged, ready to implement |
| #193 + #203 island steps assume a Pfam HMM is given | fix not started; user never answered whether to start |
| After relaunched full runs | #187 outgroup gene-positions check (Afu); publish Afu; replace the Cocci preview (needs OK) |
| #132 | P1 cross-taxon re-test (the Afumigatus test sets fit); P3/P4 on more graphs; file issues for small-n report warnings, Leiden stability sweep, small-panel island limit (0 FDR pairs at 13-24 strains) |
| Fusarium studies | FOL (partly), FSSC_solani, FSSC_ambrosia, fusarium_FOXY_vs_FSSC untracked in NII; need provenance records before commit; none published |
| Bd_pangenome off HPCC | blocked on user: target system (Docker/Apptainer?) and genome source. No Bd dataset exists yet. Code gap: `lib/ni_pangenome.py` is SLURM-only (sbatch, squeue, `module load java`); famsa AVX2 pinning (#148, #181) |
| Visualization | decide whether the heatmap / network prototype goes into the pangenome report |

## 5. Cleanup done 2026-09-28

Deleted: NII root `S1.gbk`, `S2.gbk`, `groups.csv`; 12 `docs/fungi/*_dmnd_*` test folders
(sources kept in `results/*_dmnd_*`); NovInvenio worktrees `fix-island-pfam-merged`,
`outgroup-evidence-spec`, `island-locus-view`, `locus-view-b`, `locus-view-plan`; 15 merged
local branches. Kept: 5 detached run worktrees, `NovInvenio-worktrees/docs/test`,
`NovInvenio-worktrees/mmseqs-simd-test`, the P3 sweep outputs (271 MB,
`coccidioides_pangenome/results/p3_pair_class_sweep/`), Nextflow clone 5b83064 (v0.6.0).

## 6. Environment notes

- `/usr/bin/python3.12` user site now has scikit-learn 1.9.1, scipy 1.18.1, splitspy.
- SplitsPy NeighborNet split weights are too slow at 529 taxa (> 11 min at 100); the cycle
  alone takes ~20 s.
- The group SLURM memory limit (`AssocGrpMemLimit`) blocked a 240 GB request; 64 GB ran.
- Nextflow GitHub auth is in `~/.nextflow/scm` (set up by the parallel session).
