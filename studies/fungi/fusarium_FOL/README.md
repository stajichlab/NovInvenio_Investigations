# Fusarium oxysporum f. sp. lycopersici (FOL) pangenome

Single-direction `pangenome.nf` study: 21 FOL strains (IN) and 4 FOSC
outgroup strains (OUT), all from BFD funannotate annotations (one annotation
method). Replaces the dropped FOXY vs FSSC cross-complex design
(`../fusarium_FOXY_vs_FSSC/`, left as is).

Two runs:

| run | params file | launch dir | outdir |
|---|---|---|---|
| 1, rescue OFF | `params_norescue.yaml` | `.nf_launch/mmseqs_norescue` | `results/mmseqs_norescue` |
| 2, rescue ON | `params_rescue.yaml` | `.nf_launch/mmseqs_rescue` | `results/mmseqs_rescue` |

Start or resume: `sbatch .nf_launch/mmseqs_<mode>/submit_nextflow_head.sh`
(calls `run_pangenome.sh <mode>`).

## Pipeline

Detached worktree `/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/fusarium-runs-af6fd68`,
NovInvenio commit af6fd68b69311f75417e06a962c906cfe736af41 (same as the FOXY
vs FSSC study; origin/main was a72b4bc on 2026-09-25). Nextflow 26.04.6.
mmseqs tier-1 0.9 / 0.8. Island Pfam HMM:
`/bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm`.

Params files (not CLI flags; a CLI boolean `false` is the truthy String
"false", NovInvenio issue #191):
- `pangenome_rescue_enable` false / true.
- `pangenome_pfam_chunk_size: 2000` (default 3000; see queue config).
- `pangenome_mash_threshold: 0.0005` (default 0.001). The study's clone rule
  is mash < 0.0005. At 0.001 the pipeline would also merge
  {Fol002, Fol004, Fol026, Fol073}, {Fol016, Fol072, Fol075} and
  {VPRI44100, VPRI44274, VPRI44277} (BFD mash 0.0006-0.00098).

`stajichlab_queue.config`: copy of the FOXY vs FSSC one plus
FAMILY_PFAM_SCAN on epyc with a 4 h first attempt (17 of 62 chunks hit the
2 h short-queue limit there).

## Strain selection (`bin/select_fol_strains.py`)

Source: the FOXY vs FSSC inventory and all-vs-all mash
(`../fusarium_FOXY_vs_FSSC/selection/`), BUSCO genome from BFD.

Ingroup:
1. SPECIES_IN contains "f. sp. lycopersici", not radicis-lycopersici: 28.
2. BUSCO genome C >= 90%: 28 pass (98.6-99.6%).
3. Protein floor 11,000: 28 pass (13,800-16,388). FOXY median 14,956,
   MAD 328; median - 3 MAD = 13,972 would drop Fol59 (13,800 proteins, N50
   3.0 Mb) and VPRI44276, so the 3-MAD floor was not applied.
4. Clone dedup, single linkage at mash < 0.0005, keep highest N50: 21 kept.
   Dropped: Fol074 (0.000498 to Fol038); VPRI44099, VPRI44107, VPRI44272,
   VPRI44275, VPRI44276 (to VPRI44100; 44276 is at distance 0, and 44107 =
   44272 at distance 0, same N50 and length, i.e. duplicate assemblies);
   VPRI44273 (to VPRI44274). Log: `selection/fol_selection_log.tsv`.

Outgroup rule: FOSC-named, not lycopersici, BUSCO >= 90%, proteins >=
11,000, N50 >= 1 Mb, N50 not within 0.1% of Fol4287's 4,589,937 (several
BFD strains share it: likely scaffolded onto Fol4287), and mash >= 0.010 to
every one of the 28 FOL strains. Take the 3 nearest (by minimum distance to
FOL) with distinct named formae speciales, plus the highest-N50
*F. odoratissimum*:

| strain | f. sp. / species | min mash to FOL | N50 |
|---|---|---|---|
| GL1080 | fragariae | 0.0107 | 4.27 Mb |
| AJ275 | narcissi | 0.0110 | 4.08 Mb |
| LA1E | vasinfectum | 0.0113 | 3.33 Mb |
| 36102 | *F. odoratissimum* | 0.0243 | 4.86 Mb |

Note: within-FOL mash reaches 0.025 (RR23I-1.2 is 0.024 from every other
FOL strain), so the outgroup is outside the FOL clades, not outside the
whole FOL span. f. sp. lycopersici is polyphyletic; see the notes.

Build chain:
```
/usr/bin/python3.12 bin/select_fol_strains.py
/usr/bin/python3.12 bin/build_species_csv.py
<NI pixi python> bin/build_study_config_symlinks.py --study-dir $PWD --skip-fetch
```
`data_dir/{pep,dna,gff3}` are symlinks into BFD `predict_results/`.
TaxonGroup is empty; the pipeline fills ingroup clades from mash.

## Pre-launch checks (`.nf_launch/precheck_cluster/`)

Same mmseqs parameters, run by hand (job 29099749): 42,465 families.
Co-occurrence eligible (shell+cloud, >= 5 of 21 ingroup strains): 6,601
families, 21.8 M pairs (Coccidioides immitis-in: 25.1 M).
Rescue estimate (2.0e-6 cpu-h per query x genome Mb, from the Coccidioides
traces): mean 27,561 absent families per strain, 72.5 cpu-h. Under the
150 cpu-h limit, so run 2 was launched.

## Status (2026-09-25)

Both runs finished: run 1 (job 29099753, 39 min, 4.6 used cpu-h) and
run 2 (job 29099764, 2 h 01 min, 72.8 used cpu-h). Results and comparison:
`notes/2026-09-25-first-runs.md`.
