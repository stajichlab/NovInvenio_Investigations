# Fusarium oxysporum SC (FOXY) vs F. solani SC (FSSC) — reciprocal pangenome

Two reciprocal `pangenome.nf` runs over the same 328 strains:

| direction | IN | OUT | samplesheet |
|---|---|---|---|
| `FOXY_in_FSSC_out` | 209 FOXY | 119 FSSC | `config_FOXY_in_FSSC_out.csv` |
| `FSSC_in_FOXY_out` | 119 FSSC | 209 FOXY | `config_FSSC_in_FOXY_out.csv` |

Purpose:
1. Compare two species complexes that are more divergent than *C. immitis* vs
   *C. posadasii* (mash distance FOXY to FSSC anchor 0.22-0.31, saturated;
   see `selection/nearest_complex_anchor.tsv`).
2. Re-test the tier-1 clustering result of NovInvenio issue #132 on a
   recombining, accessory-chromosome-rich taxon (caveat (b) there). These
   runs are the 0.9/0.8 baseline only; no clustering sweep was run.

Why a new directory: `studies/fungi/fusarium_FOXY/` holds a 23-row
NCBI-annotation species.csv (2026-09-14) for a different design. This study
uses one annotation method (BFD funannotate) for every strain, so it is a
sibling directory. `fusarium_FOXY/species.csv` is unchanged.

## Data source

BFD: `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD`
- `genome_annotation/<dir>/predict_results/<dir>.proteins.fa`, `.gff3`,
  `.scaffolds.fa` (the genome the GFF3 refers to).
- Directory -> samples.csv row join is by LOCUSTAG (first protein ID prefix),
  not by directory name. Strain-less rows get a bare species directory name
  (812 "bare-fallback" Fusarium_oxysporum rows in
  `results/_ambiguous_basenames_report.csv.gz`), so the name is not a safe key.
- Assembly stats: `results/genome_stats/asm_stats`; BUSCO (genome mode,
  fungi_odb12): `results/genome_stats/BUSCO_genome`; mash sketches (k=21,
  s=10000): `results/mash_sketch/<ASMID>.fa.msh`.

## Strain selection (`bin/select_strains.py`, seed 20260925)

Inventory: `bin/inventory_bfd_fusarium.py` -> `selection/bfd_fusarium_inventory.tsv`.
2,070 Fusarium_*/Neocosmospora_* directories; 1,964 have proteins + gff3 +
genome and a unique LOCUSTAG match. All-vs-all mash over those 1,964:
`selection/bfd_fusarium_all_vs_all.mash_dist.tsv.zst`.

The mash matrix is NOT tracked in git (42 MB, regenerable). `bin/select_strains.py`
reads it from `selection/bfd_fusarium_all_vs_all.mash_dist.tsv.zst` and needs it to run.
To rebuild it:
- Inputs: the BFD mash sketches `results/mash_sketch/<ASMID>.fa.msh` (k=21, s=10000)
  of the 1,964 inventory rows with `status=ok` in `selection/bfd_fusarium_inventory.tsv`.
- Output: `mash dist` over all pairs, compressed with `zstd`, saved under the name above.
- The exact command that made the committed-era file (dated 2026-09-25) was not
  recorded in this repo. An untested outline: `mash dist -p <threads> <all .msh files>
  <all .msh files> | zstd -T0 > selection/bfd_fusarium_all_vs_all.mash_dist.tsv.zst`.
  Check the column layout against what `select_strains.py` parses (around line 139)
  before relying on it.
- Record the `mash --version`, the sketch paths and the date in a note when it is rebuilt.

Complex membership needs BOTH a taxonomic name AND mash placement (nearest of
nine complex anchors, `ANCHORS` in the script):
- FOXY names: *F. oxysporum* + FOSC segregate species present in BFD
  (Lombard et al. 2019; Maryani et al. 2019): *odoratissimum,
  tardichlamydosporum, purpurascens, phialophorum, nirenbergiae,
  veterinarium*. Nearest anchor Fol4287 and d <= 0.035.
  Excluded by name even though mash places them next to Fol4287:
  *F. foetens*, *F. inflexum*, one *F. fujikuroi* (FUS01, likely mislabelled),
  *Fusarium* sp. TH15. *F. commune* (d 0.052) and *F. redolens* are outside.
- FSSC names: the `FSSC_NAMES` list in the script (Neocosmospora sensu
  Sandoval-Denis et al. 2019 / FSSC of O'Donnell et al. 2020, Geiser et al.
  2021) plus genus *Neocosmospora*. Nearest anchor *F. vanettenii* 77-13-4
  and d <= 0.17. Two names (*F. phalaenopsidis*, *F. drepaniforme*) are in
  the list on mash placement only; their literature placement was not checked.
  Excluded by name even though mash places them in FSSC: 12 unnamed
  *Fusarium* sp. strains and one *F. lateritium* (389579, likely mislabelled).
  Excluded as outside FSSC by mash (d 0.19-0.26 to the FSSC anchor, nearest
  anchor FSSC only because all anchors are saturated): *F. nematophilum,
  decemcellulare, albosuccineum, buxicola, staphyleae, dimerum, setosum*, etc.
  (they are not in the name list either).

Filters, counts (FOXY / FSSC):

| step | FOXY | FSSC |
|---|---|---|
| name match, files OK | 830 | 172 |
| BUSCO genome C >= 90% | 773 | 155 |
| mash placement | 762 | 155 |
| near-clone dedup (mash < 0.001, keep highest N50) | 499 | 119 |
| FOXY stratified subsample (target 200) + species floor | 209 | 119 |

- 11 *F. oxysporum*-labelled strains fail placement: 8 are nearest another
  anchor (6 FFSC, 1 FRSC, 1 FSSC), 3 are nearest Fol4287 but at d 0.053-0.054
  (09-002, F38, Tf1208).
- 61 of the 74 BUSCO drops (both complexes) are small, fragmented Ethiopian
  isolates (EtdFoc/EthFoc, PRJNA353637; e.g. 11-37 Mb, N50 3-5 kb).
  See `selection/selection_log.tsv` for every drop and its reason.
- The pipeline's own dedup (`pangenome_dereplicate`, mash < 0.001) only
  restricts the frequency denominator. It does not reduce clustering or
  rescue compute, so near-clones were removed before the run at the same
  threshold.
- FOXY strata: UPGMA on mash distances of the 499 dedup'd strains, cut at
  0.01 -> 89 strata (`selection/foxy_strata.tsv`). Allocation
  max(1, floor(200 * n_k / N)) plus largest remainders -> 207. Inside a
  stratum the draw takes one strain per distinct label first (forma
  specialis for *F. oxysporum*, species name otherwise), then random. Then
  *F. odoratissimum* and *F. purpurascens* (not drawn) were added -> 209.
  40 formae speciales are represented; 119 FOXY strains have no f. sp. label.
- FSSC is kept whole after dedup: 119 strains, 42 names.

Files: `selection/selected_strains.tsv` (one row per selected strain, with
N50, BUSCO, mash distances, source paths).

## Build chain

```
/usr/bin/python3.12 bin/inventory_bfd_fusarium.py           # ~16 min on GPFS
<NI pixi python> bin/select_strains.py                     # needs scipy
/usr/bin/python3.12 bin/build_species_csv.py                # species.csv
<NI pixi python> bin/build_study_config_symlinks.py --study-dir $PWD --skip-fetch
<NI pixi python> bin/build_reciprocal_configs.py
```

`bin/build_study_config_symlinks.py` runs NII `bin/build_study_config.py`
unchanged except that its `shutil.copyfile` is replaced with `os.symlink`,
so `data_dir/{pep,dna,gff3}` are symlinks into BFD (no 30 GB copy).
DATA_MANIFEST.yaml checksums are of the link targets; the records' text
still says "copied".

`TaxonGroup` is empty in species.csv. pangenome.nf fills ingroup rows with
its own mash clade labels. Coccidioides used the species name (one ingroup
clade), so pair-classification clade gating differs between the studies.

## Running

```
sbatch .nf_launch/mmseqs_norescue_FOXY_in_FSSC_out/submit_nextflow_head.sh
sbatch .nf_launch/mmseqs_norescue_FSSC_in_FOXY_out/submit_nextflow_head.sh
```

Both call `run_pangenome.sh <direction>`. Pipeline: detached worktree
`/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/fusarium-runs-af6fd68`
(NovInvenio commit af6fd68b69311f75417e06a962c906cfe736af41, main on
2026-09-25), with its own pixi env. Do not move or delete it while a run can
still be resumed. Resource overrides: `stajichlab_queue.config` (copied from
Coccidioides, plus attempt-scaled memory for other aggregate steps).
Island Pfam: `--pangenome_island_pfam_hmm
/bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm` (as in
Coccidioides). Outputs: `results/mmseqs_norescue_<direction>/`.

`RUN_RESCUE=true run_pangenome.sh <direction>` gives a separate
`mmseqs_rescue_<direction>` launch dir and outdir.

## Compute estimate (why rescue is off in the first runs)

Measured on the Coccidioides reciprocal runs (NovInvenio 91e3157, traces in
`../coccidioides_pangenome/results/rescue_freqpol_*/trace.txt`):
TBLASTN_PER_STRAIN used 1,147.5 and 1,139.8 cpu-h (realtime x %cpu), 529
tasks each, 33 min mean. Every other process together used < 12 cpu-h per run.
That run had 47,745 families, 39,235 absent families (tblastn queries) per
strain on average, and a 27.8 Mb mean genome. So the rescue cost is
about 2.0e-6 cpu-h per (query x genome Mb).

This study: 328 strains, mean genome 52.1 Mb (FOXY) / 51.3 Mb (FSSC).
Even at a low assumption of 20,000 absent queries per strain, one direction
costs 328 x 2.0e-6 x 20,000 x 51.8 = ~680 cpu-h, both ~1,360 cpu-h, over
the 800 cpu-h limit. The two directions run identical tblastn work (same
strains, same families), but the task hashes differ (the PRESENCE_MATRIX
command line carries the group labels), so nothing is shared.
The first runs therefore run with rescue disabled.

Measured after tier-1 clustering (FOXY_in_FSSC_out, 2026-09-25): 399,749
families (Coccidioides: 47,745), mean 14,485 present and 385,264 absent
families per strain. With each strain's real genome size the full rescue
estimate is 13,347 cpu-h per direction, about 26,700 cpu-h for both.

How rescue is disabled: `params_norescue.yaml` via `-params-file`.
Do not use `--pangenome_rescue_enable false` on the command line. Under
Nextflow 26.04.6 a CLI `false` arrives as the String "false", which is
truthy in `if (params.pangenome_rescue_enable)`. Tested with a two-line
script: CLI gives class String / TRUTHY, a YAML params file gives Boolean /
FALSY. The first launch (jobs 29084981/29084982) hit this and started
EXTRACT_ABSENT_QUERIES; it was cancelled before any TBLASTN task was
submitted.

## Status (2026-09-25)

Both directions ran to COOCCURRENCE and were stopped there. That step does
not scale to 399,749 families (792 M / 1,029 M candidate pairs). See
`notes/2026-09-25-first-runs.md` for the measured numbers, the blocker and
the options. Job IDs: 29084981/29084982 (attempt 1, cancelled),
29086546/29086547 (attempt 2, cancelled). A resume (`sbatch
.nf_launch/<dir>/submit_nextflow_head.sh`) reuses every finished task.
