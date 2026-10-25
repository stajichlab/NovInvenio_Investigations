# yeast_convergence: run log

What was run, in order, with job IDs and results. Dates are 2026-10-23 to 2026-10-25.
Design and caveats: `README.md` and `PLAN.md`. Comparison provenance:
`analysis/convergence/PROVENANCE.yaml`. Open items: PR #43 and issue #44.

## 1. Sources for the 106 species

- `pixi run python bin/ni resolve` filled rows with blank sources (UniProt first, NCBI second).
  Plain `python` has no `datasets` CLI and fails every NCBI lookup.
- Rows resolve could not fill were set by hand to BFD annotations
  (`/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation/<strain>/predict_results/`),
  with `local_faa`, `local_genome`, `local_gff3`.
- Five species were replaced (table in README, "Accessions").
- Taxon IDs from BFD `.gbk` headers (7) were checked against NCBI with `datasets`; all match.

## 2. BUSCO of all proteomes

- `bin/busco_manifest.py` then `bin/run_busco_proteomes.sh` (BUSCO 5.8.0, protein mode, `fungi_odb10`, 758 groups).
- Jobs 30326454, 30326455, 30326456 (106 proteomes); 30335621, 30335622 (6 candidate genomes).
  Two earlier jobs failed (30324081: `module load busco` calls `conda activate --stack`; 30324462: CRLF in the manifest).
- Result: `analysis/busco/summary.tsv`, with `summary.tsv.provenance.yaml`.
  Median 96.2% complete. Six proteomes under 85% are kept on purpose (README).

## 3. Config build

- `bin/build_study_config.py --study-dir ...` for all 16 contrasts (`config.csv`, `DATA_MANIFEST.yaml`, `data_dir/`, `annotations/`).

## 4. Pipeline runs (pairwise diamond)

- `bin/run_all.sh <ids>`; one head job per contrast, `NII_ANNOTATE_PFAM=1`.
  Pipeline: `nf_NovInvenio` commit d03f5a2 (`/rhome/jstajich/projects/NovInvenio-worktrees/main-d03f5a2`), Nextflow 26.04.6.
- Smoke runs: 30342586 (`yc_saccharomycotina`), 30342587 (`yc_pseudo_dermatophyte`). Pfam: 30342588, 30342589.
- Other 13 contrasts: 30420551 to 30420563 (about 1 h 20 min each). All completed, no failed task.
- `yc_pseudo_lichen` was not run (proteomes not checked for algal or bacterial contamination).
- Candidate counts that look unusual are explained in README (Neolecta 0, Ustilaginomycotina 10, AMF 31,206).

## 5. Enrichment

- `sbatch bin/run_enrichment_all.sh`. It needs `--constraint=milan`: the pixi `mmseqs` fails on nodes c01 to c03 (no AVX2).
- Job 30439653 ran on c01 and every run failed (the script prints a warning and still exits 0).
  Job 30444374 skipped all runs because the failed runs had left `proteome_clustered/` folders.
  I deleted those 30 folders (each held only `coverage.tsv`).
- Job 30452604 (milan, 34 min 29 s) completed all 15 contrasts, novelty and loss.

## 6. Comparison

- `/usr/bin/python3.12 bin/convergence_compare.py` (defaults), then `bin/write_comparison_provenance.py`.
- Output in `analysis/convergence/`. 9 origin, 4 pseudo (the fifth, lichen, has no run) and 2 reversal contrasts have tables.
- Result: 0 yeast-convergent terms. 4 streamlining-generic terms (F-box x3 in novelty; RNase_H clan in loss).
  This does not show that convergence is absent (README and PLAN list the limits).
- Shared outgroups of `yc_saccharomycotina` and `yc_taphrinomycotina`: neither run makes a recurrent call. No leave-shared-out rerun (PLAN section 6, A).

## 7. Publish

- Release assets, per contrast: `bin/publish_alignment_release.sh fungi/<id>/pairwise-diamond` and
  `bin/publish_report_release.sh fungi/<id>/pairwise-diamond`, both with `NII_SKIP_DEPLOY_TRIGGER=1`.
  30 releases (`alignments-fungi-<id>--pairwise-diamond`, `reports-fungi-<id>--pairwise-diamond`). Log: `logs/yc_publish/publish.log`.
  I did not use `bin/publish_all_studies.sh`: it publishes every study that has a `run.json`.
- Landing pages (class 2: `report.html`, `index.html`, `study.json`, `run.json`, `alignment.html`) were committed to `main` (4351fdf) from a worktree of `origin/main`.
- Gallery: `bin/generate_docs.py` in a clean worktree from `origin/main` with this folder's `species.csv` files added; pushed to `main` (2aa6043). The fungi page lists 32 studies.
- A push to `main` starts the Pages deploy. I did not check that the deployed site shows the 15 contrasts.

## 8. Not done

- `yc_pseudo_lichen`: no run, no enrichment, not in the comparison.
- Growth-form codings are not in `config_support/traits/` (issue #44). Three codings conflict with `traits.csv`.
- 32 species are not in the NovInvenio master pool (not added).
- Merge of PR #43. `studies/` and `bin/` changes are on the branch `yeast-convergence-contrasts`, not on `main`.

## 9. After the first comparison (2026-10-25)

- Enrichment pages, figures and the comparison page were published to `main` (docs only).
- Review fixes: `run_enrichment_all.sh` redoes a folder that has no `PROVENANCE.yaml`; control
  figure titles name the control type and focal species; the comparison page takes its counts
  from the files.
- Shared-code PRs: #45 (nested studies in `generate_enrichment_page.py`), #46 (TE flags in
  figures and pages, plain notice for 0 candidates), and this follow-up (related links and a
  gallery page card for the comparison page).
- Lichen screen: job 30474764 was cancelled (it ran on 1 thread). Job 30478145 uses 32 threads.
- `yc_rev_neolecta` has 0 novelty candidates because a one-species ingroup cannot yield any
  (see README).

