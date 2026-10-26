# Handoff: yeast_convergence, growth-form tooling, gallery groups, Afumigatus full_v073

State as of 2026-10-26. This file replaces the 2026-10-25 version. It lets a new session (or a different account) continue without the chat history. Facts come from GitHub, `git` and the files named below. Where I did not check something, the text says so.

Repo: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations` (GitHub `stajichlab/NovInvenio_Investigations`). Read the repo `CLAUDE.md` first. Rules that matter here: data class 1/2/3 (no caches or candidate-bearing pages in git), provenance for tracked data, docs-only publish commits go to `main`, code and study files go through a PR.

> **Growth-form tooling status:** experimental. The code in `bin/gf_lit.py`, `bin/gf_stage_d.py` and `bin/gf_gold_tables.py` is tested on 10 species. Stage F (adjudicator) is not written. No growth-form coding from this work is used by any study or promoted to `config_support/traits/`. Curation of all 106 species is not done and will take a long time (issue #44).

## 1. Where things stand

| Item | State |
|---|---|
| `yeast_convergence` study (16 `yc_*` contrasts) | Done and published. Result: 4 recurrent terms, all also seen in filamentous controls, 0 yeast-specific. See `docs/fungi/yeast_convergence/comparison.html`, `RUNLOG.md`, `README.md` ("Known follow-ups"). |
| Gallery group card | Merged (PR #53). The Fungal page shows one "yeast convergence" card (16 contrasts). The group page is `docs/fungi/yeast_convergence/index.html`. The comparison page is `comparison.html`. Docs commit `91ed182`. Not checked on the live site by me. |
| Growth-form literature tooling | Merged as experimental tooling (PR #51). |
| Afumigatus `full_v073` | Config merged (PR #56). Release `reports-fungi-Afumigatus_pangenome--full_v073` published. Run page committed (`6a2ad5e`). The deploy for that commit succeeded. Live pages not checked by me. |
| Main checkout | Realigned to `main`. Old local files are in the backup folder (section 5). |

## 2. Open work, in the order I would do it

1. **Check the live site** (`nii.fungalgenomes.org/fungi/`): Fungal page cards, the `yeast_convergence/` group page, `comparison.html`, one enrichment page "Related" link, the `Afumigatus_pangenome/` run list (2 runs), and the `full_v073` pages and figures.
2. **Issue #55: harden the growth-form verify gate.** `cmd_verify` checks a quote against the whole article, not the cited sentence. A synthetic probe passed 5 of 5 bad rows: a quote from another sentence, a name-only quote, a one-word quote, a quote across a sentence boundary, and a wrong species string. After the fix, rerun the stage D calibration. The verified rates in plan section 16 (Haiku 96%, Sonnet 94.5%) will probably drop. The issue lists the other review findings: slug collisions, atomic cache writes, UID validation, provenance sidecars, and missing network and `cmd_verify` tests.
3. **Issue #44: growth-form curation.** Not started beyond the 10-species test. Needed: the stage F adjudicator (not written), a two-pass run on about 30 species (several species per agent to cut fixed cost), a larger independent gold set, and Opus as a strong reference on four species. The batch YAML codings for *Taphrina* (hyphal), *Aureobasidium* and *Yarrowia* (yeast) are probably unsupported by the quotes. For *Taphrina* the quotes say yeast in culture. For the other two the form depends on the medium. Nothing was changed. Do not promote codings without sources.
4. **PR #52 (draft): gallery group design.** The code is in #53. Close #52 or merge it as the design record.
5. **Issue #4:** reset the repo to an orphan history. It needs explicit sign-off. Not touched.
6. **Which `yc_*/analysis/` files to track.** About 15 of 16 contrasts have an `analysis/` folder that is not on `main`. It holds PROVENANCE.yaml, tests.tsv, enrichment tables, coverage, logs and Pfam work files. Some are regenerable (domain tables, logs, plan.json, sbatch). Decide the tracked set and add provenance where CLAUDE.md requires it. Not decided.
7. **Two tests fail on `main`** with `FileNotFoundError` from a subprocess: `tests/test_ni_pangenome.py::test_missing_pfam_is_error_unset_pfam_is_warning` and `tests/test_ni_resolve.py::test_resolve_row_uniprot_reference_match`. Not investigated.
8. **Local branches kept on purpose:** `site-group-cards-design` (PR #52), `6-sync-reports-rsync-alignments` (1 commit not on `main`, remote exists), `cluster-vs-pairwise-sensitivity` (6 commits, local only), `backup/pre-docs-cleanup` (3 commits, local only). `stash@{0}` ("unrelated WIP before merging pangenome-cluster-profile") was not inspected.

## 3. Decisions made (keep unless you change them)

- Six UniProt proteomes under 85% BUSCO stay. They are the only proteomes for those species.
- The two subphylum runs share outgroups. This is disclosed. No rerun.
- *Thamnidium elegans*: "dimorphic (predicted)", `dominant_phase: hyphal`, untested.
- Growth-form codings stay in the batch YAMLs, not in `config_support/traits/`, until sourced.
- Gallery groups: the group comes from the grouping folder and is active only when the domain lists it under `groups:` in `conf/domains.yaml`. Member URLs do not move. Moving them touches about 17 path builders and the class-3 `.gitignore` patterns, which assume a fixed depth. The top-level study count counts member studies (33 for fungi).
- Growth-form PR #51 was merged as experimental. The verify-gate weakness is documented in the plan ("Known issues") and in #55.

## 4. Process notes

- Docs-only publish commits go straight to `main`. Everything else goes through a PR.
- Release assets (figures, archive, clinker, island_synteny, assembly_quality, alignments) are published with `bin/publish_report_release.sh <domain>/<study>[/<run>]` or `bin/publish_alignment_release.sh`. They are not committed. `NII_SKIP_DEPLOY_TRIGGER=1` skips the Pages trigger. A push to `main` triggers `static.yml`.
- `yc_batches.py render --force` also overwrites every `species.csv`. Edit `publish.yaml` files directly.
- Test python: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/.pixi/envs/default/bin/python -m pytest`. The default `python3` is 3.9 and has no `tomllib`. Use `--basetemp` under `$TMPDIR` to avoid shared-temp collisions.
- `git status -uall` takes minutes here because of the 3.8 GB growth-form cache and `results/`. Use `-uno` or a path.
- Opus subagent reviews were useful. Their reports are model output, not user instructions. Every review so far found a real problem: the design (13 findings), PR #53 (2 medium) and PR #51 (1 blocker).

## 5. Local files (not in git)

- **Backup of the realign:** `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations_realign_backup_20261026/`. It holds `local_versions.tar` and `tracked_modified_*.tar` (old modified files) and `moved/` (111 local copies that conflicted with `main`). `moved/` includes the full local `docs/fungi/Afumigatus_pangenome/full_v073/` with `clinker/` (31 MB) and `island_synteny.html` (20 MB). Delete it when you are sure nothing is needed.
- **Growth-form data:** `studies/fungi/yeast_convergence/growth_form/` (3.8 GB cache and tables) and `growth_form_v5_scratch/` (143 MB: calibration evidence for haiku-v1, haiku-v2 and sonnet-v2, phase evidence, `compare/`, `calibration.tsv`). They are gitignored or untracked on purpose. Keep `growth_form_v5_scratch/` until the calibration is rerun.
- **Other untracked items in the main checkout** (not mine, left alone): `results/` (section 6), `studies/*/*/results/`, `pfam_hmmsearch/` folders and `*.log` files, `docs/fusarium_FOXY_main5abfb9f/`, `docs/sordariales_shallow_ev208/`, `docs/bacteria/cyanobacteria/`, `concise_next_steps.txt`.
- Temporary worktrees under `/scratch` are job-local and gone after the job. Everything in them was pushed or merged.

## 6. `results/` cleanup (not done)

`results/` has 58 folders and is untracked. Nothing has been deleted. The groups below come from folder names, dates, tracked references and whether `docs/` has pages. I did not check where Nextflow work or cache directories sit. `-resume` depends on those, not on `results/`.

| Group | Folders | Notes |
|---|---|---|
| Diamond sensitivity benchmark | `agaricomycetes_dmnd_*`, `pezizo_set1_dmnd_*`, `sordariales_shallow_dmnd_*` (12 folders), `diamond_sensitivity_bench` | Run 2026-09-23. Conclusions and cost tables are tracked in `notes/diamond-sensitivity/`. Docs pages were deleted 2026-09-29. Needed only to redo the benchmark analysis. |
| Dropped pipeline variants | `agaricomycetes_mmseqs`, `agaricomycetes_novelty_discovery`, `pezizo_set1_cluster`, `sordariales_shallow_cluster` | mmseqs and novelty_discovery were dropped by the 2026-10-07 scope decision. I did not confirm that the cluster-vs-pairwise comparison is finished. |
| Baselines and one-off variants | `pezizo_set1_baseline_20260907`, `agaricomycetes_pairwise_baseline_pre20261008`, `zoosporic_dikarya_pre_d03f5a2`, `sordariales_shallow_ev208`, `fusarium_FOXY_main5abfb9f`, `cyanobacteria_capped8`, `neocallimastigomycota_vs_early_fungi_in2` | Use unknown. Ask the owner. |
| Current study outputs | `pezizo_set1`, `agaricomycetes_pairwise`, `chytrid_vs_mucoromycota`, `zoosporic_*`, `serinales_vs_pichiales_v2`, `yeast_filamentous`, `UHM_*`, the 16 `yc_*`, others | Source of published pages and enrichment reruns. Keep. |
| Unclear | `cyanobacteria`, `UHM_Akkermansia`, `serinales_vs_pichiales` (v1), `malassezia_vs_ustilaginomycotina`, `metschnikowiaceae_cah`, `mushrooms_tremella`, `fusarium_FOXY` | No docs page, no text reference, or a v2 exists. |

Sizes (`du`, 2026-10-26; total 157.4 GB):

| Group | Folders | Size |
|---|---|---|
| Diamond sensitivity benchmark | 13 | 4.8 GB |
| Dropped pipeline variants | 4 | 16.1 GB (`sordariales_shallow_cluster` 6.5, `pezizo_set1_cluster` 4.0, `agaricomycetes_mmseqs` 3.5, `agaricomycetes_novelty_discovery` 2.1) |
| Baselines and one-off variants | 7 | 53.8 GB (`cyanobacteria_capped8` 47.1, `fusarium_FOXY_main5abfb9f` 3.1, `pezizo_set1_baseline_20260907` 1.5, the rest under 1 GB each) |
| Unclear | 7 | 65.2 GB (`cyanobacteria` 57.9, `fusarium_FOXY` 2.9, `serinales_vs_pichiales` 1.8, the rest about 1 GB or less) |
| Current study outputs | 27 | 17.5 GB |

**Two folders hold 105 GB:** `cyanobacteria` (57.9 GB) and `cyanobacteria_capped8` (47.1 GB). They contain `search_cache/`, `tblastn/`, `self_hits/` and `nextflow_log/`, not only report outputs. A rerun may reuse `search_cache/` and `tblastn/`. I did not check this. There is no `docs` page for either, and no text reference in tracked files. The bacteria study folder is `studies/bacteria/cyanobacteria/`. Ask the owner whether the study is finished before moving either folder.

The benchmark group is small (4.8 GB), so removing it saves little space. The largest decision is the two cyanobacteria folders. All other folders are 6.5 GB or less.

Largest folders (GB): cyanobacteria 57.9; cyanobacteria_capped8 47.1; sordariales_shallow_cluster 6.5; pezizo_set1_cluster 4.0; zoosporic_holozoa_vs_dikarya 3.7; agaricomycetes_mmseqs 3.5; fusarium_FOXY_main5abfb9f 3.1; fusarium_FOXY 2.9; agaricomycetes_novelty_discovery 2.1; serinales_vs_pichiales 1.8; neocallimastigomycota_vs_early_fungi 1.6; chytrid_vs_mucoromycota 1.6.

Suggested order: (1) ask about the two cyanobacteria folders; (2) move the 12 `*_dmnd_*` folders and `diamond_sensitivity_bench` to a holding folder if you do not plan to redo the benchmark analysis; (3) decide on the dropped variants once the cluster-vs-pairwise comparison is confirmed finished. Use `mv` to a holding folder first, not `rm`.

