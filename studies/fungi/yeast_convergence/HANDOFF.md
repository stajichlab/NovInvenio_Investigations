# Handoff: yeast_convergence and the growth-form literature work

State as of 2026-10-25. Counts come from the run logs and PR #51 tests. Where a figure is a first-pass estimate, the text says so.

## 1. Yeast_convergence study: done

All 16 `yc_*` contrasts ran, including the lichen control. Enrichment, figures, pages, the comparison page, the gallery card and the release assets are published. PRs #43, #45, #46, #47, #48, #49 and #50 are merged.

- **Result:** 4 recurrent terms across the 9 yeast origins (F-box in 3 Pfam/clan forms in novelty; RNase_H clan in loss). All 4 are also called in at least one of 5 filamentous controls. 0 yeast-specific convergent terms. This does not show that convergence is absent. The caveats are on the comparison page.
- **Read it:** `docs/fungi/yeast_convergence/index.html` (site: `nii.fungalgenomes.org/fungi/yeast_convergence/`). Run history: `RUNLOG.md`. Open follow-ups: `README.md`, "Known follow-ups".
- **Decisions:**
  - Six UniProt proteomes under 85% BUSCO stay, because they are the only proteomes for those species.
  - The two subphylum runs share outgroups. This is disclosed. No rerun, because neither run contributes a recurrent term.
  - *Thamnidium elegans* is coded "dimorphic (predicted)", `dominant_phase: hyphal`, untested.
  - Growth-form codings stay in the batch YAMLs, not in `config_support/traits/`, until sourced (issue #44).

## 2. Growth-form literature pipeline (issue #44): built, not finished

Goal: source the growth form of all 106 species from literature, with a verifiable quote for each.

| Stage | What | Status |
|---|---|---|
| A names | synonyms from the local NCBI `names.dmp` | done |
| B search | PubMed and PMC queries: dimorph, review, growth, name | done |
| B2 selection | merge duplicates (normalized PMCID); hybrid of score order and abstract-screen order | done |
| C sentences | species, abbreviation, genus and article-level sentences | done |
| D extraction | cheap model reads chunks of 200 sentences, returns quoted rows | tested on 10 species |
| E verify | script checks each quote is verbatim in the source | done |
| phase pass | second prompt on verified rows: normal phase in culture | tested on 5 species |
| `compare` | agreement of passes per species; flags what needs a person | done |
| F adjudicator | final decision per species | not written |

Code (PR #51): `bin/gf_lit.py`, `bin/gf_stage_d.py`, `bin/gf_gold_tables.py`, 59 tests, prompts, two gold sets. Plan: `GROWTH_FORM_LITERATURE_PLAN.md` (multi-agent protocol: section 14; gold sets: 15; calibration: 16). The first `gf_lit.py` (about 2,100 lines) came from another agent. PR #51 puts it in git. It needs a full code review.

## 3. What the tests showed

- **Recall.** Pilot gold set (19 rows, built from my own pilot sources, so a regression test): 17 of 19 selected. Independent gold set (82 rows, 29 species): "about" papers selected 93% of the time outside the three most studied species. For those three (*C. albicans*, *C. neoformans*, *S. cerevisiae*) selection is weak. Per species, at least one good source is selected for 17 of 20.
- **Extraction.** Haiku: 96% of 789 rows verified. Sonnet: 94.5% of 543. 367 of 368 first-pass quotes were verbatim. Gold sources found: Haiku 12 of 13, Sonnet 10 of 13. Class call matches the reference on 9 of 10 species for both. I changed the class rule after seeing results, so this figure is not independent.
- **Phase.** Haiku leaves the stated phase blank. The separate phase prompt fixes this. The two models agree on 4 of 5 species. Only *Mucor* (hyphal) has a firm reference.
- **Possible effect on #44.** If the quotes hold, three batch `dominant_phase` codings are unsupported: *Taphrina* (quotes say yeast in culture), *Aureobasidium* and *Yarrowia* (depend on the medium). No file was changed.
- Row precision is measured only on a 12-row hand sample. Ten species is too few to choose models for 106.

## 4. Mistakes found and fixed

- "20 quotes" checked was 17 (corrected in #44 and the plan).
- PMC ID merge bug ("PMC123" vs "123"). It made the first two-step run look worse.
- Two gate bugs rejected valid rows (article/genus rows; abbreviation sentences).
- Agent counts and attribution of one hand-back were wrong. Corrected.
- Study folder was pushed to `main` once by mistake. Reverted (`f390435`).

## 5. Not done

- Stage F adjudicator.
- Opus as a strong-model reference.
- Any run on more than 10 species.
- A larger independent gold set.
- Table-row extraction in the sentence stage.
- Three species have no candidate sentence (*Cladonia rubrotincta*, *Melampsora larici-populina*, *Umbilicaria pustulata*). *Trichophyton rubrum* has no PMC full text.
- Nothing is promoted to `config_support/traits/`.

## 6. Where things live

- **In git:** this branch (`gf-lit-selection-gold`, PR #51), issue #44.
- **On `/bigdata`, untracked:** `growth_form/` (tables; a 3.7 GB cache, gitignored) and `growth_form_v2_scratch` to `growth_form_v5_scratch`. v5 holds the Stage D input, evidence passes (haiku-v1, haiku-v2, sonnet-v2), phase evidence, `compare/` and `calibration.tsv`. These results are not in git.
- **Temporary:** my checkout under `$TMPDIR` (`/scratch`). It is gone after the job. Everything in it is pushed.

## 7. Next steps (recommended)

1. Review and merge #51, with a close read of `gf_lit.py`.
2. Run the two-pass pipeline on about 30 species, several per agent to cut fixed cost.
3. Write the adjudicator, driven by `compare` flags (about 3 of 10 species needed one).
4. Have a person review proposed phase corrections before any file changes.
