# Diamond sensitivity benchmark for Tier P

Date: 2026-09-23. Status: run and analysed 2026-09-23 (results below).

Follows on from `notes/cluster-vs-pairwise/README.md` and the nf_NovInvenio todo `diamond-very-sensitive-main-search.md`.

## Question

Tier P (pairwise diamond) was the cheapest and most accurate tier on the curated controls. It still runs diamond in default (fast) mode. The cluster-vs-pairwise comparison found that 82–98% of the Tier C+H candidates that Tier P rejects are proteins where default diamond found few or no seed-group homologs. This benchmark measures what `--sensitive` and `--very-sensitive` change:

1. **Cost:** search wall-hours and CPU-hours per mode.
2. **Candidates:** how many gain and loss candidates each mode gains or loses against default mode.
3. **Agreement with Tier C+H:**
   - `ch_extras_recovered`: of the C+H candidates that default P rejects, the fraction this mode calls.
   - `ch_misses_resolved`: of the default-P candidates that C+H rejects, the fraction this mode drops.
4. **Controls:** positive-control recall and BUSCO negative false positives per mode.
5. **Risk:** a more sensitive search can find new weak or partial cross-hits (paralogs, shared domains) that remove real novelties. For each candidate a mode loses, the tool bins the best new other-group hit by E-value and query coverage.

## Design

| | |
|---|---|
| Clades | pezizo_set1 (5 IN / 6 OUT), agaricomycetes (4 / 3), sordariales_shallow (4 / 8) |
| Modes | default, `--diamond_sensitivity sensitive`, `--diamond_sensitivity very-sensitive` |
| Runs | 9 studies: `studies/fungi/<clade>_dmnd_{default,sensitive,very_sensitive}` |
| Inputs | each study symlinks its base study's `config.csv`, `data_dir`, `species.csv`, `annotations` |
| Pipeline | nf_NovInvenio worktree `NovInvenio-worktrees/bench-diamond-sensitivity`, detached at origin/main `a8b68df` |
| Fixed params | `--run_tool diamond --cluster_tool pairwise`; all other params default; no Pfam/SwissProt/modelorgs annotation |
| Hardware | workers on `preempt`; DIAMOND_SEARCH/SELF/MAKEDB pinned to `-C milan` (`conf/diamond_sensitivity_bench.config`) |
| Comparator | the existing Tier C+H runs (`results/<clade>_cluster` or `agaricomycetes_mmseqs`) |

**Why the default mode is rerun.** The existing Tier P results were made with older pipeline code (2026-09-10 and 2026-09-22). The pezizo_set1 search cost is also missing from every trace. Rerunning default mode with the same commit and hardware makes all three modes comparable.

**What changes between modes.** Only `DIAMOND_SEARCH` (the main pairwise search). `DIAMOND_SELF`, the paralog self-search, is hard-coded to `--very-sensitive` in every mode.

**Cost measurement.** `DIAMOND_SEARCH` no longer uses storeDir on this commit, so each run's trace holds its search rows. Only COMPLETED task realtime is counted. Preempted attempts are retried and not counted. The tool reads the newest trace in each run's `nextflow_log/`. If a head job is restarted, check that this trace is the complete one.

## How to run

```bash
# 1. once: the pinned worktree needs its own pixi environment
cd /bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/bench-diamond-sensitivity && pixi install

# 2. submit all 9 head jobs (head on stajichlab, workers on preempt)
bin/run_diamond_sensitivity_bench.sh            # or name individual studies

# 3. after all 9 finish
bin/analyze_diamond_sensitivity_bench.sh
```

## Outputs

In `notes/diamond-sensitivity/` (class 2, tracked):

- `<clade>.cost.tsv`: per mode, the search, self-search and Tier P wall/cpu hours, with the trace name.
- `<clade>.concordance.tsv`: per mode and direction, the counts and metrics listed under "Question".
- `<clade>.lost_evidence.tsv`: per non-default mode and direction, the lost candidates binned by best new other-group hit.
- `<clade>.<mode>.controls_summary.tsv`: `score_controls.py` summary.

In `results/diamond_sensitivity_bench/` (class 3, gitignored):

- `<clade>.lost_detail.tsv`: every lost candidate with its new hit.
- `<clade>.<mode>.controls_scored.tsv`: per-control results.

Tools:

- `bin/diamond_sensitivity_report.py`, with tests in `tests/test_diamond_sensitivity_report.py`.
- `bin/run_diamond_sensitivity_bench.sh`.
- `bin/analyze_diamond_sensitivity_bench.sh`.

## Decision criteria

State them before the results arrive:

- Adopt a mode as the Tier P default if both of these hold:
  - its controls results are no worse than default mode;
  - its lost candidates are mostly supported by full-length hits (qcov ≥ 60 and E < 1e-20), not by short partial hits.
- Cost matters less here. Default Tier P used 0.65–2.24 search cpu-h at these sizes, so even a 10× increase is small against Tier C+H (878–3305 cpu-h).
- If many lost candidates rest on low-coverage hits, the mode is removing real novelties through shared domains. In that case, pair it with the other-group coverage floor (`--other_coverage_floor_qcov 15`, issue #158) and re-measure.

## Worked example: Q7RXD0_NEUCR (sordariales_shallow)

This candidate shows the failure the benchmark targets. It was reported by the user as a likely false novelty: TBLASTN hits in every outgroup genome, but no protein hits.

- **What it is.** Q7RXD0 is 1997 aa. Residues 1 to about 1750 are a low-complexity region (60–69% S/T/G/P/A per 200-aa window). The last ~230 aa are a conserved domain. Its Chaetomium ortholog Q2HI71 is annotated "DUF7371 domain-containing protein".
- **Tier C+H called it novel.** Its mmseqs family has 2 members (Neurospora Q7RXD0, Sordaria F7W4B3). The 2111-column family HMM is mostly low-complexity sequence from two close relatives. It detects the C-terminal domain in the ingroup (Chaetomium, HMM columns 1980–2110, E 7e-29). It did not detect it in any outgroup protein. In Colletotrichum H1V9Q8, which default diamond found at E 1e-31, the HMM reported no hit.
- **Tier P did not call it novel.** Default diamond found one outgroup hit (Collhigg H1V9Q8, E 1e-31, 11% query coverage). Default diamond also missed the Chaetomium ortholog.
- **TBLASTN** hits all 8 outgroup genomes (E 1e-13 to 1e-42), always on query residues ~1750–1997.
- **Diamond sensitivity check** (full protein and C-terminal 1750–1997 against the 8 outgroup proteomes, E ≤ 0.01):

  | Mode | Outgroup proteomes with a protein hit |
  |---|---|
  | default | 1 (Collhigg) |
  | `--very-sensitive` | 5 (Clonrose, Collhigg, Cordmili, Metaanis, Vertdahl), E 1e-17 to 1e-34, all on the C-terminal domain |

  Conilign, Cryppara and Eutylata have TBLASTN hits but no protein hit in either mode. Their gene models may be missing or split there. This was not checked.

This is a false novelty in Tier C+H. Tier P rejected it, but on a single outgroup hit. The `tblastn_outgroup_hits` column already flags it in the report, because `MAKE_NOVELTIES` runs with `--skip_tblastn_filter`. The benchmark should show whether `--very-sensitive` catches this class of protein, conserved-domain proteins with long low-complexity regions, across whole proteomes.


## Results (2026-09-23)

All 9 runs completed on pipeline commit `a8b68df`, with the diamond steps on `-C milan`. The two sordariales sensitive-mode runs were first OOM-killed in the loss-direction `BUILD_PRESENCE_MATRIX` step (exit 137; default mode peaked at 2.7 GB). They were resumed with 16 GB for that step. Their searches had already completed, so costs are read from the first trace (`--trace`, chosen by the most COMPLETED `DIAMOND_SEARCH` rows).

**Cost** (search = `DIAMOND_SEARCH`+`MAKEDB`; `<clade>.cost.tsv`):

| Clade | default search cpu-h | sensitive | very-sensitive |
|---|---|---|---|
| pezizo_set1 | 1.54 | 1.56 | 1.54 |
| agaricomycetes | 0.57 | 0.70 | 0.73 |
| sordariales_shallow | 1.94 | 2.71 | 3.25 |

The search cost rises at most 1.7×. `--very-sensitive` returns about 2.4× more raw hits (Amega vs Afum: 17,248 → 41,705 lines), but each 32-thread query-vs-all task takes about 25 s instead of 18 s. Even the largest value is about 1000× below Tier C+H (878–3305 cpu-h).

**Candidates** (`<clade>.concordance.tsv`; gained/lost are against default mode):

| Clade | Dir | default | sensitive (gained/lost) | very-sensitive (gained/lost) |
|---|---|---|---|---|
| pezizo_set1 | gain | 3515 | 3222 (+1348/−1641) | 3193 (+1438/−1760) |
| pezizo_set1 | loss | 101 | 122 (+93/−72) | 141 (+113/−73) |
| agaricomycetes | gain | 9441 | 9933 (+3095/−2603) | 9961 (+3408/−2888) |
| agaricomycetes | loss | 215 | 242 (+148/−121) | 241 (+149/−123) |
| sordariales_shallow | gain | 568 | 531 (+227/−264) | 509 (+229/−288) |
| sordariales_shallow | loss | 1620 | 1363 (+658/−915) | 1330 (+682/−972) |

The candidate sets change a lot. A more sensitive search removes candidates (new other-group hits) and also adds them (new seed-group hits lift proteins over the 75% seed-group threshold).

**Agreement with Tier C+H.** The more sensitive modes move Tier P toward Tier C+H in every clade and direction:
- `ch_extras_recovered`: of the C+H candidates that default P rejects, sensitive modes now call 15–54%. For losses in pezizo_set1 this is 15%. Every other clade/direction is 41–54%.
- `ch_misses_resolved`: of the default-P candidates that C+H rejects, the sensitive modes also drop 38–83%.

**Controls** (`<clade>.<mode>.controls_summary.tsv`): identical in all three modes.
- pezizo_set1: 6/6 positives; negatives not resolved.
- agaricomycetes: positives not resolved in Tier P; 0/5 negatives false positive.
- sordariales_shallow: 1/1 positive; 0/8 false positive.

**Evidence behind lost candidates** (`<clade>.lost_evidence.tsv`): the best new other-group hit of each candidate a mode loses.

| Clade | Dir | Mode | Lost | qcov ≥ 60 and E < 1e-20 | qcov < 30 |
|---|---|---|---|---|---|
| pezizo_set1 | gain | sensitive | 1641 | 30% | 20% |
| pezizo_set1 | gain | very-sensitive | 1760 | 29% | 19% |
| agaricomycetes | gain | sensitive | 2603 | 27% | 20% |
| agaricomycetes | gain | very-sensitive | 2888 | 25% | 20% |
| sordariales_shallow | gain | sensitive | 264 | 19% | 25% |
| sordariales_shallow | gain | very-sensitive | 288 | 18% | 25% |
| sordariales_shallow | loss | very-sensitive | 972 | 47% | 9% |

The loss rows for pezizo_set1 and agaricomycetes are 36–37% full-length and 12–16% qcov < 30.

**Decision, by the rule stated before the runs:**
- The controls are no worse, so the first condition holds.
- The lost candidates are not mostly supported by full-length hits. Only 18–49% are, and 9–25% rest on hits covering less than 30% of the query. The second condition fails.
- The rule's action for this outcome is to pair the sensitive mode with the other-group coverage floor (`--other_coverage_floor_qcov 15`, issue #158) and re-measure.
- The step to run next: `--diamond_sensitivity very-sensitive --other_coverage_floor_qcov 15`, on the same three clades.
- Adoption as the default is not decided yet.


## Follow-up: very-sensitive + coverage floor qcov 15 (2026-09-23)

Studies `<clade>_dmnd_very_sensitive_cov15`, run as `--diamond_sensitivity very-sensitive --other_coverage_floor_qcov 15`, on the same pinned pipeline and hardware. All three completed (14–15 min each). The tables above were regenerated with this as a fourth mode.

**Effect of the floor.** The floor rejected 8,947–45,248 other-group hits per direction (`*.coverage_floor_rejections.tsv`). Because low-coverage hits no longer count as presence, it adds candidates:

| Clade | Dir | default | very-sensitive | very-sensitive + qcov 15 |
|---|---|---|---|---|
| pezizo_set1 | gain | 3515 | 3193 | 3625 |
| pezizo_set1 | loss | 101 | 141 | 281 |
| agaricomycetes | gain | 9441 | 9961 | 10452 |
| agaricomycetes | loss | 215 | 241 | 394 |
| sordariales_shallow | gain | 568 | 509 | 598 |
| sordariales_shallow | loss | 1620 | 1330 | 1501 |

- Controls: unchanged (pezizo_set1 6/6 positives; agaricomycetes 0/5 negatives false positive; sordariales 1/1 positive, 0/8 negatives false positive).
- Search cost: 0.86–3.37 cpu-h.

**What the dropped candidates rest on** (candidates in default mode that very-sensitive + qcov 15 drops; best other-group raw hit):

| Clade | Dir | Dropped | E < 1e-20, qcov ≥ 60 | E < 1e-20, qcov 30–60 | E < 1e-20, qcov < 30 | 1e-20 ≤ E < 1e-5 |
|---|---|---|---|---|---|---|
| pezizo_set1 | gain | 1647 | 31% | 9% | 1% | 58% |
| pezizo_set1 | loss | 71 | 38% | 17% | 1% | 42% |
| agaricomycetes | gain | 2701 | 26% | 7% | 2% | 62% |
| agaricomycetes | loss | 117 | 39% | 13% | 2% | 44% |
| sordariales_shallow | gain | 275 | 19% | 11% | 2% | 67% |
| sordariales_shallow | loss | 955 | 48% | 7% | 0% | 44% |

The pre-stated full-length condition still fails: 19–48% are full-length.

The earlier reading was that the drops come from low-coverage (shared-domain) hits. That reading was wrong. Strong hits with qcov < 30 are at most 2% of drops. The largest group, 42–67%, is **weak E-value** hits. A coverage floor cannot address those, which is why it changed little.

**Independent check: TBLASTN.** The default runs' `tblastn_summary.tsv` / `loss_tblastn_summary.tsv` hold, for every default-mode candidate, whether TBLASTN found it in any other-group genome. For the dropped candidates:

| Clade | Dir | strong-hit drops with a TBLASTN hit | weak-hit drops with a TBLASTN hit |
|---|---|---|---|
| pezizo_set1 | gain | 644/678 (95%) | 630/962 (65%) |
| pezizo_set1 | loss | 40/40 (100%) | 23/30 (77%) |
| agaricomycetes | gain | 897/960 (93%) | 988/1687 (59%) |
| agaricomycetes | loss | 55/63 (87%) | 29/52 (56%) |
| sordariales_shallow | gain | 86/86 (100%) | 144/184 (78%) |
| sordariales_shallow | loss | 527/527 (100%) | 307/422 (73%) |

Overall, 71–90% of the candidates very-sensitive drops have genome-level evidence that the gene is present in the other group. So most drops remove false novelties/losses. TBLASTN is not ground truth, since it can also hit shared domains, but it is a separate method on separate data. The part without support is 22–44% of the weak-hit drops. That is 10–26% of all drops (for example 332 of 1640 for pezizo_set1 gains).

**Where this leaves the decision (not made here):**
- The pre-stated rule (a full-length majority) is not met.
- The rule's proxy looks wrong: the drops are mostly weak-E-value hits, and those are mostly corroborated by TBLASTN.
- `--very-sensitive` changes results substantially and mostly in the direction independent evidence supports. The coverage floor at qcov 15 adds little on top.
- Two follow-ups would narrow the rest:
  - a stricter other-group E-value for presence;
  - spot-checking the TBLASTN-unsupported weak-hit drops.


## Follow-up: stricter other-group E-value (simulated, 2026-09-23)

`bin/simulate_other_evalue.py` recomputes each finished run's candidates as if other-group presence required E < T (T = 1e-10, 1e-20, 1e-30). Seed-group cells are left as the run called them. It uses the run's `presence_matrix.evalues.tsv`, which records the E-value of every hit that counted as presence. At the run's own 1e-5 it reproduces every run's `candidates.txt` / `loss_candidates.txt` exactly (12/12 checked). Output: `sim_other_evalue.summary.tsv`.

The judge is independent genome evidence: the fraction of candidates that have a TBLASTN hit in any other-group genome. The TBLASTN data are the default and very-sensitive runs' own summaries; candidates in neither are counted as unknown, and the fraction is taken over known candidates. For a true novelty/loss this fraction should be low.

| Clade | Dir | default | very-sensitive | + E < 1e-10 | + E < 1e-20 | + E < 1e-30 |
|---|---|---|---|---|---|---|
| pezizo_set1 | gain | 45% | 14% | 19% | 27% | 32% |
| pezizo_set1 | loss | 69% | 33% | 33% | 40% | 44% |
| agaricomycetes | gain | 29% | 9% | 13% | 17% | 20% |
| agaricomycetes | loss | 45% | 10% | 14% | 18% | 24% |
| sordariales_shallow | gain | 54% | 25% | 30% | 39% | 43% |
| sordariales_shallow | loss | 64% | 33% | 37% | 43% | 50% |

(The E-value columns apply to the very-sensitive run.)

**Results:**
- A stricter other-group E-value makes the candidate sets worse by this measure in every clade and direction, and the candidate counts grow sharply (at 1e-30: 1.7× for agaricomycetes gains up to 27× for pezizo_set1 losses, 141 → 3751). The weak other-group hits that very-sensitive finds are mostly real: removing them brings back candidates that TBLASTN contradicts. **The stricter-E-value follow-up does not warrant a change.**
- `--very-sensitive` at the standard 1e-5 cuts the TBLASTN-contradicted fraction by 2–4× against default mode (29–69% → 9–33%).
- The coverage floor (qcov 15) adds little on top.
- Default-mode candidate lists: 29–69% of candidates have an other-group TBLASTN hit. The pipeline reports but does not filter on TBLASTN (`--skip_tblastn_filter`). This measures how many default-mode candidates are likely false.
- Not done for the thresholded sets: controls scoring. The TBLASTN trend was the same everywhere.

**Recommendation from these data:** make `--diamond_sensitivity very-sensitive` the Tier P default, and keep `--evalue 1e-5`. This is a user decision. It changes results for every study, so it needs a `.living/decisions.md` entry and reruns (nf_NovInvenio todo `diamond-very-sensitive-main-search.md`).
