# Diamond sensitivity benchmark for Tier P

Date: 2026-09-23. Status: constructed, not yet run.

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
