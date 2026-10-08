# Pairwise pipeline cost: what the existing traces show

Date: 2026-10-07. Status: retrospective note. Nothing here was measured on purpose for a cost model. Empirical follow-up is deferred (see "Not measured").

## Context

On 2026-10-07 the maintainer decided to drop the family/HMM pathways (`--cluster_tool mmseqs`, `novelty_discovery`; NovInvenio #137, closed) and keep (a) the pangenome pipeline(s) and (b) the pairwise pipeline, with an upper bound on the number of taxa to avoid a compute explosion. A starting bound of **5,000 CPU-hours per run** was proposed. It is provisional and will be revisited.

## Data

`runs.tsv`: 19 pairwise runs under `results/*/nextflow_log/*trace.txt` that contain `DIAMOND_SEARCH`, joined to `studies/*/<run>/config.csv` (IN/OUT counts) and the size of each `data_dir/pep` file. Completed tasks only; CPU-hours = realtime x %cpu. When a run has several traces, all are summed. The sensitivity-benchmark runs (`*_dmnd_*`) are separate rows of the same three clades. Collected with a throwaway script (not kept).

## Findings

1. **Search cost is about proportional to W = (ingroup pep MB) x (all pep MB), so about quadratic in taxa.** CPU-hours per 1,000 W: agaricomycetes 0.49, pezizo_set1 0.74, sordariales_shallow 0.74, fusarium_FOXY 1.00 (19 IN / 4 OUT), cyanobacteria 1.22 (484 / 60). Small bacterial sets (9-14 taxa) show 6-30 because fixed overhead dominates. A log-log fit over the 9 default-mode runs gives exponent 0.78, R^2 0.75; the exponent is below 1 because of that overhead.
2. **Search dominates at large N, not at small N.** Search is 82% of total CPU-hours in fusarium_FOXY and 99% in cyanobacteria. For the 9-17 taxa bacterial sets total CPU-hours are 8-17 for 0.6-2 of search.
3. **Where 5,000 CPU-hours falls (search only, 80% ingroup, rates from fusarium and cyanobacteria):** about 205-230 taxa at 11 MB per proteome (fungal); about 985 taxa at 2.3 MB (bacterial). Other stages and diamond mode lower this. A working guess is 150-200 fungal proteomes; that figure was not computed.
4. **Scheduler tasks grow faster than compute.** cyanobacteria: 62,004 completed tasks, 58,478 of them `PARSE_HITS` at about 3.6 s each (about 2.4 CPU-hours), 206 failed. UHM_Koxytoca (34 taxa): 6,242 tasks for 26 search CPU-hours. No task-count formula was derived. This is probably the first limit hit.
5. **Diamond mode changes the rate by 1x to 2.8x.** Search CPU-hours default / sensitive / very-sensitive: sordariales_shallow 1.9 / 3.7 / 5.4; pezizo_set1 1.5 / 1.5 / 1.5; agaricomycetes 0.55 / 0.69 / 0.72. The shipped default is very-sensitive (NovInvenio #171). The mode used for fusarium_FOXY and cyanobacteria was not checked.

## Caveats

- Largest fungal run is 23 taxa. The 200+ taxa figure extrapolates the fungal rate by about 100x in W. The cyanobacteria run (bacterial, W 1.3 million) is consistent with it but is a different kind of proteome.
- W uses protein-file bytes (including headers), not residues.
- Hardware differs between runs. Failed and preempted attempts are not counted. Whether the cyanobacteria run completed was not checked.

## Not measured (left for later)

- A fungal scaling series (for example a Fusarium BFD subsample at N = 20, 40, 80, 160, very-sensitive, one commit, one node type, with one repeat for noise).
- TBLASTN, annotation and report costs at large N; candidate counts; `novelties.html` and alignment-shard sizes.
- Task count and scheduler overhead as a function of N; queue wait.
- What budget besides CPU-hours should bound a run (scheduler tasks, disk, wall clock).
