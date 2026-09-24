# Cluster-vs-pairwise sensitivity: summary across three clades

Date: 2026-09-23. Design: `notes/superpowers/specs/2026-09-13-cluster-vs-pairwise-sensitivity-design.md`.

This note compares Tier P (pairwise diamond, `--cluster_tool pairwise`) with Tier C+H (mmseqs families + family HMMs, `--cluster_tool mmseqs`) on three fungal clades. It also reports Tier C (cluster membership only) and Tier R (Tier C plus a paralog split of ambiguous families) where measured.

## Clades

| Clade | Ingroup | Outgroup | P run | C+H run |
|---|---|---|---|---|
| pezizo_set1 | 5 | 6 | `results/pezizo_set1` | `results/pezizo_set1_cluster` |
| agaricomycetes | 4 | 3 | `results/agaricomycetes_pairwise` | `results/agaricomycetes_mmseqs` |
| sordariales_shallow | 4 | 8 | `results/sordariales_shallow` | `results/sordariales_shallow_cluster` |

All runs use diamond (default sensitivity), `--ingroup_min_frac 0.75`, `--outgroup_min_frac 0.75` and `--loss_ingroup_max_frac 0.0`. All C+H runs use `--hmm_presence_cov 0.3 --hmm_presence_min_residues 100`. agaricomycetes_mmseqs also uses `--hmm_presence_domain_evalue 0.01`.

## Take-home

1. **At these scales, Tier P is much cheaper than Tier C+H.** Tier C+H used about 465x (agaricomycetes) and 1215x (sordariales_shallow) the CPU-hours of Tier P. For pezizo_set1, the Tier P search cost is not in any trace, so only a lower bound exists.
2. **Tier P was the most accurate tier on the curated controls.** It recovered every resolved positive control. It had 0 false positives in the one clade where negatives resolved. Tier C+H missed 3 of 5 resolved pezizo_set1 positives.
3. **Tier C and Tier R cannot be used alone.** Families are clustered from the seed group only, so no outgroup protein is ever in a family. Every family looks outgroup-absent. sordariales_shallow shows this directly: 8/8 BUSCO negatives were called novel.
4. **Tier P and Tier C+H agree poorly genome-wide** (Jaccard 0.08–0.21). The disagreement has two sides:
   - C+H misses: 23–57% are clustering-level. The protein's family was never profiled, usually because it is a singleton. The other 43–77% are HMM-level: the family HMM calls the family present in the other group.
   - C+H extras: 82–98% are proteins that Tier P finds in fewer than 75% of the seed group, or with no hit at all. The family HMM detects seed-group homologs that default-mode diamond does not.
5. **Tier P is the reference, but it is not ground truth.** The data here cannot say which method is right when the two disagree, except on the few curated controls. The C+H extras show that default-mode diamond misses homologs. This matches the open todo `diamond-very-sensitive-main-search.md` in nf_NovInvenio. Tier P with `--very-sensitive` would still cost far less than Tier C+H, but that has not been measured.

For comparisons of this size (7–12 proteomes), the data support using Tier P. They do not show where Tier C+H becomes cheaper, if it ever does (see "Cost", below).

## Analysis 1: curated controls

| Clade | Positives (P / C+H / C / R hits of resolved) | Negatives, false positives (P / C+H / C / R) |
|---|---|---|
| pezizo_set1 | 6/6, 2/5, 4/5, 4/5 | not measured (0/10 negatives resolved) |
| agaricomycetes | P unresolved (0/2); C+H 2/2, C 2/2, R 2/2 | no negatives |
| sordariales_shallow | 1/1, 1/1, 1/1, 1/1 | 0/8, 0/8, 8/8, 8/8 |

There are 9 positive controls in total. That is too few for a recall estimate, and the numbers are not pooled here. Per-control tables: `studies/fungi/pezizo_set1_cluster/tier_comparison/` and `studies/fungi/sordariales_shallow_cluster/tier_comparison/`.

## Analysis 2: genome-wide concordance and miss causes

Source: `cross_clade_concordance.tsv`, built from the per-clade `*.miss_summary.tsv` files (`bin/classify_concordance_misses.py`). All three clades were recomputed from the current result files.

| Clade | Dir | P | C+H | Shared | Precision vs P | Recall vs P | Misses: clustering / HMM | HMM misses with co-member paralog |
|---|---|---|---|---|---|---|---|---|
| pezizo_set1 | gain | 3393 | 1400 | 731 | 0.52 | 0.22 | 936 / 1726 | 235 |
| pezizo_set1 | loss | 127 | 39 | 13 | 0.33 | 0.10 | 46 / 68 | 33 |
| agaricomycetes | gain | 8553 | 5344 | 2425 | 0.45 | 0.28 | 2463 / 3665 | 1379 |
| agaricomycetes | loss | 234 | 88 | 25 | 0.28 | 0.11 | 119 / 90 | 31 |
| sordariales_shallow | gain | 568 | 208 | 110 | 0.53 | 0.19 | 231 / 227 | 21 |
| sordariales_shallow | loss | 1620 | 749 | 343 | 0.46 | 0.21 | 291 / 986 | 411 |

Miss classes, in the order the classifier checks them:

- **Clustering-level.**
  - `unprofiled_family`: the protein's family is below the size floor, usually a singleton. This is almost all of the clustering-level misses.
  - `oversized_family`: the family was flagged as oversized and not profiled.
  - `cluster_fragmentation`: family membership covers less than 75% of the seed group. This is rare, at most 19 per clade and direction.
- **HMM-level.**
  - `hmm_outgroup_paralog`: the HMM calls the family present in the other group, and the protein's registered paralog is in the same family. This is the HEX1-like case.
  - `hmm_outgroup_other`: the same outgroup call, with no co-member paralog. This is the ADA1/HAM5-like case.
  - `hmm_query_undercall`: membership covers the seed group but the HMM call does not. This was 0 in every clade and direction.

No miss or extra was left `unexplained` in any clade or direction.

**Raw diamond evidence for the HMM-outgroup misses** (`hmm_outgroup_miss_raw_diamond_evidence.tsv`). This is the best raw diamond E-value, from `search_cache/`, of the protein against any other-group proteome:

- 51–81% had no raw diamond hit at all. For these, the HMM finds something that diamond does not. The data here cannot tell whether that is real remote homology or a partial/domain-level HMM match.
- 10–43% had a raw hit with E < 1e-5. Tier P saw a significant hit, and its paralog-competition filter removed it. Tier P is making a paralogy call here, and Tier C+H is not.
- 6–12% had a weak raw hit (1e-5 ≤ E < 1e-2), below Tier P's cutoff.

The earlier committed pezizo_set1 concordance (`pezizo_set1.gain.concordance.tsv`, 2026-09-13 15:15: 1044 C+H gain candidates) is from an output set that the `pezizo_set1_cluster` rerun of 2026-09-13 20:40 replaced (now 1400). The figures above use the current files. Tier C and Tier R were not recomputed here. Their outgroup-blindness (take-home 3) makes their genome-wide precision uninformative.

## Analysis 3: compute cost

Source: `<clade>.cost.tsv` (`bin/trace_cost_report.py`, fixed 2026-09-23 to measure the pairwise search whenever its trace rows exist). Wall-h is the summed task realtime. cpu-h is realtime × %cpu.

| Clade | P wall-h | P cpu-h | C+H wall-h | C+H cpu-h | C+H / P (cpu-h) | P trace | C+H trace |
|---|---|---|---|---|---|---|---|
| pezizo_set1 | ≥0.116 | ≥1.02 | 89.7 | 2408 | ≤2372 (P search not in trace) | `20260910_170400` | `20260913_204012` |
| agaricomycetes | 0.120 | 1.89 | 37.4 | 878 | 465 | `20260910_194438` | `20260911_114647` |
| sordariales_shallow | 0.192 | 2.72 | 127.5 | 3305 | 1215 | `20260922_101317` | `20260923_103407` |

- The P totals include the self-search, TBLASTN and matrix building. The C+H totals include only the `PROFILE_SEARCH:*`/`PROFILE_LOSS_SEARCH:*` stages, so the comparison is slightly in C+H's favour.
- In sordariales_shallow, BUILD_CHUNK (famsa + hmmbuild) is 84% of C+H wall-h. The cost is in building the profiles, not in searching with them.
- C+H failed attempts are not counted (preemptions on the preempt partition; 3.96 task-h in sordariales_shallow).
- The P and C+H runs did not use identical hardware.

**Scaling (not measured).** The pairwise search grows with the number of proteome pairs. At N = 7 and N = 12 it cost 0.65 and 2.24 cpu-h. Whether Tier P stays cheaper at larger N depends on how Tier C+H's family-building cost grows. No run here measures that. A run at a larger N, or a Tier P run with `--very-sensitive`, would be needed to state a crossover.

## Open items

- A third method is needed to decide the HMM-outgroup misses with no raw diamond hit (for example TBLASTN evidence per cell, or a check of the HMM's aligned span).
- The pezizo_set1 Tier P search cost is still unmeasured. A timed rerun of its diamond all-vs-all would close this.
- Loss-direction positive controls still do not exist in any clade.
- `refine_ambiguous_families.py --query-group` (nf_NovInvenio issue #164) removes the IN/OUT-swapped config workaround for the loss-side Tier R.

## Files

- `README.md`: this note.
- `cross_clade_concordance.tsv`: the concordance and miss-class table above.
- `<clade>.<gain|loss>.miss_summary.tsv`: per-class counts.
- `hmm_outgroup_miss_raw_diamond_evidence.tsv`: the raw diamond E-value bins.
- `<clade>.cost.tsv`: per-tier costs.
- Per-protein miss detail (`*.miss_detail.tsv`) is class 3 (it lists candidate IDs). It is kept, gitignored, in each clade's `tier_comparison/` directory.
