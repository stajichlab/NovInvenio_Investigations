# Cluster-vs-pairwise sensitivity: sordariales_shallow

## Analysis 1: controls recall/FP by tier

| tier | n_controls | n_positive | n_positive_resolved | positive_hits | recall | n_negative | n_negative_resolved | negative_fp | fp_rate | n_unresolved | n_placeholder_skipped |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 9 | 1 | 1 | 1 | 1.0 | 8 | 8 | 0 | 0.0 | 0 | 0 |
| C+H | 9 | 1 | 1 | 1 | 1.0 | 8 | 8 | 0 | 0.0 | 0 | 0 |
| C | 9 | 1 | 1 | 1 | 1.0 | 8 | 8 | 8 | 1.0 | 0 | 0 |
| R | 9 | 1 | 1 | 1 | 1.0 | 8 | 8 | 8 | 1.0 | 0 | 0 |


### Per-control outcomes

| control_id | class | expected_call | anchor_type | resolved_family | actual_call | outcome | note | tier |
|---|---|---|---|---|---|---|---|---|
| POS_EAS | positive | novel | protein_id | sp|Q04571|RODL_NEUCR | novel | hit |  | P |
| NEG_BUSCO01 | negative | core | busco | tr|V5IP41|V5IP41_NEUCR | not_novel | tn |  | P |
| NEG_BUSCO02 | negative | core | busco | sp|Q7SC06|RS25_NEUCR | not_novel | tn |  | P |
| NEG_BUSCO03 | negative | core | busco | tr|Q7RY29|Q7RY29_NEUCR | not_novel | tn |  | P |
| NEG_BUSCO04 | negative | core | busco | tr|Q7RWT5|Q7RWT5_NEUCR | not_novel | tn |  | P |
| NEG_BUSCO05 | negative | core | busco | tr|Q7SCW2|Q7SCW2_NEUCR | not_novel | tn |  | P |
| NEG_BUSCO06 | negative | core | busco | tr|Q7SI31|Q7SI31_NEUCR | not_novel | tn |  | P |
| NEG_BUSCO07 | negative | core | busco | tr|Q7SH99|Q7SH99_NEUCR | not_novel | tn |  | P |
| NEG_BUSCO08 | negative | core | busco | tr|Q1K5V9|Q1K5V9_NEUCR | not_novel | tn |  | P |
| POS_EAS | positive | novel | protein_id | tr|B2B0H0|B2B0H0_PODAN | novel | hit |  | C+H |
| NEG_BUSCO01 | negative | core | busco | tr|Q2H025|Q2H025_CHAGB | not_novel | tn |  | C+H |
| NEG_BUSCO02 | negative | core | busco | sp|Q7SC06|RS25_NEUCR | not_novel | tn |  | C+H |
| NEG_BUSCO03 | negative | core | busco | tr|B2AEZ4|B2AEZ4_PODAN | not_novel | tn |  | C+H |
| NEG_BUSCO04 | negative | core | busco | tr|F7VS92|F7VS92_SORMK | not_novel | tn |  | C+H |
| NEG_BUSCO05 | negative | core | busco | tr|F7W1A4|F7W1A4_SORMK | not_novel | tn |  | C+H |
| NEG_BUSCO06 | negative | core | busco | tr|B2AT65|B2AT65_PODAN | not_novel | tn |  | C+H |
| NEG_BUSCO07 | negative | core | busco | tr|Q7SH99|Q7SH99_NEUCR | not_novel | tn |  | C+H |
| NEG_BUSCO08 | negative | core | busco | tr|Q1K5V9|Q1K5V9_NEUCR | not_novel | tn |  | C+H |
| POS_EAS | positive | novel | protein_id | tr|B2B0H0|B2B0H0_PODAN | novel | hit |  | C |
| NEG_BUSCO01 | negative | core | busco | tr|Q2H025|Q2H025_CHAGB | novel | fp |  | C |
| NEG_BUSCO02 | negative | core | busco | sp|Q7SC06|RS25_NEUCR | novel | fp |  | C |
| NEG_BUSCO03 | negative | core | busco | tr|B2AEZ4|B2AEZ4_PODAN | novel | fp |  | C |
| NEG_BUSCO04 | negative | core | busco | tr|F7VS92|F7VS92_SORMK | novel | fp |  | C |
| NEG_BUSCO05 | negative | core | busco | tr|F7W1A4|F7W1A4_SORMK | novel | fp |  | C |
| NEG_BUSCO06 | negative | core | busco | tr|B2AT65|B2AT65_PODAN | novel | fp |  | C |
| NEG_BUSCO07 | negative | core | busco | tr|Q7SH99|Q7SH99_NEUCR | novel | fp |  | C |
| NEG_BUSCO08 | negative | core | busco | tr|Q1K5V9|Q1K5V9_NEUCR | novel | fp |  | C |
| POS_EAS | positive | novel | protein_id | tr|B2B0H0|B2B0H0_PODAN | novel | hit |  | R |
| NEG_BUSCO01 | negative | core | busco | tr|Q2H025|Q2H025_CHAGB | novel | fp |  | R |
| NEG_BUSCO02 | negative | core | busco | sp|Q7SC06|RS25_NEUCR | novel | fp |  | R |
| NEG_BUSCO03 | negative | core | busco | tr|B2AEZ4|B2AEZ4_PODAN | novel | fp |  | R |
| NEG_BUSCO04 | negative | core | busco | tr|F7VS92|F7VS92_SORMK | novel | fp |  | R |
| NEG_BUSCO05 | negative | core | busco | tr|F7W1A4|F7W1A4_SORMK | novel | fp |  | R |
| NEG_BUSCO06 | negative | core | busco | tr|B2AT65|B2AT65_PODAN | novel | fp |  | R |
| NEG_BUSCO07 | negative | core | busco | tr|Q7SH99|Q7SH99_NEUCR | novel | fp |  | R |
| NEG_BUSCO08 | negative | core | busco | tr|Q1K5V9|Q1K5V9_NEUCR | novel | fp |  | R |


## Analysis 2: genome-wide concordance (gains)

| baseline_count_vs_C | delta_pct_vs_C | delta_vs_C | jaccard_vs_P | n_candidates | precision_vs_P | recall_vs_P | refined_count_vs_C | tier |
|---|---|---|---|---|---|---|---|---|
| 23925 | -99.1 | -23717 | 0.1652 | 208 | 0.5288 | 0.1937 | 208 | C+H |
|  |  |  | 0.0074 | 23925 | 0.0076 | 0.3187 |  | C |
| 23925 | -0.0 | -5 | 0.0074 | 23920 | 0.0076 | 0.3187 | 23920 | R |


## Analysis 2: genome-wide concordance (losses)

| baseline_count_vs_C | delta_pct_vs_C | delta_vs_C | jaccard_vs_P | n_candidates | precision_vs_P | recall_vs_P | refined_count_vs_C | tier |
|---|---|---|---|---|---|---|---|---|
| 47765 | -98.4 | -47016 | 0.1693 | 749 | 0.4579 | 0.2117 | 749 | C+H |
|  |  |  | 0.0153 | 47765 | 0.0156 | 0.4593 |  | C |
| 47765 | -0.1 | -66 | 0.0153 | 47699 | 0.0156 | 0.4593 | 47699 | R |


## Analysis 3: real compute cost (wall-hours)

| tier | wall_hours | note |
|---|---|---|
| P | 0.101 | PARTIAL/lower-bound only -- true all-vs-all DIAMOND_SEARCH cost is unmeasurable: storeDir cache hit means Nextflow never logs a trace row for it (unlike ordinary -resume CACHED rows). This total covers only DIAMOND_SELF, PARSE_HITS, PARSE_SELF_HITS, TBLASTN(+MAKEDB), BUILD_PRESENCE_MATRIX -- do NOT read this as Tier P's full cost. |
| C+H | 31.287 | PROFILE_SEARCH:* (gain-side family-profile pathway), measured directly. |
| C+H_loss | 96.206 | PROFILE_LOSS_SEARCH:* (loss-side equivalent), measured directly. |
| C | 0.0 | Free byproduct of Tier C+H's own clustering step. |
| R | 0.011 | Measured directly with `time` around refine_ambiguous_families.py's diamond step -- not read from any trace file (plain subprocess call, no trace row). |


### Analysis 3 correction (hand-written, 2026-09-23)

The table above is `bin/trace_cost_report.py` output. Two of its rows are wrong for this clade:

- **Tier P is not a lower bound here.** The tool hard-codes the "storeDir cache hit, unmeasurable" note. For sordariales_shallow the all-vs-all diamond search ran fresh in trace `sordariales_shallow-20260922_101317`, so its trace rows exist. Measured from that trace (realtime, completed tasks):
  - gain `SEARCH:DIAMOND_SEARCH` + `DIAMOND_MAKEDB`: 16 tasks, 0.027 h, 0.642 cpu-h
  - loss `LOSS_SEARCH:DIAMOND_SEARCH` + `DIAMOND_MAKEDB`: 20 tasks, 0.064 h, 1.601 cpu-h
  - **total search: 0.091 task-h, 2.243 cpu-h.** Tier P overall is about 0.19 task-h (search plus the 0.101 h of other P steps in the table).
  - The resumed run `20260922_220343` re-executed the 8 loss-side DIAMOND_SEARCH tasks (0.063 h, 1.215 cpu-h). This is counted once, not summed.
- **Tier R omits the loss side.** Gain refine: 38.2 s (0.011 h). Loss refine: 365.9 s (0.102 h), with 2926/11895 families ambiguous.

Tier C+H by process (trace `sordariales_shallow_cluster-20260923_103407`, completed tasks, realtime):

| process | tasks | task-h | cpu-h |
|---|---|---|---|
| PROFILE_SEARCH BUILD_CHUNK | 164 | 21.63 | 610.6 |
| PROFILE_LOSS_SEARCH BUILD_CHUNK | 238 | 85.94 | 2634.2 |
| PROFILE_SEARCH HMMSEARCH_CHUNK | 24 | 9.58 | 30.4 |
| PROFILE_LOSS_SEARCH HMMSEARCH_CHUNK | 36 | 10.20 | 29.2 |

Failed attempts (26 preemptions, exit 143; 1 OOM, exit 137) added 3.96 task-h more. That time is not counted above.

**Result:** at this scale (4 ingroup, 8 outgroup proteomes), the Tier P search cost 0.091 task-h and 2.24 cpu-h. The Tier C+H family-profile stages cost about 127.5 task-h and 3305 cpu-h. The difference is about three orders of magnitude (about 1400x task-h, about 1470x cpu-h). Limits:
- This is one data point at N=12. Pairwise cost grows with N^2, so the result does not tell where a crossover (if any) occurs.
- The two runs did not use identical hardware. C+H workers ran on the mixed-CPU preempt partition.
- BUILD_CHUNK (famsa + hmmbuild) is 84% of C+H task-h. The cost is in building the profiles, not in searching with them.

Full per-row values are in `sordariales_shallow.cost_supplement.tsv`.

## Analysis 4: findings for this clade (hand-written, 2026-09-23)

1. **Controls (Analysis 1).** All four tiers recover the one positive control (eas/bli-7/ccg-2 hydrophobin). With n=1, recall does not separate the tiers.
2. **This is the first clade with a measured fp_rate.** The BUSCO map resolved all 8 negative controls. P and C+H: 0/8 false positives. C and R: 8/8 false positives. The C/R result is structural, not noise. Families are clustered from the seed (ingroup) group only, so cluster membership never contains outgroup proteins. Every family therefore looks outgroup-absent. Tier C and Tier R cannot make an absence call by themselves.
3. **Genome-wide concordance (Analysis 2).** C+H against P: gain precision 0.53, recall 0.19 (208 vs 568 candidates); loss precision 0.46, recall 0.21 (749 vs 1620). C+H misses about 80% of the P candidates in both directions. The per-miss clustering-vs-HMM classification was not done here either (same open item (d) as pezizo_set1).
4. **Tier R effect.** 838/8191 gain families (about 10%) and 2926/11895 loss families (about 25%) were ambiguous. R changed the candidate count by 5 (gain) and 66 (loss) against C, which is under 0.2%. The loss side is near the "about a third" threshold named in the pezizo_set1 report, but the effect is still negligible, because C/R cannot make the outgroup call at all (finding 2).
5. **Cost (Analysis 3).** At 12 proteomes, Tier P is both the reference and about 1000x cheaper than Tier C+H. For clades of this size, the data do not support using C+H in place of P.
6. **Method note.** The loss-side Tier R used a config with IN and OUT swapped, because `refine_ambiguous_families.py` has no `--query-group` option.
