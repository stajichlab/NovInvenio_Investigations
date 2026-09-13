# Cluster-vs-pairwise sensitivity: pezizo_set1

## Analysis 1: controls recall/FP by tier

| tier | n_controls | n_positive | n_positive_resolved | positive_hits | recall | n_negative | n_negative_resolved | negative_fp | fp_rate | n_unresolved | n_placeholder_skipped |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 16 | 6 | 6 | 6 | 1.0 | 10 | 0 | 0 |  | 10 | 1 |
| C+H | 16 | 6 | 5 | 2 | 0.4 | 10 | 0 | 0 |  | 11 | 1 |
| C | 16 | 6 | 5 | 4 | 0.8 | 10 | 0 | 0 |  | 11 | 1 |
| R | 16 | 6 | 5 | 4 | 0.8 | 10 | 0 | 0 |  | 11 | 1 |


### Per-control outcomes

| control_id | class | expected_call | anchor_type | resolved_family | actual_call | outcome | note | tier |
|---|---|---|---|---|---|---|---|---|
| POS_HEX1 | positive | novel | protein_id | sp|P87252|HEX1_NEUCR | novel | hit |  | P |
| POS_LAH | positive | novel | protein_id | tr|V5IRA6|V5IRA6_NEUCR | novel | hit |  | P |
| POS_ADA1 | positive | novel | protein_id | tr|Q7SE74|Q7SE74_NEUCR | novel | hit |  | P |
| POS_HAM5 | positive | novel | protein_id | tr|V5IN79|V5IN79_NEUCR | novel | hit |  | P |
| POS_HAM8 | positive | novel | protein_id | tr|Q7SEK3|Q7SEK3_NEUCR | novel | hit |  | P |
| POS_SPA1 | positive | novel | protein_id | tr|Q7SI25|Q7SI25_NEUCR | novel | hit |  | P |
| NEG_BUSCO01 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO02 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO03 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO04 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO05 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO06 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO07 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO08 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO09 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| NEG_BUSCO10 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | P |
| POS_HEX1 | positive | novel | protein_id | tr|A0ACF5BXR3|A0ACF5BXR3_9PEZI | not_novel | miss |  | C+H |
| POS_LAH | positive | novel | protein_id |  | unresolved | unresolved | protein not in any profiled family | C+H |
| POS_ADA1 | positive | novel | protein_id | tr|A0ACF5CBF7|A0ACF5CBF7_9PEZI | not_novel | miss |  | C+H |
| POS_HAM5 | positive | novel | protein_id | tr|V5IN79|V5IN79_NEUCR | not_novel | miss |  | C+H |
| POS_HAM8 | positive | novel | protein_id | tr|J3K3I3|J3K3I3_COCIM | novel | hit |  | C+H |
| POS_SPA1 | positive | novel | protein_id | tr|J3K8R1|J3K8R1_COCIM | novel | hit |  | C+H |
| NEG_BUSCO01 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO02 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO03 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO04 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO05 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO06 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO07 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO08 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO09 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| NEG_BUSCO10 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C+H |
| POS_HEX1 | positive | novel | protein_id | tr|A0ACF5BXR3|A0ACF5BXR3_9PEZI | novel | hit |  | C |
| POS_LAH | positive | novel | protein_id |  | unresolved | unresolved | protein not in any profiled family | C |
| POS_ADA1 | positive | novel | protein_id | tr|A0ACF5CBF7|A0ACF5CBF7_9PEZI | novel | hit |  | C |
| POS_HAM5 | positive | novel | protein_id | tr|V5IN79|V5IN79_NEUCR | novel | hit |  | C |
| POS_HAM8 | positive | novel | protein_id | tr|J3K3I3|J3K3I3_COCIM | novel | hit |  | C |
| POS_SPA1 | positive | novel | protein_id | tr|J3K8R1|J3K8R1_COCIM | not_novel | miss |  | C |
| NEG_BUSCO01 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO02 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO03 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO04 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO05 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO06 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO07 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO08 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO09 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| NEG_BUSCO10 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | C |
| POS_HEX1 | positive | novel | protein_id | tr|A0ACF5BXR3|A0ACF5BXR3_9PEZI | novel | hit |  | R |
| POS_LAH | positive | novel | protein_id |  | unresolved | unresolved | protein not in any profiled family | R |
| POS_ADA1 | positive | novel | protein_id | tr|A0ACF5CBF7|A0ACF5CBF7_9PEZI | novel | hit |  | R |
| POS_HAM5 | positive | novel | protein_id | tr|V5IN79|V5IN79_NEUCR | novel | hit |  | R |
| POS_HAM8 | positive | novel | protein_id | tr|J3K3I3|J3K3I3_COCIM | novel | hit |  | R |
| POS_SPA1 | positive | novel | protein_id | tr|J3K8R1|J3K8R1_COCIM | not_novel | miss |  | R |
| NEG_BUSCO01 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO02 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO03 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO04 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO05 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO06 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO07 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO08 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO09 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |
| NEG_BUSCO10 | negative | core | busco |  | unresolved | unresolved | busco protein not in any profiled family | R |


## Analysis 2: genome-wide concordance (gains)

| baseline_count_vs_C | delta_pct_vs_C | delta_vs_C | jaccard_vs_P | n_candidates | precision_vs_P | recall_vs_P | refined_count_vs_C | tier |
|---|---|---|---|---|---|---|---|---|
| 21258 | -95.1 | -20214 | 0.125 | 1044 | 0.4722 | 0.1453 | 1044 | C+H |
|  |  |  | 0.06 | 21258 | 0.0656 | 0.4111 |  | C |
| 21258 | -0.0 | -7 | 0.06 | 21251 | 0.0656 | 0.4111 | 21251 | R |


## Analysis 2: genome-wide concordance (losses)

| baseline_count_vs_C | delta_pct_vs_C | delta_vs_C | jaccard_vs_P | n_candidates | precision_vs_P | recall_vs_P | refined_count_vs_C | tier |
|---|---|---|---|---|---|---|---|---|
| 9901 | -99.8 | -9877 | 0.0709 | 24 | 0.4167 | 0.0787 | 24 | C+H |
|  |  |  | 0.0029 | 9901 | 0.0029 | 0.2283 |  | C |
| 9901 | -0.2 | -15 | 0.0029 | 9886 | 0.0029 | 0.2283 | 9886 | R |


## Analysis 3: real compute cost (CPU-hours)

| tier | cpu_hours | note |
|---|---|---|
| P | 0.116 | PARTIAL/lower-bound only -- true all-vs-all DIAMOND_SEARCH cost is unmeasurable: storeDir cache hit means Nextflow never logs a trace row for it (unlike ordinary -resume CACHED rows). This total covers only DIAMOND_SELF, PARSE_HITS, PARSE_SELF_HITS, TBLASTN(+MAKEDB), BUILD_PRESENCE_MATRIX -- do NOT read this as Tier P's full cost. |
| C+H | 42.707 | PROFILE_SEARCH:* (gain-side family-profile pathway), measured directly. |
| C+H_loss | 35.394 | PROFILE_LOSS_SEARCH:* (loss-side equivalent), measured directly. |
| C | 0.0 | Free byproduct of Tier C+H's own clustering step. |
| R | 0.032 | Measured directly with `time` around refine_ambiguous_families.py's diamond step -- not read from any trace file (plain subprocess call, no trace row). |


## Analysis 4: decision framework

_Fill in by hand after reading the tables above — see the spec's Analysis 4 for the exact framework (escalate Tier R near-misses to Tier P; treat candidate-count inflation and HEX1-style clustering-inherent misses as the two conditions that would invalidate the framework)._
