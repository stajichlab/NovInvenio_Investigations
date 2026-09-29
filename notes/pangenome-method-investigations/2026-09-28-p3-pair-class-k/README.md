# #132 priority 3: pair-class window k (10) and linkage thresholds (0.5 / 0.05)

2026-09-28. `run_p3_sweep.sh` (SLURM job 29196151, 32 min, 16 CPUs, 4 settings at a time)
reran `bin/pangenome_pair_classification.py` alone, at the pinned commit a394ace, on the exact
inputs of two v0.7.0 PAIR_CLASSIFICATION tasks (work dirs in the script). `p3_summarize.py`
wrote `p3_summary.tsv`. Per-setting outputs (~1-18 MB each) are in
`studies/fungi/coccidioides_pangenome/results/p3_pair_class_sweep/` (not committed).
Each Afumigatus run took 385-651 s and 1.6 GB.

| dataset | run | FDR pairs | linkage_fraction exactly 0 |
|---|---|---|---|
| Cocci_Ci | `coccidioides_pangenome/.nf_launch/immitis_in_posadasii_out_v2` | 77,059 | 68,983 (89.5%) |
| Afumigatus | `Afumigatus_pangenome/.nf_launch/full_v070` | 1,363,573 | 1,347,005 (98.8%) |

Pairs near the thresholds at k = 10: within 0.05 of 0.5, 345 (Cocci) and 715 (Afu); within
0.02 of 0.05, 456 and 1,003.

## Results: k (thresholds 0.5 / 0.05)

| k | Afu physical | Afu ambiguous | Afu trans | Afu same label as k=10 | Cocci physical | Cocci ambiguous | Cocci same label |
|---|---|---|---|---|---|---|---|
| 3 | 5,401 | 1,548 | 571,252 | 0.9934 | 2,810 | 788 | 0.9491 |
| 5 | 7,249 | 2,435 | 568,869 | 0.9954 | 3,657 | 1,162 | 0.9644 |
| 7 | 8,728 | 3,384 | 566,792 | 0.9972 | 4,327 | 1,487 | 0.9784 |
| **10** | 10,541 | 4,624 | 564,245 | 1 | 5,059 | 1,986 | 1 |
| 15 | 12,721 | 6,446 | 561,004 | 0.9965 | 5,940 | 2,740 | 0.9728 |
| 20 | 14,321 | 7,886 | 558,771 | 0.9942 | 6,469 | 3,312 | 0.9577 |
| 30 | 16,727 | 10,175 | 555,557 | 0.9907 | 6,980 | 4,201 | 0.9393 |
| 50 | 20,162 | 13,846 | 551,208 | 0.9855 | 7,236 | 5,278 | 0.9219 |

"physical" = `unexplained_physical` (no captain-gene HMM was given, so `starship_explained` is 0).
Cocci has 0 `trans` at every setting (all low-linkage pairs are `trans_unconfirmed`, the
single-species gate); Cocci `trans_unconfirmed` goes 71,066 (k=5) to 66,104 (k=20).

- **The low-linkage side is stable.** Over k = 5-20, Afumigatus `trans` changes by -1.8%
  (568,869 to 558,771), and 99.4-99.7% of all pairs keep their k = 10 label. For Cocci,
  95.8-97.8% keep it.
- **The physical side is not stable.** `unexplained_physical` roughly doubles from k = 5 to
  k = 20 (Afu 7,249 to 14,321, x1.98; Cocci 3,657 to 6,469, x1.77) and keeps rising to 50.
  There is no window size where it levels off. Moving pairs go low -> ambiguous -> physical.

## Results: thresholds (k = 10)

| setting | Afu pairs that change label | Cocci pairs that change label |
|---|---|---|
| physical 0.3 | 0.11% | 1.03% |
| physical 0.7 | 0.11% | 0.88% |
| physical 0.9 | 0.26% | 2.29% |
| trans 0.01 | 0.09% | 0.84% |
| trans 0.1 | 0.07% | 0.48% |
| trans 0.2 | 0.18% | 1.10% |

The thresholds matter little because almost all pairs sit at linkage fraction 0 or near 1.
They do set the size of the small physical class (Afu 12,014 at 0.3, 7,032 at 0.9).

## Recommendation and evidence class

| parameter | recommendation | evidence class |
|---|---|---|
| k = 10 | keep for the trans call; report physical-pair counts as k-dependent | empirical: trans stable within 2% over k 5-20; physical count doubles over the same range |
| physical 0.5 / trans 0.05 | keep | empirical: moving either changes < 0.3% (Afu) / < 2.3% (Cocci) of labels |

The issue's test ("stable across k = 5-20 is defensible") is passed by the trans /
low-linkage classification and failed by the physical-linkage count. A physical count
reported in a paper should give k, or be shown as a curve over k.

## Limits

- Two runs. No captain-gene HMM, so the starship split is untested.
- The trans class in Cocci is empty for a different reason (the clade gate), so the Cocci
  test covers only the physical / ambiguous / low split.
