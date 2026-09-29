# #132 priority 4: Leiden resolution (default 1.0) for trans-pair modules

2026-09-28. Scripts: `p4_leiden_sweep.py` (pipeline functions from
`bin/pangenome_detect_trans_modules.py` at a394ace: RBConfigurationVertexPartition,
edge weight -log10 fdr_q clipped at 300), `p4_ami.py` (scikit-learn 1.9.1 AMI from the
saved partitions), `run_p4.sh` (SLURM job 29196421, 4 min 28 s for both graphs).
Outputs: `p4_leiden.tsv` (module counts, largest-module share, seed ARI/NMI),
`p4_ami.tsv`. Partitions (`*.labels.npz`) are in
`studies/fungi/coccidioides_pangenome/results/p3_pair_class_sweep/` (not committed).

## Graphs

| dataset | pair table | families | trans edges |
|---|---|---|---|
| Afumigatus | `Afumigatus_pangenome/.nf_launch/full_v070` PAIR_CLASSIFICATION task (v0.7.0) | 8,299 | 564,245 |
| Cocci_genus_0.9_0.8 | `coccidioides_pangenome/results/tier1_sweep/coccidioides_sweep_id0.9_cov0.8` (530 proteomes incl. U. reesii, af6fd68) | 14,602 | 543,867 |

The single-species Cocci runs (`immitis_in_posadasii_out_v2`) have **0 trans pairs**: all
68,840 low-linkage pairs are `trans_unconfirmed` (fail the permutation or >= 2 clades gate),
so no modules exist there. That is the documented behaviour for a one-species ingroup.

## Results

Resolutions 0.25-16, 10 seeds each. `seed AMI` = mean pairwise AMI between the 10 seeds.

| resolution | Afu modules | Afu largest share | Afu seed AMI | Cocci modules | Cocci largest share | Cocci seed AMI |
|---|---|---|---|---|---|---|
| 0.25 | 6.6 | 0.673 | 0.869 | 2.0 | 0.568 | 1.000 |
| 0.5 | 7.9 | 0.357 | 0.878 | 4.0 | 0.510 | 0.997 |
| **1.0** | 11.0 | 0.341 | **0.936** | 6.6 | 0.483 | 0.919 |
| 1.5 | 15.6 | 0.290 | 0.828 | 10.9 | 0.246 | 0.845 |
| 2 | 136.1 | 0.217 | 0.882 | 15.7 | 0.196 | 0.834 |
| 3 | 335.7 | 0.180 | 0.862 | 26.1 | 0.143 | 0.795 |
| 4 | 495.2 | 0.179 | 0.896 | 45.5 | 0.114 | 0.764 |
| 6 | 662.0 | 0.128 | 0.913 | 92.8 | 0.076 | 0.755 |
| 8 | 752.4 | 0.122 | 0.926 | 158.5 | 0.054 | 0.743 |
| 12 | 863.6 | 0.107 | 0.866 | 324.1 | 0.036 | 0.703 |
| 16 | 942.5 | 0.094 | 0.832 | 510.6 | 0.027 | 0.690 |

(Afumigatus module counts include many small modules above resolution 2: modules with
>= 10 families are 26.5 at 2, 51.8 at 6, 53.8 at 8.)

- **Cocci: stability falls monotonically with resolution.** The most stable partitions
  are the coarsest (2-4 modules, one holding more than half the families). "Pick the most
  stable plateau" therefore picks a partition with no detail. At 1.0 the largest module
  holds 48% of families.
- **Afumigatus: two stability peaks**, at 1.0 (11 modules, largest 34%, AMI 0.936) and at
  6-8 (AMI 0.913-0.926, largest 12-13%, about 52-54 modules of >= 10 families).
- AMI between neighbouring resolutions is 0.62-0.87: partitions change substantially at
  every step. There is no resolution range where the partition stays the same.
- Not comparable with the 2026-09-19 standalone note ("2.0 gives 66 at AMI 0.79"): that
  used a different pair set and Leiden call. Here 2.0 gives 15.7 modules for Cocci.

## Recommendation and evidence class

- **Do not treat 1.0 as a validated default.** Evidence class: empirical, negative. For
  Afumigatus it is a stable point; for pooled Cocci it gives one module with 48% of
  families.
- A stability-only selection rule fails on Cocci (it selects the trivial partition).
  A rule needs a second criterion, for example the most stable resolution among those
  where the largest module holds < 25% of families. That gives 8 for Afumigatus
  (AMI 0.926) and 2 for Cocci (AMI 0.834). The 25% threshold is itself unvalidated,
  so this is a proposal, not a result.
- The sweep is cheap (about 2 min per graph for 110 Leiden runs), so the pipeline can run it
  and report the stability curve with every study instead of one fixed resolution.
  Consensus clustering over seeds is the standard alternative to picking one seed.

## Limits

- Two graphs. The Cocci graph includes U. reesii and both species; module structure there
  may largely follow the species split.
- AMI across seeds measures reproducibility, not biological correctness.
