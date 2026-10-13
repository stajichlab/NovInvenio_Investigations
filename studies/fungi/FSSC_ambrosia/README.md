# Ambrosia Fusarium Clade (AFC) — pangenome (ingroup) vs nearest non-AFC FSSC

One `pangenome.nf` run. IN = 28 AFC strains (beetle symbionts). OUT = 4
strains of the nearest named non-AFC FSSC species. Rescue is off in the
first run.

Data, inventory and mash matrix are reused read-only from
`../fusarium_FOXY_vs_FSSC/selection/` (BFD funannotate output; mash k=21,
s=10000). Nothing there was changed.

## Ingroup rule (`bin/select_fssc_subclade.py --study ambrosia`)

1. Pool: 194 BFD strains whose nearest complex anchor is *F. vanettenii*
   77-13-4, at mash d <= 0.20.
2. Candidates (72 strains): the AFC names *euwallaceae, ambrosium,
   floridanum, obliquiseptatum, tuaranense, kuroshium, oligoseptatum,
   pseudensiforme, protoensiforme, rekanum, awan, akasia, mekan, variasi,
   warna*; the conditional names *papillatum, duplospermum*; the unnamed
   "Neocosmospora sp. B6999-1..5"; and the 14 "Fusarium cf. solani" strains
   that are not EtdFoc/EthFoc (UCR*, TW*, NRRL62590/62591/22643/66088).
3. Placement test, the same for every candidate label (one species name; each
   cf. solani strain is its own label):
   keep if min d to AFC strains of *other* labels < min d to any strain of a
   named non-AFC FSSC species.
   - "Named non-AFC FSSC" excludes "sp."/"cf." labels, one mislabelled
     "F. oxysporum" (14003, nearest *F. pseudensiforme* 0.036), and
     *F. drepaniforme* (not on the AFC list, but mash nests it in the AFC:
     0.024 to *F. tuaranense*, 0.089 to any other named species). So
     *F. drepaniforme* is neither ingroup, comparator nor outgroup.
   - Kept: *papillatum* (0.027 AFC vs 0.089), *duplospermum* (0.054 vs 0.086),
     NRRL22643 (0.055 vs 0.088) and NRRL66088 (0.055 vs 0.088). NRRL22643/66088
     then fail the protein floor (below).
   - Dropped: *F. protoensiforme* NRRL 22178 (0.088 to AFC, *F. awan*; 0.069
     to *F. piperis*) and *F. awan* (0.088 to AFC; 0.081 to
     *F. neocosmosporiellum*). Both are on the brief's AFC list; mash places
     them closer to non-AFC species. *F. awan* would also fail the protein floor.
4. QC:
   - BUSCO genome C >= 90%. Dropped 1: *F. papillatum* NRRL 62944 (84.0; the
     other *papillatum* strain NRRL62944 passes).
   - Proteins >= max(11,000, median - 3 x scaled MAD). Median 14,579, scaled
     MAD 2,380, so median - 3 MAD = 7,440 and the 11,000 floor applies.
     Dropped 23 (7,233-9,656 proteins, except one at 416):
     all 14 cf. solani (UCR4674/4676/4677/4678/4679/5499/5508/6395, TW44, TW56,
     NRRL62590, NRRL62591, NRRL22643, NRRL66088), *mekan* x2, *obliquiseptatum*
     x3, *pseudensiforme* x2, *warna* x1, *tuaranense* NRRL 46518
     (GCA_013363205.1, 416 proteins). Their genomes are 45-54 Mb with genome
     BUSCO 99.2-99.8%, so the low counts look like annotation shortfall, not
     gene loss. I did not test this. This QC removes 5 named species
     (*mekan, obliquiseptatum, pseudensiforme, warna*, and all cf. solani).
   - Same strain deposited twice: *F. ambrosium* NRRL 20438 (GCA_003947045.1,
     N50 122 kb) dropped for NRRL20438 (GCA_057676125.1, N50 196 kb); d 0.0010.
   - Clones (single linkage, mash d < 0.0005, keep highest N50): 8 groups,
     19 strains dropped. B6999-1 7->1 (CBS153107), B6999-4 7->1 (CBS153114),
     B6999-2 2->1, euwallaceae CBS135854/135859 ->1, floridanum 3->1
     (NRRL62628), rekanum 2->1, tuaranense 2->1 (NRRL46518), variasi 2->1.
     Table: `selection/clone_groups.tsv`. Pairs just above the threshold stay:
     euwallaceae MAFF243812/243813 (0.0005), ambrosium CBS571.94/MAFF246289
     (0.0006).
5. Ingroup: 28 strains, 17 labels: *euwallaceae* 5, *ambrosium* 4,
   *kuroshium* 2, *oligoseptatum* 2, *floridanum* 2, *akasia* 2,
   *duplospermum* 2, *tuaranense, rekanum, variasi, papillatum* 1 each,
   B6999-1..5 1 each. Proteins 11,983-17,233.

Mash-AFC strains not used (not on the name list): *F. drepaniforme* x2,
"Fusarium sp." AF-6 (NRRL62590), NRRL 62610, NRRL 62957, NRRL 66088
(step `afc_by_mash_not_by_name` in `selection/selection_log.tsv`).

## Outgroup rule

Named non-AFC FSSC species (the comparator set above). *F. solani* identity
is by mash: members of the *F. solani* s. str. cluster (`../FSSC_solani`)
count as *F. solani*; "F. solani"-labelled strains outside it are not used.
Strains must have BUSCO >= 90%, >= 11,000 proteins, N50 >= 100 kb. Species
ranked by min mash d to the ingroup; the 4 nearest each give their
highest-N50 strain.

| species | strain | N50 | proteins | species min d | chosen strain d (min / median to IN) |
|---|---|---|---|---|---|
| *F. vanettenii* | Fs6.1 | 3.36 Mb | 19,538 | 0.0829 | 0.0932 / 0.0958 |
| *F. keratoplasticum* | Fu6.1 | 4.21 Mb | 14,285 | 0.0832 | 0.0832 / 0.0905 |
| *F. solani* s. str. | F421 | 4.03 Mb | 17,704 | 0.0839 | 0.0857 / 0.0903 |
| *F. phalaenopsidis* | FUZ10S | 3.66 Mb | 14,463 | 0.0857 | 0.0857 / 0.0912 |

All named non-AFC species are 0.083-0.093 from the AFC. The ranking spans
only 0.003 and is near mash saturation, so "nearest" is weak here.
*F. falciforme* (0.0864) is rank 5.

## Mash structure of the ingroup (`bin/mash_structure.py`)

Pairwise mash 0.0005-0.0707 (median 0.0293). UPGMA clusters:
- cut 0.04: 4 clusters. One of 24 strains; *duplospermum* (2);
  *variasi* (1); B6999-3 (1).
- cut 0.02: 11 clusters. The largest (12) = *euwallaceae, kuroshium,
  floridanum*, B6999-1/4/5. Then *ambrosium* (4), *oligoseptatum*,
  *akasia*, *duplospermum* (2 each), and singletons *tuaranense, rekanum,
  papillatum, variasi*, B6999-2, B6999-3.
- *variasi* and *duplospermum* are 0.054-0.071 from the rest. They are the
  deep lineages of this ingroup.
Per-strain table: `selection/ingroup_mash_clades.tsv`.

## Build chain, running

Same as `../FSSC_solani/README.md`, with `--study ambrosia`.

```
sbatch .nf_launch/mmseqs_norescue/submit_nextflow_head.sh
```

Pipeline: detached worktree
`/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/fssc-runs-a72b4bc`
(NovInvenio a72b4bc082a5915bb1e80dc6c7f0d935ad0fc1e8). Rescue off via
`-params-file params_norescue.yaml`. `stajichlab_queue.config` adds
FAMILY_PFAM_SCAN `time = 8.h`. Outputs: `results/mmseqs_norescue/`.

## Status (2026-09-25)

First run finished (job 29101565, 2 h 07 min, 10.3 used / ~109 allocated
cpu-h). 82,758 families; 5,535 core. Co-occurrence: 56.6 M candidate pairs,
0 FDR pairs, so no islands, trans pairs or modules. At 24 representative
strains no pair can reach the BH threshold (min Fisher p 3.7e-7 vs 8.8e-10
needed). See `notes/2026-09-25-first-run.md`.
