# Fusarium solani sensu stricto — pangenome (ingroup) vs nearest FSSC species

One `pangenome.nf` run. IN = *F. solani* s. str. (13 strains). OUT = 4 strains
of the nearest named FSSC species. Rescue is off in the first run.

Data, inventory and mash matrix are reused read-only from
`../fusarium_FOXY_vs_FSSC/selection/` (BFD funannotate output for every
strain; mash k=21, s=10000; see that README). Nothing there was changed.

## Ingroup rule (`bin/select_fssc_subclade.py --study solani`)

1. Pool: 194 BFD strains whose nearest complex anchor is *F. vanettenii*
   77-13-4, at mash d <= 0.20.
2. Candidates: the single-linkage cluster at mash d < 0.03 that contains
   *F. solani* CBS 140079 (GCA_045529305.1). Membership is by mash, not label.
   Recomputed 2026-09-25: 23 strains, max pairwise d 0.0176.
   - 11 "Fusarium solani", 10 "Fusarium cf. solani" EtdFoc/EthFoc
     (Ethiopia, PRJNA353637), and 2 more that the brief did not list:
     "Fusarium haematococcum" LQ1 (d 0.0106 to CBS 140079) and
     "Fusarium sp. R1" (d 0.0109). Both sit inside the cluster
     (nearest neighbours Karbala-1 / NAPSME6.1 at 0.006-0.007), so they are
     kept as *F. solani* s. str. by the mash rule. Their labels are likely wrong.
3. QC:
   - BUSCO genome (fungi_odb12) C >= 90%. Dropped 5: EtdFoc-178 (84.0),
     EtdFoc-209 (87.8), EtdFoc-29 (73.0), EtdFoc-35 (58.9), EtdFoc-63 (89.4).
   - Proteins >= max(11,000, median - 3 x scaled MAD). Median 16,580, scaled
     MAD 1,597 (1.4826 x raw MAD), so the floor used was 11,790.
     Dropped 5: EtdFoc-155 (9,952), EtdFoc-208 (10,177), EtdFoc-214 (9,955),
     EthFoc-38 (10,205), EthFoc-8 (9,965). All 10 Ethiopian strains are gone.
   - Same strain deposited twice: none. Clones (mash d < 0.0005): none.
4. Ingroup: 13 strains (11 *F. solani*, LQ1, R1). Proteins 14,493-18,084.

## Other "F. solani"-labelled lineages (not in the ingroup)

14 BFD strains labelled *F. solani* are 0.046-0.076 from the cluster. With
single linkage at d < 0.03 over the whole pool they fall with:

| lineage (d<0.03 group) | "F. solani" strains | min d to s. str. cluster |
|---|---|---|
| with *F. falciforme* / *F. paranaense* | CC_50, CR12, CSH_1, CSH_2, CSH_7, S2_22 | 0.049-0.051 |
| with *F. keratoplasticum* Fu6.1 | BAT 10, NK-NH1, SC02 | 0.047-0.048 |
| with *F. metavorans* FSSC_6 | 2021000008 | 0.072 |
| with *F. vanettenii* | GU-3 | 0.068 |
| with *F. liriodendri* | KMZW-1 | 0.075 |
| singletons | JS-169, IISc-1 | 0.074, 0.076 |

The 14 non-Ethiopian "Fusarium cf. solani" (UCR*, TW*, NRRL*) are
0.088-0.093 away; they are in `../FSSC_ambrosia` candidates.
Full list: `selection/selection_log.tsv` (step `not_in_solani_cluster`).

## Outgroup rule

From named FSSC species other than *F. solani* (so not the mislabelled
lineages above). Strains must have BUSCO >= 90%, >= 11,000 proteins and
N50 >= 100 kb. Species are ranked by min mash d to the 13-strain ingroup.
The 4 nearest species each give their highest-N50 strain.

| species | strain | N50 | proteins | species min d | chosen strain d (min / median to IN) |
|---|---|---|---|---|---|
| *F. phalaenopsidis* | FUZ10S | 3.66 Mb | 14,463 | 0.0423 | 0.0423 / 0.0431 |
| *F. vanettenii* | Fs6.1 | 3.36 Mb | 19,538 | 0.0464 | 0.0683 / 0.0706 |
| *F. falciforme* | Progeny4 | 4.10 Mb | 16,187 | 0.0483 | 0.0503 / 0.0512 |
| *F. keratoplasticum* | Fu6.1 | 4.21 Mb | 14,285 | 0.0488 | 0.0488 / 0.0503 |

Note: the highest-N50 *F. vanettenii* strain (Fs6.1) is not the nearest one
(0.068 vs 0.046). *F. phalaenopsidis* is on the FSSC list on mash placement
only (see the FOXY_vs_FSSC README). *F. solani-melongenae* is 0.0588 (rank 6).

## Mash structure of the ingroup (`bin/mash_structure.py`)

Pairwise mash 0.0060-0.0119 (median 0.0104). UPGMA: one cluster at 0.02;
5 clusters at 0.01 (sizes 6/3/2/1/1); no pair below 0.005. So the ingroup is
one species with shallow, fairly even structure. Per-strain table:
`selection/ingroup_mash_clades.tsv`. Matrix of the 17 selected strains:
`selection/selected_mash_matrix.tsv`.

## Build chain

```
/usr/bin/python3.12 bin/select_fssc_subclade.py --study solani --study-dir $PWD   # ~5 min (grep counts on GPFS)
<pixi python> bin/mash_structure.py --study-dir $PWD                              # scipy
/usr/bin/python3.12 bin/build_species_csv.py --study-dir $PWD
<pixi python> bin/build_study_config_symlinks.py --study-dir studies/fungi/FSSC_solani --skip-fetch   # from NII root
```

`<pixi python>` = `/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/fssc-runs-a72b4bc/.pixi/envs/default/bin/python`.
species.csv Group is IN/OUT directly, so `config.csv` is the samplesheet.
`data_dir/{pep,dna,gff3}` are symlinks into BFD. TaxonGroup is empty;
pangenome.nf fills ingroup clade labels itself.

## Running

```
sbatch .nf_launch/mmseqs_norescue/submit_nextflow_head.sh
```

Calls `run_pangenome.sh`. Pipeline: detached worktree
`/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/fssc-runs-a72b4bc`
(NovInvenio a72b4bc082a5915bb1e80dc6c7f0d935ad0fc1e8, origin/main on
2026-09-25), own pixi env. Do not remove it while a resume may be needed.
Rescue off via `-params-file params_norescue.yaml` (not the CLI boolean).
`stajichlab_queue.config` = the FOXY_vs_FSSC overrides plus
FAMILY_PFAM_SCAN `time = 8.h` on the first attempt.
Outputs: `results/mmseqs_norescue/`. Run notes: `notes/`.

## Status (2026-09-25)

First run finished (job 29101564, 57 min, 3.9 used / ~27 allocated cpu-h).
50,251 families; 8,626 core. Co-occurrence found 0 FDR pairs, so there are
no islands, trans pairs or modules. At 13 strains no pair can reach the BH
threshold (min Fisher p 5.8e-4 vs 3.0e-9 needed). See
`notes/2026-09-25-first-run.md`.

## Scope decision (2026-09-25)

Use this study for family and frequency-bin patterns only (family counts, core /
shell / singleton bins, pangenome openness). Do not interpret co-occurrence,
islands, trans pairs or modules from it.

Reason: with 13 ingroup strains, the smallest possible one-sided Fisher p is
1/C(13,6) = 5.8e-4. BH runs over all 16.5 M pairs, so a block of about 190,000
pairs tied at that p would be needed before any pair passes. The first run found
0 FDR-significant pairs. By the same arithmetic FOL (21 strains) needs about
1,300, and found 820,493.
