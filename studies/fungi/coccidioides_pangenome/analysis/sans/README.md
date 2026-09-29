# SANS test: splits and trees for the 529 Coccidioides pangenome strains

Started 2026-09-28. Exploratory. Nothing here is wired into a pipeline run.

## Question

Can SANS (https://github.com/gi-bielefeld/sans, v2.8_2, commit 3d3cdd9) show the
pangenome's gene presence/absence patterns on or against a phylogeny? How do its
trees compare to the SNP trees?

## Tool

- SANS has no gene-content mode. It builds splits from shared k-mers
  (genomes, or proteins with `-a`).
- `SANS -s` loads a splits file and runs only the filters (`-f strict` gives a
  tree, `-f weakly` a network). This gives a gene-content mode:
  `gene_content_splits.py` writes one split per family carrier set, weight =
  number of families with that split.
- Build: `/bigdata/stajichlab/jstajich/software/sans`, makefile changed to
  `-march=x86-64-v2 -DmaxN=600` (default maxN=100 is too small) and
  `-static-libstdc++ -static-libgcc` (the system libstdc++ is too old for gcc 12).
- Name lists must hold the bare strain names that the splits file uses
  (`out/strains.txt`). For k-mer runs, the kmtricks form `ID : path` sets the names
  (`genomes.kmt`, `proteins.kmt`).

## Inputs

- Presence matrix: `../../results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome/presence_matrix.rescued.tsv`
  (47,745 families x 529 strains). Present = `present` or `genome_only`.
- Genomes and proteins: `../../data_dir/dna`, `../../data_dir/pep` (symlinked in `genomes/`, `proteins/`).
- Reference trees:
  - SNP, IQ-TREE consensus (`GTR+F+ASC`, 446 tips, 323,660 columns) and FastTree:
    `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/strain_tree/unpruned/Cocciall_v1.All.SNP.mfa.contree`, `Cocciall_v1.All.SNP.fasttree.tre`.
    Tip rename: `Coahuilla_2=Coahuila_2`. `Ref_CimmRS` is not renamed: the SNP trees also hold `RS`, and the Cp tree uses `Ref_CimmRS` as a Ci outgroup. 412 SNP tips are pangenome strains. (Before 2026-09-28 20:10 a `Ref_CimmRS=RS` rename made a duplicate RS tip; all numbers below are from the corrected rerun.)
  - BUSCO CDS species tree with 51 mash grafts: `../species_tree/coccidioides_coccidioides_only.full.nwk`.

## Files

| File | What |
|---|---|
| `gene_content_splits.py` | presence matrix -> SANS splits file |
| `run_sans.sh` | SLURM job: SANS k-mer splits, `dna` (k=31) or `aa` (`-a`, k=10); `-t 10n`, then strict tree + weakly network |
| `compare_trees.py` | split precision/recall/RF between trees on shared tips; family fit (clade/conflict); `--parsimony` changes per family; `--shuffle` label-permuted control |
| `out/` | splits, trees, networks, comparison text |

## Results

SANS k-mer runs (`run_sans.sh`, jobs 29181645/29181646, 2026-09-28): DNA k=31 4 min, 15.6 GB peak;
protein k=10 38 s, 1.1 GB. Strict trees: `out/dna.strict.nwk`, `out/aa.strict.nwk`; networks `out/*.weakly.nex`.
Full comparison tables: `out/compare_all.txt`, `out/compare_Ci.txt` (Ci SNP IQ-TREE), `out/compare_Cp.txt` (Cp SNP FastTree; no Cp IQ-TREE contree exists).

### Gene-content splits

`out/gene_content.splits`: 32,006 non-trivial families give 28,794 distinct splits.
Heaviest: a clonal `UTAH_20380X*` group (49 families). Ci | Cp: 42 families.
Strict filter keeps 284 splits (a binary tree of 529 tips has 526).

Split agreement (`out/compare_all.txt`, trees restricted to shared tips):

| query | reference | shared tips | shared splits / query splits | shared / reference splits |
|---|---|---|---|---|
| SNP FastTree | SNP IQ-TREE | 446 | 0.553 | 0.553 (0.691 of UFboot >= 95) |
| BUSCO CDS | SNP IQ-TREE | 412 | 0.186 | 0.186 (0.264 of UFboot >= 95) |
| gene content | SNP IQ-TREE | 412 | 0.285 | 0.134 (0.187 of UFboot >= 95) |

Family fit on the SNP IQ-TREE (412 matrix strains): of 28,830 non-trivial families,
716 (2.5%) are exactly one clade; the rest conflict with at least one split.

Parsimony (412 strains common to all trees, 28,830 families with minor side >= 2):

| tree | total changes | mean per family | minor 2 | 3-5 | 6-20 | 21-100 | >100 |
|---|---|---|---|---|---|---|---|
| SNP IQ-TREE | 549,451 | 19.06 | 1.91 | 3.50 | 9.70 | 36.84 | 48.35 |
| SNP FastTree | 550,751 | 19.10 | | | | | 48.44 |
| BUSCO CDS | 563,902 | 19.56 | 1.92 | 3.53 | 9.85 | 37.94 | 49.62 |
| gene content (SANS strict) | 805,735 | 27.95 | 1.87 | 3.47 | 9.98 | 41.53 | 101.15 |
| SNP IQ-TREE, labels shuffled | 915,999 | 31.77 | | | | | 116.81 |
| families per band | | | 5,224 | 5,889 | 6,701 | 6,803 | 4,213 |

Caveat: the gene-content tree has many polytomies. The rule used at a polytomy
(Hartigan, `k - max` per node) scores the star, not the best resolution, so it
overstates that tree's length. Only the three binary trees compare fairly.

### All six trees (`out/compare_all.txt`)

Split agreement with the SNP IQ-TREE (412 shared tips):

| tree | splits | in SNP tree | SNP UFboot >= 95 splits recovered |
|---|---|---|---|
| SANS DNA | 198 | 0.465 | 0.308 of 273 |
| SANS protein | 195 | 0.472 | 0.319 of 273 |
| SANS gene content | 193 | 0.285 | 0.187 of 273 |
| BUSCO CDS | 409 | 0.186 | 0.264 of 273 |
| shuffled control | 410 | 0.005 | 0.004 of 273 |

Parsimony, mean changes per informative family (pooled 412, Ci 159, Cp 253 strains):

| tree | pooled (412) | Ci only (159) | Cp only (254) |
|---|---|---|---|
| SNP tree (Ci: IQ-TREE; Cp: FastTree; pooled: IQ-TREE) | 19.06 | 10.88 | 17.65 |
| BUSCO CDS | 19.56 | 11.19 | 17.97 |
| SANS DNA | 20.53 | 11.79 | 18.97 |
| SANS protein | 20.48 | 11.71 | 18.96 |
| SANS gene content | 27.95 | 13.64 | 19.90 |
| shuffled SNP tree | 31.77 | 14.34 | 20.74 |

Within a species, the best tree is only 24% (Ci) and 15% (Cp) shorter than
the shuffled control. So within-species gene presence/absence has little
tree structure. These data do not separate recombination from presence-call
noise. The SANS trees are partly unresolved, so their lengths are upper bounds.

## Provenance

Derived analysis; no external data is committed. Inputs (sha256, 2026-09-28):

| input | sha256 |
|---|---|
| `presence_matrix.rescued.tsv` (run `rescue_freqpol_immitis_in_posadasii_out`) | `9823da8b2b15b43b6f24f0b0d3647db2125669d9771f9e8ee139b5a9973b0469` |
| `Cocciall_v1.All.SNP.mfa.contree` | `bafeea8d30520d8dcf5ebdd3721c49e109e8d05bbc6aa71de4184b12325870b1` |
| `Cocciall_v1.All.SNP.fasttree.tre` | `4cfbc53daca62f399ffa858444cf04b78f85dae49f3a969b41b14c22cd53609c` |
| `Cocciall_v1.Cimmitis.SNP.mfa.contree` | `886fb1c0b56a263f4f96650f5691f0bab090934e9ea27ed6d98a6d34ba106b73` |
| `Cocciall_v1.Cposadasii.SNP.fasttree.tre` | `90e8b8ea6672439cab1126d3b8f80ef333ee26493b13b427ccbf5266f2ade850` |
| `../species_tree/coccidioides_coccidioides_only.full.nwk` | `bd53141a8ec9494dce9f97b9ffbe072ac7c51f684b14029ad0b13e5f4ba9e9b5` |

Derived by: `gene_content_splits.py` (default `--present present,genome_only`),
`run_sans.sh dna|aa` (SANS 2.8_2, commit 3d3cdd9, GPL-3.0), SANS `-s ... -f strict|weakly`
on `out/gene_content.splits`, and `compare_trees.py` (commands in this README).
Not committed (regenerable, see `.gitignore`): `genomes/`, `proteins/`, `logs/`, `out/gene_content.splits`.

## Reproduce

```bash
S=/bigdata/stajichlab/jstajich/software/sans/SANS
M=../../results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome/presence_matrix.rescued.tsv
D=/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/strain_tree/unpruned
python3.12 gene_content_splits.py --matrix $M --out out/gene_content.splits --strains-out out/strains.txt
$S -i out/strains.txt -s out/gene_content.splits -f strict -N out/gene_content.strict.nwk -o out/gene_content.strict.splits
$S -i out/strains.txt -s out/gene_content.splits -f weakly -o out/gene_content.weakly.splits -X out/gene_content.weakly.nex
mkdir -p genomes proteins logs   # symlink <Short>.fa to ../../data_dir/{dna,pep}/<Short>.{dna,pep}.fa
sbatch --job-name sans-dna run_sans.sh dna; sbatch --job-name sans-aa run_sans.sh aa
python3.12 compare_trees.py --rename Coahuilla_2=Coahuila_2 \
  --families $M --parsimony --shuffle snp_iqtree \
  --families-ref snp_iqtree \
  --tree snp_iqtree=$D/Cocciall_v1.All.SNP.mfa.contree --tree snp_fasttree=$D/Cocciall_v1.All.SNP.fasttree.tre \
  --tree busco_cds=../species_tree/coccidioides_coccidioides_only.full.nwk \
  --tree sans_dna=out/dna.strict.nwk --tree sans_aa=out/aa.strict.nwk \
  --tree sans_gene_content=out/gene_content.strict.nwk > out/compare_all.txt
# Ci / Cp: the same with the Cimmitis contree / Cposadasii fasttree as the SNP tree -> out/compare_Ci.txt, out/compare_Cp.txt
```

## Prototype page (2026-09-28)

`prototype/`: tree-ordered heatmap (option 1) and split networks (option 2) in one page.
Published privately: https://claude.ai/artifact/Wz53NNHTojwx4gCbfdAy8y

- `build_network.py`: NeighborNet circular ordering (SplitsPy 0.x, `pip install --user splitspy`)
  on a distance matrix, then SANS splits that are contiguous on it, laid out with the outline
  algorithm. SplitsPy's least-squares weights were too slow for 529 taxa (>11 min for 100), so
  SANS weights are used instead. Gene content: Jaccard over the 32,006 informative families,
  104 of 566 SANS weakly splits circular (33.1% of weight). Genomes: mash distance
  (`results/preserved_from_work/genus_vs_ureesii/all_strains.mash_dist.tsv.gz`), 663 of 1,028
  splits (97.2%).
- `build_viz_data.py`: per tree, the best-matching clade per family (Jaccard on the minority side),
  Fitch changes, packed presence bits; Pfam names from the run's `pfam.domtblout` (island families only).
- `template.html`: the page; `__DATA__` is replaced by `viz_data.json`.

```bash
cd prototype
python3.12 build_network.py --splits ../out/gene_content.weakly.splits --dist jacc.npy --labels ../out/strains.txt --out net_gene_content.json
python3.12 build_network.py --splits ../out/dna.weakly.splits --dist ../../../results/preserved_from_work/genus_vs_ureesii/all_strains.mash_dist.tsv.gz --labels ../out/strains.txt --out net_dna.json
python3.12 build_viz_data.py --matrix <run>/presence_matrix.rescued.tsv --config ../../../config.csv --pfam <run>/pfam.domtblout \
  --tree snp=<SNP contree> --tree busco=../../species_tree/coccidioides_coccidioides_only.full.nwk \
  --rename Coahuilla_2=Coahuila_2 --network dna=net_dna.json --network gene_content=net_gene_content.json --out viz_data.json
```
`jacc.npy` is 1 - Jaccard between strains over families with minority side >= 2 (not yet a script; see the commit message).
