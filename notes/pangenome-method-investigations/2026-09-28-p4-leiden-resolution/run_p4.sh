#!/usr/bin/env bash
#SBATCH -p epyc
#SBATCH -c 2
#SBATCH --mem 16G
#SBATCH -t 08:00:00
#SBATCH -J p4-leiden
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/notes/pangenome-method-investigations/2026-09-28-p4-leiden-resolution/p4_%j.log
set -euo pipefail
C=/rhome/jstajich/.nextflow/assets/.repos/stajichlab/NovInvenio/clones/a394ace8d288b19591ea148e5a12d55aeac668b1
D=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/notes/pangenome-method-investigations/2026-09-28-p4-leiden-resolution
S=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi
$C/.pixi/envs/default/bin/python $D/p4_leiden_sweep.py --code $C \
  --pairs Afumigatus=$S/Afumigatus_pangenome/.nf_launch/full_v070/work/ff/9f2719203cd38a0920c83231edddb5/pair_classification.tsv.zst \
  --pairs Cocci_genus_0.9_0.8=$S/coccidioides_pangenome/results/tier1_sweep/coccidioides_sweep_id0.9_cov0.8/pangenome/pair_classification.tsv.zst \
  --out $D/p4_leiden.tsv --labels_dir /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/results/p3_pair_class_sweep
