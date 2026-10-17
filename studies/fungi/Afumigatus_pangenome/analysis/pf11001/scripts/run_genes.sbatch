#!/bin/bash -l
#SBATCH -p epyc --constraint=milan -c 2 --mem 14G -t 1:00:00
#SBATCH -J pf11001_genes
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/analysis/pf11001/logs/genes.%j.out
# Stage 2 gene-level tables (needs ~10-30 GB RAM: loads gene_positions of three runs).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
D=$NII/studies/fungi/Afumigatus_pangenome/analysis/pf11001
"$NII/.pixi/envs/default/bin/python3" "$D/scripts/pf11001_genes.py"
echo done
