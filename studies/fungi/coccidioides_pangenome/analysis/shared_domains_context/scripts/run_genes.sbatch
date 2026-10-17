#!/bin/bash -l
#SBATCH -p short --qos short -c 2 --mem 48G -t 1:00:00
#SBATCH -J sdc_genes
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/shared_domains_context/logs/genes.%j.out
# Stage 1 gene-level tables (loads gene_positions of three runs).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
D=$NII/studies/fungi/coccidioides_pangenome/analysis/shared_domains_context
/usr/bin/time -v "$NII/.pixi/envs/default/bin/python3" "$D/scripts/sdc_genes.py" 2>&1 | grep -v -E '^\s+(Average|Exit|Page|Swaps|File|Socket|Signals|Voluntary|Involuntary|Major|Minor|Percent|System|User|Command)'
echo done
