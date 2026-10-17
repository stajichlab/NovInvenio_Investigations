#!/bin/bash -l
#SBATCH -p epyc --constraint=milan -c 2 --mem 12G -t 1:00:00
#SBATCH -J af_sweep_metric
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/analysis/starship_island_overlap/logs/sweep_metric.%j.out
# Reproduce the memory-note numbers (21.3% of placed island genes in Af293 Starships vs
# 5.2% of non-core Af293 genes; full_v070) with the original scorer, and apply the
# same scorer to full_v072. Original scorer: analysis/island_genome_count/scripts/score_starships.py.
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
AF=$NII/studies/fungi/Afumigatus_pangenome
D=$AF/analysis/starship_island_overlap
mkdir -p "$D/reproduce_sweep_metric"
for RUN in full_v070 full_v072; do
  P=$AF/results/$RUN/config/pangenome
  "$NII/.pixi/envs/default/bin/python3" "$AF/analysis/island_genome_count/scripts/score_starships.py" \
    --islands "$P/significant_islands.tsv" --frequency_table "$P/frequency_table.tsv" \
    --presence_matrix "$P/presence_matrix.rescued.tsv" --cluster_tsv "$P/cluster/tier1_cluster.tsv" \
    --gene_positions "$P/gene_positions.tsv.zst" \
    --starships "$AF/analysis/island_genome_count/inputs/af293_starships_table_s7.tsv" \
    --presence_truth "$AF/results/id_crosswalk/ground_truth_starships_by_short.tsv" \
    --subset_strains "$D/inputs/strains_$RUN.txt" \
    --out_prefix "$D/reproduce_sweep_metric/$RUN"
done
