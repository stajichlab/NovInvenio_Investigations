#!/bin/bash
# Gap check (2026-10-08): why N=121 gave 8,226 islands vs 11,775 in full_v070. Reruns PAIR_CLASSIFICATION and
# BUILD_ISLANDS on full_v070's co-occurrence output with all-strain vs sweep-subset family_positions. See ISLAND_GENOME_COUNT.md.
#SBATCH --job-name=gapcheck
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=3:00:00
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/gapcheck_%j.out
set -euo pipefail
STUDY=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome
PG=$STUDY/results/full_v070/config/pangenome
AN=$STUDY/analysis/island_genome_count
B=/rhome/jstajich/.nextflow/assets/.repos/stajichlab/NovInvenio/clones/f87fd1e9073b0d0d8927baca8a93f65ae3d86913/bin
OUT=$AN/outputs/gap_check
mkdir -p $OUT
export PATH="$HOME/.pixi/bin:$PATH"
eval "$(pixi shell-hook -s bash --frozen --manifest-path /bigdata/stajichlab/jstajich/projects/NovInvenio/pixi.toml)"
W=${SCRATCH:?}/gapcheck; mkdir -p $W; cd $W
: > empty_evalues.tsv
cp $PG/cluster/tier1_cluster.tsv .
# subset family_positions exactly as make_subset.py did: rows of the 123 kept strains
zstd -dc $PG/family_positions.tsv.zst | awk -F'\t' 'NR==FNR{k[$1]=1;next} FNR==1||($1 in k)' $AN/outputs/runs/n121_s0/subset_strains.txt - | zstd -q -o fp_subset.tsv.zst
for V in full subset; do
  fp=$PG/family_positions.tsv.zst; [ $V = subset ] && fp=$W/fp_subset.tsv.zst
  python $B/pangenome_pair_classification.py --cooccurring_pairs $PG/cooccurring_pairs.tsv.zst \
    --family_positions $fp --cluster_tsv tier1_cluster.tsv --captain_tblout empty_evalues.tsv --k 10 \
    --physical_threshold 0.5 --trans_threshold 0.05 --min_co_carrying 5 --perm_alpha 0.05 --min_clades 2 \
    --id_sep '|' --output pc_$V.tsv.zst
  python $B/pangenome_build_islands.py --family_positions $fp --frequency_table $PG/frequency_table.tsv \
    --pair_classification pc_$V.tsv.zst --cluster_tsv tier1_cluster.tsv --min_island_size 2 --id_sep '|' \
    --output islands_$V.tsv
  echo -e "$V\t$(( $(wc -l < islands_$V.tsv) - 1 ))" >> $OUT/gap_check_islands.tsv
  zstd -dc pc_$V.tsv.zst | awk -F'\t' 'NR>1{c[$3]++}END{for(k in c)print "'$V'\t"k"\t"c[k]}' >> $OUT/gap_check_classes.tsv
done
echo done
