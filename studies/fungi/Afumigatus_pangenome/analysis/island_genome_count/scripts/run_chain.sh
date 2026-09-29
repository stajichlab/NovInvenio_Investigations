#!/bin/bash
#SBATCH --job-name=island_ncount
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=8:00:00
# Rerun the pangenome island steps (FREQUENCY_BINS -> COOCCURRENCE ->
# PAIR_CLASSIFICATION -> BUILD_ISLANDS) on genome subsets of full_v070, with
# the step commands and parameters copied from the pipeline's own
# .command.sh files (NovInvenio f87fd1e, run in10_v1), then score the islands
# against known Starships.
# Usage: sbatch run_chain.sh <seed> <N> [<N> ...]
# Absolute paths only (no BASH_SOURCE: it breaks under SLURM).
set -euo pipefail

SEED=$1; shift
NS=("$@")

STUDY=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome
AN=$STUDY/analysis/island_genome_count
PG=$STUDY/results/full_v070/config/pangenome
PIPE=/rhome/jstajich/.nextflow/assets/.repos/stajichlab/NovInvenio/clones/f87fd1e9073b0d0d8927baca8a93f65ae3d86913
B=$PIPE/bin

export PATH="$HOME/.pixi/bin:$PATH"
eval "$(pixi shell-hook -s bash --frozen --manifest-path /bigdata/stajichlab/jstajich/projects/NovInvenio/pixi.toml)"

WORK=${SCRATCH:?SCRATCH not set}/island_ncount_s${SEED}
mkdir -p "$WORK"
cp "$PG/cluster/tier1_cluster.tsv" "$WORK/"
: > "$WORK/empty_evalues.tsv"   # the pipeline passes an empty captain tblout too

step() {  # step <name> <cmd...>: run, append wall seconds to timings.tsv
    local name=$1; shift
    local t0; t0=$(date +%s)
    "$@"
    echo -e "${name}\t$(( $(date +%s) - t0 ))" >> timings.tsv
}

for N in "${NS[@]}"; do
    TAG=n${N}_s${SEED}
    D=$WORK/$TAG
    OUT=$AN/outputs/runs/$TAG
    mkdir -p "$D" "$OUT"
    cd "$D"
    : > timings.tsv
    ln -sf "$WORK/tier1_cluster.tsv" "$WORK/empty_evalues.tsv" .

    step make_subset python "$AN/scripts/make_subset.py" --pangenome_dir "$PG" \
        --n "$N" --seed "$SEED" --out_dir "$D"

    step frequency_bins python "$B/pangenome_frequency_bins.py" \
        --matrix presence_matrix.rescued.tsv --config samplesheet.with_clades.csv \
        --ingroup_label IN --outgroup_label OUT --outgroup_min_bin_strains 3 \
        --inventory strain_inventory.tsv --core_cutoff 0.95 --softcore_cutoff 0.90 \
        --shell_cutoff 0.15 --output frequency_table.tsv

    step cooccurrence python "$B/pangenome_cooccurrence.py" \
        --matrix presence_matrix.rescued.tsv --frequency_table frequency_table.tsv \
        --config samplesheet.with_clades.csv --ingroup_label IN --outgroup_label OUT \
        --inventory strain_inventory.tsv --min_strain_count 5 --fdr_alpha 0.05 \
        --screen_alpha 0.2 --polarity_loss_min_frac 0.9 --polarity_gain_max_frac 0.1 \
        --output cooccurring_pairs.tsv.zst > cooccurrence.log 2>&1

    step pair_classification python "$B/pangenome_pair_classification.py" \
        --cooccurring_pairs cooccurring_pairs.tsv.zst --family_positions family_positions.tsv.zst \
        --cluster_tsv tier1_cluster.tsv --captain_tblout empty_evalues.tsv --k 10 \
        --physical_threshold 0.5 --trans_threshold 0.05 --min_co_carrying 5 \
        --perm_alpha 0.05 --min_clades 2 --id_sep '|' --output pair_classification.tsv.zst

    step build_islands python "$B/pangenome_build_islands.py" \
        --family_positions family_positions.tsv.zst --frequency_table frequency_table.tsv \
        --pair_classification pair_classification.tsv.zst --cluster_tsv tier1_cluster.tsv \
        --min_island_size 2 --id_sep '|' --output significant_islands.tsv

    step score_starships python "$AN/scripts/score_starships.py" \
        --islands significant_islands.tsv --frequency_table frequency_table.tsv \
        --presence_matrix presence_matrix.rescued.tsv --cluster_tsv tier1_cluster.tsv \
        --gene_positions "$PG/gene_positions.tsv.zst" \
        --starships "$AN/inputs/af293_starships_table_s7.tsv" \
        --presence_truth "$STUDY/results/id_crosswalk/ground_truth_starships_by_short.tsv" \
        --subset_strains subset_strains.txt --out_prefix "$D/starship"

    # Small per-run tables go back to /bigdata; the large pair tables stay in $SCRATCH.
    zstd -dc cooccurring_pairs.tsv.zst | tail -n +2 | wc -l > n_cooccurring_pairs.txt
    zstd -dc pair_classification.tsv.zst | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="classification")c=i;next}{n[$c]++}END{for(k in n)print k"\t"n[k]}' \
        > classification_counts.tsv
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="bin")c=i;next}{n[$c]++}END{for(k in n)print k"\t"n[k]}' \
        frequency_table.tsv > bin_counts.tsv
    zstd -q -f significant_islands.tsv -o significant_islands.tsv.zst
    cp subset_strains.txt timings.tsv cooccurrence.log n_cooccurring_pairs.txt \
        classification_counts.tsv bin_counts.tsv significant_islands.tsv.zst \
        starship.summary.json starship.starships.tsv "$OUT/"
    cd "$WORK"
    rm -rf "$D"
done
