#!/bin/bash
#SBATCH -J strat_loci
#SBATCH -p exfab
#SBATCH -A exfab
#SBATCH -c 16
#SBATCH --mem 24G
#SBATCH -t 06:00:00
# Offline repeat of ISLAND_DNA_TARGETS -> ISLAND_DNA_CHECK -> ISLAND_LOCI with frequency-stratified
# candidates (PR branch locus-frequency-strata). Inputs are the full_v072 task directories.
# Usage: sbatch run_stratified_loci.sh   (env: STRATA, NCAND, OUT)
set -euo pipefail
ENVBIN=/rhome/jstajich/.nextflow/assets/.repos/stajichlab/NovInvenio/clones/edfb88399327879b59c2f341a91ad8d30825be37/.pixi/envs/default/bin
export PATH="$ENVBIN:$PATH"
PIPE=/rhome/jstajich/projects/NovInvenio-worktrees/locusviz
# Defaults are the Afumigatus full_v072 run. Override with env vars for another study:
#   STUDY, T (ISLAND_DNA_TARGETS task dir), CFG (samplesheet name in T), REFS, PROJECT, TAG
STUDY=${STUDY:-/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome}
T=${T:-$STUDY/.nf_launch/full_v072/work/c3/c5563d3d30598d1f67db1b1e5e9a0f}
CFG=${CFG:-config.csv}
REFS=${REFS:-Asfu_Af293,Asfu_A1163}
PROJECT=${PROJECT:-config}
STRATA=${STRATA:-0.03,0.10,0.30,0.70}
NCAND=${NCAND:-250}
OUT=${OUT:-$STUDY/analysis/locus_summary/strat}
mkdir -p "$OUT"; cd "$OUT"
COMMON=(--islands_with_domains $T/islands_with_domains.tsv --presence_matrix $T/presence_matrix.rescued.tsv
  --family_positions $T/family_positions.tsv.zst --frequency_table $T/frequency_table.tsv
  --assembly_quality $T/assembly_quality_vs_content.tsv --config $T/$CFG --domtblout $T/pfam.domtblout
  --domain_evalue 0.001 --gene_positions $T/gene_positions.tsv.zst --cluster_tsv $T/tier1_cluster.tsv
  --id_sep '|' --flank 5 --flank_min 3 --k 10 --empty_frac 0.8 --min_column_strains 2 --model_diff_min_frac 0.5
  --island_reference_strains "$REFS" --containment 0.5 --rank_by presence
  --top_loci "$NCAND" --per_rank 20 --poly_min_strains 20 --poly_min_frac 0.05 --poly_max_frac 0.95
  --fixed_diff 0.95 --candidates "$NCAND" --candidate_strata "$STRATA" --rescue_positions $T/rescue_positions.tsv
  --dna_min_id 90 --dna_min_cov 80 --min_strains 2 --clinker_max_strains 0 --project "$PROJECT")
python3 "$PIPE/bin/pangenome_island_loci.py" "${COMMON[@]}" \
  --dna_targets_dir dna_targets --dna_batch 50 --rescue_tblastn $T/rescue_tblastn/*.tsv.zst \
  --output island_loci.pre_dna.json
mkdir -p calls
for b in dna_targets/batch_*.tsv; do
  n=$(basename "$b" .tsv)
  python3 "$PIPE/bin/pangenome_island_dna_check.py" --targets "$b" --config $T/$CFG \
    --data_dir $STUDY/data_dir --min_id 90 --min_cov 80 --cpus "${SLURM_CPUS_PER_TASK:-16}" \
    --output "calls/dna_calls_${n#batch_}.tsv"
done
python3 "$PIPE/bin/pangenome_island_loci.py" "${COMMON[@]}" --dna_check true --dna_calls calls/dna_calls_*.tsv \
  --regions_out island_regions.tsv --output island_loci.json
