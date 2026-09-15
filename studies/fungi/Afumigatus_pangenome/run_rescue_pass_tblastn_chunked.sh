#!/usr/bin/bash
#SBATCH -p stajichlab -c 8 --mem 16gb --time=12:00:00 --array=0-9 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/rescue_pass_tblastn_chunk_%a.log

# Chunked, faster version of run_rescue_pass_tblastn.sh -- for next time, not
# a replacement for a run already in progress. A single 32-thread tblastn job
# is capped at ONE node's cores; this splits the 47,983 tier-1 family-rep
# queries into 10 independent SLURM array tasks (each still searching the
# SAME shared genome database), so SLURM can schedule them across the
# `stajichlab` partition's 5-node pool concurrently -- real wall-clock
# parallelism a single job can never reach, not just a smaller per-task
# thread count. -c 8 per task (not 32): tblastn's threading benefit
# saturates well below 32 for typical workloads, and 10 x 8 = 80 requested
# cores spread across nodes beats 1 x 32 confined to one.
#
# Output is zstd-compressed per chunk (per this session's new general
# storage-compression convention) -- rescue_pass.py's parse_tblastn_hits
# needs updating to read .zst input (or decompress via `zstd -dc` into a
# pipe) before this chunked output can feed it directly; not yet done, see
# the study notes' open items.
#
# Prerequisite: the genome blast database must already exist (built once by
# run_rescue_pass_tblastn.sh's first stage, or by a small standalone
# `makeblastdb` step run ahead of this array job) -- this script does NOT
# rebuild it, since the database is query-independent and rebuilding it once
# per array task would be 10x redundant work for zero benefit.
#
# Query chunks must already exist (bin/split_fasta_chunks.py --n_chunks 10),
# e.g.:
#   pixi run python3 bin/split_fasta_chunks.py \
#       --input results/full_293run/tier1_rep_seq.fasta --n_chunks 10 \
#       --out_prefix results/rescue_pass/chunks/query_chunk

set -euo pipefail

STUDY="/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
RESC="$STUDY/results/rescue_pass"

: "${SLURM_ARRAY_TASK_ID:?run this as a SLURM array job (--array=0-9), not interactively}"
CHUNK_ID=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")
QUERY_CHUNK="$RESC/chunks/query_chunk_${CHUNK_ID}.fa"

if [[ ! -f "$QUERY_CHUNK" ]]; then
    echo "ERROR: $QUERY_CHUNK not found -- run bin/split_fasta_chunks.py first" >&2
    exit 1
fi
if [[ ! -f "$RESC/all_genomes_db.nin" && ! -f "$RESC/all_genomes_db.00.nin" ]]; then
    echo "ERROR: $RESC/all_genomes_db not found -- build it once before submitting this array job" >&2
    exit 1
fi

cd "$NII_ROOT"
echo "Task $SLURM_ARRAY_TASK_ID: tblastn $QUERY_CHUNK vs $RESC/all_genomes_db"
pixi run tblastn \
    -query "$QUERY_CHUNK" \
    -db "$RESC/all_genomes_db" \
    -outfmt "6 std qcovs" \
    -evalue 1e-10 \
    -max_target_seqs 5 \
    -num_threads 8 \
    | pixi run zstd -T8 -o "$RESC/chunks/tblastn_chunk_${CHUNK_ID}.tsv.zst"

echo "Task $SLURM_ARRAY_TASK_ID done: $RESC/chunks/tblastn_chunk_${CHUNK_ID}.tsv.zst"
