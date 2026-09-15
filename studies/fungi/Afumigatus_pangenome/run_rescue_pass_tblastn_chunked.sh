#!/usr/bin/bash
#SBATCH -p stajichlab -c 8 --mem 16gb --time=12:00:00 --array=0-15 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/rescue_pass_tblastn_chunk_%a.log

# Chunked, faster version of run_rescue_pass_tblastn.sh. A single 32-thread
# tblastn job is capped at ONE node's cores; this splits the 47,983 tier-1
# family-rep queries into 16 independent SLURM array tasks (each still
# searching the SAME shared genome database), so SLURM can schedule them
# across the `stajichlab` partition's 5-node pool concurrently -- real
# wall-clock parallelism a single job can never reach, not just a smaller
# per-task thread count. -c 8 per task (not 32): tblastn's threading benefit
# saturates well below 32 for typical workloads, and 16 x 8 = 128 requested
# cores spread across nodes beats 1 x 32 confined to one.
#
# 2026-09-15 CORRECTNESS FIX + rerun: -max_target_seqs was 5, which is a cap
# on distinct SUBJECT SEQUENCES per query across the WHOLE combined
# 295-strain database, not 5 per strain -- verified against the completed
# run's own output: 21,223 of 23,747 queries with any hit (89.4%) hit
# exactly 5 distinct strains, i.e. were silently truncated for any family
# genuinely present in more than a handful of strains (most shell/core
# families). That run's "100,952 ABSENT -> GENOME_ONLY" result is very
# likely an undercount and is being redone from scratch with
# -max_target_seqs 1000 (>3x the 295-strain panel, generous headroom for a
# strain occasionally contributing >1 hit contig). Also bumped from 10 to
# 16 array tasks for more real parallelism (still ~1-1.5h/chunk target per
# the HPCC scatter-gather sizing guidance, now scaled for ~3k queries/chunk
# instead of ~4.8k).
#
# Output is written to $SCRATCH (node-local, not the shared /bigdata
# filesystem) and zstd-compressed AS IT'S PRODUCED (streamed straight from
# tblastn's stdout into zstd, never touching disk uncompressed) -- keeps
# heavy per-task write I/O off the shared filesystem while 16 tasks run
# concurrently, per this session's HPCC $SCRATCH convention. The finished
# .tsv.zst is copied back to $RESC/chunks/ (shared storage) once complete,
# not written there directly during the run.
#
# Prerequisite: the genome blast database must already exist (built once by
# run_rescue_pass_tblastn.sh's first stage) -- this script does NOT rebuild
# it, since the database is query-independent and rebuilding it once per
# array task would be 16x redundant work for zero benefit.
#
# Query chunks must already exist (bin/split_fasta_chunks.py --n_chunks 16):
#   pixi run python3 bin/split_fasta_chunks.py \
#       --input results/full_293run/tier1_rep_seq.fasta --n_chunks 16 \
#       --out_prefix results/rescue_pass/chunks/query_chunk

set -euo pipefail

STUDY="/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
RESC="$STUDY/results/rescue_pass"
: "${SCRATCH:?no \$SCRATCH set -- must run as a SLURM job on a node with node-local scratch}"

: "${SLURM_ARRAY_TASK_ID:?run this as a SLURM array job (--array=0-15), not interactively}"
CHUNK_ID=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")
QUERY_CHUNK="$RESC/chunks/query_chunk_${CHUNK_ID}.fa"
SCRATCH_OUT="$SCRATCH/tblastn_chunk_${CHUNK_ID}.tsv.zst"

if [[ ! -f "$QUERY_CHUNK" ]]; then
    echo "ERROR: $QUERY_CHUNK not found -- run bin/split_fasta_chunks.py first" >&2
    exit 1
fi
if [[ ! -f "$RESC/all_genomes_db.nin" && ! -f "$RESC/all_genomes_db.00.nin" ]]; then
    echo "ERROR: $RESC/all_genomes_db not found -- build it once before submitting this array job" >&2
    exit 1
fi

cd "$NII_ROOT"
echo "Task $SLURM_ARRAY_TASK_ID: tblastn $QUERY_CHUNK vs $RESC/all_genomes_db -> $SCRATCH_OUT"
pixi run tblastn \
    -query "$QUERY_CHUNK" \
    -db "$RESC/all_genomes_db" \
    -outfmt "6 std qcovs" \
    -evalue 1e-10 \
    -max_target_seqs 1000 \
    -num_threads 8 \
    | pixi run zstd -T8 -o "$SCRATCH_OUT"

mkdir -p "$RESC/chunks"
cp "$SCRATCH_OUT" "$RESC/chunks/tblastn_chunk_${CHUNK_ID}.tsv.zst"
echo "Task $SLURM_ARRAY_TASK_ID done: $RESC/chunks/tblastn_chunk_${CHUNK_ID}.tsv.zst"
