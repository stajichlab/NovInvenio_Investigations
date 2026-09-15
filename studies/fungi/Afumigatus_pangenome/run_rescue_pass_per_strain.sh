#!/usr/bin/bash
#SBATCH -p stajichlab -c 8 --mem 16gb --time=2-00:00:00 --array=0-23 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/rescue_pass_per_strain_%a.log

# Genome-level tblastn rescue pass, REDESIGNED 2026-09-15 from a single
# combined-genome-database search to a per-strain search (Fable/independent
# pipeline review + verified against the prior run's own output).
#
# WHY the redesign, not just a bigger -max_target_seqs: the original design
# (run_rescue_pass_tblastn.sh / run_rescue_pass_tblastn_chunked.sh) searched
# all 47,983 tier-1 family-rep queries against ONE combined database of all
# 295 strains' genomes. -max_target_seqs 5 there caps DISTINCT SUBJECT
# SEQUENCES **across the whole combined database, per query** -- not per
# strain. Verified against that run's real output: 21,223 of 23,747 hit
# queries (89.4%) landed at exactly 5 distinct strains hit, meaning any
# family genuinely present in more than a handful of strains (most
# shell/core families) had its true rescue candidates silently dropped once
# 5 strains' hits filled the budget. A per-strain database has no such
# cross-strain competition at all -- max_target_seqs only has to be large
# enough to catch real intra-genome duplication (e.g. this study's own
# DUF3435 Starship-captain gene, median 6-8 copies/strain), which is cheap
# to set generously (200, well above any observed copy number) since each
# per-strain database is tiny.
#
# This also cuts real search volume, not just fixes correctness: each
# strain's query set is restricted to only the families CURRENTLY ABSENT in
# that strain (bin/extract_absent_family_queries.py), skipping the ~20% of
# (family, strain) cells already known PRESENT/GENOME_ONLY -- the original
# design's "restricting the query set barely reduces work" reasoning was
# specific to a single shared combined database (search cost there scales
# with the FULL 8GB db regardless of query count); it doesn't hold once
# each strain gets its own ~1/295th-sized database.
#
# Load balance: 295 strains split across N_TASKS array tasks, NOT evenly by
# strain count -- per-strain query-set size (n_absent_families) varies, so
# strains are sorted by descending query count and dealt round-robin across
# tasks (a simple longest-processing-time-first greedy schedule) to keep
# task runtimes comparable rather than letting one task draw all the
# expensive strains.
#
# $SCRATCH staging: the per-strain genome FASTA (Short-prefixed), its BLAST
# database, and the tblastn run itself all happen on $SCRATCH (node-local),
# not the shared /bigdata filesystem -- per this session's HPCC $SCRATCH
# convention, this keeps 24 concurrent tasks' database-build + search I/O
# off shared storage. Output is compressed AS PRODUCED (tblastn's stdout
# piped straight into zstd, never touching disk uncompressed) and only the
# final small .tsv.zst is copied back to shared storage.
#
# Runtime estimate (NOT yet measured for this redesign): the original
# combined-DB 10-way query-chunked run measured ~4h20m/task (8 cores,
# ~4,798 queries against the full 8GB db, capped output). This redesign's
# per-task cost is roughly (queries_per_strain x that strain's own genome
# size, summed over the ~12 strains/task at 24 tasks) -- smaller per-pair
# database but a comparable total query count, so expect a similar order of
# magnitude; --time is set generously (2 days) until a real first run gives
# a measured number to revise this comment and the array task count by.
#
# Prerequisite: bin/extract_absent_family_queries.py must already have been
# run once (cheap, ~1 min, no SLURM job needed) to produce
# results/rescue_pass/per_strain_queries/<Short>.absent.fa + manifest.tsv:
#   pixi run python3 bin/extract_absent_family_queries.py \
#       --matrix results/full_293run/presence_matrix.tsv \
#       --rep_fasta results/full_293run/tier1_rep_seq.fasta \
#       --out_dir results/rescue_pass/per_strain_queries

set -euo pipefail

STUDY="/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
RESC="$STUDY/results/rescue_pass"
QDIR="$RESC/per_strain_queries"
OUTDIR="$RESC/per_strain_chunks"
CONFIG="$STUDY/config.csv"
N_TASKS=24

: "${SCRATCH:?no \$SCRATCH set -- must run as a SLURM job on a node with node-local scratch}"
: "${SLURM_ARRAY_TASK_ID:?run this as a SLURM array job (--array=0-23), not interactively}"

if [[ ! -f "$QDIR/manifest.tsv" ]]; then
    echo "ERROR: $QDIR/manifest.tsv not found -- run bin/extract_absent_family_queries.py first" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
cd "$NII_ROOT"

declare -A DNA_OF
while IFS=, read -r GROUP SPECIES STRAIN PROTEIN DNA GFF3 SHORT TAXON; do
    [[ "$GROUP" == "GROUP" ]] && continue
    DNA_OF["$SHORT"]="$DNA"
done < "$CONFIG"

mapfile -t MY_STRAINS < <(tail -n +2 "$QDIR/manifest.tsv" | sort -t$'\t' -k2,2nr | \
    awk -F'\t' -v tid="$SLURM_ARRAY_TASK_ID" -v n="$N_TASKS" 'NR % n == tid {print $1}')

echo "Task $SLURM_ARRAY_TASK_ID: ${#MY_STRAINS[@]} strains: ${MY_STRAINS[*]}"

for SHORT in "${MY_STRAINS[@]}"; do
    DNA_FILE="${DNA_OF[$SHORT]:-}"
    if [[ -z "$DNA_FILE" ]]; then
        echo "ERROR: no DNA file for $SHORT in $CONFIG" >&2
        exit 1
    fi
    GENOME_SRC="$STUDY/data_dir/dna/$DNA_FILE"
    QUERY_SRC="$QDIR/${SHORT}.absent.fa"
    if [[ ! -f "$QUERY_SRC" ]]; then
        echo "  $SHORT: no absent-family queries, skipping" >&2
        continue
    fi

    WORK="$SCRATCH/$SHORT"
    mkdir -p "$WORK"
    awk -v s="$SHORT" '/^>/{sub(/^>/, ">" s "|")} {print}' "$GENOME_SRC" > "$WORK/genome.fa"
    pixi run makeblastdb -in "$WORK/genome.fa" -dbtype nucl -out "$WORK/genome_db" > /dev/null

    N_QUERIES=$(grep -c "^>" "$QUERY_SRC")
    echo "  $SHORT: tblastn ($N_QUERIES absent-family queries) vs own genome"
    pixi run tblastn \
        -query "$QUERY_SRC" \
        -db "$WORK/genome_db" \
        -outfmt "6 std qcovs" \
        -evalue 1e-10 \
        -max_target_seqs 200 \
        -num_threads 8 \
        | pixi run zstd -T8 -o "$WORK/${SHORT}.tblastn.tsv.zst"

    cp "$WORK/${SHORT}.tblastn.tsv.zst" "$OUTDIR/${SHORT}.tblastn.tsv.zst"
    rm -rf "$WORK"
done

echo "Task $SLURM_ARRAY_TASK_ID done."
