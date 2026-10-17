#!/usr/bin/bash
#SBATCH -p stajichlab -N 1 -n 8 --mem 16G -t 03:00:00 --constraint=milan
#SBATCH -J pfam_reps
# Pfam-A hmmsearch (HMMER 3.4, --cut_ga) of a range of family representative proteins
# (mmseqs tier1_rep_seq.fasta of one nf_NovInvenio pangenome run). Same settings as
# bin/annotate_proteomes_pfam.py chunk_script(). --constraint=milan: AVX2 nodes only.
# Absolute paths only (no BASH_SOURCE). Raw files stay in node-local $SCRATCH.
# Usage: sbatch -o <log> pfam_reps_chunk.sbatch <rep_fasta> <start(0-based)> <end> <out.domtbl.zst>
set -euo pipefail
FASTA=$1; START=$2; END=$3; OUT=$4
PFAM_HMM=/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm
export PATH=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/.pixi/envs/default/bin:"$PATH"
W="${SCRATCH:?SCRATCH must be set (run under SLURM)}/pfam_reps_${SLURM_JOB_ID}"
mkdir -p "$W"
N=$(( END - START ))
awk -v s="$START" -v e="$END" '/^>/{n++} n>s && n<=e' "$FASTA" > "$W/in.fa"
got=$(grep -c '^>' "$W/in.fa" || true)
echo "proteins: $got (expected $N) from $FASTA"
if [ "$got" -ne "$N" ]; then echo "ERROR: chunk has $got proteins, expected $N" >&2; exit 1; fi
hmmsearch -h | sed -n 2p
echo "host: $(hostname)"
t0=$(date +%s)
hmmsearch --cpu "${SLURM_CPUS_ON_NODE:-8}" --cut_ga --noali -o /dev/null \
    --domtblout "$W/out.domtbl" "$PFAM_HMM" "$W/in.fa"
echo "hmmsearch seconds: $(( $(date +%s) - t0 ))"
tail -n 1 "$W/out.domtbl" | grep -q '^# \[ok\]' || { echo "ERROR: domain table has no [ok] line" >&2; exit 1; }
zstd -T0 -q -f "$W/out.domtbl" -o "$OUT.tmp"
mv -f "$OUT.tmp" "$OUT"
rm -rf "$W"
echo "wrote $OUT"
