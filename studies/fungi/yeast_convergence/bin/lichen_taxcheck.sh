#!/usr/bin/bash
#SBATCH -p stajichlab -N 1 -n 1 -c 32 --mem 128G -t 12:00:00 --constraint=milan
#SBATCH --job-name yc-lichen-tax
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/yc_lichen_tax_%j.out
#
# Taxonomic screen of lichen proteomes: DIAMOND blastp against NCBI nr (clustered,
# taxonomy built in), best hit per protein. Letharia lupina WasteWater1 is the test;
# Ascocoryne sarcoides (non-lichen fungus, same BFD pipeline) is the baseline.
# The two proteomes go into one search so the 195 GB database is read once.
# Summary: bin/lichen_taxcheck_summary.py  ->  analysis/lichen_check/summary.tsv
# Absolute paths; no BASH_SOURCE (SLURM).
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/diamond/2.1.12/bin:$PATH   # `module load diamond` did not put it on PATH in a test

OUT=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/yeast_convergence/analysis/lichen_check
BFD=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation
DB=/srv/projects/db/ncbi/diamond/current/nr_cluster_seq.dmnd
W=${SCRATCH:?}
mkdir -p "$OUT"

: > "$W/q.faa"
for pair in "Lethlupi:Letharia_lupina_WasteWater1" "Ascosarc:Ascocoryne_sarcoides_NRRL_50072"; do
    short=${pair%%:*}; dir=${pair#*:}
    zcat "$BFD/$dir/predict_results/$dir.proteins.fa.gz" | sed "s/^>/>${short}|/" >> "$W/q.faa"
done
grep '^>' "$W/q.faa" | sed 's/^>//; s/ .*//' > "$OUT/queries.txt"

diamond blastp -d "$DB" -q "$W/q.faa" -o "$W/hits.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore staxids skingdoms sphylums sscinames \
    --max-target-seqs 5 --evalue 1e-5 --fast -b 4 -c 1 --threads "${SLURM_CPUS_PER_TASK:?set by -c}" --tmpdir "$W" \
    2> "$OUT/diamond.log"
zstd -q -T0 -o "$OUT/hits.tsv.zst" "$W/hits.tsv"
diamond version > "$OUT/diamond_version.txt"
ls -l "$DB" >> "$OUT/diamond_version.txt"
