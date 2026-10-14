#!/usr/bin/bash
#SBATCH -p stajichlab -N 1 -n 8 --mem 16G -t 03:00:00 --constraint=milan
#SBATCH --job-name hmmsearch-pfam
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/hmmsearch_%x_%j.out
# Pfam-A hmmsearch (HMMER 3.4, --cut_ga) of one NCBI/JGI Neocallimastigomycota proteome.
# hmmsearch replaces hmmscan: on 500 PiroE2sp proteins (4 cores) hmmsearch took 81 s, hmmscan 362 s,
# and both gave the same 525 (protein, Pfam accession) pairs (job 29943954).
# Usage: sbatch run_hmmsearch.sh <Short>
# Output: analysis/hmmscan/<Short>.pfam.hmmsearch.domtbl.zst (hmmsearch layout: col 1 = protein, col 4 = Pfam name, col 5 = Pfam accession).
set -euo pipefail

SHORT="${1:?usage: sbatch run_hmmsearch.sh <Short>}"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
STUDY="$NII_ROOT/studies/fungi/neocallimastigomycota_vs_early_fungi"
PFAM="/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm"
OUT="$STUDY/analysis/hmmscan"
export PATH="$NII_ROOT/.pixi/envs/default/bin:$PATH"

mkdir -p "$OUT"
cd "${SCRATCH:?}"
cp "$STUDY/data_dir/pep/${SHORT}.pep.fa" in.fa
echo "proteins: $(grep -c '^>' in.fa)"
hmmsearch -h | sed -n 2p
start=$(date +%s)
hmmsearch --cpu "${SLURM_CPUS_ON_NODE:-8}" --cut_ga --noali -o /dev/null --domtblout "${SHORT}.pfam.hmmsearch.domtbl" "$PFAM" in.fa
end=$(date +%s)
echo "hmmsearch seconds: $((end - start))"
zstd -T0 -q "${SHORT}.pfam.hmmsearch.domtbl" -o "$OUT/${SHORT}.pfam.hmmsearch.domtbl.zst"
echo "wrote $OUT/${SHORT}.pfam.hmmsearch.domtbl.zst"
