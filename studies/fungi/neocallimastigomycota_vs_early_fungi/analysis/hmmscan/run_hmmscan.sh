#!/usr/bin/bash
#SBATCH -p stajichlab
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem 16G
#SBATCH -t 04:00:00
#SBATCH --job-name hmmscan-pilot
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/hmmscan_%x_%j.out
#SBATCH -e /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/hmmscan_%x_%j.err

# Pfam-A hmmscan of one NCBI/JGI Neocallimastigomycota proteome (annotation gap: these
# proteomes have no UniProt annotation). Pilot to measure runtime before sizing the other three.
# Usage: sbatch run_hmmscan.sh <Short>      (default PiroE2sp)
# Output (compressed): <outdir>/<Short>.pfam.domtbl.gz, plus a timing line in the slurm log.
set -euo pipefail

SHORT="${1:-PiroE2sp}"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
STUDY="$NII_ROOT/studies/fungi/neocallimastigomycota_vs_early_fungi"
PFAM="/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm"
OUT="$STUDY/analysis/hmmscan"
export PATH="$NII_ROOT/.pixi/envs/default/bin:$PATH"

mkdir -p "$OUT"
cd "${SCRATCH:?}"
cp "$STUDY/data_dir/pep/${SHORT}.pep.fa" in.fa
echo "proteins: $(grep -c '^>' in.fa)"
hmmscan -h | sed -n 2p
start=$(date +%s)
hmmscan --cpu "${SLURM_CPUS_ON_NODE:-16}" --cut_ga --noali -o /dev/null --domtblout "${SHORT}.pfam.domtbl" "$PFAM" in.fa
end=$(date +%s)
echo "hmmscan seconds: $((end - start))"
zstd -T0 -q "${SHORT}.pfam.domtbl" -o "$OUT/${SHORT}.pfam.domtbl.zst"
echo "wrote $OUT/${SHORT}.pfam.domtbl.zst"
