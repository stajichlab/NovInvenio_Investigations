#!/usr/bin/env bash
# BUSCO (protein mode, fungi_odb12) against all 12 sordariales_shallow proteomes,
# for identifying real single-copy-complete-in-all-proteomes negative controls
# (same method/lineage as pezizo_set1's/agaricomycetes' controls files).

#SBATCH -p stajichlab
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 16G
#SBATCH -t 1-00:00:00
#SBATCH --job-name busco-sordariales
#SBATCH -o logs/slurm/busco_sordariales_%j.out
#SBATCH -e logs/slurm/busco_sordariales_%j.err

set -euo pipefail

REPO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
STUDY_DIR="$REPO_ROOT/studies/fungi/sordariales_shallow"
PEP_DIR="$STUDY_DIR/data_dir/pep"
OUT_DIR="$STUDY_DIR/busco"
mkdir -p "$REPO_ROOT/logs/slurm" "$OUT_DIR"

source /etc/profile.d/modules.sh 2>/dev/null || true
# Compute nodes (unlike the login node) pre-load a baseline miniconda3 module
# that conflicts with busco/6.0.0's own internal `module load miniconda3` --
# the busco module load silently fails to reach PATH unless the baseline one
# is unloaded first (confirmed via a direct srun test).
module unload miniconda3 2>/dev/null || true
module load busco/6.0.0

cd "$OUT_DIR"
for pep in "$PEP_DIR"/*.pep.fa; do
    short="$(basename "$pep" .pep.fa)"
    busco -i "$pep" -o "$short" -m protein -l fungi_odb12 -c 8 -f
done
