#!/usr/bin/env bash
# sbatch launcher: sordariales_shallow pairwise (diamond) re-run on NovInvenio PR #216
# (issue #208: other-group evidence + outgroup signal class) to check the new sidecars on
# real data. Writes to a NEW project, sordariales_shallow_ev208, so the results behind the
# published pairwise-diamond page are not touched. Not resumed: the old work dir is gone.
#
# Submit with: sbatch studies/fungi/sordariales_shallow/run_ev208.sh

#SBATCH -p stajichlab
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 3-00:00:00
#SBATCH --job-name nf-sordariales-ev208
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/nf_sordariales_ev208_%j.out
#SBATCH -e /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/nf_sordariales_ev208_%j.err

set -euo pipefail

REPO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
STUDY_DIR="$REPO_ROOT/studies/fungi/sordariales_shallow"
PIPE="/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/i208"   # i208-outgroup-evidence, PR #216
mkdir -p "$REPO_ROOT/logs/slurm"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow

LAUNCH_DIR="$REPO_ROOT/.nf_launch/fungi/sordariales_shallow_ev208"
mkdir -p "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

nextflow run "$PIPE/main.nf" \
    --config "$STUDY_DIR/config.csv" \
    --data_dir "$STUDY_DIR/data_dir" \
    --project sordariales_shallow_ev208 \
    --outdir "$REPO_ROOT/results" \
    -profile slurm \
    -c "$PIPE/conf/ucr_hpcc_slurm.config" \
    -c "$STUDY_DIR/conf_preempt.config" \
    --run_tool diamond
