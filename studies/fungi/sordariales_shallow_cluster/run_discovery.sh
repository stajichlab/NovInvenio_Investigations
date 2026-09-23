#!/usr/bin/env bash
# sbatch launcher for sordariales_shallow_cluster's Tier C+H (--cluster_tool
# mmseqs) run -- compared against sordariales_shallow's own Tier P (pairwise,
# default) run for the cluster-vs-pairwise-sensitivity investigation's third
# clade. --run_tool diamond / --cluster_tool mmseqs / --hmm_presence_cov 0.3 /
# --hmm_presence_min_residues 100 all come from this study's own
# run_params.txt (read automatically by bin/run_study.sh), not repeated here --
# same convention as pezizo_set1_cluster/cl_pezizo_set1.sh and
# agaricomycetes_mmseqs's own launcher.
#
# Head process on 'stajichlab' (lab-owned, low contention), worker tasks on
# 'preempt' via a study-specific conf_preempt.config override -- same split
# validated working for sordariales_shallow's own Tier P run.

#SBATCH -p stajichlab
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 3-00:00:00
#SBATCH --job-name nf-sordariales-cluster
#SBATCH -o logs/slurm/nf_sordariales_cluster_%j.out
#SBATCH -e logs/slurm/nf_sordariales_cluster_%j.err

set -euo pipefail

REPO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
STUDY_DIR="$REPO_ROOT/studies/fungi/sordariales_shallow_cluster"
mkdir -p "$REPO_ROOT/logs/slurm"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow

export NII_PIPELINE="/bigdata/stajichlab/jstajich/projects/NovInvenio/main.nf"

cd "$REPO_ROOT"
bin/run_study.sh fungi/sordariales_shallow_cluster \
    -profile slurm \
    -c /bigdata/stajichlab/jstajich/projects/NovInvenio/conf/ucr_hpcc_slurm.config \
    -c /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/sordariales_shallow/conf_preempt.config \
    -resume
