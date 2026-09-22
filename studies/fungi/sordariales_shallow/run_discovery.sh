#!/usr/bin/env bash
# sbatch launcher for sordariales_shallow's Tier P (pairwise) discovery pass --
# the third clade's exploratory run for building a real controls.csv (positive
# controls sourced from literature cross-referencing of candidates.txt, not
# guessed). See notes/superpowers/specs/2026-09-13-cluster-vs-pairwise-sensitivity-design.md
# for why Tier P is the sensitivity reference used for this kind of discovery.
#
# Submit with: sbatch studies/fungi/sordariales_shallow/run_discovery.sh
#
# Head process runs on 'epyc' (long time limit, light footprint) -- NOT
# 'preempt' -- since a preempted head would kill the whole run's coordination.
# Per-task worker processes route to 'preempt' via conf_preempt.config (a
# study-specific override, not a change to the shared
# conf/ucr_hpcc_slurm.config) -- this is a cheap, retry-tolerant exploratory
# run, not a production study, so preemption risk on individual tasks is
# acceptable.

#SBATCH -p epyc
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 3-00:00:00
#SBATCH --job-name nf-sordariales-discovery
#SBATCH -o logs/slurm/nf_sordariales_discovery_%j.out
#SBATCH -e logs/slurm/nf_sordariales_discovery_%j.err

set -euo pipefail

REPO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
STUDY_DIR="$REPO_ROOT/studies/fungi/sordariales_shallow"
mkdir -p "$REPO_ROOT/logs/slurm"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow

# Local pipeline checkout (enables bin/sync_reports.sh after a successful run,
# and lets this run pick up any uncommitted pipeline-side fixes without
# waiting on a release).
export NII_PIPELINE="/bigdata/stajichlab/jstajich/projects/NovInvenio/main.nf"

cd "$REPO_ROOT"
bin/run_study.sh fungi/sordariales_shallow \
    -profile slurm \
    -c /bigdata/stajichlab/jstajich/projects/NovInvenio/conf/ucr_hpcc_slurm.config \
    -c "$STUDY_DIR/conf_preempt.config" \
    -resume \
    --run_tool diamond
