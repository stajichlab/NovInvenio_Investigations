#!/usr/bin/env bash
# Submit the diamond-sensitivity benchmark: one sbatch head job per study
# (<clade>_dmnd_{default,sensitive,very_sensitive}). See
# notes/diamond-sensitivity/README.md.
#
# Usage: bin/run_diamond_sensitivity_bench.sh [study ...]
#   no arguments = all 9 studies. Prints the submitted job IDs.
#
# The pipeline is a pinned worktree of nf_NovInvenio origin/main (a8b68df), not
# the shared checkout, so no other session's git operations can change the code
# during the runs (see nf_NovInvenio CLAUDE.md, "live-checkout race"). Paths are
# absolute on purpose: BASH_SOURCE does not resolve under SLURM.

set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
PIPE_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/bench-diamond-sensitivity"
PIPE_COMMIT="a8b68df"

if [ "$(git -C "$PIPE_ROOT" rev-parse --short HEAD)" != "$PIPE_COMMIT" ]; then
    echo "ERROR: $PIPE_ROOT is not at $PIPE_COMMIT" >&2
    exit 1
fi
if [ ! -d "$PIPE_ROOT/.pixi/envs/default" ]; then
    echo "ERROR: run 'pixi install' in $PIPE_ROOT first" >&2
    exit 1
fi

if [ "$#" -gt 0 ]; then
    STUDIES=("$@")
else
    STUDIES=()
    for c in pezizo_set1 agaricomycetes sordariales_shallow; do
        for m in default sensitive very_sensitive; do STUDIES+=("${c}_dmnd_${m}"); done
    done
fi

mkdir -p "$NII_ROOT/logs/slurm"
for s in "${STUDIES[@]}"; do
    [ -f "$NII_ROOT/studies/fungi/$s/run_params.txt" ] || { echo "ERROR: no study $s" >&2; exit 1; }
    sbatch --parsable -p stajichlab -N 1 -n 2 --mem 8G -t 2-00:00:00 \
        --job-name "nf-$s" \
        -o "$NII_ROOT/logs/slurm/nf_${s}_%j.out" -e "$NII_ROOT/logs/slurm/nf_${s}_%j.err" \
        --wrap "source /etc/profile.d/modules.sh 2>/dev/null || true; module load nextflow; \
cd $NII_ROOT && NII_PIPELINE=$PIPE_ROOT/main.nf bin/run_study.sh fungi/$s \
-profile slurm -c $PIPE_ROOT/conf/ucr_hpcc_slurm.config -c $NII_ROOT/conf/diamond_sensitivity_bench.config -resume" \
        | sed "s/^/$s\t/"
done
