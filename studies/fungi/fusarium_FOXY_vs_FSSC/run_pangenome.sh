#!/usr/bin/bash
# Run nf_NovInvenio's pangenome.nf for one direction of the reciprocal
# Fusarium oxysporum SC (FOXY) vs F. solani SC (FSSC) comparison.
#
# Usage:
#   run_pangenome.sh <FOXY_in_FSSC_out|FSSC_in_FOXY_out> [extra nextflow args...]
# Normally started by the per-direction SLURM head script
#   .nf_launch/<direction>/submit_nextflow_head.sh
#
# All paths are absolute (no BASH_SOURCE). Each direction has its own
# launch dir (.nf_launch/<direction>) so the two runs never share Nextflow
# session state, and its own --outdir under this study's results/.
#
# Pipeline: a pinned, detached worktree of NovInvenio at af6fd68 (not the
# live checkout, which other sessions change). Override with
# NII_PIPELINE_DIR. The commit is printed at start.
#
# Rescue: RUN_RESCUE=false (the default for the first runs) passes
# -params-file params_norescue.yaml (pangenome_rescue_enable: false). The CLI
# form "--pangenome_rescue_enable false" did NOT disable rescue under
# Nextflow 26.04.6. See README.md "Compute estimate".
set -euo pipefail

DIRECTION="${1:?Usage: run_pangenome.sh <FOXY_in_FSSC_out|FSSC_in_FOXY_out> [extra nextflow args...]}"
shift || true

STUDY_DIR="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/fusarium_FOXY_vs_FSSC"
PIPELINE_DIR="${NII_PIPELINE_DIR:-/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/fusarium-runs-af6fd68}"
PFAM_HMM="/bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm"
RUN_RESCUE="${RUN_RESCUE:-false}"

case "$DIRECTION" in
    FOXY_in_FSSC_out|FSSC_in_FOXY_out) SAMPLESHEET="$STUDY_DIR/config_${DIRECTION}.csv" ;;
    *) echo "ERROR: direction must be FOXY_in_FSSC_out or FSSC_in_FOXY_out" >&2; exit 1 ;;
esac
[ -f "$SAMPLESHEET" ] || { echo "ERROR: $SAMPLESHEET not found (run bin/build_reciprocal_configs.py)" >&2; exit 1; }

if [ "$RUN_RESCUE" = "true" ]; then
    TAG="mmseqs_rescue_${DIRECTION}"
    RESCUE_ARGS=()
else
    TAG="mmseqs_norescue_${DIRECTION}"
    RESCUE_ARGS=(-params-file "$STUDY_DIR/params_norescue.yaml")
fi

OUTDIR="$STUDY_DIR/results/$TAG"
LAUNCH_DIR="$STUDY_DIR/.nf_launch/$TAG"
mkdir -p "$OUTDIR" "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

echo "== pipeline dir: $PIPELINE_DIR =="
echo "== pipeline commit: $(git -C "$PIPELINE_DIR" rev-parse HEAD) =="
echo "== samplesheet: $SAMPLESHEET  rescue: $RUN_RESCUE =="

export PATH="$HOME/.pixi/bin:$PATH"
nextflow run "$PIPELINE_DIR/pangenome.nf" \
    --pangenome_samplesheet "$SAMPLESHEET" \
    --pangenome_data_dir "$STUDY_DIR/data_dir" \
    --pangenome_cluster_backend mmseqs \
    --pangenome_project "fusarium_${TAG}" \
    --outdir "$OUTDIR" \
    --pangenome_island_pfam_hmm "$PFAM_HMM" \
    "${RESCUE_ARGS[@]}" \
    -profile slurm \
    -c "$PIPELINE_DIR/conf/ucr_hpcc_slurm.config" \
    -c "$STUDY_DIR/stajichlab_queue.config" \
    -with-trace "$OUTDIR/trace.txt" \
    -resume \
    "$@"
