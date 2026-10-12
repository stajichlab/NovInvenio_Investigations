#!/usr/bin/bash
# Run nf_NovInvenio's pangenome.nf for the FSSC_ambrosia study (IN/OUT from config.csv).
#
# Usage: run_pangenome.sh [extra nextflow args...]
# Normally started by .nf_launch/<tag>/submit_nextflow_head.sh.
#
# Adapted from studies/fungi/fusarium_FOXY_vs_FSSC/run_pangenome.sh.
# All paths are absolute (no BASH_SOURCE). Own launch dir (.nf_launch/<tag>)
# so the run never shares Nextflow session state with another run.
#
# Pipeline: a pinned, detached worktree of NovInvenio at a72b4bc (origin/main
# on 2026-09-25), with its own pixi env. Override with NII_PIPELINE_DIR.
#
# Rescue: RUN_RESCUE=false (default) passes -params-file params_norescue.yaml
# (pangenome_rescue_enable: false). Do not use "--pangenome_rescue_enable
# false" on the command line: under Nextflow 26.04.6 it arrives as the String
# "false", which is truthy (see ../fusarium_FOXY_vs_FSSC/README.md).
set -euo pipefail

STUDY_DIR="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/FSSC_ambrosia"
PIPELINE_DIR="${NII_PIPELINE_DIR:-/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/fssc-runs-a72b4bc}"
PFAM_HMM="/bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm"
RUN_RESCUE="${RUN_RESCUE:-false}"
SAMPLESHEET="$STUDY_DIR/config.csv"

if [ "$RUN_RESCUE" = "true" ]; then
    TAG="mmseqs_rescue"
    RESCUE_ARGS=()
else
    TAG="mmseqs_norescue"
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
    --pangenome_project "FSSC_ambrosia_${TAG}" \
    --outdir "$OUTDIR" \
    --pangenome_island_pfam_hmm "$PFAM_HMM" \
    "${RESCUE_ARGS[@]}" \
    -profile slurm \
    -c "$PIPELINE_DIR/conf/ucr_hpcc_slurm.config" \
    -c "$STUDY_DIR/stajichlab_queue.config" \
    -with-trace "$OUTDIR/trace.txt" \
    -resume \
    "$@"
