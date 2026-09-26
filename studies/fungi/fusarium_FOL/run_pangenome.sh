#!/usr/bin/bash
# Run nf_NovInvenio's pangenome.nf for the F. oxysporum f. sp. lycopersici
# (FOL) study: 21 FOL ingroup strains (IN), 4 FOSC outgroup strains (OUT).
#
# Usage:
#   run_pangenome.sh <norescue|rescue> [extra nextflow args...]
# Normally started by .nf_launch/mmseqs_<mode>/submit_nextflow_head.sh
#
# All paths are absolute (no BASH_SOURCE). Each mode has its own launch dir
# and its own --outdir under this study's results/.
#
# Rescue is set by a params file (params_norescue.yaml / params_rescue.yaml),
# never by --pangenome_rescue_enable on the command line: a CLI "false" is
# the String "false", which is truthy (NovInvenio issue #191).
#
# Pipeline: pinned detached worktree of NovInvenio at af6fd68 (the same one
# as ../fusarium_FOXY_vs_FSSC). Override with NII_PIPELINE_DIR.
set -euo pipefail

MODE="${1:?Usage: run_pangenome.sh <norescue|rescue> [extra nextflow args...]}"
shift || true

STUDY_DIR="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/fusarium_FOL"
PIPELINE_DIR="${NII_PIPELINE_DIR:-/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/fusarium-runs-af6fd68}"
PFAM_HMM="/bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm"

case "$MODE" in
    norescue|rescue) PARAMS="$STUDY_DIR/params_${MODE}.yaml" ;;
    *) echo "ERROR: mode must be norescue or rescue" >&2; exit 1 ;;
esac

TAG="mmseqs_${MODE}"
OUTDIR="$STUDY_DIR/results/$TAG"
LAUNCH_DIR="$STUDY_DIR/.nf_launch/$TAG"
mkdir -p "$OUTDIR" "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

echo "== pipeline dir: $PIPELINE_DIR =="
echo "== pipeline commit: $(git -C "$PIPELINE_DIR" rev-parse HEAD) =="
echo "== params file: $PARAMS =="
cat "$PARAMS"

export PATH="$HOME/.pixi/bin:$PATH"
nextflow run "$PIPELINE_DIR/pangenome.nf" \
    --pangenome_samplesheet "$STUDY_DIR/config.csv" \
    --pangenome_data_dir "$STUDY_DIR/data_dir" \
    --pangenome_cluster_backend mmseqs \
    --pangenome_project "fusarium_FOL_${TAG}" \
    --outdir "$OUTDIR" \
    --pangenome_island_pfam_hmm "$PFAM_HMM" \
    -params-file "$PARAMS" \
    -profile slurm \
    -c "$PIPELINE_DIR/conf/ucr_hpcc_slurm.config" \
    -c "$STUDY_DIR/stajichlab_queue.config" \
    -with-trace "$OUTDIR/trace.txt" \
    -resume \
    "$@"
