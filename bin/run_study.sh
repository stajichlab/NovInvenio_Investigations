#!/usr/bin/bash
# Build a study's config+data_dir (if not already built) and run nf_NovInvenio against it.
#
# Usage: bin/run_study.sh <domain>/<set_name> [extra nextflow args...]
# Example: bin/run_study.sh fungal/pezizo_set1 --run_tool diamond --pfam_hmm /path/to/Pfam-A.hmm
#
# PIPELINE resolution: the NovInvenio -> nf_NovInvenio rename is deliberately
# deferred (DESIGN.md Sec 2/9 -- no purge/rename yet), so this defaults to the
# pipeline's *current* GitHub name, which is runnable today via Nextflow's own
# git-fetch support. Override for a local checkout during pipeline development:
#   NII_PIPELINE=/bigdata/stajichlab/jstajich/projects/NovInvenio/main.nf bin/run_study.sh ...
# Once the rename lands, change PIPELINE_DEFAULT below (one line) -- no other
# changes needed here or in any study.

set -euo pipefail

PIPELINE_DEFAULT="stajichlab/NovInvenio"   # TODO: -> stajichlab/nf_NovInvenio after rename

STUDY="${1:?Usage: bin/run_study.sh <domain>/<set_name> [extra nextflow args...]}"
shift || true

REPO_ROOT="$(cd "$(dirname "${0}")/.." && pwd)"
STUDY_DIR="$REPO_ROOT/studies/$STUDY"
PIPELINE="${NII_PIPELINE:-$PIPELINE_DEFAULT}"

if [ ! -f "$STUDY_DIR/species.csv" ]; then
    echo "ERROR: $STUDY_DIR/species.csv not found" >&2
    exit 1
fi

if [ ! -f "$STUDY_DIR/config.csv" ] || [ ! -d "$STUDY_DIR/data_dir" ]; then
    echo "== building config + data_dir for $STUDY ==" >&2
    python3 "$REPO_ROOT/bin/build_study_config.py" --study-dir "$STUDY_DIR"
else
    echo "== $STUDY_DIR/config.csv + data_dir already present, skipping fetch (delete to rebuild) ==" >&2
fi

LAUNCH_DIR="$REPO_ROOT/.nf_launch/$STUDY"
mkdir -p "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

echo "== running $PIPELINE against $STUDY ==" >&2
nextflow run "$PIPELINE" \
    --config "$STUDY_DIR/config.csv" \
    --data_dir "$STUDY_DIR/data_dir" \
    --project "$(basename "$STUDY")" \
    --outdir "$REPO_ROOT/results" \
    "$@"
