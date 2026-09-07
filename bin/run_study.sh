#!/usr/bin/bash
# Build a study's config+data_dir (if not already built) and run nf_NovInvenio against it.
#
# Usage: bin/run_study.sh <domain>/<set_name> [extra nextflow args...]
# Example: bin/run_study.sh fungi/pezizo_set1 --pfam_hmm /path/to/Pfam-A.hmm
#
# A study's own studies/<domain>/<set>/run_params.txt (if present) supplies this
# study's committed nextflow params (e.g. --run_tool diamond) -- applied first, so
# any matching flag given on the command line still wins (nextflow/most CLIs take
# the last occurrence of a repeated flag).
#
# PIPELINE resolution: the NovInvenio -> nf_NovInvenio rename is deliberately
# deferred (DESIGN.md Sec 2/9 -- no purge/rename yet), so this defaults to the
# pipeline's *current* GitHub name, which is runnable today via Nextflow's own
# git-fetch support. Override for a local checkout during pipeline development:
#   NII_PIPELINE=/bigdata/stajichlab/jstajich/projects/NovInvenio/main.nf bin/run_study.sh ...
# Once the rename lands, change PIPELINE_DEFAULT below (one line) -- no other
# changes needed here or in any study.
#
# After a successful run, bin/sync_reports.sh merges this study's UniProt-derived
# annotation (gene name/description/GO/Pfam/InterPro -- studies/<domain>/<set>/
# annotations/, built by bin/build_study_config.py) into the presence matrices and
# regenerates novelties/core/losses.html, syncing them into view/ and docs/. This
# needs a local pipeline checkout (same NII_PIPELINE requirement as above) -- when
# NII_PIPELINE isn't a local main.nf path, report sync is skipped with a warning,
# not a hard failure, since the nextflow run itself still succeeded.

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

STUDY_PARAMS=()
if [ -f "$STUDY_DIR/run_params.txt" ]; then
    while IFS= read -r line; do
        line="${line%%#*}"                 # strip comments
        [ -n "${line// }" ] || continue    # skip blank/comment-only lines
        read -ra words <<< "$line"
        STUDY_PARAMS+=("${words[@]}")
    done < "$STUDY_DIR/run_params.txt"
    echo "== applying ${STUDY_DIR}/run_params.txt: ${STUDY_PARAMS[*]} ==" >&2
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
    "${STUDY_PARAMS[@]}" \
    "$@"

if [[ "$PIPELINE" == */main.nf ]]; then
    NOVINVENIO_ROOT="${NOVINVENIO_ROOT:-$(dirname "$PIPELINE")}" "$REPO_ROOT/bin/sync_reports.sh" "$STUDY" \
        || echo "== WARNING: bin/sync_reports.sh failed -- nextflow run itself succeeded, reports just weren't UniProt-annotation-synced ==" >&2
else
    echo "== PIPELINE ($PIPELINE) is not a local checkout -- skipping bin/sync_reports.sh; set NII_PIPELINE=/path/to/local/NovInvenio/main.nf to enable it ==" >&2
fi
