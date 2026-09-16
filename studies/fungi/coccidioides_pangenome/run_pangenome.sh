#!/usr/bin/bash
# Run nf_NovInvenio's pangenome.nf against this study, for one of three
# variants (whole-set / immitis / posadasii). bin/run_study.sh isn't usable
# here -- it's wired to main.nf's --config/--data_dir params, not
# pangenome.nf's --pangenome_samplesheet/--pangenome_data_dir.
#
# Usage: studies/fungi/coccidioides_pangenome/run_pangenome.sh <wholeset|immitis|posadasii> [extra nextflow args...]
#
# NII_PIPELINE_DIR must point at a local nf_NovInvenio checkout (pangenome.nf
# isn't runnable via a bare `stajichlab/nf_NovInvenio` git-fetch the way
# main.nf is -- confirm this is still true before relying on it, or set
# NII_PIPELINE_DIR explicitly either way):
#   NII_PIPELINE_DIR=/bigdata/stajichlab/jstajich/projects/NovInvenio \
#       studies/fungi/coccidioides_pangenome/run_pangenome.sh wholeset

set -euo pipefail

VARIANT="${1:?Usage: run_pangenome.sh <wholeset|immitis|posadasii> [extra nextflow args...]}"
shift || true

STUDY_DIR="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome"
PIPELINE_DIR="${NII_PIPELINE_DIR:-/bigdata/stajichlab/jstajich/projects/NovInvenio}"

case "$VARIANT" in
    wholeset)  SAMPLESHEET="$STUDY_DIR/config.csv" ;;
    immitis)   SAMPLESHEET="$STUDY_DIR/config_immitis.csv" ;;
    posadasii) SAMPLESHEET="$STUDY_DIR/config_posadasii.csv" ;;
    *) echo "ERROR: variant must be wholeset, immitis, or posadasii" >&2; exit 1 ;;
esac

if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: $SAMPLESHEET not found -- run build_coccidioides_species_csv.py, build_study_config.py, and filter_config_by_taxon.py first" >&2
    exit 1
fi

# Named per README's ask: keep mmseqs results distinguishable from a future
# diamond rerun (diamond backend is currently hard-disabled in pangenome.nf).
OUTDIR="$STUDY_DIR/results/mmseqs_${VARIANT}"
mkdir -p "$OUTDIR"

LAUNCH_DIR="$STUDY_DIR/.nf_launch/${VARIANT}"
mkdir -p "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

# PIPELINE_DIR is a live branch another agent may still be committing to --
# record exactly which commit this run used.
echo "== pipeline commit: $(git -C "$PIPELINE_DIR" rev-parse HEAD) ==" >&2

nextflow run "$PIPELINE_DIR/pangenome.nf" \
    --pangenome_samplesheet "$SAMPLESHEET" \
    --pangenome_data_dir "$STUDY_DIR/data_dir" \
    --pangenome_project "coccidioides_${VARIANT}" \
    --outdir "$OUTDIR" \
    "$@"
