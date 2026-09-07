#!/usr/bin/env bash
# One-off driver: pull genome DNA (for TBLASTN validation) for the UHM_Koxytoca
# study's 16 NCBI RefSeq outgroup accessions via the existing
# bin/fetch_genome_assembly.py recipe (ephemeral pull -> data/ncbi/, gitignored,
# provenance sidecar per accession -- see DESIGN.md Sec 4).
#
# Only the outgroup needs DNA: nf_NovInvenio's VALIDATE workflow only ever runs
# TBLASTN against outgroup_dna_ch when cluster_tool=pairwise (main.nf:246-249) --
# the ingroup UHM MAG proteins have no DNA requirement here.
#
# Re-run any time to re-populate data/ncbi/ from scratch (e.g. after clearing the
# cache); studies/bacteria/UHM_Koxytoca/bin/build_koxytoca_config.py then
# materializes these into studies/bacteria/UHM_Koxytoca/data_dir/dna/.
#
# This is a study-specific script (lives under the study's own bin/, not NII's
# shared bin/ -- see CLAUDE.md's "Where new code goes"), so paths are hardcoded
# absolute rather than derived via $(dirname "${BASH_SOURCE[0]}"), which breaks
# under SLURM (see NII's own CLAUDE.md and ~/.claude/CLAUDE.md).
#
# Usage: studies/bacteria/UHM_Koxytoca/bin/fetch_koxytoca_outgroup_dna.sh [--outdir data/ncbi]
set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
OUTDIR="$NII_ROOT/data/ncbi"
if [[ "${1:-}" == "--outdir" ]]; then
    OUTDIR="$2"
fi

if ! command -v datasets >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load ncbi_datasets/18.30.1 2>/dev/null || true
fi
command -v datasets >/dev/null 2>&1 || {
    echo "ERROR: \`datasets\` (ncbi-datasets-cli) not found on PATH -- module load ncbi_datasets or install via pixi." >&2
    exit 1
}

BIN_DIR="$NII_ROOT/bin"

# GCF_* accession, Short code -- keyed off studies/bacteria/UHM_Koxytoca/species.csv's
# own OUT rows (Accession, Short columns), duplicated here as a literal list so this
# script has no dependency on species.csv's column order/parsing.
ACCESSIONS=(
    "GCF_050916275.1 KoxO_050916275"
    "GCF_017310465.1 KoxO_017310465"
    "GCF_929608445.1 KoxO_929608445"
    "GCF_040561445.1 KoxO_040561445"
    "GCF_900083795.1 KoxO_900083795"
    "GCF_024918735.1 KoxO_024918735"
    "GCF_902162875.1 KoxO_902162875"
    "GCF_053276915.1 KoxO_053276915"
    "GCF_030343105.1 KoxO_030343105"
    "GCF_030343375.1 KoxO_030343375"
    "GCF_016734995.1 KoxO_016734995"
    "GCF_030343565.1 KoxO_030343565"
    "GCF_900083995.1 KoxO_900083995"
    "GCF_004360035.1 KoxO_004360035"
    "GCF_030344075.1 KoxO_030344075"
    "GCF_051951865.1 KoxO_051951865"
)

for entry in "${ACCESSIONS[@]}"; do
    read -r acc short <<< "$entry"
    python3 "$BIN_DIR/fetch_genome_assembly.py" \
        --accession "$acc" --outdir "$OUTDIR" --short "$short"
done

echo "Done. Now run: studies/bacteria/UHM_Koxytoca/bin/build_koxytoca_config.py --study-dir studies/bacteria/UHM_Koxytoca" >&2
