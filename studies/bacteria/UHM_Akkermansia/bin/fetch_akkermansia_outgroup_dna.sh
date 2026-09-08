#!/usr/bin/env bash
# One-off driver: pull genome DNA (for TBLASTN validation) for the 5 of 8
# UHM_Akkermansia outgroup clades whose representative proteome is an NCBI WGS
# genome assembly (clade02/03/04/07/08 -- headers like "JAUNER010000020.1_1")
# rather than an internal metashot MAG (clade01/05/06 -- "k141_<contig>_<gene>",
# no public genome to fetch, skipped here).
#
# The 6-letter WGS project prefix embedded in each proteome's contig ids
# (JAUNER, JAFWZG, CAJGCV, CANPYO, JAUNQN) isn't itself a genome assembly
# accession -- resolved once via NCBI's esearch/esummary (db=assembly, free-text
# search on the prefix, cross-checked against the returned record's own "wgs"
# field) to the GCA_* accession each is hardcoded below; see
# studies/bacteria/UHM_Akkermansia/species.csv's own NCBI_Accession column for
# the same mapping. Pulls via the existing bin/fetch_genome_assembly.py recipe
# (ephemeral pull -> data/ncbi/, gitignored, provenance sidecar per accession --
# see DESIGN.md Sec 4).
#
# Re-run any time to re-populate data/ncbi/ from scratch (e.g. after clearing the
# cache); studies/bacteria/UHM_Akkermansia/bin/build_akkermansia_config.py then
# materializes these into studies/bacteria/UHM_Akkermansia/data_dir/dna/.
#
# This is a study-specific script (lives under the study's own bin/, not NII's
# shared bin/ -- see CLAUDE.md's "Where new code goes"), so paths are hardcoded
# absolute rather than derived via $(dirname "${BASH_SOURCE[0]}"), which breaks
# under SLURM (see NII's own CLAUDE.md and ~/.claude/CLAUDE.md).
#
# Usage: studies/bacteria/UHM_Akkermansia/bin/fetch_akkermansia_outgroup_dna.sh [--outdir data/ncbi]
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
    echo "ERROR: \`datasets\` (ncbi-datasets-cli) not found on PATH -- module load ncbi_datasets, or run this via \`pixi run\` (NII's pixi.toml pins it)." >&2
    exit 1
}

BIN_DIR="$NII_ROOT/bin"

# GCA_* accession, Short code -- keyed off studies/bacteria/UHM_Akkermansia/
# species.csv's own OUT rows (NCBI_Accession, Short columns), duplicated here as
# a literal list so this script has no dependency on species.csv's column order.
ACCESSIONS=(
    "GCA_030525455.1 clade02"
    "GCA_017517385.1 clade03"
    "GCA_904501915.1 clade04"
    "GCA_947588645.1 clade07"
    "GCA_030536145.1 clade08"
)

for entry in "${ACCESSIONS[@]}"; do
    read -r acc short <<< "$entry"
    python3 "$BIN_DIR/fetch_genome_assembly.py" \
        --accession "$acc" --outdir "$OUTDIR" --short "$short"
done

echo "Done. Now run: studies/bacteria/UHM_Akkermansia/bin/build_akkermansia_config.py --study-dir studies/bacteria/UHM_Akkermansia" >&2
