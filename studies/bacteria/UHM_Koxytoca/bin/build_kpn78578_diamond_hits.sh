#!/usr/bin/env bash
# One-off: build the diamond_fasta lookup hits for the Kpn78578 model organism
# (studies/bacteria/UHM_Koxytoca/modelorgs.yaml). Ingroup proteins are
# prodigal-called locus tags with no direct relationship to Kpn78578's UniProt
# accessions, so lib/model_organisms.py's `id_transform: diamond_fasta` needs a
# precomputed best-hit TSV (query_id -> ref_protein_id) it can just look up --
# it does NOT run diamond itself (see NovInvenio/lib/model_organisms.py).
#
# Run this AFTER fetch_kpn78578_modelorg.py (needs
# config_support/modelorgs/Kpn78578_protein.faa) and re-run whenever the ingroup
# .faa set changes.
#
# This is a study-specific script (lives under this study's own bin/, not NII's
# shared bin/ -- see CLAUDE.md's "Where new code goes"), so NII_ROOT is hardcoded
# below rather than derived via $(dirname "${BASH_SOURCE[0]}").
#
# Usage: studies/bacteria/UHM_Koxytoca/bin/build_kpn78578_diamond_hits.sh [--ingroup-faa-dir DIR]
set -euo pipefail

INGROUP_FAA_DIR="/bigdata/stajichlab/jpere468/klebsiella_story/koxytoca_ingroup_faa"
if [[ "${1:-}" == "--ingroup-faa-dir" ]]; then
    INGROUP_FAA_DIR="$2"
fi

if ! command -v diamond >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load diamond/2.1.12 2>/dev/null || true
fi
command -v diamond >/dev/null 2>&1 || {
    echo "ERROR: \`diamond\` not found on PATH -- module load diamond." >&2
    exit 1
}

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
MODELORG_DIR="$NII_ROOT/config_support/modelorgs"
REF_FASTA="$MODELORG_DIR/Kpn78578_protein.faa"
[[ -f "$REF_FASTA" ]] || { echo "ERROR: $REF_FASTA not found -- run studies/bacteria/UHM_Koxytoca/bin/fetch_kpn78578_modelorg.py first." >&2; exit 1; }

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

diamond makedb --in "$REF_FASTA" --db "$WORKDIR/Kpn78578_ref" --quiet

QUERY_FASTA="$WORKDIR/ingroup_query.faa"
cat "$INGROUP_FAA_DIR"/*.faa > "$QUERY_FASTA"

OUT_TSV="$MODELORG_DIR/Kpn78578_vs_UHM_ingroup.diamond.tsv"
diamond blastp \
    --query "$QUERY_FASTA" \
    --db "$WORKDIR/Kpn78578_ref" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --max-target-seqs 1 --evalue 1e-5 --threads 2 --block-size 0.4 --index-chunks 4 \
    --out "$OUT_TSV" \
    --quiet

echo "Wrote $OUT_TSV ($(wc -l < "$OUT_TSV") hits)" >&2
