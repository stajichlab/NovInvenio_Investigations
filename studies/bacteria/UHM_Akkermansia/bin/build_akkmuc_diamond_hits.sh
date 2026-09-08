#!/usr/bin/env bash
# One-off: build the diamond_fasta lookup hits for the AkkMuc model organism
# (studies/bacteria/UHM_Akkermansia/modelorgs.yaml). Ingroup proteins are
# prodigal-called MAG locus tags with no direct relationship to AkkMuc's UniProt
# accessions, so lib/model_organisms.py's `id_transform: diamond_fasta` needs a
# precomputed best-hit TSV (query_id -> ref_protein_id) it can just look up -- it
# does NOT run diamond itself (see nf_NovInvenio/lib/model_organisms.py).
#
# Unlike UHM_Koxytoca's build_kpn78578_diamond_hits.sh (which pulls its ingroup
# query proteins from an external staging dir), this study's 6 candidate-genus
# ingroup proteomes already live locally, header-fixed, under this study's own
# data_dir/pep/ (built by bin/build_akkermansia_config.py) -- so the query set is
# just those 6 files, named directly rather than resolved from a --ingroup-faa-dir.
#
# Run this AFTER fetch_akkmuc_modelorg.py (needs
# config_support/modelorgs/AkkMuc_protein.faa) and re-run whenever the ingroup
# Short set (species.csv's Source=ingroup rows) changes.
#
# This is a study-specific script (lives under this study's own bin/, not NII's
# shared bin/ -- see CLAUDE.md's "Where new code goes"), so NII_ROOT/STUDY_DIR are
# hardcoded below rather than derived via $(dirname "${BASH_SOURCE[0]}").
#
# Usage: studies/bacteria/UHM_Akkermansia/bin/build_akkmuc_diamond_hits.sh
set -euo pipefail

if ! command -v diamond >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load diamond 2>/dev/null || true
fi
command -v diamond >/dev/null 2>&1 || {
    echo "ERROR: \`diamond\` not found on PATH -- module load diamond." >&2
    exit 1
}

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
STUDY_DIR="$NII_ROOT/studies/bacteria/UHM_Akkermansia"
MODELORG_DIR="$NII_ROOT/config_support/modelorgs"
REF_FASTA="$MODELORG_DIR/AkkMuc_protein.faa"
[[ -f "$REF_FASTA" ]] || { echo "ERROR: $REF_FASTA not found -- run studies/bacteria/UHM_Akkermansia/bin/fetch_akkmuc_modelorg.py first." >&2; exit 1; }

# The 6 candidate-genus ingroup Shorts (species.csv Source=ingroup rows) -- these
# are the only ones modelorgs.yaml has an entry for; see that file's header
# comment for why (a model_organisms entry only ever fires for candidates whose
# OWN source_proteome Short matches its `short:`, same gotcha UHM_Koxytoca hit).
INGROUP_SHORTS=(C286 C287 C288 C289 C294 C298)

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

diamond makedb --in "$REF_FASTA" --db "$WORKDIR/AkkMuc_ref" --quiet

QUERY_FASTA="$WORKDIR/ingroup_query.faa"
: > "$QUERY_FASTA"
for short in "${INGROUP_SHORTS[@]}"; do
    pep="$STUDY_DIR/data_dir/pep/${short}.pep.fa"
    [[ -f "$pep" ]] || { echo "ERROR: expected $pep not found" >&2; exit 1; }
    cat "$pep" >> "$QUERY_FASTA"
done

OUT_TSV="$MODELORG_DIR/AkkMuc_vs_UHM_ingroup.diamond.tsv"
diamond blastp \
    --query "$QUERY_FASTA" \
    --db "$WORKDIR/AkkMuc_ref" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --max-target-seqs 1 --evalue 1e-5 --threads 2 --block-size 0.4 --index-chunks 4 \
    --out "$OUT_TSV" \
    --quiet

echo "Wrote $OUT_TSV ($(wc -l < "$OUT_TSV") hits)" >&2
