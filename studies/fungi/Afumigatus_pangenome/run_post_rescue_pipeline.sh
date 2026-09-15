#!/usr/bin/bash
#SBATCH -p stajichlab -c 1 --mem 32gb --time=6:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/post_rescue_pipeline.log

# Fold the corrected per-strain tblastn rescue results (job 28399806, see
# run_rescue_pass_per_strain.sh) back into the presence matrix and rerun
# every downstream step that depends on it, in one submission -- written
# ahead of that job finishing so there's no dead time between "rescue done"
# and "corrected figures exist".
#
# Everything upstream of the presence matrix (tier1 clustering, gene/family
# positions, the DUF3435 captain-gene hmmsearch) is UNCHANGED by the rescue
# fix -- rescue only flips some ABSENT calls to GENOME_ONLY, it doesn't
# touch protein clustering or genomic coordinates -- so those inputs are
# reused as-is, not regenerated.
#
# --mem 32gb / --time 6h, not the old rescued-matrix run's 64gb/48h: that
# estimate was sized for cooccurrence.py's Monte Carlo permutation null
# (2026-09-15 commit f7726e4 replaced it with an exact closed-form test,
# benchmarked ~460x faster -- projected well under an hour for the
# eligible-family scale this study reaches). 6h/32gb is generous headroom
# over that projection, not a tight fit; revise if the real run disagrees.

set -euo pipefail

STUDY="/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome"
NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
RUN="$STUDY/results/full_293run"
RESC="$STUDY/results/rescue_pass"

export NOVINVENIO_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio"
cd "$NII_ROOT"

N_CHUNKS=$(ls "$RESC/per_strain_chunks"/*.tblastn.tsv.zst 2>/dev/null | wc -l)
echo "Found $N_CHUNKS per-strain tblastn chunks in $RESC/per_strain_chunks"
if [[ "$N_CHUNKS" -lt 295 ]]; then
    echo "ERROR: expected 295 (one per strain) -- job 28399806 may not be finished yet." >&2
    echo "Check: squeue -j 28399806" >&2
    exit 1
fi

echo "Step 1/5: folding rescue hits into the presence matrix..."
TBLASTN_ARGS=()
for f in "$RESC/per_strain_chunks"/*.tblastn.tsv.zst; do
    TBLASTN_ARGS+=(--tblastn_tsv "$f")
done
pixi run python3 "$STUDY/bin/rescue_pass.py" \
    --matrix "$RUN/presence_matrix.tsv" \
    "${TBLASTN_ARGS[@]}" \
    --output "$RUN/presence_matrix.rescued.tsv"

echo "Step 2/5: frequency binning on the rescued matrix..."
pixi run python3 "$STUDY/bin/frequency_bins.py" \
    --matrix "$RUN/presence_matrix.rescued.tsv" \
    --config "$STUDY/config.csv" \
    --inventory "$RUN/strain_inventory.tsv" \
    --output "$RUN/frequency_table.rescued.tsv"

echo "Step 3/5: co-occurrence (exact stratified test, not Monte Carlo)..."
pixi run python3 "$STUDY/bin/cooccurrence.py" \
    --matrix "$RUN/presence_matrix.rescued.tsv" \
    --frequency_table "$RUN/frequency_table.rescued.tsv" \
    --config "$STUDY/config.csv" \
    --inventory "$RUN/strain_inventory.tsv" \
    --output "$RUN/cooccurring_pairs.rescued.tsv"

echo "Step 4/5: pair classification (physical linkage / Starship mechanism)..."
pixi run python3 "$STUDY/bin/pair_classification.py" \
    --cooccurring_pairs "$RUN/cooccurring_pairs.rescued.tsv" \
    --family_positions "$RUN/family_positions.tsv" \
    --cluster_tsv "$RUN/tier1_cluster.tsv" \
    --captain_tblout "$STUDY/results/captain_gene/DUF3435_vs_study.tblout" \
    --output "$RUN/pair_classification.rescued.tsv"

echo "Step 5/5: regenerating summary figures from the corrected data..."
pixi run python3 "$STUDY/bin/plot_pangenome_summary.py" \
    --frequency_table "$RUN/frequency_table.rescued.tsv" \
    --matrix "$RUN/presence_matrix.rescued.tsv" \
    --pair_classification "$RUN/pair_classification.rescued.tsv" \
    --out_dir "$STUDY/results/figures_rescued"

echo "Done. Corrected outputs:"
echo "  $RUN/presence_matrix.rescued.tsv"
echo "  $RUN/frequency_table.rescued.tsv"
echo "  $RUN/cooccurring_pairs.rescued.tsv"
echo "  $RUN/pair_classification.rescued.tsv"
echo "  $STUDY/results/figures_rescued/"
