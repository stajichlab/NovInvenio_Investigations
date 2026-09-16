#!/usr/bin/bash
#SBATCH -p stajichlab -c 8 --mem 32gb --time=4:00:00 --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/post_rescue_pipeline.log

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
# --mem 32gb / --time 4h / -c 8, not the old rescued-matrix run's 64gb/48h:
# that estimate was sized for cooccurrence.py's Monte Carlo permutation null
# (2026-09-15 commit f7726e4 replaced it with an exact closed-form test,
# benchmarked ~460x faster). Real measured timing from the first full run
# (2026-09-15): steps 1-3 (fold-back through co-occurrence) ~48 min total;
# step 4 (extract_rescue_positions.py, re-parsing all 295 per-strain tblastn
# outputs) ~68 min; step 6 (pair_classification.py against the ~3x-larger
# merged family_positions) ~40 min.
#
# Step 4's 68 minutes was NOT actually I/O or file-parsing cost (raw zstd
# decompression of all 295 files takes seconds) -- profiling found the real
# cost was `family not in matrix.families`, an O(47,983) linear scan on a
# LIST executed once per entry in the hit-position map (hundreds of
# thousands to millions of times). Fixed by converting to a set inside
# extract_rescue_positions.py (~7.7x speedup measured on a 20-file subset:
# 290s -> 38s, identical output) -- ALWAYS worth checking for this pattern
# before reaching for multiprocessing, since it fixes the actual bottleneck
# instead of parallelizing something that wasn't the real cost. Added
# `--processes` multiprocessing on top of that (parses independent files
# concurrently) for a further, smaller, genuine win once the real bottleneck
# was gone (~1.2x further on the same subset) -- `-c 8` here matches the
# `--processes 8` passed to extract_rescue_positions.py below. 4h is
# generous headroom over the (now much shorter) expected total; revise once
# a second real run at full 295-strain scale confirms a new number.
#
# Steps 4/5 (added after the first real run of this script, 2026-09-15):
# pair_classification.py can't resolve physical linkage for a rescue-pass
# (GENOME_ONLY) presence call using family_positions.tsv alone -- that file
# only has GFF3-annotated positions, and a GENOME_ONLY call is a raw tblastn
# genomic hit with no annotated protein. Verified on the real first run:
# insufficient_data jumped to 96.2% of FDR-significant pairs before this fix,
# down to 0.3% after. extract_rescue_positions.py + build_family_positions.py
# --rescue_positions close that gap by deriving a genomic position directly
# from each GENOME_ONLY call's own tblastn hit coordinates -- see both
# scripts' module docstrings for the full detail.

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

TBLASTN_ARGS=()
for f in "$RESC/per_strain_chunks"/*.tblastn.tsv.zst; do
    TBLASTN_ARGS+=(--tblastn_tsv "$f")
done

echo "Step 1/7: folding rescue hits into the presence matrix..."
pixi run python3 "$STUDY/bin/rescue_pass.py" \
    --matrix "$RUN/presence_matrix.tsv" \
    "${TBLASTN_ARGS[@]}" \
    --output "$RUN/presence_matrix.rescued.tsv"

echo "Step 2/7: frequency binning on the rescued matrix..."
pixi run python3 "$STUDY/bin/frequency_bins.py" \
    --matrix "$RUN/presence_matrix.rescued.tsv" \
    --config "$STUDY/config.csv" \
    --inventory "$RUN/strain_inventory.tsv" \
    --output "$RUN/frequency_table.rescued.tsv"

echo "Step 3/7: co-occurrence (exact stratified test, not Monte Carlo)..."
pixi run python3 "$STUDY/bin/cooccurrence.py" \
    --matrix "$RUN/presence_matrix.rescued.tsv" \
    --frequency_table "$RUN/frequency_table.rescued.tsv" \
    --config "$STUDY/config.csv" \
    --inventory "$RUN/strain_inventory.tsv" \
    --output "$RUN/cooccurring_pairs.rescued.tsv"

echo "Step 4/7: extracting genomic positions for GENOME_ONLY (rescue) calls..."
pixi run python3 "$STUDY/bin/extract_rescue_positions.py" \
    --matrix "$RUN/presence_matrix.rescued.tsv" \
    "${TBLASTN_ARGS[@]}" \
    --processes 8 \
    --output "$RUN/rescue_positions.tsv"

echo "Step 5/7: rebuilding family_positions.tsv with rescue positions merged in..."
pixi run python3 "$STUDY/bin/build_family_positions.py" \
    --gene_positions "$RUN/gene_positions.tsv" \
    --cluster_tsv "$RUN/tier1_cluster.tsv" \
    --rescue_positions "$RUN/rescue_positions.tsv" \
    --output "$RUN/family_positions.rescued.tsv"
rm -f "$RUN/rescue_positions.tsv"  # fully consumed above, regenerable from the tblastn chunks

echo "Step 6/7: pair classification (physical linkage / Starship mechanism)..."
pixi run python3 "$STUDY/bin/pair_classification.py" \
    --cooccurring_pairs "$RUN/cooccurring_pairs.rescued.tsv" \
    --family_positions "$RUN/family_positions.rescued.tsv" \
    --cluster_tsv "$RUN/tier1_cluster.tsv" \
    --captain_tblout "$STUDY/results/captain_gene/DUF3435_vs_study.tblout" \
    --output "$RUN/pair_classification.rescued.tsv"

echo "Step 7/7: regenerating summary figures from the corrected data..."
pixi run python3 "$STUDY/bin/plot_pangenome_summary.py" \
    --frequency_table "$RUN/frequency_table.rescued.tsv" \
    --matrix "$RUN/presence_matrix.rescued.tsv" \
    --pair_classification "$RUN/pair_classification.rescued.tsv" \
    --out_dir "$STUDY/results/figures"

echo "Done. Corrected outputs:"
echo "  $RUN/presence_matrix.rescued.tsv"
echo "  $RUN/frequency_table.rescued.tsv"
echo "  $RUN/cooccurring_pairs.rescued.tsv"
echo "  $RUN/family_positions.rescued.tsv"
echo "  $RUN/pair_classification.rescued.tsv"
echo "  $STUDY/results/figures/"
