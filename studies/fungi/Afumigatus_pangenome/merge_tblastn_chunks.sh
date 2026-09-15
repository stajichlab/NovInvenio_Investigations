#!/usr/bin/bash
# Merge run_rescue_pass_tblastn_chunked.sh's per-chunk .zst outputs into one
# compressed tblastn_all_vs_all.tsv.zst -- run after all array tasks finish
# (`squeue` shows none of them running/pending). Decompression happens
# streaming (zstd -dc), never materializing the uncompressed concatenation
# on shared storage.
set -euo pipefail

STUDY="/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome"
RESC="$STUDY/results/rescue_pass"

shopt -s nullglob
chunks=("$RESC"/chunks/tblastn_chunk_*.tsv.zst)
if [[ ${#chunks[@]} -eq 0 ]]; then
    echo "ERROR: no chunk outputs found in $RESC/chunks/" >&2
    exit 1
fi

echo "Merging ${#chunks[@]} chunks -> $RESC/tblastn_all_vs_all.tsv.zst"
pixi run zstd -dc "${chunks[@]}" | pixi run zstd -T0 -o "$RESC/tblastn_all_vs_all.tsv.zst"
echo "Done."
