#!/usr/bin/env bash
# Analyse the diamond-sensitivity benchmark once all 9 runs are finished.
# Writes notes/diamond-sensitivity/<clade>.{cost,concordance,lost_evidence}.tsv
# (class 2) and a gitignored per-candidate lost-detail file per clade (class 3),
# plus controls scoring per mode. See notes/diamond-sensitivity/README.md.
set -euo pipefail

NII_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations"
PIPE_ROOT="/bigdata/stajichlab/jstajich/projects/NovInvenio-worktrees/bench-diamond-sensitivity"
OUT="$NII_ROOT/notes/diamond-sensitivity"
DETAIL="$NII_ROOT/results/diamond_sensitivity_bench"
mkdir -p "$OUT" "$DETAIL"
cd "$NII_ROOT"

# clade label | config study | Tier C+H results dir | controls csv | busco map
# sordariales_shallow.controls.csv is not on nf_NovInvenio origin yet (local commit
# eae0c97 only), so its NII copy under legacy/ is used.
CLADES=(
  "pezizo_set1|pezizo_set1|results/pezizo_set1_cluster|$PIPE_ROOT/configs/controls/pezizo_set1.controls.csv|studies/fungi/pezizo_set1_cluster/tier_comparison/pezizo_set1.busco_map.tsv"
  "agaricomycetes|agaricomycetes_pairwise|results/agaricomycetes_mmseqs|$PIPE_ROOT/configs/controls/Agaricales.controls.csv|"
  "sordariales_shallow|sordariales_shallow|results/sordariales_shallow_cluster|legacy/novinvenio_configs/controls/sordariales_shallow.controls.csv|studies/fungi/sordariales_shallow/busco_map.tsv"
)

for spec in "${CLADES[@]}"; do
    IFS='|' read -r label cstudy chdir controls busco <<< "$spec"
    cfg="studies/fungi/$cstudy/config.csv"
    python3 bin/diamond_sensitivity_report.py --label "$label" --config "$cfg" \
        --mode default "results/${label}_dmnd_default" \
        --mode sensitive "results/${label}_dmnd_sensitive" \
        --mode very_sensitive "results/${label}_dmnd_very_sensitive" \
        --ch-dir "$chdir" \
        --output-prefix "$OUT/$label" \
        --output-lost-detail "$DETAIL/$label.lost_detail.tsv"
    for m in default sensitive very_sensitive; do
        bm=(); [ -n "$busco" ] && bm=(--busco-map "$busco")
        pixi run -q --manifest-path "$PIPE_ROOT/pixi.toml" python "$PIPE_ROOT/bin/score_controls.py" \
            --controls "$controls" \
            --matrix "results/${label}_dmnd_${m}/presence_matrix.tsv" \
            --config "$cfg" "${bm[@]}" \
            --output "$DETAIL/$label.$m.controls_scored.tsv" \
            --summary "$OUT/$label.$m.controls_summary.tsv" 2>&1 | grep -v WARN | tail -2
    done
done
