#!/bin/bash -l
#SBATCH -p epyc --constraint=milan -c 2 --mem 8G -t 1:00:00
#SBATCH -J cocci_starship_island
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/starship_island_overlap/logs/overlap.%j.out
# Coccidioides: no Starship element calls exist. Captain (DUF3435) families vs islands,
# and antiSMASH SM regions vs captain genes, for both runs.
# Needs scripts/run_captain_antismash.sh outputs in inputs/.
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
CP=$NII/studies/fungi/coccidioides_pangenome
D=$CP/analysis/starship_island_overlap
PY=$NII/.pixi/envs/default/bin/python3
declare -A RUNDIR=(
  [immitis_in_posadasii_out_v5]=$CP/results/immitis_in_posadasii_out_v5/config_immitis_in_posadasii_out/pangenome
  [posadasii_in_immitis_out_v2]=$CP/results/posadasii_in_immitis_out_v2/config_posadasii_in_RS_WA211_out/pangenome
)
declare -A CFG=(
  [immitis_in_posadasii_out_v5]=$CP/config_immitis_in_posadasii_out.csv
  [posadasii_in_immitis_out_v2]=$CP/config_posadasii_in_RS_WA211_out.csv
)
for RUN in immitis_in_posadasii_out_v5 posadasii_in_immitis_out_v2; do
  "$PY" "$NII/bin/starship_island_overlap.py" \
    --run-dir "${RUNDIR[$RUN]}" \
    --captain-tblout "$D/inputs/DUF3435_vs_study.tblout.zst" \
    --gain-loss "$CP/analysis/gain_loss/$RUN/family_gain_loss.tsv.gz" \
    --run-label "$RUN" --out-dir "$D/$RUN"
  "$PY" "$D/scripts/sm_captain_overlap.py" \
    --antismash-dir "$D/inputs/antismash" \
    --captain-tblout "$D/inputs/DUF3435_vs_antismash_annotation.tblout.zst" \
    --config "${CFG[$RUN]}" --out-dir "$D/$RUN"
done
