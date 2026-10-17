#!/bin/bash -l
#SBATCH -p epyc --constraint=milan -c 2 --mem 8G -t 1:00:00
#SBATCH -J af_starship_island
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/analysis/starship_island_overlap/logs/overlap.%j.out
# Starship (Table S7, Af293) and DUF3435 captain families vs pangenome islands, full_v072.
# Two variants: (1) the run's own gene_positions (default; 8,686 of 9,647 Af293 members
# placed), (2) plus coordinates for the unplaced Af293 members (Swiss-Prot entries)
# recovered by scripts/af293_missing_coords.py -> full_v072/with_unplaced_af293_coords/.
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
AF=$NII/studies/fungi/Afumigatus_pangenome
D=$AF/analysis/starship_island_overlap
P=$AF/results/full_v072/config/pangenome
PY=$NII/.pixi/envs/default/bin/python3
EXTRA=$D/inputs/af293_unplaced_swissprot_coords.tsv
for V in default with_unplaced_af293_coords; do
  OUT=$D/full_v072; XARG=()
  if [ "$V" != default ]; then OUT=$D/full_v072/$V; XARG=(--extra-gene-positions "$EXTRA"); fi
  "$PY" "$NII/bin/starship_island_overlap.py" --run-dir "$P" \
    --starship-annotation "$D/inputs/af293_starships_table_s7.elements.tsv" \
    --captain-tblout "$AF/results/captain_gene/DUF3435_vs_study.tblout" \
    --gain-loss "$AF/analysis/gain_loss/full_v072/family_gain_loss.tsv.gz" \
    --run-label "full_v072${V/default/}" --out-dir "$OUT" "${XARG[@]}"
  "$PY" "$D/scripts/sm_backbone_starship.py" --run-dir "$P" \
    --starship-annotation "$D/inputs/af293_starships_table_s7.elements.tsv" \
    --sm-tblout "$AF/results/sm_backbone/SM_backbone_vs_study.tblout" \
    --captain-tblout "$AF/results/captain_gene/DUF3435_vs_study.tblout" \
    --out-dir "$OUT" "${XARG[@]}"
done
