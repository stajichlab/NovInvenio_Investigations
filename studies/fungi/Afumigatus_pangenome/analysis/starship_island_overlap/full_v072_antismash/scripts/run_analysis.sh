#!/bin/bash
# Parse the Af293 antiSMASH JSON and compute SM cluster / Starship / island overlap (login node, ~15 s).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
AF=$NII/studies/fungi/Afumigatus_pangenome
D=$AF/analysis/starship_island_overlap
O=$D/full_v072_antismash
PY=$NII/.pixi/envs/default/bin/python3
"$PY" "$O/scripts/antismash_sm_starship.py" \
  --antismash-json "$O/antismash_run/Af293.json.gz" \
  --gff3 "$AF/data_dir/gff3/UP000002530_330879.gff3" \
  --pep "$AF/data_dir/pep/UP000002530_330879.pep.fa" \
  --genome "$AF/data_dir/dna/UP000002530_330879.dna.fa" \
  --run-dir "$AF/results/full_v072/config/pangenome" \
  --extra-gene-positions "$D/inputs/af293_unplaced_swissprot_coords.tsv" \
  --starship-annotation "$D/inputs/af293_starships_table_s7.elements.tsv" \
  --captain-tblout "$AF/results/captain_gene/DUF3435_vs_study.tblout" \
  --gain-loss "$AF/analysis/gain_loss/full_v072/family_gain_loss.tsv.gz" \
  --proxy-genes "$D/full_v072/with_unplaced_af293_coords/sm_backbone_af293_genes.tsv.gz" \
  --out-dir "$O"
"$PY" "$O/scripts/make_sm_figure.py" --dir "$O" --out-prefix "$O/figures/af293_antismash_sm_starship"
