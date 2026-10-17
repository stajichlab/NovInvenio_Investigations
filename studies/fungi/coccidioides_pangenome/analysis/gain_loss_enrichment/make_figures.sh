#!/usr/bin/bash
# Slide figures for the three gain/loss Pfam enrichment runs (primary scope: shell + soft-core).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
PY=$NII/.pixi/envs/default/bin/python3
F=$NII/bin/pangenome_gain_loss_enrichment_figures.py
AF=$NII/studies/fungi/Afumigatus_pangenome/analysis
CP=$NII/studies/fungi/coccidioides_pangenome/analysis
$PY $F --name "A. fumigatus vs A. lentulus + A. fischeri" \
    --enrichment-dir $AF/gain_loss_enrichment/full_v072/shell_soft_core \
    --gain-loss-dir $AF/gain_loss/full_v072 --out-dir $AF/gain_loss_enrichment/figures --stem full_v072_pfam_top
$PY $F --name "C. immitis (in) vs C. posadasii (out)" \
    --enrichment-dir $CP/gain_loss_enrichment/immitis_in_posadasii_out_v5/shell_soft_core \
    --gain-loss-dir $CP/gain_loss/immitis_in_posadasii_out_v5 --out-dir $CP/gain_loss_enrichment/figures \
    --stem immitis_in_posadasii_out_v5_pfam_top
$PY $F --name "C. posadasii (in) vs C. immitis RS, WA_211 (out)" \
    --enrichment-dir $CP/gain_loss_enrichment/posadasii_in_immitis_out_v2/shell_soft_core \
    --gain-loss-dir $CP/gain_loss/posadasii_in_immitis_out_v2 --out-dir $CP/gain_loss_enrichment/figures \
    --stem posadasii_in_immitis_out_v2_pfam_top
