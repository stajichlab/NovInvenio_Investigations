#!/usr/bin/bash
#SBATCH -p stajichlab -n 8 --mem 32G -t 04:00:00 -J nii_enrich_zoholo
#SBATCH --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/eukaryotes/zoosporic_holozoa_vs_dikarya/analysis/enrichment/run_enrichment.%j.log
# Primary enrichment run, same settings as studies/fungi/zoosporic_dikarya/analysis/enrichment.
# Absolute paths only (no BASH_SOURCE; see ~/.claude/CLAUDE.md).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
TOOL=$NII/studies/fungi/zoosporic_dikarya/analysis/enrichment
PY=$NII/.pixi/envs/default/bin/python3
NAME=zoosporic_holozoa_vs_dikarya
S=$NII/studies/eukaryotes/$NAME
E=$S/analysis/enrichment
O=$E/outputs
export SCRATCH="${SCRATCH:?SCRATCH must be set (run under SLURM)}"
mkdir -p "$O" "$E/figures"
$PY "$TOOL/enrich_novelties_v2.py" --results "$NII/results/$NAME" --config "$S/config.csv" \
    --annotations "$S/annotations" --study "$NAME" --cache-dir "$O/cache" --threads "${SLURM_NTASKS:-8}" \
    --pfam-dat /srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat \
    --interpro-entries /srv/projects/db/interproscan/6.0.0/110.0/interpro/110.0/entries.json \
    --controls "$TOOL/actin_controls.tsv" --ingroup-min-frac 0.75 \
    --universe proteome --cluster-background \
    --out-dir "$O/proteome_clustered" > "$O/run_proteome_clustered.log" 2>&1
$PY "$TOOL/make_summary_figure.py" --enrich "$O/proteome_clustered" --config "$S/config.csv" \
    --title "Zoosporic fungi + Holozoa (25 IN) vs Dikarya (7 OUT): enrichment of novelty candidates" \
    --footnote "Annotation: UniProt reference proteomes (Pfam, InterPro, GO). Universe: all proteins of the IN proteomes. Family units: candidate families from the pipeline; background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units." \
    --out "$E/figures/${NAME}_summary" > "$O/figure.log" 2>&1
