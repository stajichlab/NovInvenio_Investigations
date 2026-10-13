#!/usr/bin/bash
#SBATCH -p stajichlab -n 8 --mem 32G -t 04:00:00 -J nii_enrich_lowann --constraint=milan
#SBATCH --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/enrichment_overview/run_enrichment_lowann.%j.log
# Primary enrichment run (same settings as studies/fungi/zoosporic_dikarya/analysis/enrichment/
# run_enrichment.sh, primary run only) for the studies that pass the coverage rule in
# coverage_all_studies.tsv (>=2 IN species annotated and >=50% candidates annotated).
# Absolute paths only (no BASH_SOURCE; see ~/.claude/CLAUDE.md).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
TOOL=$NII/studies/fungi/zoosporic_dikarya/analysis/enrichment
PY=$NII/.pixi/envs/default/bin/python3
THREADS="${SLURM_CPUS_ON_NODE:-${SLURM_NTASKS:-4}}"
export SCRATCH="${SCRATCH:?SCRATCH must be set (run under SLURM)}"

# study|ingroup_min_frac|title|run label
STUDIES=(
  "fungi/malassezia_vs_ustilaginomycotina|0.75|Malassezia (10 IN) vs Ustilaginomycotina|diamond pairwise, nf_NovInvenio run 2026-10-11"
  "fungi/metschnikowiaceae_cah|0.75|Metschnikowiaceae incl. C. auris (9 IN) vs OUT|diamond pairwise, nf_NovInvenio run 2026-10-11"
  "fungi/sordariales_shallow|0.75|Sordariales shallow (4 IN) vs OUT|diamond pairwise, nf_NovInvenio run"
)
for entry in "${STUDIES[@]}"; do
  IFS='|' read -r ST FRAC TITLE RUNLABEL <<< "$entry"
  NAME=$(basename "$ST")
  S=$NII/studies/$ST
  E=$S/analysis/enrichment
  O=$E/outputs
  mkdir -p "$O" "$E/figures"
  echo "[$(date +%T)] $NAME start"
  $PY "$TOOL/enrich_novelties_v2.py" --results "$NII/results/$NAME" --config "$S/config.csv" \
      --annotations "$S/annotations" --study "$NAME" --cache-dir "$O/cache" --threads "$THREADS" \
      --pfam-dat /srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat \
      --interpro-entries /srv/projects/db/interproscan/6.0.0/110.0/interpro/110.0/entries.json \
      --controls "$TOOL/actin_controls.tsv" --ingroup-min-frac "$FRAC" \
      --universe proteome --cluster-background \
      --out-dir "$O/proteome_clustered" > "$O/run_proteome_clustered.log" 2>&1
  $PY "$TOOL/make_summary_figure.py" --enrich "$O/proteome_clustered" --config "$S/config.csv" \
      --title "$TITLE: enrichment of novelty candidates ($RUNLABEL)" \
      --footnote "Annotation: UniProt reference proteomes (Pfam, InterPro, GO); IN species without annotation contribute no annotated units. Universe: all proteins of the IN proteomes. Family units: candidate families from the pipeline; background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Labels: k = candidate units with term / all units with term." \
      --out "$E/figures/${NAME}_summary" > "$O/figure.log" 2>&1
  echo "[$(date +%T)] $NAME done"
done
