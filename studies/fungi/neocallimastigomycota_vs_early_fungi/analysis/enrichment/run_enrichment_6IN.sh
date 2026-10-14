#!/usr/bin/bash
#SBATCH -p stajichlab -n 8 --mem 32G -t 04:00:00 -J neocall_enrich_6IN --constraint=milan
#SBATCH --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/neocallimastigomycota_vs_early_fungi/analysis/enrichment/outputs/run_enrichment_6IN.%j.log
# --constraint=milan: the pixi mmseqs binary needs AVX2 (nodes c01-c03 have none; SIGILL).
# Enrichment for neocallimastigomycota_vs_early_fungi after the 4 NCBI/JGI IN proteomes got
# Pfam-only annotation files (bin/hmm_domtbl_to_annotations.py; analysis/hmmscan/README.md).
# Settings are the primary-run settings of studies/enrichment_overview/run_enrichment_all.sh
# (novelty) and run_enrichment_loss_all.sh (loss): --universe proteome --cluster-background,
# --ingroup-min-frac 0.75 (pipeline ingroup_min_frac / outgroup_min_frac 0.75 in the
# nextflow report of results/neocallimastigomycota_vs_early_fungi), loss --other-max-frac 0.0.
# Steps:
#   1. novelty, baseline: only the 15 UniProt annotation files (2 of 6 IN species annotated),
#      same as before the Pfam files existed -> outputs/proteome_clustered_2annotated_IN/
#   2. loss: copy the existing outputs/loss/proteome_clustered/ (2026-10-12/13 run) to
#      outputs/loss/proteome_clustered_2annotated_IN/ (only if that copy does not exist yet)
#   3. novelty, all annotation files (6 of 6 IN species have Pfam; 2 have GO/InterPro)
#      -> outputs/proteome_clustered/ + figures/neocallimastigomycota_vs_early_fungi_summary.*
#   4. loss rerun -> outputs/loss/proteome_clustered/ + figures/..._loss_summary.*
# Absolute paths only (no BASH_SOURCE).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
NAME=neocallimastigomycota_vs_early_fungi
S=$NII/studies/fungi/$NAME
E=$S/analysis/enrichment
O=$E/outputs
TOOL=$NII/studies/fungi/zoosporic_dikarya/analysis/enrichment
PY="$NII/.pixi/envs/default/bin/python3 -I"
THREADS="${SLURM_NTASKS:-8}"
export SCRATCH="${SCRATCH:?SCRATCH must be set (run under SLURM)}"
COMMON=(--results "$NII/results/$NAME" --config "$S/config.csv" --study "$NAME" --threads "$THREADS"
        --pfam-dat /srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat
        --interpro-entries /srv/projects/db/interproscan/6.0.0/110.0/interpro/110.0/entries.json
        --controls "$TOOL/actin_controls.tsv" --ingroup-min-frac 0.75
        --universe proteome --cluster-background)
mkdir -p "$O/loss" "$E/figures"

# 1. novelty baseline with the UniProt files only
echo "[$(date +%T)] novelty baseline (UniProt annotation files only) start"
UNI="$SCRATCH/annotations_uniprot_only"
mkdir -p "$UNI"
for f in "$S"/annotations/UP*.tsv.gz; do ln -sf "$f" "$UNI/"; done
T0=$(date +%s)
$PY "$TOOL/enrich_novelties_v2.py" "${COMMON[@]}" --annotations "$UNI" --cache-dir "$O/cache" \
    --out-dir "$O/proteome_clustered_2annotated_IN" > "$O/run_proteome_clustered_2annotated_IN.log" 2>&1
printf 'note: "Baseline for comparison. Annotation dir = symlinks to the 15 UniProt files only (%s, removed after the job), so the 4 NCBI/JGI IN species have no annotation, as before 2026-10-14. Made by %s."\nwall_s_enrichment: %s\n' \
    "$UNI" "$E/run_enrichment_6IN.sh" "$(( $(date +%s) - T0 ))" >> "$O/proteome_clustered_2annotated_IN/PROVENANCE.yaml"
echo "[$(date +%T)] novelty baseline done ($(( $(date +%s) - T0 )) s)"

# 2. keep the existing loss output
if [ ! -e "$O/loss/proteome_clustered_2annotated_IN" ]; then
  cp -a "$O/loss/proteome_clustered" "$O/loss/proteome_clustered_2annotated_IN"
  cp -a "$O/loss/run_proteome_clustered.log" "$O/loss/run_proteome_clustered_2annotated_IN.log"
  mkdir -p "$E/figures/2annotated_IN"
  for ext in png svg pdf; do
    cp -a "$E/figures/${NAME}_loss_summary.$ext" "$E/figures/2annotated_IN/" 2>/dev/null || true
  done
  echo "[$(date +%T)] kept old loss output in $O/loss/proteome_clustered_2annotated_IN"
fi

# 3. novelty, all annotation files
echo "[$(date +%T)] novelty (6 IN with Pfam) start"
T0=$(date +%s)
$PY "$TOOL/enrich_novelties_v2.py" "${COMMON[@]}" --annotations "$S/annotations" --cache-dir "$O/cache" \
    --out-dir "$O/proteome_clustered" > "$O/run_proteome_clustered.log" 2>&1
T1=$(date +%s)
printf 'wall_s_enrichment: %s\n' "$((T1 - T0))" >> "$O/proteome_clustered/PROVENANCE.yaml"
NCAND=$(grep -c . "$NII/results/$NAME/candidates.txt")
$PY "$TOOL/make_summary_figure.py" --enrich "$O/proteome_clustered" --config "$S/config.csv" \
    --annot-label "GO/Pfam/InterPro (UniProt) or Pfam (HMMER)" \
    --title "Neocallimastigomycota (6 IN) vs early fungi (13 OUT): enrichment of novelty candidates (present in >=0.75 of IN, absent from all OUT; diamond pairwise, nf_NovInvenio main d03f5a2 run 2026-10-11)" \
    --footnote "Annotation: UniProt reference proteomes (Pfam, InterPro, GO) for Piromyces finnis and Anaeromyces robustus. Pfam only (Pfam-A 38.2, HMMER 3.4 hmmsearch or hmmscan, --cut_ga; no GO, no InterPro) for the 4 NCBI/JGI proteomes: Neocallimastix californiae, Neocallimastix lanati, Neocallimastix sp. 'constans', Piromyces sp. E2. Pfam and Pfam-clan tests cover 6 of 6 IN species; GO and InterPro tests cover 2 of 6. Universe: all proteins of the 6 IN proteomes. Family units: candidate families from the pipeline; background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Novelty candidates: $NCAND." \
    --out "$E/figures/${NAME}_summary" > "$O/figure.log" 2>&1
echo "[$(date +%T)] novelty done (enrichment $((T1 - T0)) s)"

# 4. loss rerun (settings and PROVENANCE addition as in run_enrichment_loss_all.sh)
echo "[$(date +%T)] loss start"
T0=$(date +%s)
$PY "$TOOL/enrich_novelties_v2.py" "${COMMON[@]}" --annotations "$S/annotations" --cache-dir "$O/cache" \
    --direction loss --other-max-frac 0.0 \
    --out-dir "$O/loss/proteome_clustered" > "$O/loss/run_proteome_clustered.log" 2>&1
T1=$(date +%s)
FRAC=$($PY -c "import csv,sys; r=[x for x in csv.DictReader(open(sys.argv[1]),delimiter='\t') if x['species']=='ALL_OUT' and x['subset']=='candidates']; print(r[0]['frac_any'])" "$O/loss/proteome_clustered/coverage.tsv")
LOW=$($PY -c "import sys; print('true' if float(sys.argv[1]) < 0.5 else 'false')" "$FRAC")
printf 'candidate_annotation:\n  frac_loss_candidates_any_annotation: %s\n  low_coverage: %s\n  rule: "low_coverage = fraction < 0.5 (flag only; no run was excluded on this value)"\n  wall_s_enrichment: %s\n' \
    "$FRAC" "$LOW" "$((T1 - T0))" >> "$O/loss/proteome_clustered/PROVENANCE.yaml"
NCAND=$(grep -c . "$NII/results/$NAME/loss_candidates.txt")
$PY "$TOOL/make_summary_figure.py" --enrich "$O/loss/proteome_clustered" --config "$S/config.csv" \
    --group OUT --candidate-label "loss candidate proteins" \
    --title "early fungi (13 OUT) vs Neocallimastigomycota (6 IN): enrichment of loss-direction candidates (present in >=0.75 of OUT, absent from all IN; diamond pairwise)" \
    --footnote "Annotation: UniProt reference proteomes (Pfam, InterPro, GO). Universe: all proteins of the 13 OUT proteomes. Family units: loss-candidate families from the pipeline (mmseqs on loss_candidates.fa); background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Loss candidates: $NCAND; annotated fraction $FRAC." \
    --out "$E/figures/${NAME}_loss_summary" > "$O/loss/figure.log" 2>&1
echo "[$(date +%T)] loss done (enrichment $((T1 - T0)) s)"
