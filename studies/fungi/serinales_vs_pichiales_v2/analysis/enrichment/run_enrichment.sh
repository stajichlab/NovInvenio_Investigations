#!/usr/bin/bash
#SBATCH -p stajichlab -n 8 --mem 32G -t 01:30:00 -J serin2_enrich --constraint=milan
#SBATCH --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/serinales_vs_pichiales_v2/analysis/enrichment/run_enrichment.%j.log
# serinales_vs_pichiales_v2 primary enrichment (novelty + loss), 2026-10-16.
# Same settings as studies/enrichment_overview/run_enrichment_serinales.sh (v1 novelty) and
# run_enrichment_loss_all.sh (v1 loss): proteome universe, clustered background,
# ingroup_min_frac 0.75 (novelty); --direction loss with 0.75 / --other-max-frac 0.0 (loss).
# Parameters confirmed from the v2 nextflow BUILD_PRESENCE_MATRIX .command.sh
# (.nf_launch/fungi/serinales_vs_pichiales_v2/work/31/c3e269..., work/41/97c071...).
# Annotation: all 26 species are UniProt reference proteomes (annotations/ hard-linked from v1).
# --constraint=milan: mmseqs needs AVX2 (c01-c03 have none).
# Absolute paths only (no BASH_SOURCE; see ~/.claude/CLAUDE.md).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
TOOL=$NII/studies/fungi/zoosporic_dikarya/analysis/enrichment
PY="$NII/.pixi/envs/default/bin/python3 -I"
THREADS="${SLURM_NTASKS:-8}"
export SCRATCH="${SCRATCH:?SCRATCH must be set (run under SLURM)}"
NAME=serinales_vs_pichiales_v2
S=$NII/studies/fungi/$NAME
E=$S/analysis/enrichment
R=$NII/results/$NAME
PFAMDAT=/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat
IPR=/srv/projects/db/interproscan/6.0.0/110.0/interpro/110.0/entries.json
RUNLABEL="diamond pairwise, nf_NovInvenio main d03f5a2 run 2026-10-13"
NIN=$(awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["GROUP"]=="IN"' "$S/config.csv" | wc -l)
NOUT=$(awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["GROUP"]=="OUT"' "$S/config.csv" | wc -l)
mkdir -p "$E/outputs/loss" "$E/figures"

echo "[$(date +%T)] novelty start (IN $NIN, OUT $NOUT, candidates $(grep -c . "$R/candidates.txt"))"
T0=$(date +%s)
$PY "$TOOL/enrich_novelties_v2.py" --results "$R" --config "$S/config.csv" \
    --annotations "$S/annotations" --study "$NAME" --cache-dir "$E/outputs/cache" --threads "$THREADS" \
    --pfam-dat "$PFAMDAT" --interpro-entries "$IPR" \
    --controls "$TOOL/actin_controls.tsv" --ingroup-min-frac 0.75 \
    --universe proteome --cluster-background \
    --out-dir "$E/outputs/proteome_clustered" > "$E/outputs/run_proteome_clustered.log" 2>&1
echo "[$(date +%T)] novelty enrichment done ($(( $(date +%s) - T0 )) s)"
$PY "$TOOL/make_summary_figure.py" --enrich "$E/outputs/proteome_clustered" --config "$S/config.csv" \
    --title "Serinales v2 ($NIN IN) vs Pichiales ($NOUT OUT): enrichment of novelty candidates ($RUNLABEL)" \
    --footnote "Annotation: UniProt reference proteomes (Pfam, InterPro, GO); IN species without annotation contribute no annotated units. Universe: all proteins of the IN proteomes. Family units: candidate families from the pipeline; background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Labels: k = candidate units with term / all units with term." \
    --out "$E/figures/${NAME}_summary" > "$E/outputs/figure.log" 2>&1

echo "[$(date +%T)] loss start"
T0=$(date +%s)
NCAND=$(grep -c . "$R/loss_candidates.txt")
$PY "$TOOL/enrich_novelties_v2.py" --results "$R" --config "$S/config.csv" \
    --annotations "$S/annotations" --study "$NAME" --cache-dir "$E/outputs/cache" --threads "$THREADS" \
    --pfam-dat "$PFAMDAT" --interpro-entries "$IPR" \
    --controls "$TOOL/actin_controls.tsv" \
    --direction loss --ingroup-min-frac 0.75 --other-max-frac 0.0 \
    --universe proteome --cluster-background \
    --out-dir "$E/outputs/loss/proteome_clustered" > "$E/outputs/loss/run_proteome_clustered.log" 2>&1
T1=$(date +%s)
FRAC=$($PY -c "import csv,sys; r=[x for x in csv.DictReader(open(sys.argv[1]),delimiter='\t') if x['species']=='ALL_OUT' and x['subset']=='candidates']; print(r[0]['frac_any'])" "$E/outputs/loss/proteome_clustered/coverage.tsv")
LOW=$($PY -c "import sys; print('true' if float(sys.argv[1]) < 0.5 else 'false')" "$FRAC")
printf 'candidate_annotation:\n  frac_loss_candidates_any_annotation: %s\n  low_coverage: %s\n  rule: "low_coverage = fraction < 0.5 (flag only; no run was excluded on this value)"\n  wall_s_enrichment: %s\n' \
    "$FRAC" "$LOW" "$((T1 - T0))" >> "$E/outputs/loss/proteome_clustered/PROVENANCE.yaml"
LOWNOTE=""
if [ "$LOW" = true ]; then
  LOWNOTE=" LOW COVERAGE: only $FRAC of loss candidates carry any UniProt annotation; read results with care."
fi
echo "[$(date +%T)] loss enrichment done ($((T1 - T0)) s)"
$PY "$TOOL/make_summary_figure.py" --enrich "$E/outputs/loss/proteome_clustered" --config "$S/config.csv" \
    --group OUT --candidate-label "loss candidate proteins" \
    --title "Pichiales ($NOUT OUT) vs Serinales v2 ($NIN IN): enrichment of loss-direction candidates (present in >=0.75 of OUT, absent from all IN; $RUNLABEL)" \
    --footnote "Annotation: UniProt reference proteomes (Pfam, InterPro, GO). Universe: all proteins of the $NOUT OUT proteomes. Family units: loss-candidate families from the pipeline (mmseqs on loss_candidates.fa); background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Loss candidates: $NCAND; annotated fraction $FRAC.$LOWNOTE" \
    --out "$E/figures/${NAME}_loss_summary" > "$E/outputs/loss/figure.log" 2>&1
echo "[$(date +%T)] all done"
