#!/usr/bin/bash
#SBATCH -p stajichlab -n 8 --mem 32G -t 04:00:00 -J nii_enrich_zoholo_loss
#SBATCH --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/eukaryotes/zoosporic_holozoa_vs_dikarya/analysis/enrichment/run_enrichment_loss.%j.log
# Loss-direction enrichment (genes present in the Dikarya OUT group, absent from the 25 IN
# species). Same settings as the novelty run in run_enrichment.sh, with --direction loss.
#   1. loss_near_miss/loss_loose.py: rebuild loss candidates from the search cache
#      (K = 0 must equal loss_candidates.txt) and the looser K sets with mmseqs families.
#   2. primary: pipeline loss_candidates.txt (332) -> outputs/loss/proteome_clustered
#   3. loose:   K = 2 set (<= 2 IN species with a hit) -> outputs/loss/loss_loose_max2
#   4. summary figure for the primary run.
# Absolute paths only (no BASH_SOURCE; see ~/.claude/CLAUDE.md).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
TOOL=$NII/studies/fungi/zoosporic_dikarya/analysis/enrichment
PY="$NII/.pixi/envs/default/bin/python3 -I"
NAME=zoosporic_holozoa_vs_dikarya
S=$NII/studies/eukaryotes/$NAME
R=$NII/results/$NAME
E=$S/analysis/enrichment
O=$E/outputs/loss
NM=$S/analysis/loss_near_miss
export SCRATCH="${SCRATCH:?SCRATCH must be set (run under SLURM)}"
THREADS="${SLURM_NTASKS:-8}"
mkdir -p "$O" "$E/figures" "$NM/outputs"

$PY "$NM/loss_loose.py" --results "$R" --config "$S/config.csv" --out-dir "$NM/outputs" \
    --outgroup-min-frac 0.75 --ks 0,1,2,3 --in-max-frac 0.15 --threads "$THREADS" \
    > "$NM/outputs/loss_loose.log" 2>&1 || {
    echo "loss_loose.py: K=0 does not reproduce loss_candidates.txt (or it failed); see $NM/outputs/loss_loose.log" >&2
    RUN_LOOSE=0
}
RUN_LOOSE=${RUN_LOOSE:-1}

COMMON=(--results "$R" --config "$S/config.csv" --annotations "$S/annotations" --study "$NAME"
        --cache-dir "$E/outputs/cache" --threads "$THREADS"
        --pfam-dat /srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat
        --interpro-entries /srv/projects/db/interproscan/6.0.0/110.0/interpro/110.0/entries.json
        --controls "$TOOL/actin_controls.tsv"
        --direction loss --ingroup-min-frac 0.75 --other-max-frac 0.0
        --universe proteome --cluster-background)

$PY "$TOOL/enrich_novelties_v2.py" "${COMMON[@]}" --out-dir "$O/proteome_clustered" \
    > "$O/run_proteome_clustered.log" 2>&1

if [ "$RUN_LOOSE" = 1 ]; then
    # --other-max-frac 0.0 here only labels the matrix rule check; the K=2 candidates have
    # up to 2 IN hits by design, so the check reports them as differences (expected).
    $PY "$TOOL/enrich_novelties_v2.py" "${COMMON[@]}" \
        --candidates "$NM/outputs/candidates_K2.txt" --families "$NM/outputs/families_K2.tsv" \
        --out-dir "$O/loss_loose_max2" > "$O/run_loss_loose_max2.log" 2>&1
fi

$PY "$TOOL/make_summary_figure.py" --enrich "$O/proteome_clustered" --config "$S/config.csv" \
    --group OUT --candidate-label "loss candidate proteins" \
    --title "Dikarya (7 OUT) vs zoosporic fungi + Holozoa (25 IN): enrichment of loss-direction candidates (present in >=0.75 of Dikarya, absent from all IN)" \
    --footnote "Annotation: UniProt reference proteomes (Pfam, InterPro, GO). Universe: all proteins of the 7 Dikarya proteomes. Family units: loss-candidate families from the pipeline (mmseqs on loss_candidates.fa); background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Small candidate set (332): low power." \
    --out "$E/figures/${NAME}_loss_summary" > "$O/figure.log" 2>&1
