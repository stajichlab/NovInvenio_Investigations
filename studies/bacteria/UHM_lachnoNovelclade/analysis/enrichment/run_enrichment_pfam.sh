#!/usr/bin/bash
#SBATCH -p stajichlab -n 8 --mem 32G -t 02:00:00 -J lachno_enrich_pfam --constraint=milan
#SBATCH --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/bacteria/UHM_lachnoNovelclade/analysis/enrichment/run_enrichment_pfam.%j.log
# UHM_lachnoNovelclade Pfam-only enrichment (novelty + loss), 2026-10-14.
# Annotation: Pfam-A hmmsearch for all 9 proteomes (bin/annotate_proteomes_pfam.py);
# GO and InterPro are empty for every species, so only the Pfam and Pfam-clan tests have data.
# Settings: as studies/fungi/agaricomycetes_pairwise/analysis/enrichment/run_enrichment_pfam.sh.
#   ingroup_min_frac 0.75, outgroup_min_frac 0.75, loss_ingroup_max_frac 0.0: from the
#     nf_NovInvenio run report results/UHM_lachnoNovelclade/nextflow_log/
#     UHM_lachnoNovelclade-20260911_153453-report.html.
#   --presence-matrix results/UHM_lachnoNovelclade/presence_matrix.tsv (novelty): no
#     presence_matrix.uniprot.tsv in this run; the matrix is only used for counts.
# --config analysis/enrichment/config.fasta_stems.csv: config.csv Short differs from the
#   FASTA stem (= matrix column, species.csv Short) for 6 of 9 species. The derived copy is
#   made by studies/bacteria/UHM_lachnoNovelclade/bin/make_fasta_stem_config.py (sidecar
#   config.fasta_stems.csv.provenance.yaml); this script reruns it. --species-csv and
#   --proteins-dir are given explicitly because the derived config is not in the study dir.
# --constraint=milan: mmseqs needs AVX2 (c01-c03 have none).
# Absolute paths only (no BASH_SOURCE; see ~/.claude/CLAUDE.md).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
TOOL=$NII/studies/fungi/zoosporic_dikarya/analysis/enrichment
PY="$NII/.pixi/envs/default/bin/python3 -I"
THREADS="${SLURM_NTASKS:-8}"
export SCRATCH="${SCRATCH:?SCRATCH must be set (run under SLURM)}"
NAME=UHM_lachnoNovelclade
S=$NII/studies/bacteria/$NAME
E=$S/analysis/enrichment
R=$NII/results/$NAME
CFG=$E/config.fasta_stems.csv
PFAMDAT=/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat
IPR=/srv/projects/db/interproscan/6.0.0/110.0/interpro/110.0/entries.json
RUNLABEL="diamond pairwise, NovInvenio run 2026-09-11"
ANNOT="Annotation: Pfam-A hmmsearch (HMMER 3.4, --cut_ga, Pfam 38.2) of all 9 proteomes (prodigal MAG proteins and RefSeq proteins); GO and InterPro not available."
COMMON=(--config "$CFG" --species-csv "$S/species.csv" --proteins-dir "$S/data_dir/pep"
        --annotations "$S/annotations" --study "$NAME" --cache-dir "$E/outputs/cache" --threads "$THREADS"
        --pfam-dat "$PFAMDAT" --interpro-entries "$IPR" --controls "$TOOL/actin_controls.tsv"
        --universe proteome --cluster-background)
mkdir -p "$E/outputs/loss" "$E/figures"

$PY "$S/bin/make_fasta_stem_config.py" --out "$CFG"

echo "[$(date +%T)] novelty start"
T0=$(date +%s)
$PY "$TOOL/enrich_novelties_v2.py" --results "$R" "${COMMON[@]}" \
    --presence-matrix "$R/presence_matrix.tsv" --ingroup-min-frac 0.75 \
    --out-dir "$E/outputs/proteome_clustered" > "$E/outputs/run_proteome_clustered.log" 2>&1
T1=$(date +%s)
printf 'wall_s_enrichment: %s\n' "$((T1 - T0))" >> "$E/outputs/proteome_clustered/PROVENANCE.yaml"
echo "[$(date +%T)] novelty enrichment done ($((T1 - T0)) s)"
$PY "$TOOL/make_summary_figure.py" --enrich "$E/outputs/proteome_clustered" --config "$CFG" \
    --annot-label "Pfam (hmmsearch)" \
    --title "Lachnospiraceae IN set (5 IN: 1 UHM MAG + 4 RefSeq) vs UHM MAG clade (4 OUT): enrichment of novelty candidates ($RUNLABEL)" \
    --footnote "$ANNOT Universe: all proteins of the IN proteomes. Family units: candidate families from the pipeline; background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Labels: k = candidate units with term / all units with term." \
    --out "$E/figures/${NAME}_summary" > "$E/outputs/figure.log" 2>&1

echo "[$(date +%T)] loss start"
T0=$(date +%s)
NCAND=$(grep -c . "$R/loss_candidates.txt")
$PY "$TOOL/enrich_novelties_v2.py" --results "$R" "${COMMON[@]}" \
    --direction loss --ingroup-min-frac 0.75 --other-max-frac 0.0 \
    --out-dir "$E/outputs/loss/proteome_clustered" > "$E/outputs/loss/run_proteome_clustered.log" 2>&1
T1=$(date +%s)
FRAC=$($PY -c "import csv,sys; r=[x for x in csv.DictReader(open(sys.argv[1]),delimiter='\t') if x['species']=='ALL_OUT' and x['subset']=='candidates']; print(r[0]['frac_any'])" "$E/outputs/loss/proteome_clustered/coverage.tsv")
LOW=$($PY -c "import sys; print('true' if float(sys.argv[1]) < 0.5 else 'false')" "$FRAC")
printf 'candidate_annotation:\n  frac_loss_candidates_any_annotation: %s\n  low_coverage: %s\n  rule: "low_coverage = fraction < 0.5 (flag only; no run was excluded on this value)"\n  wall_s_enrichment: %s\n' \
    "$FRAC" "$LOW" "$((T1 - T0))" >> "$E/outputs/loss/proteome_clustered/PROVENANCE.yaml"
LOWNOTE=""
if [ "$LOW" = true ]; then
  LOWNOTE=" LOW COVERAGE: only $FRAC of loss candidates carry a Pfam annotation; read results with care."
fi
echo "[$(date +%T)] loss enrichment done ($((T1 - T0)) s)"
$PY "$TOOL/make_summary_figure.py" --enrich "$E/outputs/loss/proteome_clustered" --config "$CFG" \
    --group OUT --candidate-label "loss candidate proteins" --annot-label "Pfam (hmmsearch)" \
    --title "UHM MAG clade (4 OUT) vs Lachnospiraceae IN set (5 IN): enrichment of loss-direction candidates (present in >=0.75 of OUT, absent from all IN; $RUNLABEL)" \
    --footnote "$ANNOT Universe: all proteins of the 4 OUT proteomes. Family units: loss-candidate families from the pipeline (mmseqs on loss_candidates.fa); background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Loss candidates: $NCAND; annotated fraction $FRAC.$LOWNOTE" \
    --out "$E/figures/${NAME}_loss_summary" > "$E/outputs/loss/figure.log" 2>&1
echo "[$(date +%T)] all done"
