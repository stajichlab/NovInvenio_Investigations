#!/usr/bin/bash
#SBATCH -p stajichlab -n 8 --mem 32G -t 02:00:00 -J pezizo_enrich_pfam --constraint=milan
#SBATCH --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/pezizo_set1/analysis/enrichment/run_enrichment_pfam.%j.log
# pezizo_set1 enrichment rerun (novelty + loss), 2026-10-14, after Amega got Pfam annotation.
# Change from the earlier run (studies/enrichment_overview/run_enrichment_all.sh and
# run_enrichment_loss_all.sh, pezizo_set1 entries): annotations/UP001658139.tsv.gz (Amega)
# now has Pfam from a Pfam-A hmmsearch (bin/annotate_proteomes_pfam.py --fill-empty-uniprot;
# the UniProt entries had no GO/Pfam/InterPro). The original file is in annotations_orig/.
# GO and InterPro still cover only the 4 other IN species; Pfam and Pfam clan cover all 5.
# The earlier outputs were moved (not deleted) to
#   outputs/proteome_clustered_Amega_unannotated
#   outputs/loss/proteome_clustered_Amega_unannotated
# Settings are the same as the earlier run: ingroup_min_frac 0.75 (novelty);
# outgroup_min_frac 0.75 and loss_ingroup_max_frac 0.0 (loss).
# --constraint=milan: mmseqs needs AVX2 (c01-c03 have none).
# Absolute paths only (no BASH_SOURCE; see ~/.claude/CLAUDE.md).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
TOOL=$NII/studies/fungi/zoosporic_dikarya/analysis/enrichment
PY="$NII/.pixi/envs/default/bin/python3 -I"
THREADS="${SLURM_NTASKS:-8}"
export SCRATCH="${SCRATCH:?SCRATCH must be set (run under SLURM)}"
NAME=pezizo_set1
S=$NII/studies/fungi/$NAME
E=$S/analysis/enrichment
R=$NII/results/$NAME
PFAMDAT=/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat
IPR=/srv/projects/db/interproscan/6.0.0/110.0/interpro/110.0/entries.json
RUNLABEL="diamond pairwise, nf_NovInvenio main05fcf3a run 2026-10-08"
ANNOT="Annotation: UniProt reference proteomes (Pfam, InterPro, GO); Amega (UP001658139, UniProt entries without cross-references): Pfam-A hmmsearch (HMMER 3.4, --cut_ga, Pfam 38.2), no GO/InterPro."
for d in "$E/outputs/proteome_clustered" "$E/outputs/loss/proteome_clustered"; do
  if [ ! -d "${d}_Amega_unannotated" ]; then
    echo "missing ${d}_Amega_unannotated: move the earlier output there first" >&2; exit 2
  fi
done
mkdir -p "$E/outputs/loss" "$E/figures"

echo "[$(date +%T)] novelty start"
T0=$(date +%s)
$PY "$TOOL/enrich_novelties_v2.py" --results "$R" --config "$S/config.csv" \
    --annotations "$S/annotations" --study "$NAME" --cache-dir "$E/outputs/cache" --threads "$THREADS" \
    --pfam-dat "$PFAMDAT" --interpro-entries "$IPR" \
    --controls "$TOOL/actin_controls.tsv" --ingroup-min-frac 0.75 \
    --universe proteome --cluster-background \
    --out-dir "$E/outputs/proteome_clustered" > "$E/outputs/run_proteome_clustered.log" 2>&1
echo "[$(date +%T)] novelty enrichment done ($(( $(date +%s) - T0 )) s)"
$PY "$TOOL/make_summary_figure.py" --enrich "$E/outputs/proteome_clustered" --config "$S/config.csv" \
    --annot-label "GO/Pfam/InterPro (UniProt) or Pfam (hmmsearch)" \
    --title "Pezizomycotina set1 (5 IN) vs OUT (6): enrichment of novelty candidates ($RUNLABEL)" \
    --footnote "$ANNOT Universe: all proteins of the IN proteomes. Family units: candidate families from the pipeline; background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Labels: k = candidate units with term / all units with term." \
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
  LOWNOTE=" LOW COVERAGE: only $FRAC of loss candidates carry any annotation; read results with care."
fi
echo "[$(date +%T)] loss enrichment done ($((T1 - T0)) s)"
$PY "$TOOL/make_summary_figure.py" --enrich "$E/outputs/loss/proteome_clustered" --config "$S/config.csv" \
    --group OUT --candidate-label "loss candidate proteins" \
    --title "OUT (6 OUT) vs Pezizomycotina set1 (5 IN): enrichment of loss-direction candidates (present in >=0.75 of OUT, absent from all IN; $RUNLABEL)" \
    --footnote "Annotation: UniProt reference proteomes (Pfam, InterPro, GO). Universe: all proteins of the 6 OUT proteomes. Family units: loss-candidate families from the pipeline (mmseqs on loss_candidates.fa); background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Loss candidates: $NCAND; annotated fraction $FRAC.$LOWNOTE" \
    --out "$E/figures/${NAME}_loss_summary" > "$E/outputs/loss/figure.log" 2>&1
echo "[$(date +%T)] all done"
