#!/usr/bin/bash
#SBATCH -p stajichlab -n 8 --mem 32G -t 04:00:00 -J nii_enrich_loss --constraint=milan
# --constraint=milan: the pixi mmseqs binary needs AVX2. stajichlab nodes c01-c03 (AMD
# Opteron 6376, feature abu_dhabi) have no AVX2; mmseqs dies there with SIGILL (rc 132) and an
# empty stderr (seen in jobs 29893156 and 29893409, both on c01).
#SBATCH --out /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/enrichment_overview/run_enrichment_loss_all.%j.log
# Primary loss-direction enrichment (same settings as
# studies/eukaryotes/zoosporic_holozoa_vs_dikarya/analysis/enrichment/run_enrichment_loss.sh,
# primary run only: no loss_loose K sets) for the studies that support it:
# loss_candidates.txt non-empty, loss_presence_matrix.tsv present, >= 2 OUT species with a
# UniProt annotation file, not a sensitivity/baseline copy.
# Usage: sbatch run_enrichment_loss_all.sh STUDY [STUDY ...]   (STUDY = folder name in the table)
# No minimum candidate-annotation fraction is applied. The annotated fraction of loss
# candidates (coverage.tsv ALL_OUT candidates frac_any) is appended to PROVENANCE.yaml, and a
# value < 0.5 adds a low-coverage note to the figure footnote.
# Absolute paths only (no BASH_SOURCE; see ~/.claude/CLAUDE.md).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
TOOL=$NII/studies/fungi/zoosporic_dikarya/analysis/enrichment
PY="$NII/.pixi/envs/default/bin/python3 -I"
THREADS="${SLURM_NTASKS:-8}"
export SCRATCH="${SCRATCH:?SCRATCH must be set (run under SLURM)}"

# folder|domain/folder|OUT label|IN label|run label
TABLE=(
  "chytrid_vs_mucoromycota|fungi/chytrid_vs_mucoromycota|Mucoromycota|chytrids|diamond pairwise, nf_NovInvenio run 2026-10-11"
  "malassezia_vs_ustilaginomycotina|fungi/malassezia_vs_ustilaginomycotina|Ustilaginomycotina|Malassezia|diamond pairwise, nf_NovInvenio run 2026-10-11"
  "metschnikowiaceae_cah|fungi/metschnikowiaceae_cah|other Saccharomycotina|Metschnikowiaceae incl. C. auris|diamond pairwise, nf_NovInvenio run 2026-10-11"
  "mushrooms_tremella|fungi/mushrooms_tremella|Tremella|mushrooms|diamond pairwise, nf_NovInvenio run 2026-09-10"
  "neocallimastigomycota_vs_early_fungi|fungi/neocallimastigomycota_vs_early_fungi|early fungi|Neocallimastigomycota|diamond pairwise"
  "pezizo_set1|fungi/pezizo_set1|OUT|Pezizomycotina set1|diamond pairwise, nf_NovInvenio main05fcf3a run 2026-10-08"
  "serinales_vs_pichiales|fungi/serinales_vs_pichiales|Pichiales|Serinales|diamond pairwise"
  "sordariales_shallow|fungi/sordariales_shallow|OUT|Sordariales shallow|diamond pairwise"
  "yeast_filamentous|fungi/yeast_filamentous|filamentous Ascomycota|yeasts|diamond pairwise, nf_NovInvenio run 2026-09-10"
  "zoosporic_dikarya|fungi/zoosporic_dikarya|Dikarya|zoosporic fungi|diamond pairwise"
)

for want in "$@"; do
  hit=""
  for entry in "${TABLE[@]}"; do
    IFS='|' read -r NAME ST OUTL INL RUNLABEL <<< "$entry"
    if [ "$NAME" = "$want" ]; then hit=1; break; fi
  done
  if [ -z "$hit" ]; then echo "unknown study $want" >&2; exit 2; fi
  S=$NII/studies/$ST
  E=$S/analysis/enrichment
  O=$E/outputs/loss
  mkdir -p "$O" "$E/figures"
  NOUT=$(awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["GROUP"]=="OUT"' "$S/config.csv" | wc -l)
  NIN=$(awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["GROUP"]=="IN"' "$S/config.csv" | wc -l)
  NCAND=$(grep -c . "$NII/results/$NAME/loss_candidates.txt")
  echo "[$(date +%T)] $NAME start (OUT $NOUT, IN $NIN, loss candidates $NCAND)"
  T0=$(date +%s)
  $PY "$TOOL/enrich_novelties_v2.py" --results "$NII/results/$NAME" --config "$S/config.csv" \
      --annotations "$S/annotations" --study "$NAME" --cache-dir "$E/outputs/cache" --threads "$THREADS" \
      --pfam-dat /srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat \
      --interpro-entries /srv/projects/db/interproscan/6.0.0/110.0/interpro/110.0/entries.json \
      --controls "$TOOL/actin_controls.tsv" \
      --direction loss --ingroup-min-frac 0.75 --other-max-frac 0.0 \
      --universe proteome --cluster-background \
      --out-dir "$O/proteome_clustered" > "$O/run_proteome_clustered.log" 2>&1 || {
    echo "[$(date +%T)] $NAME FAILED enrichment; see $O/run_proteome_clustered.log"; continue; }
  T1=$(date +%s)
  # Annotated fraction of loss candidates (any of Pfam/InterPro/GO) -> PROVENANCE.yaml
  FRAC=$($PY -c "import csv,sys; r=[x for x in csv.DictReader(open(sys.argv[1]),delimiter='\t') if x['species']=='ALL_OUT' and x['subset']=='candidates']; print(r[0]['frac_any'])" "$O/proteome_clustered/coverage.tsv")
  LOW=$($PY -c "import sys; print('true' if float(sys.argv[1]) < 0.5 else 'false')" "$FRAC")
  printf 'candidate_annotation:\n  frac_loss_candidates_any_annotation: %s\n  low_coverage: %s\n  rule: "low_coverage = fraction < 0.5 (flag only; no run was excluded on this value)"\n  wall_s_enrichment: %s\n' \
      "$FRAC" "$LOW" "$((T1 - T0))" >> "$O/proteome_clustered/PROVENANCE.yaml"
  LOWNOTE=""
  if [ "$LOW" = true ]; then
    LOWNOTE=" LOW COVERAGE: only $FRAC of loss candidates carry any UniProt annotation; read results with care."
  fi
  $PY "$TOOL/make_summary_figure.py" --enrich "$O/proteome_clustered" --config "$S/config.csv" \
      --group OUT --candidate-label "loss candidate proteins" \
      --title "$OUTL ($NOUT OUT) vs $INL ($NIN IN): enrichment of loss-direction candidates (present in >=0.75 of OUT, absent from all IN; $RUNLABEL)" \
      --footnote "Annotation: UniProt reference proteomes (Pfam, InterPro, GO). Universe: all proteins of the $NOUT OUT proteomes. Family units: loss-candidate families from the pipeline (mmseqs on loss_candidates.fa); background clustered with mmseqs (0.3 id, 0.8 cov, cov-mode 0). One-sided hypergeometric; BH within each vocabulary. Bars: FDR<0.05 and k>=3 candidate units. Loss candidates: $NCAND; annotated fraction $FRAC.$LOWNOTE" \
      --out "$E/figures/${NAME}_loss_summary" > "$O/figure.log" 2>&1 || echo "[$(date +%T)] $NAME FAILED figure; see $O/figure.log"
  echo "[$(date +%T)] $NAME done (enrichment ${T1}-${T0} = $((T1 - T0)) s, total $(( $(date +%s) - T0 )) s)"
done
