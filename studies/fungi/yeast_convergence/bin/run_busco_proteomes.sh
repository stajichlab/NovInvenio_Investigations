#!/usr/bin/env bash
#SBATCH --job-name=yc_busco
#SBATCH -p stajichlab
#SBATCH --constraint=milan
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/logs/slurm/yc_busco_%j.out
# BUSCO (protein mode, fungi_odb10) for every proteome in a manifest from busco_manifest.py.
# Usage: sbatch run_busco_proteomes.sh MANIFEST.tsv
# Rows with an empty fasta column are skipped. Results: analysis/busco/summary.tsv
set -euo pipefail
MANIFEST=${1:?manifest tsv}
OUT=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/yeast_convergence/analysis/busco
LINEAGE=/bigdata/stajichlab/shared/lib/funannotate_db/fungi_odb10
WORK=${SCRATCH:?SCRATCH not set}/busco
mkdir -p "$OUT/runs" "$WORK"
BUSCO_ENV=/opt/linux/rocky/8.x/x86_64/pkgs/busco/5.8.0   # `module load busco` needs `conda activate --stack`, which fails in batch jobs
export PATH="$BUSCO_ENV/bin:$PATH"
SUMMARY="$OUT/summary.tsv"
[ -s "$SUMMARY" ] || printf 'name\tsource\taccession\tcomplete_pct\tsingle_pct\tdup_pct\tfragmented_pct\tmissing_pct\tn_groups\tlineage\n' > "$SUMMARY"
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r name source acc fasta; do
  [ -n "$fasta" ] || { echo "skip $name (no FASTA)"; continue; }
  # Skip only when this name AND accession are already scored: a Short re-pointed to a new FASTA is rerun.
  awk -F'\t' -v n="$name" -v a="$acc" '$1==n && $3==a {f=1} END{exit !f}' "$SUMMARY" && { echo "done $name"; continue; }
  d="$WORK/$name"; rm -rf "$d"; mkdir -p "$d"
  case "$fasta" in
    *.gz) zcat "$fasta" > "$d/in.faa" ;;
    *) cp "$fasta" "$d/in.faa" ;;
  esac
  (cd "$d" && busco -i in.faa -m proteins -l "$LINEAGE" -c "$SLURM_CPUS_PER_TASK" -o run -f --offline > busco.out 2>&1) || { echo "BUSCO failed $name"; continue; }
  line=$(grep -h "C:.*n:" "$d"/run/short_summary*.txt | head -1 | tr -d ' \t')
  # C:90.4%[S:90.4%,D:0.0%],F:0.0%,M:9.6%,n:1312
  c=$(sed -E 's/.*C:([0-9.]+)%.*/\1/' <<<"$line"); s=$(sed -E 's/.*S:([0-9.]+)%.*/\1/' <<<"$line")
  dd=$(sed -E 's/.*D:([0-9.]+)%.*/\1/' <<<"$line"); f=$(sed -E 's/.*F:([0-9.]+)%.*/\1/' <<<"$line")
  m=$(sed -E 's/.*M:([0-9.]+)%.*/\1/' <<<"$line"); n=$(sed -E 's/.*n:([0-9]+).*/\1/' <<<"$line")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tfungi_odb10\n' "$name" "$source" "$acc" "$c" "$s" "$dd" "$f" "$m" "$n" >> "$SUMMARY"
  cp "$d"/run/short_summary*.txt "$OUT/runs/${name}.short_summary.txt"
  rm -rf "$d"
done
echo finished
