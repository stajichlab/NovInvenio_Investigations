#!/bin/bash -l
#SBATCH -p epyc --constraint=milan -c 16 --mem 16G -t 2:00:00
#SBATCH -J pf11001_hmmsearch
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/analysis/pf11001/logs/hmmsearch.%j.out
# hmmsearch PF11001 (single profile, --cut_ga) against every proteome in manifest.tsv.
# One domtblout per proteome, written to $SCRATCH, then combined (zstd) into hmmsearch/.
# Also records the number of sequences per proteome (seq_counts.tsv).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
D=$NII/studies/fungi/Afumigatus_pangenome/analysis/pf11001
B=$NII/.pixi/envs/default/bin
HMM=$D/hmm/PF11001.hmm
MAN=$D/hmmsearch/manifest.tsv
W=${SCRATCH:?}/pf11001
mkdir -p "$W/out"
tail -n +2 "$MAN" | awk -F'\t' '{print NR"\t"$1"\t"$5"\t"$7}' > "$W/jobs.tsv"
run_one() {
  idx=$1; set_=$2; short=$3; pep=$4
  o="$W/out/${idx}.domtbl"
  "$B/hmmsearch" --cut_ga --noali -o /dev/null --cpu 1 --domtblout "$o" "$HMM" "$pep"
  n=$(grep -c '^>' "$pep")
  printf '%s\t%s\t%s\t%s\n' "$idx" "$set_" "$short" "$n" > "$W/out/${idx}.n"
}
export -f run_one
export B HMM W
cut -f1-4 "$W/jobs.tsv" | xargs -P 16 -L 1 bash -c 'run_one "$@"' _
# combine: prefix every data line with idx, set, short
{
  printf 'idx\tset\tshort\tn_seqs\n'
  cat "$W"/out/*.n | sort -n
} > "$D/hmmsearch/seq_counts.tsv"
for f in "$W"/out/*.domtbl; do
  idx=$(basename "$f" .domtbl)
  grep -v '^#' "$f" | awk -v i="$idx" 'NF{print i"\t"$0}' || true
done | zstd -q -T0 -o "$D/hmmsearch/all_proteomes.PF11001.domtbl.idx.zst" -f
grep '^#' "$W/out/1.domtbl" | head -3 > "$D/hmmsearch/domtbl_header.txt"
grep '^# Option\|^# Program\|^# Version\|^# Date' "$W/out/1.domtbl" >> "$D/hmmsearch/domtbl_header.txt" || true
echo done
