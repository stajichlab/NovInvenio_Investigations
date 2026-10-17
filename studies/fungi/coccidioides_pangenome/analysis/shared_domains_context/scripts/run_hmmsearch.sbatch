#!/bin/bash -l
#SBATCH -p epyc --constraint=milan -c 16 --mem 24G -t 3:00:00
#SBATCH -J sdc_hmmsearch
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/shared_domains_context/logs/hmmsearch.%j.out
# hmmsearch --cut_ga of the 13 shared-domain Pfam models (hmm/shared_domains.hmm,
# hmmfetch from Pfam 38.2) against every proteome in the PF11001 manifest
# (295 A. fumigatus-study, 529 Coccidioides, 144 survey proteomes).
# Each query model is scored independently (same as one single-profile run per model).
# One domtblout per proteome in $SCRATCH, then combined (idx-prefixed, zstd) into hmmsearch/.
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
D=$NII/studies/fungi/coccidioides_pangenome/analysis/shared_domains_context
B=$NII/.pixi/envs/default/bin
HMM=$D/hmm/shared_domains.hmm
MAN=$D/hmmsearch/manifest.tsv
W=${SCRATCH:?}/sdc
mkdir -p "$W/out"
tail -n +2 "$MAN" | awk -F'\t' '{print NR"\t"$1"\t"$5"\t"$7}' > "$W/jobs.tsv"
run_one() {
  idx=$1; pep=$4
  "$B/hmmsearch" --cut_ga --noali -o /dev/null --cpu 1 --domtblout "$W/out/${idx}.domtbl" "$HMM" "$pep"
}
export -f run_one
export B HMM W
cut -f1-4 "$W/jobs.tsv" | xargs -P 16 -L 1 bash -c 'run_one "$@"' _
n=$(ls "$W"/out/*.domtbl | wc -l)
echo "domtbl files: $n"
for f in "$W"/out/*.domtbl; do
  idx=$(basename "$f" .domtbl)
  grep -v '^#' "$f" | awk -v i="$idx" 'NF{print i"\t"$0}' || true
done | zstd -q -T0 -f -o "$D/hmmsearch/all_proteomes.shared_domains.domtbl.idx.zst"
grep '^#' "$W/out/1.domtbl" > "$D/hmmsearch/domtbl_header.txt" || true
echo done
