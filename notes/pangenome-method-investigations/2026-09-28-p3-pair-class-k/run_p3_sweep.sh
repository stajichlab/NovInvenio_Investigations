#!/usr/bin/env bash
# #132 priority 3: rerun PAIR_CLASSIFICATION alone over k and the two linkage
# thresholds, on the exact inputs of two v0.7.0 runs (pinned code a394ace).
# Outputs (not committed, ~1 MB each): $OUT/<dataset>/k<k>_p<phys>_t<trans>.tsv.zst
#SBATCH -p epyc
#SBATCH -c 16
#SBATCH --mem 24G
#SBATCH -t 06:00:00
#SBATCH -J p3-pairk
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/results/p3_pair_class_sweep/logs/%x_%j.out
set -euo pipefail
C=/rhome/jstajich/.nextflow/assets/.repos/stajichlab/NovInvenio/clones/a394ace8d288b19591ea148e5a12d55aeac668b1
PY=$C/.pixi/envs/default/bin/python
OUT=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/results/p3_pair_class_sweep
declare -A W=(
  [Cocci_Ci]=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/.nf_launch/immitis_in_posadasii_out_v2/work/c7/15058e841ae6857ceb15dab867381a
  [Afumigatus]=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/.nf_launch/full_v070/work/a3/6396d4c0491625651ee210dbe913a0
)
SETTINGS="3,0.5,0.05 5,0.5,0.05 7,0.5,0.05 10,0.5,0.05 15,0.5,0.05 20,0.5,0.05 30,0.5,0.05 50,0.5,0.05
10,0.3,0.05 10,0.7,0.05 10,0.9,0.05 10,0.5,0.01 10,0.5,0.1 10,0.5,0.2"
run() {  # dataset k phys trans
  local d=$1 k=$2 p=$3 t=$4 w=${W[$1]}
  mkdir -p $OUT/$d
  local o=$OUT/$d/k${k}_p${p}_t${t}.tsv.zst
  [ -s $o ] && return 0
  cd "$w"
  args=$(sed -n 's/.*pangenome_pair_classification.py//p' .command.sh | sed "s/--k [0-9]*/--k $k/; s/--physical_threshold [0-9.]*/--physical_threshold $p/; s/--trans_threshold [0-9.]*/--trans_threshold $t/; s#--output [^ ]*#--output $o#")
  /usr/bin/time -f "$d k=$k p=$p t=$t %e s %M KB" $PY $C/bin/pangenome_pair_classification.py $args
}
export -f run; export OUT PY C
for d in Cocci_Ci Afumigatus; do
  for s in $SETTINGS; do IFS=, read k p t <<< "$s"; echo "$d $k $p $t"; done
done | xargs -P 4 -L 1 bash -c 'declare -A W=([Cocci_Ci]='"${W[Cocci_Ci]}"' [Afumigatus]='"${W[Afumigatus]}"'); '"$(declare -f run)"'; run "$@"' _
echo DONE
