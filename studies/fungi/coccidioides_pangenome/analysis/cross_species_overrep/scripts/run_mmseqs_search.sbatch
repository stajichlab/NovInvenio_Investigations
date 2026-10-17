#!/bin/bash -l
#SBATCH -p stajichlab --constraint=milan -c 32 --mem 48G -t 4:00:00
#SBATCH -J xsp_mmseqs
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/cross_species_overrep/logs/mmseqs.%j.out
# All-vs-all family-representative searches between A. fumigatus full_v072 and each
# Coccidioides run (v5, v2), both directions, for best hits and reciprocal best hits.
# Parameters: mmseqs easy-search -e 1e-5 -c 0.5 --cov-mode 0 (query AND target coverage
# >= 50%), default sensitivity (-s 5.7), default --max-seqs 300.
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
B=$NII/.pixi/envs/default/bin
O=$NII/studies/fungi/coccidioides_pangenome/analysis/cross_species_overrep/mmseqs
AF=$NII/studies/fungi/Afumigatus_pangenome/results/full_v072/config/pangenome/cluster/tier1_rep_seq.fasta
V5=$NII/studies/fungi/coccidioides_pangenome/results/immitis_in_posadasii_out_v5/config_immitis_in_posadasii_out/pangenome/cluster/tier1_rep_seq.fasta
V2=$NII/studies/fungi/coccidioides_pangenome/results/posadasii_in_immitis_out_v2/config_posadasii_in_RS_WA211_out/pangenome/cluster/tier1_rep_seq.fasta
W=${SCRATCH:?}/xsp_mmseqs
mkdir -p "$W" "$O"
"$B/mmseqs" version
FMT=query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov
run() {  # name query target
  "$B/mmseqs" easy-search "$2" "$3" "$W/$1.m8" "$W/tmp_$1" \
    -e 1e-5 -c 0.5 --cov-mode 0 --threads 32 --format-output "$FMT" > "$W/$1.log" 2>&1
  zstd -q -f -T0 -o "$O/$1.m8.zst" "$W/$1.m8"
  tail -3 "$W/$1.log"
  rm -rf "$W/tmp_$1"
}
run af_vs_v5 "$AF" "$V5"
run v5_vs_af "$V5" "$AF"
run af_vs_v2 "$AF" "$V2"
run v2_vs_af "$V2" "$AF"
echo done
