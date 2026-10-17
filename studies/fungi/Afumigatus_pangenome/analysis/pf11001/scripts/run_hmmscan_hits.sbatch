#!/bin/bash -l
#SBATCH -p epyc --constraint=milan -c 16 --mem 6G -t 2:00:00
#SBATCH -J pf11001_hmmscan
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/analysis/pf11001/logs/hmmscan.%j.out
# Full Pfam-A (--cut_ga) scan of every PF11001 hit protein, to list the other
# Pfam domains on the same proteins.
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
D=$NII/studies/fungi/Afumigatus_pangenome/analysis/pf11001
B=$NII/.pixi/envs/default/bin
W=${SCRATCH:?}/pf11001_scan
mkdir -p "$W"
"$B/python3" "$D/scripts/extract_hit_proteins.py" "$D/hmmsearch" "$W/hits.fa"
cp "$W/hits.fa.lengths.tsv" "$D/hmmsearch/hit_protein_lengths.tsv"
zstd -q -f -o "$D/hmmsearch/hit_proteins.fa.zst" "$W/hits.fa"
cp /srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm* "$W/"
"$B/hmmscan" --cut_ga --noali -o /dev/null --cpu 16 --domtblout "$W/hits_vs_pfam.domtbl" \
  "$W/Pfam-A.hmm" "$W/hits.fa"
zstd -q -f -o "$D/hmmsearch/hit_proteins_vs_PfamA.cutga.domtbl.zst" "$W/hits_vs_pfam.domtbl"
echo done
