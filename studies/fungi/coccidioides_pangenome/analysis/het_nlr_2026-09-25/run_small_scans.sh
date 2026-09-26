#!/bin/bash
# Small scans on sequences extracted by het_extract_seqs.py into $SCRATCH/het_nlr.
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
A=$NII/studies/fungi/coccidioides_pangenome/analysis/het_nlr_2026-09-25
PF=/bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm
S=${SCRATCH:?}/het_nlr
PX="pixi run --manifest-path $NII/pixi.toml"
# 1. relaxed P-loop NTPase check on all HET/NLR-class proteins (no GA; E <= 1e-3)
printf "NACHT\nNB-ARC\nAAA_16\nAAA_22\n" > $S/nod.txt
$PX hmmfetch -f $PF $S/nod.txt > $S/nod.hmm
$PX hmmsearch -E 1e-3 --domE 1e-3 --cpu $(nproc) --domtblout $S/nod_relaxed.domtblout -o /dev/null $S/nod.hmm $S/het_prots.fa
gzip -c $S/nod_relaxed.domtblout > $A/het_prots_nod_relaxed.domtblout.gz
# 2. rep pairwise identity (ortholog-splitting check)
$PX diamond makedb --in $S/het_reps.fa -d $S/het_reps --quiet
$PX diamond blastp --ultra-sensitive -q $S/het_reps.fa -d $S/het_reps -e 1e-5 --max-target-seqs 100 --quiet \
  -f 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore -o $A/het_rep_pairwise.tsv
# 3. full Pfam (GA) on HET/NLR reps and neighbour reps
cat $S/het_reps.fa $S/nb_reps.fa > $S/scan.fa
$PX hmmscan --cut_ga --cpu $(nproc) --domtblout $S/reps_pfam.domtblout -o /dev/null $PF $S/scan.fa
gzip -c $S/reps_pfam.domtblout > $A/reps_pfam.domtblout.gz
echo done > $A/small_scans.done
