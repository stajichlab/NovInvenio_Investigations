#!/bin/bash
# HET/NLR Pfam panel (het_nlr_models.hmm, 44 models, Pfam-A 38.2) vs all 529 proteomes, --cut_ga.
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
STUDY=$NII/studies/fungi/coccidioides_pangenome
A=$STUDY/analysis/het_nlr_2026-09-25
S=${SCRATCH:?}/het_nlr; mkdir -p $S
cd $STUDY
tail -n +2 results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome/samplesheet.with_clades.csv | cut -d, -f7 | while read s; do awk -v s="$s" '/^>/{sub(/^>/,">"s"|")}{print}' data_dir/pep/$s.pep.fa; done > $S/all.fa
grep -c ">" $S/all.fa > $A/allprot.nseq
pixi run --manifest-path $NII/pixi.toml hmmsearch --cut_ga --cpu $(nproc) --tblout $S/allprot.tblout --domtblout $S/allprot.domtblout -o /dev/null $A/het_nlr_models.hmm $S/all.fa
gzip -c $S/allprot.domtblout > $A/allprot.domtblout.gz
gzip -c $S/allprot.tblout > $A/allprot.tblout.gz
rm -f $S/all.fa $S/allprot.domtblout $S/allprot.tblout
echo done > $A/allprot.done
