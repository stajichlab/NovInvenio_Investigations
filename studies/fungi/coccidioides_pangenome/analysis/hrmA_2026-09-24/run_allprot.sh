#!/bin/bash
set -euo pipefail
S=${SCRATCH:?}/hrmA; mkdir -p $S
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome
tail -n +2 results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome/samplesheet.with_clades.csv | cut -d, -f7 | while read s; do awk -v s="$s" '/^>/{sub(/^>/,">"s"|")}{print}' data_dir/pep/$s.pep.fa; done > $S/all.fa
grep -c ">" $S/all.fa > /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/hrmA_2026-09-24/allprot.nseq
pixi run --manifest-path /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/pixi.toml hmmsearch --cut_ga --cpu 6 --tblout /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/hrmA_2026-09-24/allprot.tblout --domtblout /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/hrmA_2026-09-24/allprot.domtblout -o /dev/null /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/hrmA_2026-09-24/hrmA_models.hmm $S/all.fa
rm -f $S/all.fa
echo done > /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/hrmA_2026-09-24/allprot.done
