#!/bin/bash -l
#SBATCH -p epyc --constraint=milan -c 16 --mem 32G -t 2:00:00
#SBATCH -J cocci_captain_sm
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/starship_island_overlap/logs/captain_sm.%j.out
# Cocci captain (DUF3435, PF11917.14) search on the pangenome study proteins, and
# antiSMASH region collection + DUF3435 search on the antiSMASH annotation proteins.
# Same HMM and E-value cutoff (full-sequence E <= 1e-5 applied downstream) as the
# Afumigatus results/captain_gene search.
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
CP=$NII/studies/fungi/coccidioides_pangenome
OUT=$CP/analysis/starship_island_overlap/inputs
HMM=$NII/studies/fungi/Afumigatus_pangenome/results/captain_gene/DUF3435.hmm
BIN=$NII/.pixi/envs/default/bin
ANN=/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/annotation
W=${SCRATCH:?}/cocci_captain_sm
mkdir -p "$W" "$OUT"

# 1. study proteins (annotation freeze 20260112, as used by the pangenome runs)
: > "$W/study_proteins.faa"
tail -n +2 "$CP/config_immitis_in_posadasii_out.csv" | while IFS=, read -r grp sp strain prot dna gff short tg; do
  awk -v s="$short" '/^>/{sub(/^>/,""); split($0,a," "); print ">" s "|" a[1]; next} {print}' \
    "$CP/data_dir/pep/$prot" >> "$W/study_proteins.faa"
done
grep -c '^>' "$W/study_proteins.faa"
"$BIN/hmmsearch" --cpu 16 --noali --tblout "$W/DUF3435_vs_study.tblout" "$HMM" "$W/study_proteins.faa" > "$W/DUF3435_vs_study.log"

# 2. antiSMASH regions (own annotation space)
"$BIN/python3" "$CP/analysis/starship_island_overlap/scripts/parse_antismash.py" \
  --annotation-root "$ANN" --config "$CP/config_immitis_in_posadasii_out.csv" --out-dir "$W/antismash"
"$BIN/hmmsearch" --cpu 16 --noali --tblout "$W/DUF3435_vs_antismash_annotation.tblout" "$HMM" "$W/antismash/all_proteins.faa" > "$W/DUF3435_vs_antismash_annotation.log"

# copy back (small files only; protein FASTAs stay on scratch)
zstd -q -f "$W/DUF3435_vs_study.tblout" -o "$OUT/DUF3435_vs_study.tblout.zst"
zstd -q -f "$W/DUF3435_vs_antismash_annotation.tblout" -o "$OUT/DUF3435_vs_antismash_annotation.tblout.zst"
tail -n 20 "$W/DUF3435_vs_study.log" > "$OUT/DUF3435_vs_study.log"
tail -n 20 "$W/DUF3435_vs_antismash_annotation.log" > "$OUT/DUF3435_vs_antismash_annotation.log"
mkdir -p "$OUT/antismash"
cp "$W/antismash/"*.tsv "$W/antismash/"*.tsv.gz "$OUT/antismash/"
sha256sum "$W/study_proteins.faa" "$W/antismash/all_proteins.faa" "$HMM" > "$OUT/input_sha256.txt"
grep -c '^>' "$W/antismash/all_proteins.faa" >> "$OUT/input_sha256.txt"
echo done
