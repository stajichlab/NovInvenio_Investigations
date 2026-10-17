#!/bin/bash -l
#SBATCH -p stajichlab --constraint=milan -c 8 --mem 24G -t 2:00:00
#SBATCH -J af293_antismash
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/analysis/starship_island_overlap/full_v072_antismash/antismash_run/antismash.%j.log
# antiSMASH on A. fumigatus Af293 (UP000002530 study genome + NCBI GenBank GFF3 used by the
# Afumigatus pangenome, data_dir/). Gene models from the GFF3 (no gene finding).
# Run dir in node-local $SCRATCH; JSON, region GenBank files and tables are copied back.
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
AF=$NII/studies/fungi/Afumigatus_pangenome
OUT=$AF/analysis/starship_island_overlap/full_v072_antismash/antismash_run
DNA=$AF/data_dir/dna/UP000002530_330879.dna.fa
GFF=$AF/data_dir/gff3/UP000002530_330879.gff3
CPU=${SLURM_CPUS_ON_NODE:-8}
W=${SCRATCH:?SCRATCH not set}/af293_antismash
mkdir -p "$W"
module load antismash/8.0.4
echo "host: $(hostname)  cpus: $CPU  job: ${SLURM_JOB_ID:-none}"
which antismash
antismash --version
START=$(date +%s)
echo "start: $(date -Iseconds)"
CMD=(antismash --taxon fungi --genefinding-gff3 "$GFF" --output-dir "$W/run" --output-basename Af293
     --fullhmmer --clusterhmmer --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go
     -c "$CPU" "$DNA")
echo "command: ${CMD[*]}"
"${CMD[@]}"
END=$(date +%s)
echo "end: $(date -Iseconds)"
echo "runtime_seconds: $((END - START))"
du -sh "$W/run"; du -sh "$W/run"/* | sort -h | tail -20
# copy back: JSON (gz), region GenBank files (tar.gz), clusterblast/knownclusterblast/subclusterblast txt (tar.gz)
gzip -c "$W/run/Af293.json" > "$OUT/Af293.json.gz"
( cd "$W/run" && tar czf "$OUT/region_gbks.tar.gz" ./*.region*.gbk )
( cd "$W/run" && tar czf "$OUT/clusterblast_txt.tar.gz" $(ls -d clusterblast knownclusterblast subclusterblast 2>/dev/null) ) || true
SZ_GBK=$(stat -c %s "$W/run/Af293.gbk")
if [ "$SZ_GBK" -le 20000000 ]; then gzip -c "$W/run/Af293.gbk" > "$OUT/Af293.gbk.gz"; else echo "Af293.gbk $SZ_GBK bytes > 20 MB: not copied"; fi
sha256sum "$DNA" "$GFF" > "$OUT/inputs.sha256"
echo "runtime_seconds: $((END - START))" > "$OUT/runtime.txt"
echo "${CMD[*]}" > "$OUT/command.txt"
antismash --version > "$OUT/version.txt"
ls -la "$OUT"
rm -rf "$W"
