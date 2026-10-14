#!/usr/bin/bash
#SBATCH -J nf-afum-pangenome
#SBATCH -p exfab
#SBATCH -A exfab
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 30-00:00:00
#SBATCH -o /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/.nf_launch/slurm-nf-afum-%j.out
#SBATCH -e /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome/.nf_launch/slurm-nf-afum-%j.err
#
# Run nf_NovInvenio's pangenome.nf on the 295-row A. fumigatus study
# (293 IN strains + 2 OUT: A. lentulus, A. fischeri).
# Modeled on studies/fungi/coccidioides_pangenome/run_pangenome.sh and its
# .nf_launch/*/submit_nextflow_head.sh head-job scripts.
#
# Usage (head job; -resume is always on, so resubmitting continues a run):
#   mkdir -p studies/fungi/Afumigatus_pangenome/.nf_launch
#   sbatch studies/fungi/Afumigatus_pangenome/run_pangenome.sh [variant] [extra nextflow args...]
# variant defaults to nf_latest. Output: results/<variant>/.
#
# Results from results/full_293run/ came from the older bespoke study
# scripts. This script never writes there.
#
# Pipeline pinning: other sessions edit the live NovInvenio checkout, and
# Nextflow stages bin/ and lib/ live from projectDir. So the first launch
# freezes NovInvenio HEAD with `git archive` into the launch dir and runs
# that snapshot. Resubmits reuse the same snapshot (same path + mtimes, so
# -resume cache keys stay valid). The commit is recorded in
# .nf_launch/<variant>/pipeline_commit.txt. The tool env still comes from
# the live checkout's pixi env (see stajichlab_queue.config).
#
# GFF3 workaround (2026-09-25 smoke test): 3 of the 295 proteomes are
# UniProt FASTAs (the 2 outgroups and the IN reference Asfu_Af293; headers
# like `tr|A0AAN4PAS8|A0AAN4PAS8_ASPLE ... GN=IFM60648_03940`), but their
# GFF3s are NCBI/DDBJ annotation keyed by `protein_id=GFF60727.1`.
# bin/pangenome_build_gene_positions.py resolved 0/10351 (Aslen_ref) and
# 0/9647 (Asfu_Af293) IDs and hard-errors. The pipeline has no GN= <->
# locus_tag crosswalk. So this script builds a GFF3 override dir
# (--pangenome_gff3_dir): symlinks for every NCBI-FASTA strain's GFF3, and
# for each UniProt-FASTA strain (first header starts with tr| or sp|) a copy
# whose CDS protein_id= is replaced by the FASTA ID whose GN= equals the CDS
# locus_tag. Only unambiguous 1:1 locus_tag <-> GN matches are rewritten;
# others are left unchanged (they stay unresolved). Source GFF3s in
# data_dir/ are not modified.
#
# Never use BASH_SOURCE here: SLURM copies this script to a spool dir.
set -euo pipefail

VARIANT="${1:-nf_latest}"
shift || true

STUDY_DIR=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/Afumigatus_pangenome
NI_LIVE=/bigdata/stajichlab/jstajich/projects/NovInvenio
SAMPLESHEET="$STUDY_DIR/config.csv"
DATA_DIR="$STUDY_DIR/data_dir"
PFAM_HMM="$NI_LIVE/db/pfam/Pfam-A.hmm"
OUTDIR="$STUDY_DIR/results/${VARIANT}"
LAUNCH_DIR="$STUDY_DIR/.nf_launch/${VARIANT}"
GFF3_OVERRIDE="$LAUNCH_DIR/inputs/gff3"

case "$VARIANT" in
    full_293run|test_smoke_17strain) echo "ERROR: refusing to write into existing results/$VARIANT" >&2; exit 1 ;;
esac

mkdir -p "$OUTDIR" "$LAUNCH_DIR"
export PATH="$HOME/.pixi/bin:$HOME/.local/bin:$PATH"

# --- 1. freeze the pipeline commit (first launch only) ---
if [ ! -s "$LAUNCH_DIR/pipeline_commit.txt" ]; then
    git -C "$NI_LIVE" rev-parse HEAD > "$LAUNCH_DIR/pipeline_commit.txt"
fi
COMMIT="$(cat "$LAUNCH_DIR/pipeline_commit.txt")"
PIPELINE_DIR="$LAUNCH_DIR/pipeline_${COMMIT:0:7}"
if [ ! -f "$PIPELINE_DIR/pangenome.nf" ]; then
    mkdir -p "$PIPELINE_DIR"
    git -C "$NI_LIVE" archive "$COMMIT" | tar -x -C "$PIPELINE_DIR"
fi
echo "== pipeline commit: $COMMIT (snapshot $PIPELINE_DIR) ==" >&2

# --- 2. GFF3 override dir (first launch only) ---
if [ ! -f "$GFF3_OVERRIDE/.done" ]; then
    mkdir -p "$GFF3_OVERRIDE"
    /usr/bin/python3.12 - "$SAMPLESHEET" "$DATA_DIR" "$GFF3_OVERRIDE" <<'PY'
import csv, os, re, sys
from collections import defaultdict
samplesheet, data_dir, outdir = sys.argv[1:4]
for row in csv.DictReader(open(samplesheet)):
    gff = row["GFF3"].strip()
    src = os.path.join(data_dir, "gff3", gff)
    dst = os.path.join(outdir, gff)
    if os.path.lexists(dst):
        os.remove(dst)
    pep = os.path.join(data_dir, "pep", row["Protein"].strip())
    with open(pep) as fh:
        first = fh.readline()
    if not first.startswith((">tr|", ">sp|")):
        os.symlink(src, dst)
        continue
    # FASTA first token keyed by UniProt GN= (keep only unique GN values).
    gn_ids = defaultdict(list)
    with open(pep) as fh:
        for line in fh:
            if line.startswith(">"):
                m = re.search(r" GN=(\S+)", line)
                if m:
                    gn_ids[m.group(1)].append(line[1:].split()[0])
    gn_to_id = {g: ids[0] for g, ids in gn_ids.items() if len(ids) == 1}
    # locus_tag -> distinct protein_id= values on its CDS rows.
    lt_pids = defaultdict(set)
    lines = open(src).read().splitlines(keepends=True)
    for line in lines:
        f = line.split("\t")
        if len(f) >= 9 and f[2] == "CDS":
            lt = re.search(r"(?:^|;)locus_tag=([^;\n]+)", f[8])
            pid = re.search(r"(?:^|;)protein_id=([^;\n]+)", f[8])
            if lt and pid:
                lt_pids[lt.group(1)].add(pid.group(1))
    n_cds = n_rw = 0
    rewritten_lt = set()
    with open(dst, "w") as out:
        for line in lines:
            f = line.split("\t")
            if len(f) >= 9 and f[2] == "CDS":
                n_cds += 1
                lt = re.search(r"(?:^|;)locus_tag=([^;\n]+)", f[8])
                if lt and lt.group(1) in gn_to_id and len(lt_pids[lt.group(1)]) == 1:
                    new_id = gn_to_id[lt.group(1)]
                    f[8] = re.sub(r"((?:^|;)protein_id=)[^;\n]+", lambda m: m.group(1) + new_id, f[8])
                    line = "\t".join(f)
                    n_rw += 1
                    rewritten_lt.add(lt.group(1))
            out.write(line)
    print(f"gff3_override: {row['Short']}: {len(rewritten_lt)} loci / {n_rw} of {n_cds} CDS rows "
          f"rewritten to FASTA IDs; {sum(len(v) for v in gn_ids.values())} FASTA records", file=sys.stderr)
PY
    touch "$GFF3_OVERRIDE/.done"
fi

# --- 3. run ---
cd "$LAUNCH_DIR"
nextflow run "$PIPELINE_DIR/pangenome.nf" \
    --pangenome_samplesheet "$SAMPLESHEET" \
    --pangenome_data_dir "$DATA_DIR" \
    --pangenome_gff3_dir "$GFF3_OVERRIDE" \
    --pangenome_cluster_backend mmseqs \
    --pangenome_project "afumigatus_${VARIANT}" \
    --outdir "$OUTDIR" \
    --pangenome_island_pfam_hmm "$PFAM_HMM" \
    -profile slurm \
    -c "$PIPELINE_DIR/conf/ucr_hpcc_slurm.config" \
    -c "$STUDY_DIR/stajichlab_queue.config" \
    -with-trace "$OUTDIR/trace.${SLURM_JOB_ID:-local}.txt" \
    -resume \
    "$@"
