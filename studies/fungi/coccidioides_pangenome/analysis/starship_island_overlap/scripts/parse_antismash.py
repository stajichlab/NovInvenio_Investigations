#!/usr/bin/env python3
"""Collect antiSMASH 8 region calls and protein sets for Coccidioides strains.

Source: per-strain funannotate + antiSMASH runs under
/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/annotation/<strain>/
  antismash_local/*.region*.gbk     antiSMASH regions (Orig. start/end in the header)
  annotate_misc/genome.proteins.fasta  proteins of the same annotation
  annotate_results/*.gff3            gene coordinates of the same annotation

IMPORTANT: this annotation is NOT the one used by the pangenome runs (annotation
freeze 20260112). Contig names (RagTag scaffolds) and gene models differ, and the
same locus tag can name a different protein in the freeze. So everything here
stays in the antiSMASH annotation's own coordinate and gene space.

Only strains whose directory name equals a study Short ID (exact match) are used.

Outputs (in --out-dir):
  sm_regions.tsv.gz    one row per antiSMASH region
  sm_region_genes.tsv.gz  one row per CDS in a region (locus_tag, gene_kind)
  gene_counts.tsv      strain, n_proteins (genome.proteins.fasta records)
  strain_map.tsv       Short, annotation dir found (yes/no)
  all_proteins.faa     all proteins, header ">Short|protein_id" (for hmmsearch)
  gene_coords.tsv.gz   Short, protein_id(locus_tag), contig, start, end (mRNA rows of the gff3)
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
import os
import re
import sys
from pathlib import Path

from Bio import SeqIO

ORIG_START = re.compile(r"Orig\. start\s*::\s*(\d+)")
ORIG_END = re.compile(r"Orig\. end\s*::\s*(\d+)")


def read_shorts(configs: list[Path]) -> list[str]:
    shorts = []
    for c in configs:
        with c.open() as fh:
            for row in csv.DictReader(fh):
                if row["Short"] not in shorts:
                    shorts.append(row["Short"])
    return shorts


def parse_region(path: str):
    text = open(path).read()
    m1, m2 = ORIG_START.search(text), ORIG_END.search(text)
    rec = SeqIO.read(path, "genbank")
    off = int(m1.group(1)) if m1 else 0
    region = next(f for f in rec.features if f.type == "region")
    products = region.qualifiers.get("product", [])
    cats = sorted({c for f in rec.features if f.type == "protocluster" for c in f.qualifiers.get("category", [])})
    genes = []
    for f in rec.features:
        if f.type != "CDS":
            continue
        tag = f.qualifiers.get("locus_tag", [""])[0]
        kind = f.qualifiers.get("gene_kind", [""])[0]
        genes.append((tag, kind, int(f.location.start) + off, int(f.location.end) + off))
    return {
        "contig": rec.name if rec.name else rec.id,
        "record_id": rec.id,
        "start": off + int(region.location.start) + (0 if m1 else 0),
        "end": off + int(region.location.end),
        "orig_start": int(m1.group(1)) if m1 else "",
        "orig_end": int(m2.group(1)) if m2 else "",
        "products": ",".join(products),
        "categories": ",".join(cats),
        "contig_edge": region.qualifiers.get("contig_edge", [""])[0],
        "genes": genes,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation-root", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path, action="append")
    ap.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    shorts = read_shorts(a.config)

    smap = open(a.out_dir / "strain_map.tsv", "w")
    smap.write("Short\tannotation_dir\tn_region_gbk\n")
    reg_out = gzip.open(a.out_dir / "sm_regions.tsv.gz", "wt")
    reg_out.write("region_id\tShort\tcontig\tstart\tend\tproducts\tcategories\tcontig_edge\tn_cds\tn_biosynthetic\n")
    gene_out = gzip.open(a.out_dir / "sm_region_genes.tsv.gz", "wt")
    gene_out.write("region_id\tShort\tlocus_tag\tgene_kind\tstart\tend\n")
    counts = open(a.out_dir / "gene_counts.tsv", "w")
    counts.write("Short\tn_proteins\n")
    faa = open(a.out_dir / "all_proteins.faa", "w")
    coords = gzip.open(a.out_dir / "gene_coords.tsv.gz", "wt")
    coords.write("Short\tprotein_id\tcontig\tstart\tend\n")

    for short in shorts:
        d = a.annotation_root / short
        gbks = sorted(glob.glob(str(d / "antismash_local" / "*.region*.gbk")))
        prot = d / "annotate_misc" / "genome.proteins.fasta"
        gff = sorted(glob.glob(str(d / "annotate_results" / "*.gff3")))
        ok = d.is_dir() and gbks is not None and prot.exists() and len(gff) == 1 and os.path.isdir(d / "antismash_local")
        smap.write(f"{short}\t{'yes' if ok else 'no'}\t{len(gbks)}\n")
        if not ok:
            continue
        for g in gbks:
            r = parse_region(g)
            rid = f"{short}:{Path(g).name.replace('.gbk', '')}"
            nb = sum(1 for x in r["genes"] if x[1] == "biosynthetic")
            reg_out.write(f"{rid}\t{short}\t{r['contig']}\t{r['start']}\t{r['end']}\t{r['products']}\t"
                          f"{r['categories']}\t{r['contig_edge']}\t{len(r['genes'])}\t{nb}\n")
            for tag, kind, s, e in r["genes"]:
                gene_out.write(f"{rid}\t{short}\t{tag}\t{kind}\t{s}\t{e}\n")
        n = 0
        for rec in SeqIO.parse(str(prot), "fasta"):
            n += 1
            faa.write(f">{short}|{rec.id}\n{str(rec.seq)}\n")
        counts.write(f"{short}\t{n}\n")
        with open(gff[0]) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[2] != "mRNA":
                    continue
                pid = re.search(r"ID=([^;]+)", p[8]).group(1)
                coords.write(f"{short}\t{pid}\t{p[0]}\t{p[3]}\t{p[4]}\n")
    for fh in (smap, reg_out, gene_out, counts, faa, coords):
        fh.close()
    print("done", file=sys.stderr)


if __name__ == "__main__":
    main()
