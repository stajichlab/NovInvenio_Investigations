#!/usr/bin/env python3
"""Extract (a) HET/NLR family reps, (b) all HET/NLR-class proteins, (c) neighbour/marker
family reps, as FASTA into $SCRATCH/het_nlr. Read-only on pipeline outputs."""
import os, sys
from collections import Counter
from pathlib import Path
import pandas as pd
STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome")
P = STUDY / "results/rescue_freqpol_immitis_in_posadasii_out/output/pangenome"
A = STUDY / "analysis/het_nlr_2026-09-25"
S = Path(os.environ["SCRATCH"]) / "het_nlr"
c = pd.read_csv(A / "het_family_census.tsv", sep="\t")
lt = pd.read_csv(A / "het_loci.tsv", sep="\t")
nb = pd.read_csv(A / "het_neighbourhoods.tsv.gz", sep="\t")
h = pd.read_csv(A / "het_gene_locations.tsv.gz", sep="\t")
reps = set(c.family)
# neighbour families: flank >= 10% of a family's anchors, plus all locus markers
nanch = nb[nb.offset == 0].groupby("anchor_family").anchor.nunique()
fr = nb[nb.offset != 0].groupby(["anchor_family", "family"]).anchor.nunique()
fr = fr / fr.index.get_level_values(0).map(nanch).values
nbf = set(fr[fr >= 0.10].index.get_level_values(1)) - reps
for m in lt.markers.dropna():
    nbf |= set(m.split(","))
nbf -= reps
done = set()
hp = Path(STUDY / "analysis/hrmA_2026-09-24/neighbour_reps_pfam.domtblout")
for line in open(hp):
    if not line.startswith("#"):
        done.add(line.split()[3])
hr = set(l.strip() for l in open(STUDY / "analysis/hrmA_2026-09-24/neighbour_families.txt"))
todo = nbf - hr
print("neighbour families:", len(nbf), "already scanned in hrmA work:", len(nbf & hr), "to scan:", len(todo), file=sys.stderr)
(A / "neighbour_families.txt").write_text("\n".join(sorted(nbf)) + "\n")
want = reps | todo
out_r, out_n = open(S / "het_reps.fa", "w"), open(S / "nb_reps.fa", "w")
cur = None
with open(P / "cluster/tier1_rep_seq.fasta") as fh:
    for line in fh:
        if line.startswith(">"):
            i = line[1:].split()[0]
            cur = out_r if i in reps else (out_n if i in todo else None)
        if cur:
            cur.write(line)
out_r.close(); out_n.close()
# all HET/NLR proteins from per-strain proteomes
need = {}
for m in h.member:
    s, p = m.split("|", 1)
    need.setdefault(s, set()).add(p)
with open(S / "het_prots.fa", "w") as out:
    for s, ps in need.items():
        keep = False
        for line in open(STUDY / f"data_dir/pep/{s}.pep.fa"):
            if line.startswith(">"):
                pid = line[1:].split()[0]
                keep = pid in ps
                if keep:
                    out.write(f">{s}|{pid}\n"); continue
            if keep:
                out.write(line)
