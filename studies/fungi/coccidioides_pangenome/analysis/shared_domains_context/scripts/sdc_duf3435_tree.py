#!/usr/bin/env python3
"""DUF3435 (PF11917) domain tree: Af293, C. immitis RS, C. posadasii Silveira and 12 other
Dikarya survey proteomes (one proteome per species).

Steps: domain envelopes (hmmsearch --cut_ga env coordinates, tables/genes.*.tsv.gz) ->
tree/duf3435_env.fa -> hmmalign --trim --amino (hmm/DUF3435.hmm) -> esl-reformat afa ->
FastTree -lg -gamma (SH-like local support) -> midpoint root (Bio.Phylo) -> clade summary.

Writes tree/duf3435_env.fa, duf3435.trim.sto, duf3435.trim.afa, duf3435.fasttree.nwk,
duf3435_tips.tsv, duf3435_clades.tsv, duf3435_tree_summary.tsv, fasttree.log
"""
from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from collections import Counter
from pathlib import Path

from Bio import Phylo

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
D = NII / "studies/fungi/coccidioides_pangenome/analysis/shared_domains_context"
T = D / "tables"
TR = D / "tree"
B = NII / ".pixi/envs/default/bin"
GENUS = {"Asfu_Af293": "A. fumigatus Af293", "RS": "Coccidioides", "UCSF_Cp_Silveira": "Coccidioides"}
SURVEY = ["Amega", "Anid", "Chaeglob", "Clonrose", "Collhigg", "Conilign", "Cryppara", "Metaanis",
          "Vertdahl", "Ztri", "Tillcari", "Violpalu"]


def fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split(None, 1)[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main():
    TR.mkdir(exist_ok=True)
    man = {}
    with open(D / "hmmsearch/manifest.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            man.setdefault(r["short"], r)
    want = set(GENUS) | set(SURVEY)
    hits = []
    for s in ("afum", "cocci", "survey"):
        with gzip.open(T / f"genes.{s}.tsv.gz", "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["domain"] == "DUF3435" and r["short"] in want:
                    hits.append(r)
    # survey proteomes can be listed twice in the manifest under one short name; keep the first pep
    seqcache = {}
    recs = []
    tips = []
    for r in hits:
        short = r["short"]
        if short not in seqcache:
            seqcache[short] = fasta(man[short]["pep_path"])
        seq = seqcache[short][r["protein_id"]]
        for k, env in enumerate(r["env"].split(";"), 1):
            a, b = map(int, env.split("-"))
            name = f"{short}__{r['protein_id'].replace('|', '_')}__d{k}"
            recs.append((name, seq[a - 1:b]))
            tips.append(dict(tip=name, short=short, species=man[short]["species"],
                             group=GENUS.get(short, "other Dikarya"), protein_id=r["protein_id"],
                             env=env, dom_score=r["best_dom_score"]))
    with open(TR / "duf3435_env.fa", "w") as fo:
        for n, s in recs:
            fo.write(f">{n}\n{s}\n")
    print(f"DUF3435 envelopes: {len(recs)} from {len({t['short'] for t in tips})} proteomes "
          f"{dict(Counter(t['group'] for t in tips))}", file=sys.stderr)
    hmm = TR / "DUF3435.hmm"
    subprocess.run([str(B / "hmmfetch"), "-o", str(hmm), "/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm",
                    "DUF3435"], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([str(B / "hmmalign"), "--trim", "--amino", "-o", str(TR / "duf3435.trim.sto"), str(hmm),
                    str(TR / "duf3435_env.fa")], check=True)
    subprocess.run([str(B / "esl-reformat"), "-o", str(TR / "duf3435.trim.afa"), "afa", str(TR / "duf3435.trim.sto")],
                   check=True)
    with open(TR / "duf3435.fasttree.nwk", "w") as fo, open(TR / "fasttree.log", "w") as fl:
        subprocess.run([str(B / "FastTree"), "-lg", "-gamma", str(TR / "duf3435.trim.afa")], check=True,
                       stdout=fo, stderr=fl)

    tree = Phylo.read(TR / "duf3435.fasttree.nwk", "newick")
    tree.root_at_midpoint()
    Phylo.write(tree, TR / "duf3435.fasttree.midpoint.nwk", "newick")
    grp = {t["tip"]: t["group"] for t in tips}
    for cl in tree.find_clades(order="postorder"):
        cl.tipset = {cl.name} if cl.is_terminal() else set().union(*(c.tipset for c in cl.clades))
    parent = {}
    for cl in tree.find_clades():
        for c in cl.clades:
            parent[id(c)] = cl
    rows = []
    summ = []
    for g in ("A. fumigatus Af293", "Coccidioides"):
        pure = []
        for cl in tree.find_clades():
            if all(grp[t] == g for t in cl.tipset):
                p = parent.get(id(cl))
                if p is None or not all(grp[t] == g for t in p.tipset):
                    pure.append(cl)
        for i, cl in enumerate(pure, 1):
            rows.append(dict(group=g, clade=i, n_tips=len(cl.tipset),
                             support=(cl.confidence if cl.confidence is not None else ""),
                             tips=",".join(sorted(cl.tipset))))
        n_g = sum(1 for v in grp.values() if v == g)
        summ.append(dict(group=g, n_tips=n_g, n_maximal_single_group_clades=len(pure),
                         clade_sizes=",".join(str(len(c.tipset)) for c in sorted(pure, key=lambda c: -len(c.tipset)))))
    # mixed nodes: clades with Af293 AND Coccidioides tips, support >= 0.9, no larger requirement
    mixed = [cl for cl in tree.find_clades() if not cl.is_terminal()
             and {"A. fumigatus Af293", "Coccidioides"} <= {grp[t] for t in cl.tipset}]
    mixed_min = [cl for cl in mixed if not any(
        c is not cl and c.tipset < cl.tipset and {"A. fumigatus Af293", "Coccidioides"} <= {grp[t] for t in c.tipset}
        for c in mixed)]
    for cl in mixed_min:
        cnt = Counter(grp[t] for t in cl.tipset)
        rows.append(dict(group="smallest_mixed_Af293_Cocci", clade="", n_tips=len(cl.tipset),
                         support=(cl.confidence if cl.confidence is not None else ""),
                         composition=",".join(f"{k}:{v}" for k, v in cnt.items()),
                         tips=",".join(sorted(cl.tipset))))
    summ.append(dict(group="smallest clades with both Af293 and Coccidioides tips", n_tips="",
                     n_maximal_single_group_clades=len(mixed_min),
                     clade_sizes=",".join(str(len(c.tipset)) for c in mixed_min),
                     supports=",".join(str(c.confidence) for c in mixed_min)))
    # nearest neighbour tip (patristic) for Af293 and Coccidioides tips
    term = tree.get_terminals()
    for t in tips:
        if t["group"] == "other Dikarya":
            continue
        me = next(x for x in term if x.name == t["tip"])
        best = min((x for x in term if x is not me), key=lambda x: tree.distance(me, x))
        t["nearest_tip"] = best.name
        t["nearest_tip_group"] = grp[best.name]
        t["nearest_distance"] = round(tree.distance(me, best), 4)
    for path, rr in ((TR / "duf3435_tips.tsv", tips), (TR / "duf3435_clades.tsv", rows),
                     (TR / "duf3435_tree_summary.tsv", summ)):
        hdr = []
        for r in rr:
            for k in r:
                if k not in hdr:
                    hdr.append(k)
        with open(path, "w") as fo:
            fo.write("\t".join(hdr) + "\n")
            for r in rr:
                fo.write("\t".join(str(r.get(h, "")) for h in hdr) + "\n")


if __name__ == "__main__":
    main()
