#!/usr/bin/env python3
"""Annotation coverage of IN proteomes and novelty candidates, per NII study.

Uses the same lookup as studies/fungi/zoosporic_dikarya/analysis/enrichment/
enrich_novelties_v2.py: annotations/<name>.tsv.gz, <name> from lib/annotation_names.py
(species.csv rows),
matched by bare UniProt accession (tr|ACC|NAME -> ACC). "Annotated" = any of
go_ids, pfam_ids, interpro_ids non-empty.

Output: one row per (study, IN species) plus one ALL_IN row per study.
Run: .pixi/envs/default/bin/python3 -I coverage_all_studies.py --out coverage_all_studies.tsv STUDY...
STUDY = <domain>/<name>; results are read from results/<name>.
"""
import argparse
import csv
import gzip
import sys
from pathlib import Path

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII / "lib"))
from annotation_names import annotation_names, annotation_path  # noqa: E402
csv.field_size_limit(10**9)


def bare_acc(pid):
    parts = pid.split("|")
    if len(parts) >= 2 and parts[0] in ("tr", "sp"):
        return parts[1]
    return pid


def fasta_ids(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("studies", nargs="+")
    a = ap.parse_args()
    rows = []
    for st in a.studies:
        sd = NII / "studies" / st
        name = Path(st).name
        rd = NII / "results" / name
        cfg = list(csv.DictReader(open(sd / "config.csv")))
        sp_list = list(csv.DictReader(open(sd / "species.csv")))
        sp_rows = {r["Short"]: r for r in sp_list}
        ann_name = annotation_names(sp_list)
        in_sp = [r["Short"] for r in cfg if r["GROUP"] == "IN"]
        pep = {r["Short"]: r["Protein"] for r in cfg}
        cands = {l.strip() for l in open(rd / "candidates.txt") if l.strip()}
        has_matrix = (rd / "presence_matrix.uniprot.tsv").exists()
        tot = dict(n_prot=0, n_ann=0, n_cand=0, n_cand_ann=0)
        n_sp_ann = 0
        for sp in in_sp:
            r = sp_rows.get(sp, {})
            acc = r.get("Protein_Accession") or r.get("Accession") or ""
            src = r.get("Protein_Source") or r.get("Source") or ""
            af = (annotation_path(sd / "annotations", ann_name[sp]) if sp in ann_name
                  else sd / "annotations" / "NONE.tsv.gz")
            ann = set()
            if acc and af.exists():
                with gzip.open(af, "rt") as fh:
                    for row in csv.DictReader(fh, delimiter="\t"):
                        if row.get("go_ids") or row.get("pfam_ids") or row.get("interpro_ids"):
                            ann.add(row["accession"])
            fa = sd / "data_dir" / "pep" / pep[sp]
            ids = fasta_ids(fa) if fa.exists() else []
            n_ann = sum(1 for i in ids if bare_acc(i) in ann)
            sc = [c for c in cands if c.split("::", 1)[0] == sp]
            n_cann = sum(1 for c in sc if bare_acc(c.split("::", 1)[1]) in ann)
            if n_ann > 0:
                n_sp_ann += 1
            rows.append(dict(study=st, species=sp, protein_source=src, protein_accession=acc,
                             annotation_file=af.name if af.exists() else "NONE",
                             n_proteins=len(ids), n_annotated=n_ann,
                             frac_annotated=round(n_ann / len(ids), 4) if ids else 0,
                             n_candidates=len(sc), n_candidates_annotated=n_cann,
                             frac_candidates_annotated=round(n_cann / len(sc), 4) if sc else 0,
                             presence_matrix_uniprot=has_matrix, n_in_species_annotated=""))
            for k, v in (("n_prot", len(ids)), ("n_ann", n_ann), ("n_cand", len(sc)), ("n_cand_ann", n_cann)):
                tot[k] += v
        n_other = len(cands) - tot["n_cand"]
        rows.append(dict(study=st, species="ALL_IN", protein_source="", protein_accession="",
                         annotation_file="", n_proteins=tot["n_prot"], n_annotated=tot["n_ann"],
                         frac_annotated=round(tot["n_ann"] / tot["n_prot"], 4) if tot["n_prot"] else 0,
                         n_candidates=len(cands), n_candidates_annotated=tot["n_cand_ann"],
                         frac_candidates_annotated=round(tot["n_cand_ann"] / len(cands), 4) if cands else 0,
                         presence_matrix_uniprot=has_matrix,
                         n_in_species_annotated=f"{n_sp_ann}/{len(in_sp)}"
                         + (f" (candidates not matched to an IN species: {n_other})" if n_other else "")))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
