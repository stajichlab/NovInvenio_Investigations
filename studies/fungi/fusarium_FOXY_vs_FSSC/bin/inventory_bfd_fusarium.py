#!/usr/bin/env python3
"""Inventory BFD funannotate annotations for Fusarium / Neocosmospora strains.

Join rule: an annotation directory is matched to its samples.csv row by the
LOCUSTAG prefix of its first protein ID (funannotate names every protein
<LOCUSTAG>_NNNNNN-T1). This avoids the directory-name ambiguity: BFD names a
directory makeSampleTag(SPECIES, STRAIN), and a strain-less row gets the bare
species name (results/_ambiguous_basenames_report.csv.gz lists 812 such
"bare-fallback" Fusarium_oxysporum rows).

For each directory it records: proteins.fa / gff3 / scaffolds.fa presence and
size, the matched samples.csv row, whether the directory name equals the
expected makeSampleTag() value, the assembly stats (N50, total length,
contigs) from results/genome_stats/asm_stats and the BUSCO protein
completeness from results/genome_stats/BUSCO_protein.

Output: selection/bfd_fusarium_inventory.tsv (one row per annotation dir).
Hardcoded paths, per the NII "study-specific scripts" convention.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

BFD = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD")
ANN = BFD / "genome_annotation"
STATS = BFD / "results" / "genome_stats"
STUDY = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/fusarium_FOXY_vs_FSSC")
OUT = STUDY / "selection" / "bfd_fusarium_inventory.tsv"


def clean_strain(raw: str) -> str:
    s = (raw or "").strip().replace("'", "").replace('"', "").split(";")[0].strip().replace(":", " ")
    s = re.sub(r"^\s*\*+", "", s)
    s = re.sub(r"\*+\s*$", "", s)
    s = re.sub(r"\s*\*+\s*", "-", s)
    return s.strip()


def make_sample_tag(species: str, strain: str) -> str:
    sp = (species or "").strip().replace("'", "").replace('"', "")
    st = clean_strain(strain)
    return re.sub(r"[\s/#\[\]?{}]+", "_", "_".join(x for x in (sp, st) if x))


def first_locustag(pep: Path) -> str:
    with open(pep) as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].split("_")[0]
    return ""


def index_stats() -> tuple[dict, dict]:
    asm = {}
    for f in (STATS / "asm_stats").glob("*/*.stats.txt"):
        asm[f.name[: -len(".stats.txt")]] = f
    busco = {}
    for f in (STATS / "BUSCO_protein").glob("*/*.BUSCO_summary.*.txt"):
        busco[f.name.split(".BUSCO_summary")[0]] = f
    return asm, busco


def parse_asm(f: Path) -> dict:
    out = {}
    for line in f.read_text().splitlines():
        m = re.match(r"\s*([A-Z0-9% ]+?)\s*=\s*(\S+)", line)
        if m:
            out[m.group(1).strip()] = m.group(2)
    return {"n50": out.get("N50", ""), "total_len": out.get("TOTAL LENGTH", ""),
            "n_contigs": out.get("CONTIG COUNT", "")}


def parse_busco(f: Path) -> str:
    m = re.search(r"C:([\d.]+)%", f.read_text())
    return m.group(1) if m else ""


def main() -> int:
    rows_by_lt: dict[str, list[dict]] = {}
    with open(BFD / "samples.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            rows_by_lt.setdefault(r["LOCUSTAG"].strip(), []).append(r)
    asm, busco = index_stats()

    dirs = sorted(d for d in ANN.iterdir()
                  if d.is_dir() and (d.name.startswith("Fusarium_") or d.name.startswith("Neocosmospora_")))
    fields = ["dirname", "locustag", "asmid", "species", "species_in", "strain", "ncbi_taxid",
              "bioproject", "expected_tag", "tag_match", "n_samples_rows_for_locustag",
              "proteins", "gff3", "scaffolds", "proteins_bytes", "gff3_bytes", "scaffolds_bytes",
              "n50", "total_len", "n_contigs", "busco_protein_C", "status"]
    out_rows = []
    for d in dirs:
        pr = d / "predict_results"
        pep = next(iter(sorted(pr.glob("*.proteins.fa"))), None) if pr.is_dir() else None
        gff = next(iter(sorted(pr.glob("*.gff3"))), None) if pr.is_dir() else None
        dna = next(iter(sorted(pr.glob("*.scaffolds.fa"))), None) if pr.is_dir() else None
        rec = {k: "" for k in fields}
        rec["dirname"] = d.name
        for key, p in (("proteins", pep), ("gff3", gff), ("scaffolds", dna)):
            rec[key] = str(p) if p else ""
            rec[f"{key}_bytes"] = str(p.stat().st_size) if p else "0"
        status = []
        if not pep or pep.stat().st_size == 0:
            status.append("no_proteins")
        if not gff or gff.stat().st_size == 0:
            status.append("no_gff3")
        if not dna or dna.stat().st_size == 0:
            status.append("no_genome")
        if pep and pep.stat().st_size > 0:
            lt = first_locustag(pep)
            rec["locustag"] = lt
            cand = rows_by_lt.get(lt, [])
            rec["n_samples_rows_for_locustag"] = str(len(cand))
            if len(cand) == 1:
                r = cand[0]
                rec.update(asmid=r["ASMID"], species=r["SPECIES"], species_in=r["SPECIES_IN"],
                           strain=r["STRAIN"], ncbi_taxid=r["NCBI_TAXONID"], bioproject=r["BIOPROJECT"])
                rec["expected_tag"] = make_sample_tag(r["SPECIES"], r["STRAIN"])
                rec["tag_match"] = "Y" if rec["expected_tag"] == d.name else "N"
                if r["ASMID"] in asm:
                    rec.update(parse_asm(asm[r["ASMID"]]))
                if lt in busco:
                    rec["busco_protein_C"] = parse_busco(busco[lt])
            elif not cand:
                status.append("locustag_not_in_samples")
            else:
                status.append("locustag_multi_rows")
        rec["status"] = ";".join(status) or "ok"
        out_rows.append(rec)

    # a locustag seen in two directories is ambiguous (stale re-annotation)
    seen: dict[str, int] = {}
    for r in out_rows:
        if r["locustag"]:
            seen[r["locustag"]] = seen.get(r["locustag"], 0) + 1
    for r in out_rows:
        if r["locustag"] and seen[r["locustag"]] > 1:
            if r["tag_match"] != "Y":
                prev = [] if r["status"] == "ok" else [r["status"]]
                r["status"] = ";".join(prev + ["dup_locustag_nonmatching_dir"])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    n_ok = sum(1 for r in out_rows if r["status"] == "ok")
    print(f"wrote {OUT}: {len(out_rows)} dirs, {n_ok} status=ok", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
