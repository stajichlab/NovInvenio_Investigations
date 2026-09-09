#!/usr/bin/env python3
"""Merge UniProt-derived annotation (bin/extract_dat_annotations.py output) into an
nf_NovInvenio presence_matrix.tsv or novelties.<Short>.tsv, keyed by protein_id.

This replaces nf_NovInvenio's own ANNOTATE_MATRIX (Pfam hmmscan + SwissProt diamond)
for UniProt-sourced studies: candidates here already carry real UniProt accessions,
so their DR Pfam/InterPro/GO cross-references are looked up directly rather than
recomputed -- see studies/fungi/pezizo_set1/run.sh's header comment for why running
a fresh hmmscan/diamond pass on top would be redundant, and could disagree with
UniProt's own calls.

Adds new columns (prefixed uniprot_ to stay distinct from any gene_name/Pfam_Names
columns nf_NovInvenio's own ANNOTATE_MATRIX might have added, e.g. via
--modelorgs_config, which this script leaves untouched): uniprot_gene_name,
uniprot_description, uniprot_pfam_ids, uniprot_pfam_names, uniprot_interpro_ids,
uniprot_go_ids, uniprot_ec_numbers, uniprot_alphafold_id, uniprot_xrefs.
A protein_id with no matching
UniProt annotation (shouldn't happen for a study built entirely from UniProt
proteomes, but --lenient allows it) gets empty strings, not a hard error.

Input matrix must have a protein_id column -- both presence_matrix.tsv and
novelties.<Short>.tsv qualify. Writes a new file; never modifies the input in place.
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from uniprot_ids import bare_accession  # noqa: E402


def load_annotations(paths: list[Path]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for p in paths:
        with open(p, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                out[row["accession"]] = row
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True, type=Path, help="presence_matrix.tsv or novelties.<Short>.tsv (must have a protein_id column)")
    ap.add_argument("--annotations", nargs="+", required=True, type=Path, help="One or more extract_dat_annotations.py TSVs")
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--lenient", action="store_true", help="Don't error on a protein_id with no UniProt annotation (default: error, since that indicates a non-UniProt input slipped in)")
    args = ap.parse_args()

    annot = load_annotations(args.annotations)

    with open(args.matrix, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "protein_id" not in reader.fieldnames:
            sys.exit(f"ERROR: {args.matrix} has no protein_id column -- fieldnames: {reader.fieldnames}")
        rows = list(reader)
        fieldnames = list(reader.fieldnames)

    new_cols = [
        "uniprot_gene_name", "uniprot_description", "uniprot_pfam_ids", "uniprot_pfam_names",
        "uniprot_interpro_ids", "uniprot_go_ids", "uniprot_ec_numbers", "uniprot_alphafold_id",
        "uniprot_xrefs",
    ]
    for c in new_cols:
        if c not in fieldnames:
            fieldnames.append(c)

    n_missing = 0
    for row in rows:
        a = annot.get(bare_accession(row["protein_id"]))
        if a is None:
            n_missing += 1
            for c in new_cols:
                row[c] = ""
            continue
        row["uniprot_gene_name"] = a["gene_name"]
        row["uniprot_description"] = a.get("description", "")
        row["uniprot_pfam_ids"] = a["pfam_ids"]
        row["uniprot_pfam_names"] = a.get("pfam_names", "")
        row["uniprot_interpro_ids"] = a["interpro_ids"]
        row["uniprot_go_ids"] = a["go_ids"]
        row["uniprot_ec_numbers"] = a.get("ec_numbers", "")
        row["uniprot_alphafold_id"] = a.get("alphafold_id", "")
        row["uniprot_xrefs"] = a.get("xrefs", "")

    if n_missing and not args.lenient:
        sys.exit(
            f"ERROR: {n_missing}/{len(rows)} protein_id values had no matching UniProt "
            f"annotation -- pass --lenient if this is expected (e.g. a mixed-source "
            f"study), otherwise this indicates the wrong --annotations were given."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {args.output} ({len(rows)} rows, {n_missing} without UniProt annotation)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
