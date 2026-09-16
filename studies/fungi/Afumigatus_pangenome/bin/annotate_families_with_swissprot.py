#!/usr/bin/env python3
"""SwissProt-based functional annotation for eligible (shell+cloud)
family representative sequences -- a proxy for real functional
information beyond Pfam domain composition, per the report's own
methodological discussion (README/PANGENOME_CLUSTER_PROFILE_NOTES.md).

Two sources, combined:
1. **Free, no search needed**: some family rep IDs are themselves
   UniProt-formatted (`sp|ACCESSION|NAME` reviewed Swiss-Prot, or
   `tr|ACCESSION|NAME` unreviewed TrEMBL) because the strain's own
   proteome was sourced from UniProt rather than NCBI -- the ID already
   IS the UniProt record, just needs recognizing.
2. **A real diamond blastp best-hit search** against the local, already
   diamond-indexed SwissProt database
   (NovInvenio/db/uniprot/uniprot_sprot.fasta.dmnd), for every family
   that doesn't already have a source #1 self-annotation.

Usage:
  # Step 1 (already run separately, real diamond command):
  #   diamond blastp -q all_eligible_family_reps.fa \\
  #       -d .../uniprot_sprot.fasta.dmnd \\
  #       -o eligible_vs_sprot.tsv \\
  #       --outfmt 6 qseqid sseqid stitle pident length evalue bitscore \\
  #       --max-target-seqs 1 --evalue 1e-5
  annotate_families_with_swissprot.py \\
      --diamond_tsv eligible_vs_sprot.tsv \\
      --family_ids all_eligible_ids.txt \\
      --output family_swissprot_annotation.tsv
"""
from __future__ import annotations

import argparse
import re
import sys

_SELF_UNIPROT_RE = re.compile(r"\|(sp|tr)\|([A-Za-z0-9_]+)\|([A-Za-z0-9_]+)$")


def parse_self_uniprot_id(family_id: str) -> tuple[str, str, str] | None:
    """If `family_id` itself embeds a UniProt-formatted accession (the
    strain's proteome was UniProt-sourced), return (source, accession,
    name) -- source is "sp" (reviewed Swiss-Prot) or "tr" (unreviewed
    TrEMBL). Returns None for an ordinary NCBI-style ID
    (e.g. "Asfu_X|KAH123456.1")."""
    m = _SELF_UNIPROT_RE.search(family_id)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3)


def parse_diamond_stitle(stitle: str) -> tuple[str, str, str]:
    """`stitle` column from diamond's `--outfmt 6 ... stitle`, e.g.
    "sp|D4AXL1|BLML_ARTBC Beta-lactamase-like protein ARB_00930 OS=...".
    Returns (accession, name, description) -- description is everything
    between the name and the first "OS=" (organism) field."""
    parts = stitle.split("|", 2)
    accession = parts[1] if len(parts) > 1 else ""
    rest = parts[2] if len(parts) > 2 else stitle
    name, _, description = rest.partition(" ")
    description = description.split(" OS=", 1)[0].strip()
    return accession, name, description


def load_diamond_besthits(path: str) -> dict[str, dict]:
    """{family_id: {source: "blast", accession, name, description, pident,
    evalue}} -- one row per query in the diamond output (already
    --max-target-seqs 1, so no need to pick a best hit here)."""
    hits: dict[str, dict] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            qseqid, _sseqid, stitle, pident, _length, evalue, _bitscore = parts[:7]
            accession, name, description = parse_diamond_stitle(stitle)
            hits[qseqid] = {
                "source": "blast", "accession": accession, "name": name,
                "description": description, "pident": pident, "evalue": evalue,
            }
    return hits


def annotate_families(
    family_ids: list[str], diamond_hits: dict[str, dict],
) -> dict[str, dict]:
    """{family_id: annotation dict} for every ID in `family_ids`, self-
    UniProt annotation taking priority over a blast hit when both are
    available (a self-embedded ID is exact, not an inferred homolog)."""
    annotations: dict[str, dict] = {}
    for family_id in family_ids:
        self_hit = parse_self_uniprot_id(family_id)
        if self_hit:
            source, accession, name = self_hit
            annotations[family_id] = {
                "source": "self_" + source, "accession": accession, "name": name,
                "description": "", "pident": "100.0", "evalue": "0.0",
            }
        elif family_id in diamond_hits:
            annotations[family_id] = diamond_hits[family_id]
    return annotations


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diamond_tsv", required=True)
    ap.add_argument("--family_ids", required=True, help="one family ID per line")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    diamond_hits = load_diamond_besthits(args.diamond_tsv)
    with open(args.family_ids) as fh:
        family_ids = [line.strip() for line in fh if line.strip()]

    annotations = annotate_families(family_ids, diamond_hits)
    n_self = sum(1 for a in annotations.values() if a["source"].startswith("self_"))
    n_blast = sum(1 for a in annotations.values() if a["source"] == "blast")
    print(
        f"annotate_families_with_swissprot: {len(family_ids)} families -- "
        f"{n_self} self-annotated (embedded UniProt ID), {n_blast} from diamond blastp, "
        f"{len(family_ids) - len(annotations)} unannotated", file=sys.stderr,
    )

    with open(args.output, "w") as out:
        out.write("family\tsource\taccession\tname\tdescription\tpident\tevalue\n")
        for family_id in family_ids:
            a = annotations.get(family_id)
            if a is None:
                continue
            out.write(
                f"{family_id}\t{a['source']}\t{a['accession']}\t{a['name']}\t"
                f"{a['description']}\t{a['pident']}\t{a['evalue']}\n"
            )


if __name__ == "__main__":
    main()
