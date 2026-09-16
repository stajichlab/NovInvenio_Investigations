#!/usr/bin/env python3
"""Build the benchmark scorecard's ground-truth inputs
(notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md,
component 6) from `mbio.01092-25-s0002.xlsx`, in the exact shape
`bin/benchmark_scorecard.py` and `bin/id_crosswalk.py` consume:

1. `ground_truth_starships_by_short.tsv` (positive control, presence
   recovery): Table S6 (population freq) + Table S21 (presence/absence
   genotyping) re-keyed from the paper's own `isolateID`/`originalID` to
   THIS study's own `Short` strain names, via a normalized-alnum string
   match against `config.csv` -- the same method used for the 2026-09-13
   TaxonGroup fill (254/293 IN strains matched), reproduced here as real
   code (that fill's own matching logic was never checked in as a
   standalone script -- see PANGENOME_CLUSTER_PROFILE_NOTES.md's "ID
   crosswalk" open item).
2. `ground_truth_cargo_by_nameid.tsv` (positive control, cargo grouping):
   Table S13's (starshipID, geneID) rows re-grouped by `nameID` (the named
   Starship, e.g. "Gnosis-h1") rather than `starshipID` (one per strain
   instance) -- `parse_starship_supplement.cargo_gene_sets` groups by
   starshipID because that's what Table S13 keys BLAST-recovery rows on,
   but the presence/absence ground truth (Table S21) and frequency table
   (Table S6) are keyed by nameID. A "does clustering recover this named
   Starship's cargo as one family" question needs the union of every
   instance's cargo genes under that one name, not left split by strain.
3. `negative_control_conserved_genes.tsv` (negative control): Table S19's
   secondary-metabolite/virulence gene catalog (802 genes across ~30 named
   BGC/Cluster_N groups plus independently-characterized virulence loci --
   NONE of them Starship/cargo-annotated in the paper) with a real NCBI
   accession, used as a "should NOT look variably mobile" check: unlike
   Starship cargo, these are canonical core-genome secondary-metabolite
   clusters, so the expectation is high recovered presence (core/soft-core
   frequency), not a Starship-style presence/absence split.

Usage:
  build_ground_truth_tables.py --xlsx mbio.01092-25-s0002.xlsx \\
      --config config.csv --output_prefix results/id_crosswalk/ground_truth
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from parse_starship_supplement import high_confidence_starships  # noqa: E402


def normalize_id(value: object) -> str:
    """Alnum-only, lowercased -- the same normalization
    PANGENOME_CLUSTER_PROFILE_NOTES.md's TaxonGroup fill used for matching
    Table S21 isolateID/originalID against config.csv Strain/Short."""
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def build_short_lookup(config: pd.DataFrame) -> dict[str, str]:
    """{normalized(Short or Strain): Short} -- both spellings point at the
    same Short so a match on either side succeeds."""
    lookup: dict[str, str] = {}
    for _, row in config.iterrows():
        short = row["Short"]
        lookup[normalize_id(short)] = short
        if pd.notna(row.get("Strain")):
            lookup[normalize_id(row["Strain"])] = short
    return lookup


def match_isolates_to_short(
    isolate_ids: list[str], original_ids: list[str], short_lookup: dict[str, str]
) -> dict[str, str]:
    """{isolateID: Short} for every isolateID (falling back to originalID)
    that normalized-matches a Short/Strain in this study's config.csv."""
    matched: dict[str, str] = {}
    for isolate_id, original_id in zip(isolate_ids, original_ids):
        for candidate in (isolate_id, original_id):
            key = normalize_id(candidate)
            if key in short_lookup:
                matched[isolate_id] = short_lookup[key]
                break
    return matched


_AF293_XP_RE = re.compile(r"^AF293_(XP-[\d.]+)$")


def normalize_geneid_for_crosswalk(gene_id: str) -> str:
    """Table S13's AF293 rows spell a real NCBI RefSeq accession as
    `AF293_XP-<digits>.<version>` (paper-internal strain-prefix convention),
    not the plain `XP_<digits>.<version>` accession the crosswalk (built by
    fetching these same accessions from NCBI and diamond-blastp'ing them
    against this study's proteome) actually keys on. Every OTHER strain's
    geneID (e.g. `47-10_000766`) is the paper's own internal locus numbering
    with no public accession at all -- left untouched, since there is
    nothing to normalize it to (see PANGENOME_CLUSTER_PROFILE_NOTES.md's
    "ID crosswalk" section for why only the AF293-prefixed rows are
    sequence-crosswalkable in this benchmark)."""
    m = _AF293_XP_RE.match(gene_id)
    if m:
        return m.group(1).replace("-", "_", 1)
    return gene_id


def cargo_gene_sets_by_nameid(table_s13: pd.DataFrame) -> dict[str, set[str]]:
    """Union Table S13's cargo gene sets across every starshipID instance
    sharing the same nameID -- see module docstring point 2. AF293-prefixed
    geneIDs are normalized to their plain NCBI accession form so they match
    the crosswalk's keys (see `normalize_geneid_for_crosswalk`)."""
    result: dict[str, set[str]] = {}
    for name_id, gene_id in zip(table_s13["nameID"], table_s13["geneID"]):
        result.setdefault(name_id, set()).add(normalize_geneid_for_crosswalk(gene_id))
    return result


def conserved_negative_control_genes(table_s19: pd.DataFrame) -> pd.DataFrame:
    """Table S19 rows with a real accession and a named SM cluster
    (BGC*/Cluster_N) or an independently-characterized virulence locus --
    the negative-control set (see module docstring point 3). Rows with no
    accession at all (curated-name-only entries, or the two literal
    "AspGD only" placeholder cells) are dropped -- they can never be
    crosswalked to a sequence regardless of control type."""
    has_accession = table_s19["accession"].notna()
    is_real_accession = table_s19["accession"].astype(str).str.match(r"^[A-Za-z]{1,3}_?\d")
    return table_s19[has_accession & is_real_accession].copy()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--output_prefix", required=True)
    args = ap.parse_args()

    xl = pd.ExcelFile(args.xlsx)
    table_s6 = xl.parse("Table S6", header=1)
    table_s21 = xl.parse("Table S21", header=1)
    table_s13 = xl.parse("Table S13", header=1)
    table_s19 = xl.parse("Table S19", header=1)

    config = pd.read_csv(args.config)
    short_lookup = build_short_lookup(config)

    # --- 1. presence recovery ground truth, re-keyed to this study's Short ---
    starships = high_confidence_starships(table_s6, table_s21)
    isolate_ids = table_s21["isolateID"].tolist()
    original_ids = table_s21["originalID"].tolist()
    isolate_to_short = match_isolates_to_short(isolate_ids, original_ids, short_lookup)
    print(
        f"Strain-overlap match: {len(set(isolate_to_short.values()))} of "
        f"{len(config)} config.csv strains matched to a Table S21 isolateID "
        f"({len(isolate_to_short)} of {len(set(isolate_ids))} isolates matched)",
        file=sys.stderr,
    )

    starships_path = f"{args.output_prefix}_starships_by_short.tsv"
    with open(starships_path, "w") as fh:
        fh.write("nameID\tpopulation_freq\tpresence_by_short_json\n")
        for name_id, info in sorted(starships.items()):
            presence_by_short = {
                isolate_to_short[iso]: present
                for iso, present in info["presence"].items()
                if iso in isolate_to_short
            }
            fh.write(f"{name_id}\t{info['population_freq']}\t{json.dumps(presence_by_short)}\n")
    print(f"Wrote {len(starships)} Starships to {starships_path}", file=sys.stderr)

    # --- 2. cargo grouping ground truth, keyed by nameID ---
    cargo_by_name = cargo_gene_sets_by_nameid(table_s13)
    cargo_path = f"{args.output_prefix}_cargo_by_nameid.tsv"
    with open(cargo_path, "w") as fh:
        fh.write("nameID\tgeneID\n")
        for name_id, genes in sorted(cargo_by_name.items()):
            for gene in sorted(genes):
                fh.write(f"{name_id}\t{gene}\n")
    print(
        f"Wrote {sum(len(g) for g in cargo_by_name.values())} cargo-gene rows "
        f"across {len(cargo_by_name)} named Starships to {cargo_path}",
        file=sys.stderr,
    )

    # --- 3. negative control: conserved, non-Starship SM/virulence genes ---
    negative = conserved_negative_control_genes(table_s19)
    negative_path = f"{args.output_prefix}_negative_control_conserved_genes.tsv"
    negative[["gene", "locus", "accession", "cluster", "predictedFunction"]].to_csv(
        negative_path, sep="\t", index=False
    )
    print(f"Wrote {len(negative)} negative-control genes to {negative_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
