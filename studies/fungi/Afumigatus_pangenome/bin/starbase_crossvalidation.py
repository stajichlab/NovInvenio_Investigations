#!/usr/bin/env python3
"""Cross-validates this study's own DUF3435 (Starship captain gene) hmmsearch
screen against starbase's independently-curated Starship annotations for
*A. fumigatus* Af293 -- the reference strain both datasets can be compared
on at exact base-pair resolution, since this study's own Af293 genome
(GCA_000002655.1) and starbase's Af293 records both use the same real
GenBank chromosome accessions (CM000169.1-CM000176.1).

Design: starbase's own protein set for Af293 (its own `ships`/`gff` tables)
is NOT used to identify captain gene IDs directly, because this study's
DUF3435 screen was run against the UniProt reference proteome
(UP000002530), a DIFFERENT ID space than starbase's genomic gff annotation
or this study's own NCBI-GFF3-derived gene_positions.tsv. Bridging that gap
needs a real position lookup, not an ID join: this study's own DUF3435-hit
protein sequences are tblastn'd against the Af293 genome (see
Usage step 1) to get real genomic coordinates, which is then compared
directly against starbase's curated Starship boundaries (span containment,
not ID matching).

Usage (three real steps, not all inside this script):
  1. Extract this study's own DUF3435 hit protein sequences for Af293 and
     tblastn them against the Af293 genome (data_dir/dna/*330879*.dna.fa)
     with -outfmt "6 qseqid sseqid pident length evalue bitscore qstart
     qend sstart send qlen" -- see PANGENOME_CLUSTER_PROFILE_NOTES.md's
     2026-09-17 entry for the exact commands used.
  2. starbase_crossvalidation.py --tblastn af293_captain_positions.tsv \\
         --starbase_db data/starbase/starbase_v0.1.0-pre.sqlite \\
         --output starbase_crossvalidation.tsv
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
import sys

# Af293's own chromosome order (chromosome 1-8 -> real GenBank accession),
# confirmed against this study's own data_dir/dna/*330879*.dna.fa FASTA
# headers -- starbase's contigID field uses this real accession directly
# for some curation passes ("aspfum5_CM000171.1") but a human-readable
# "ChrN...Af293" label for others ("aspfum3_Chr3AfumigatusAf293"); both
# need to resolve to the same chromosome accession to be comparable.
AF293_CHR_TO_ACCESSION = {
    1: "CM000169.1", 2: "CM000170.1", 3: "CM000171.1", 4: "CM000172.1",
    5: "CM000173.1", 6: "CM000174.1", 7: "CM000175.1", 8: "CM000176.1",
}


def _normalize_af293_contig(contig_id: str) -> str | None:
    """Resolves a starbase contigID to a real Af293 chromosome accession,
    or None if it can't be resolved (e.g. an empty/malformed field)."""
    if not contig_id:
        return None
    for accession in AF293_CHR_TO_ACCESSION.values():
        if accession in contig_id:
            return accession
    import re
    m = re.search(r"Chr(\d+)", contig_id)
    if m:
        return AF293_CHR_TO_ACCESSION.get(int(m.group(1)))
    return None


def parse_tblastn_spans(lines: list[str], max_gap: int = 20_000) -> dict[str, dict]:
    """Parses tblastn -outfmt "6 qseqid sseqid pident length evalue bitscore
    qstart qend sstart send qlen" lines into {qseqid: {contig, start, end}}
    -- one genomic span per query.

    A query commonly has HSPs on more than one contig (paralogs elsewhere in
    the genome) or even multiple, spatially SEPARATE clusters of HSPs on the
    SAME contig (a real, expected pattern for Starship-associated tyrosine
    recombinases, which occur as multiple divergent copies genome-wide --
    confirmed on real data 2026-09-17: one query had a tight, >80%-identity
    HSP cluster at one locus and a second, 42-49%-identity cluster ~445kb
    away on the same contig -- two genuinely different loci, not one wide
    locus). Taking a raw min/max across ALL same-contig HSPs would silently
    merge those into one bogus 450kb span. Instead: HSPs within `max_gap` of
    each other are merged into one cluster, and only the cluster CONTAINING
    the single lowest-evalue HSP across the whole query is kept as that
    query's span -- the other cluster(s), if any, represent a different
    (unresolved, lower-confidence) locus and are dropped, not merged in."""
    hsps: dict[str, list[dict]] = {}
    for line in lines:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        qseqid, sseqid = parts[0], parts[1]
        evalue = float(parts[4])
        sstart, send = int(parts[8]), int(parts[9])
        hsps.setdefault(qseqid, []).append({
            "contig": sseqid, "evalue": evalue,
            "start": min(sstart, send), "end": max(sstart, send),
        })

    spans = {}
    for qseqid, rows in hsps.items():
        best = min(rows, key=lambda r: r["evalue"])
        same_contig = sorted(
            (r for r in rows if r["contig"] == best["contig"]), key=lambda r: r["start"]
        )
        # Merge same_contig HSPs into gap-separated clusters, then keep the
        # cluster that contains `best`.
        clusters: list[list[dict]] = []
        for r in same_contig:
            if clusters and r["start"] - clusters[-1][-1]["end"] <= max_gap:
                clusters[-1].append(r)
            else:
                clusters.append([r])
        chosen = next(c for c in clusters if best in c)
        spans[qseqid] = {
            "contig": best["contig"],
            "start": min(r["start"] for r in chosen),
            "end": max(r["end"] for r in chosen),
        }
    return spans


def load_starbase_af293_starships(db_path: str, af293_taxonomy_ids: list[int]) -> list[dict]:
    """Real, deduplicated Af293 Starship boundaries from starbase's
    starship_features table. Multiple curation-pass rows per starshipID
    (varying by a few bp) are collapsed to ONE widest span per starshipID
    -- the boundaryType='flank' row is preferred when present (it's
    starbase's own best-boundary-evidence call), else the widest span
    among all rows for that starshipID."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    placeholders = ",".join("?" * len(af293_taxonomy_ids))
    rows = conn.execute(f"""
        SELECT sf.starshipID, sf.contigID, sf.elementBegin, sf.elementEnd,
               sf.boundaryType, sf.captainID, fn.familyName
        FROM starship_features sf
        JOIN joined_ships js ON sf.ship_id = js.ship_id
        LEFT JOIN family_names fn ON js.ship_family_id = fn.id
        WHERE js.tax_id IN ({placeholders}) AND sf.elementBegin != ''
    """, af293_taxonomy_ids).fetchall()
    conn.close()

    by_starship: dict[str, list[sqlite3.Row]] = {}
    for row in rows:
        by_starship.setdefault(row["starshipID"], []).append(row)

    starships = []
    for starship_id, group in by_starship.items():
        flank_rows = [r for r in group if r["boundaryType"] == "flank"]
        chosen = flank_rows if flank_rows else group
        contig = _normalize_af293_contig(chosen[0]["contigID"])
        if contig is None:
            continue
        starships.append({
            "starshipID": starship_id,
            "contig": contig,
            "start": min(int(r["elementBegin"]) for r in chosen),
            "end": max(int(r["elementEnd"]) for r in chosen),
            "captainID": chosen[0]["captainID"],
            "family": chosen[0]["familyName"] or "",
        })
    return starships


def match_captains_to_starships(
    captain_spans: dict[str, dict], starships: list[dict]
) -> dict[str, dict]:
    """{query: {matched_starship, family}} -- a captain "matches" a starbase
    Starship when the captain's own genomic span falls entirely within (or
    overlaps substantially with) the Starship's curated boundary on the
    SAME contig. Containment, not exact-coordinate equality, since a
    protein's coding span is always narrower than the full mobile
    element's boundary."""
    results = {}
    for query, span in captain_spans.items():
        matched = None
        family = ""
        for s in starships:
            if span["contig"] != s["contig"]:
                continue
            if span["start"] >= s["start"] and span["end"] <= s["end"]:
                matched = s["starshipID"]
                family = s["family"]
                break
        results[query] = {"matched_starship": matched, "family": family}
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tblastn", required=True, help="tblastn -outfmt 6 (see module docstring for columns)")
    ap.add_argument("--starbase_db", required=True)
    ap.add_argument("--af293_taxonomy_ids", default="23,172,371",
                     help="comma-separated taxonomy.id values matching Af293 in starbase "
                          "(multiple curation-source rows exist for the same strain)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open(args.tblastn) as fh:
        lines = fh.readlines()
    captain_spans = parse_tblastn_spans(lines)
    print(f"starbase_crossvalidation: {len(captain_spans)} captain-hit queries with a genomic span",
          file=sys.stderr)

    tax_ids = [int(t) for t in args.af293_taxonomy_ids.split(",")]
    starships = load_starbase_af293_starships(args.starbase_db, tax_ids)
    print(f"starbase_crossvalidation: {len(starships)} curated Af293 Starships in starbase",
          file=sys.stderr)

    matches = match_captains_to_starships(captain_spans, starships)
    n_matched = sum(1 for m in matches.values() if m["matched_starship"])
    matched_ids = {m["matched_starship"] for m in matches.values() if m["matched_starship"]}
    n_starbase_recovered = len(matched_ids)
    print(f"starbase_crossvalidation: {n_matched}/{len(matches)} captain hits fall within a "
          f"curated starbase Starship boundary; {n_starbase_recovered}/{len(starships)} starbase "
          f"Starships recovered by this study's own screen", file=sys.stderr)

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["query", "contig", "start", "end", "matched_starbase_starship", "starbase_family"])
        for query, span in captain_spans.items():
            m = matches[query]
            writer.writerow([query, span["contig"], span["start"], span["end"],
                              m["matched_starship"] or "", m["family"]])
        writer.writerow([])
        for s in starships:
            recovered = "Y" if s["starshipID"] in matched_ids else "N"
            writer.writerow([f"# starbase_starship:{s['starshipID']}", s["contig"], s["start"],
                              s["end"], f"recovered_by_screen={recovered}", s["family"]])


if __name__ == "__main__":
    main()
