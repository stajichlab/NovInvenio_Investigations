"""UniProt/NCBI dataset resolution for bin/ni resolve.

See notes/superpowers/specs/2026-09-11-ni-dataset-resolution-design.md for the
full design and the bioinformatics rationale behind each rule enforced here.

Field paths below were verified live against rest.uniprot.org on 2026-09-11
(see task-1-report.md for the full curl transcript). Two places differ from
this task's original design sketch:

- `strain` is a TOP-LEVEL field on a proteome record (`r["strain"]`), not
  nested under `taxonomy` (`r["taxonomy"]["strain"]`) as the sketch assumed.
- The taxonomy-search endpoint's `fields` query parameter must name the field
  `id` (not `taxon_id`) to select the taxon-ID column; the response JSON still
  keys the value as `taxonId`, so the result-parsing side of the sketch was
  correct -- only the request's `fields=` value was wrong.

Everything else (proteomeCompletenessReport.buscoReport.score,
genomeAssembly.assemblyId, superkingdom, the Link-header pagination shape)
matched the sketch as written.
"""
from __future__ import annotations

import csv as _csv
import json
import re
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from fetch_uniprot_proteome import DOMAIN_FOLDERS, UNIPROT_FTP_BASE  # noqa: E402

UNIPROT_REST = "https://rest.uniprot.org"

# Keep ':', '+', '(' and ')' unescaped -- UniProt's query syntax uses them
# literally (e.g. "taxonomy_id:5334 AND proteome_type:REFERENCE") and rest.
# uniprot.org accepts them unescaped in the query string.
_QUERY_SAFE_CHARS = ":+()"


def resolve_taxon_id(species: str) -> list[dict[str, str]]:
    """Resolve a species name to NCBI taxon ID(s) via UniProt's taxonomy search.

    Empty list means no resolvable taxon -- caller reports and leaves the row
    blank. A name can resolve to more than one hit (species-level plus
    strain-level taxa, or unrelated organisms sharing part of the name);
    disambiguating among them is a later task's job, not this function's.
    """
    url = (
        f"{UNIPROT_REST}/taxonomy/search"
        f"?query={urllib.parse.quote(species)}&fields=id,scientific_name&format=json"
    )
    with urllib.request.urlopen(url) as resp:
        payload = json.load(resp)
    return [
        {"taxon_id": str(r["taxonId"]), "scientific_name": r["scientificName"]}
        for r in payload.get("results", [])
    ]


def _parse_next_link(link_header: str | None) -> str | None:
    if not link_header:
        return None
    m = re.search(r'<([^>]+)>;\s*rel="next"', link_header)
    return m.group(1) if m else None


def search_uniprot_proteomes(taxon_id: str, *, reference_only: bool) -> list[dict[str, Any]]:
    """Search UniProt proteomes for a taxon ID.

    reference_only=True restricts to proteome_type:REFERENCE (Step 1); False
    searches all proteomes for this taxon (Step 1b, strain-targeted
    widening). Paginates internally by following the response's Link header
    (rel="next") so callers never see a truncated result.
    """
    query = f"taxonomy_id:{taxon_id}"
    if reference_only:
        query += " AND proteome_type:REFERENCE"
    url = (
        f"{UNIPROT_REST}/proteomes/search"
        f"?query={urllib.parse.quote(query, safe=_QUERY_SAFE_CHARS)}&format=json&size=500"
    )

    candidates: list[dict[str, Any]] = []
    while url:
        with urllib.request.urlopen(url) as resp:
            payload = json.load(resp)
            next_url = _parse_next_link(resp.headers.get("Link"))
        for r in payload.get("results", []):
            busco = (
                r.get("proteomeCompletenessReport", {})
                .get("buscoReport", {})
                .get("score")
            )
            candidates.append({
                "proteome_id": r["id"],
                "strain": r.get("strain"),
                "busco_score": busco,
                "genome_assembly_id": r.get("genomeAssembly", {}).get("assemblyId"),
                "superkingdom": r.get("superkingdom", "Eukaryota"),
            })
        url = next_url
    return candidates


def check_ftp_available(proteome_id: str, superkingdom: str) -> bool:
    """HEAD-equivalent check that this proteome's FTP directory exists in the
    current release snapshot (the live REST API can list a proteome not yet
    mirrored to FTP -- confirmed live for Scom during this spec's review)."""
    domain_folder = DOMAIN_FOLDERS.get(superkingdom.capitalize(), "Eukaryota")
    url = f"{UNIPROT_FTP_BASE}/{domain_folder}/{proteome_id}/"
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req):
            return True
    except urllib.error.HTTPError:
        return False


def strain_matches(candidate_strain: str | None, wanted_strain: str) -> bool:
    if not candidate_strain or not wanted_strain:
        return False
    return wanted_strain.lower() in candidate_strain.lower()


def query_ncbi_assemblies(taxon_id: str, *, require_reference: bool) -> list[dict[str, Any]]:
    """Query NCBI Datasets for every assembly of this taxon. require_reference=True
    passes --reference (Step 2's default path); False queries all assemblies
    (Step 2b, strain-targeted widening).

    Field paths verified live on 2026-09-11 against
    `datasets summary genome taxon 5334 --as-json-lines` (Schizophyllum
    commune, which has a GCA/GCF reference-genome pair). Two paths differ
    from this task's brief sketch:

    - `paired_accession` is a TOP-LEVEL field on the record (`r["paired_accession"]`),
      not nested under `assembly_info`. `assembly_info.paired_accession` does
      not exist; there's a differently-shaped `assembly_info.paired_assembly`
      object instead (with its own nested `accession`/`status`), which we
      don't use here.
    - There is no `assembly_info.submission_date` field. The closest analog
      present on every record is `assembly_info.release_date` (the date the
      assembly was released to the public database), used here for the
      "submission_date" key. A `submission_date` field does exist, but only
      nested under `assembly_info.biosample.submission_date`, which
      describes the biological sample, not the assembly record.
    """
    cmd = ["datasets", "summary", "genome", "taxon", taxon_id, "--as-json-lines"]
    if require_reference:
        cmd.append("--reference")
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    candidates = []
    for line in result.stdout.strip().splitlines():
        if not line:
            continue
        r = json.loads(line)
        info = r.get("assembly_info", {})
        candidates.append({
            "accession": r["accession"],
            "paired_accession": r.get("paired_accession"),
            "refseq_category": info.get("refseq_category"),
            "assembly_level": info.get("assembly_level"),
            "assembly_status": info.get("assembly_status"),
            "strain": r.get("organism", {}).get("infraspecific_names", {}).get("strain"),
            "has_annotation": "annotation_info" in r,
            "submitter": info.get("submitter"),
            "submission_date": info.get("release_date"),
        })
    return candidates


_ASSEMBLY_LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 0, "Scaffold": 1, "Contig": 2}


def _collapse_paired(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One entry per unique underlying assembly (a GCA/GCF pair collapses to
    one, keyed by whichever of the pair is the GCA -- INSDC accessions start
    with 'GCA_', matching UniProt's own genomeAssembly.assemblyId convention).

    Both members of a GCA/GCF pair can independently carry
    refseq_category="reference genome" for what is the identical underlying
    assembly (confirmed live: GCA_000143185.2 / GCF_000143185.2, Schizophyllum
    commune) -- without this collapse, rank_ncbi_candidates would report a
    spurious tie between two records describing one genome.
    """
    by_key: dict[str, dict[str, Any]] = {}
    seen_as_pair_of: set[str] = set()
    for c in candidates:
        if c["accession"] in seen_as_pair_of:
            continue
        key = c["accession"] if c["accession"].startswith("GCA_") else c.get("paired_accession") or c["accession"]
        if key not in by_key:
            by_key[key] = c if key == c["accession"] else next(
                (x for x in candidates if x["accession"] == key), c
            )
        if c.get("paired_accession"):
            seen_as_pair_of.add(c["paired_accession"])
    return list(by_key.values())


def rank_ncbi_candidates(candidates: list[dict[str, Any]], *, require_annotation: bool) -> list[dict[str, Any]]:
    """Returns the candidates in the top surviving rank tier (0, 1, or several
    if genuinely tied). Caller: len==0 -> none found, len==1 -> accept,
    len>1 -> report as tied and leave blank."""
    pool = _collapse_paired(candidates)
    if require_annotation:
        pool = [c for c in pool if c["has_annotation"]]
    if not pool:
        return []

    def tier(c: dict[str, Any]) -> int:
        if c["refseq_category"] == "reference genome":
            return 0
        if c["refseq_category"] == "representative genome":
            return 1
        return 2 + _ASSEMBLY_LEVEL_RANK.get(c["assembly_level"], 3)

    best_tier = min(tier(c) for c in pool)
    return [c for c in pool if tier(c) == best_tier]


@dataclass
class ResolvedRow:
    short: str
    resolved: bool
    protein_source: str = ""
    protein_accession: str = ""
    taxon_id: str = ""
    genome_source: str = ""
    genome_accession: str = ""
    busco_score: float | None = None  # from UniProt, when protein_source == "uniprot"
    report_lines: list[str] = field(default_factory=list)


def resolve_row(short: str, species: str, strain: str) -> ResolvedRow:
    """Resolve one species.csv row: Steps 0 (taxon), 1 (UniProt reference),
    1b (strain-targeted UniProt widening), 2 (NCBI genome/protein fallback),
    2b (strain-targeted NCBI widening). Any ambiguity leaves the row
    unresolved with a reason in report_lines rather than guessing.
    """
    taxa = resolve_taxon_id(species)
    if len(taxa) != 1:
        reason = "not a resolvable NCBI taxon name" if not taxa else f"{len(taxa)} ambiguous taxon matches"
        return ResolvedRow(short=short, resolved=False, report_lines=[reason])
    taxon_id = taxa[0]["taxon_id"]

    protein_source = protein_accession = genome_source = genome_accession = ""
    busco_score: float | None = None
    report: list[str] = []

    # Step 1: UniProt reference-only search
    uniprot_candidates = search_uniprot_proteomes(taxon_id, reference_only=True)
    chosen_uniprot = _pick_uniprot_candidate(uniprot_candidates, strain, report)

    if chosen_uniprot is None and strain:
        # Step 1b: strain-targeted widening (non-reference proteomes)
        widened = search_uniprot_proteomes(taxon_id, reference_only=False)
        matches = [c for c in widened if strain_matches(c["strain"], strain)]
        if len(matches) == 1:
            chosen_uniprot = matches[0]
            report.append(f"strain-targeted UniProt match: {strain!r}")

    if chosen_uniprot is not None:
        if check_ftp_available(chosen_uniprot["proteome_id"], chosen_uniprot["superkingdom"]):
            protein_source = "uniprot"
            protein_accession = chosen_uniprot["proteome_id"]
            busco_score = chosen_uniprot.get("busco_score")
            if chosen_uniprot.get("strain") and not strain_matches(chosen_uniprot["strain"], strain):
                report.append(f"WARNING strain mismatch: accepted record's strain is {chosen_uniprot['strain']!r}, species.csv says {strain!r}")
            if chosen_uniprot.get("genome_assembly_id"):
                genome_source = "ncbi"
                genome_accession = chosen_uniprot["genome_assembly_id"]
                _flag_if_superseded(taxon_id, genome_accession, report)
        else:
            report.append(
                f"REFERENCE proteome {chosen_uniprot['proteome_id']} found in REST but not yet "
                "in the FTP release -- treating as absent"
            )
            chosen_uniprot = None

    if not genome_source:
        # Step 2 / Step 2b: NCBI genome (and protein, if UniProt didn't resolve one)
        need_protein_too = not protein_source
        ncbi_candidates = query_ncbi_assemblies(taxon_id, require_reference=True)
        top = rank_ncbi_candidates(ncbi_candidates, require_annotation=need_protein_too)

        if not top and strain:
            # Collapse GCA/GCF pairs before strain-matching, or a single
            # underlying assembly with a paired record shows up as a
            # spurious 2-way tie.
            widened = _collapse_paired(query_ncbi_assemblies(taxon_id, require_reference=False))
            strain_hits = [c for c in widened if strain_matches(c["strain"], strain)]
            # An unannotated assembly can't supply a protein FASTA -- same
            # rule Step 2's rank_ncbi_candidates(require_annotation=...)
            # already applies to the non-widened path.
            matches = [c for c in strain_hits if c["has_annotation"]] if need_protein_too else strain_hits
            if len(matches) == 1:
                top = matches
                report.append(f"strain-targeted NCBI match: {strain!r}")
            elif strain_hits and not matches:
                report.append(
                    f"strain-targeted NCBI match {strain!r} found but lacks annotation -- "
                    "cannot supply a protein FASTA from it"
                )

        if len(top) == 1:
            genome_source = "ncbi"
            genome_accession = top[0]["accession"]
            if top[0]["assembly_status"] != "current":
                report.append(f"NOTE: {genome_accession} is superseded")
            if need_protein_too:
                protein_source = "ncbi"
                protein_accession = genome_accession
        elif len(top) > 1:
            report.append(f"{len(top)} NCBI assembly groups tied, none disambiguated:")
            for c in top[:10]:
                report.append(f"  {c['accession']}  ({c['assembly_level']}, submitted {c['submission_date']}, {c['submitter']})")
            if len(top) > 10:
                report.append(f"  ... and {len(top) - 10} more")

    # A row is resolved once its protein source is settled -- genome_source
    # is populated alongside protein in every code path here except the
    # UniProt-widening one (a non-reference proteome with no linked genome
    # assembly), where the protein is still usable on its own.
    resolved = bool(protein_source)
    if resolved and not genome_source:
        report.append(f"protein resolved ({protein_source} {protein_accession}), but no genome found")
    if not resolved and not report:
        report.append("no UniProt reference proteome and no NCBI assembly found")

    return ResolvedRow(
        short=short, resolved=resolved,
        protein_source=protein_source if resolved else "",
        protein_accession=protein_accession if resolved else "",
        taxon_id=taxon_id if resolved else "",
        genome_source=genome_source if resolved else "",
        genome_accession=genome_accession if resolved else "",
        busco_score=busco_score if resolved else None,
        report_lines=report,
    )


def _flag_if_superseded(taxon_id: str, genome_accession: str, report: list[str]) -> None:
    """When a UniProt proteome record supplies a genome_assembly_id directly
    (Step 1's one-query-resolves-both-columns path), that accession can be
    stale in UniProt's own record even though the proteome itself is current
    (live example: UP000007431 -> GCA_000143185.1, while NCBI's current
    assembly for that taxon is GCA_000143185.2). Cross-check against NCBI
    and flag rather than silently accepting it.
    """
    candidates = query_ncbi_assemblies(taxon_id, require_reference=False)
    match = next((c for c in candidates if c["accession"] == genome_accession), None)
    if match is not None and match["assembly_status"] != "current":
        report.append(f"NOTE: {genome_accession} is superseded")


def _pick_uniprot_candidate(candidates: list[dict[str, Any]], strain: str, report: list[str]) -> dict[str, Any] | None:
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        matches = [c for c in candidates if strain_matches(c["strain"], strain)]
        if len(matches) == 1:
            return matches[0]
        report.append(f"{len(candidates)} ambiguous UniProt reference proteomes, strain did not disambiguate")
    return None


def resolve_species_csv(study_dir: Path) -> None:
    """Drive resolve_row over every blank row of study_dir/species.csv,
    rewriting the file in place. Only rows where BOTH Protein_Source and
    Genome_Source are blank are touched -- a partially-filled row (e.g. a
    local FASTA already set) is left completely untouched.

    All output goes to plain stdout (print(), not sys.stderr) -- this is a
    human-readable resolution report, not a progress log alongside piped
    data, so it does not follow bin/build_study_config.py's stdout/stderr
    split.
    """
    species_csv = study_dir / "species.csv"
    with open(species_csv, newline="") as fh:
        reader = _csv.DictReader(fh)
        fieldnames = reader.fieldnames
        rows = list(reader)

    resolved_count = 0
    unresolved: list[ResolvedRow] = []

    for row in rows:
        if row["Protein_Source"] or row["Genome_Source"]:
            continue
        result = resolve_row(row["Short"], row["Species"], row["Strain"])
        if result.resolved:
            row["Protein_Source"] = result.protein_source
            row["Protein_Accession"] = result.protein_accession
            row["Taxon_ID"] = result.taxon_id
            row["Genome_Source"] = result.genome_source
            row["Genome_Accession"] = result.genome_accession
            resolved_count += 1
            busco_note = f" (BUSCO {result.busco_score:.1f})" if result.busco_score is not None else ""
            print(f"  {result.short}  {row['Species']}  -> {result.protein_source}{busco_note} "
                  f"{result.protein_accession} (+ {result.genome_source} {result.genome_accession})")
            for line in result.report_lines:
                print(f"    {line}")
        else:
            unresolved.append(result)

    with open(species_csv, "w", newline="") as fh:
        writer = _csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    total_blank = resolved_count + len(unresolved)
    print(f"\nResolved {resolved_count} of {total_blank} species needing resolution.")
    if unresolved:
        print(f"\n{len(unresolved)} species need a decision -- species.csv left blank for these rows:")
        for r in unresolved:
            print(f"  {r.short}")
            for line in r.report_lines:
                print(f"    {line}")
