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

import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
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
