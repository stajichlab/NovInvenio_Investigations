"""NCBI species-population enumeration for bin/ni discover.

See notes/superpowers/specs/2026-09-11-ni-discover-design.md for the full
design, and its two independent bioinformatics reviews (both live-verified
against real NCBI data) for why every rule here exists -- several are
corrections of an initially-plausible-but-wrong first draft.

Field paths below were verified live against the `datasets` CLI on
2026-09-12 (see task-1-report.md for the full transcript). Deviations from
this task's original design sketch:

- `datasets summary taxonomy taxon <id1> <id2> ...` returns ONE JSON object
  (not JSON-lines, unlike the genome-summary commands), shaped
  `{"reports": [...], "total_count": N}` -- one report per queried id, in
  query order. Multi-id batching works exactly as assumed.
- Each report's taxonomy record has no top-level "name" field for the
  taxon's own name -- that lives at `taxonomy.current_scientific_name.name`.
- `classification.species` (id + name) gives the species-level ancestor,
  used for `species_name`; it is ABSENT (not merely null) on a report for a
  taxon at or above species rank (live-verified: 171631, a SPECIES_GROUP,
  has no "species" key under classification at all) -- `species_name` is
  `""` in that case.
- `parents` is the FULL ancestor lineage as a flat list of ids, ordered
  root-first with the immediate parent LAST -- not a single-element direct-
  parent field. `find_forma_specialis_group`/`find_species_complex` walk one
  step at a time using `parents[-1]` (the immediate parent), which also
  matches this test suite's single-element `parents` fixtures since a
  1-element list's last element equals its first.
- `annotation_info`'s gene-count and BUSCO field paths in the design sketch
  (`stats.gene_counts.protein_coding`, `busco.complete`) were confirmed
  exactly right by a live `--annotated` query (Fusarium oxysporum, 70
  records). `annotation_info.pipeline` (e.g. "NCBI Eukaryotic Annotation
  Propagation Pipeline") is present on RefSeq (GCF_) records and absent on
  GenBank-submitted ones -- None in that case, not a parsing bug.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).parent))
from ni_resolve import _collapse_paired, resolve_taxon_id  # noqa: F401,E402

_FSP_REGEX = re.compile(r"f\.\s*sp\.\s*(\S+)", re.IGNORECASE)


def query_species_genomes(taxon_id: str) -> list[dict[str, Any]]:
    """Every NCBI-annotated genome for a taxon (`datasets summary genome
    taxon <id> --annotated --as-json-lines`)."""
    cmd = ["datasets", "summary", "genome", "taxon", taxon_id, "--annotated", "--as-json-lines"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    genomes = []
    for line in result.stdout.strip().splitlines():
        if not line:
            continue
        r = json.loads(line)
        info = r.get("assembly_info", {})
        org = r.get("organism", {})
        infra = org.get("infraspecific_names", {})
        ann = r.get("annotation_info", {})
        genomes.append({
            "accession": r["accession"],
            "paired_accession": r.get("paired_accession"),
            "organism_name": org.get("organism_name", ""),
            "record_taxon_id": str(org.get("tax_id", "")),
            "strain": infra.get("strain"),
            "isolate": infra.get("isolate"),
            "biosample_accession": info.get("biosample", {}).get("accession"),
            "assembly_level": info.get("assembly_level"),
            "has_annotation": bool(ann),
            "annotation_pipeline": ann.get("pipeline"),
            "protein_coding_count": ann.get("stats", {}).get("gene_counts", {}).get("protein_coding"),
            "busco_score": ann.get("busco", {}).get("complete"),
        })
    return genomes


def fetch_taxonomy_ranks(taxon_ids: list[str]) -> dict[str, dict[str, Any]]:
    """One batched `datasets summary taxonomy taxon <id1> <id2> ...` call.

    Returns `{taxon_id: {"rank": str, "name": str, "parents": list[str],
    "species_name": str}}`. The real response is a single JSON object
    `{"reports": [...]}`, one report per requested id (in the order given),
    each shaped `{"query": ["<id>"], "taxonomy": {...}}` -- see the module
    docstring for the exact field paths this pulls from `taxonomy`.
    """
    if not taxon_ids:
        return {}
    cmd = ["datasets", "summary", "taxonomy", "taxon", *taxon_ids]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    payload = json.loads(result.stdout)
    ranks: dict[str, dict[str, Any]] = {}
    for report in payload.get("reports", []):
        t = report.get("taxonomy", {})
        tax_id = str(t.get("tax_id", ""))
        if not tax_id:
            continue
        species = t.get("classification", {}).get("species")
        ranks[tax_id] = {
            "rank": t.get("rank", ""),
            "name": t.get("current_scientific_name", {}).get("name", ""),
            "parents": [str(p) for p in t.get("parents", [])],
            "species_name": species["name"] if species else "",
        }
    return ranks


def find_forma_specialis_group(
    record_taxon_id: str, rank_lookup: dict[str, dict[str, Any]], organism_name: str
) -> tuple[str, bool]:
    """Walk `rank_lookup[record_taxon_id]`'s ancestor chain (one step at a
    time via each node's immediate parent, `parents[-1]`) looking for a
    FORMA_SPECIALIS-rank ancestor. Falls back to a regex on `organism_name`
    if no such ancestor is found in `rank_lookup` (e.g. the caller's second
    batched lookup didn't reach far enough up the tree). Returns
    `(group_label, used_fallback)`.
    """
    seen: set[str] = set()
    current: str | None = record_taxon_id
    while current and current not in seen:
        seen.add(current)
        node = rank_lookup.get(current)
        if node is None:
            break
        if node["rank"] == "FORMA_SPECIALIS":
            label = node["name"]
            prefix = node.get("species_name", "")
            if prefix and label.startswith(prefix):
                label = label[len(prefix):].strip()
            m = _FSP_REGEX.search(label)
            return (m.group(1) if m else label), False
        parents = node.get("parents", [])
        current = parents[-1] if parents else None

    m = _FSP_REGEX.search(organism_name)
    if m:
        return m.group(1), True
    return "no-fsp-in-name", False


def find_species_complex(taxon_id: str, rank_lookup: dict[str, dict[str, Any]]) -> dict[str, Any] | None:
    """Walk `taxon_id`'s ancestor chain (via each node's immediate parent,
    `parents[-1]`) looking for a SPECIES_GROUP-rank ancestor. Returns
    `{"taxon_id": ..., "name": ...}` or `None`."""
    seen: set[str] = set()
    current: str | None = taxon_id
    while current and current not in seen:
        seen.add(current)
        node = rank_lookup.get(current)
        if node is None:
            return None
        if node["rank"] == "SPECIES_GROUP":
            return {"taxon_id": current, "name": node["name"]}
        parents = node.get("parents", [])
        current = parents[-1] if parents else None
    return None
