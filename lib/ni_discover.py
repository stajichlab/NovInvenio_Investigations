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

import csv
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


_RACE_TOKEN_RE = re.compile(r"^(race\s*\d+|TR\d+|R\d+)$", re.IGNORECASE)


def derive_species_name(genome: dict[str, Any], rank_lookup: dict[str, dict[str, Any]]) -> str:
    """The genome's species-level name, from NCBI Taxonomy's own
    classification (`rank_lookup[...]["species_name"]`) -- never derived by
    string surgery on `organism_name`, per the design spec's explicit
    correction. Falls back to `organism_name` if the taxon is missing from
    `rank_lookup` or has no recorded species_name."""
    node = rank_lookup.get(genome["record_taxon_id"])
    if node and node.get("species_name"):
        return node["species_name"]
    return genome.get("organism_name", "")


def _sanitize(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "", s.replace(" ", "_"))


def pick_short(genome: dict[str, Any], species_abbrev: str, used_shorts: set[str]) -> str:
    """A species-prefixed, sanitized, de-duplicated short identifier for a
    genome. Prefers `isolate` over `strain` when `strain` looks like a bare
    race/pathotype token (e.g. "TR4", "race 4", "R1") and a real `isolate`
    is available -- otherwise uses `strain`, falling back to `isolate` or
    "unk". Never produces a bare-numeric Short (the species_abbrev prefix
    always precedes it). Adds a numeric suffix on collision with an
    already-used Short."""
    strain = genome.get("strain") or ""
    isolate = genome.get("isolate") or ""
    if strain and _RACE_TOKEN_RE.match(strain) and isolate:
        base = isolate
    elif strain:
        base = strain
    elif isolate:
        base = isolate
    else:
        base = "unk"
    candidate = _sanitize(f"{species_abbrev}_{base}")
    if candidate and candidate[0].isdigit():
        candidate = f"{species_abbrev}_{candidate}"
    final = candidate
    i = 2
    while final in used_shorts:
        final = f"{candidate}_{i}"
        i += 1
    return final


def _normalize_strain_key(genome: dict[str, Any]) -> str | None:
    raw = (genome.get("strain") or genome.get("isolate") or "").lower()
    raw = re.sub(r"^fo[_\s]?", "", raw)
    raw = re.sub(r"[^a-z0-9]", "", raw)
    return raw or None


def detect_duplicate_groups(genomes: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    """Partition genomes into duplicate-groups via a union-find over two
    independent criteria: matching `biosample_accession` (both non-None) OR
    matching normalized strain/isolate strings. Either criterion alone
    unions two genomes -- neither is a fallback-only path, since real NCBI
    data has cases needing each (e.g. the same strain re-registered under
    different BioSamples, or the same BioSample recorded under different
    strain spellings). A genome with no match to any other is its own
    singleton group."""
    n = len(genomes)
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj

    by_biosample: dict[str, int] = {}
    by_strain_key: dict[str, int] = {}
    for idx, g in enumerate(genomes):
        bs = g.get("biosample_accession")
        if bs:
            if bs in by_biosample:
                union(idx, by_biosample[bs])
            else:
                by_biosample[bs] = idx
        key = _normalize_strain_key(g)
        if key:
            if key in by_strain_key:
                union(idx, by_strain_key[key])
            else:
                by_strain_key[key] = idx

    groups: dict[int, list[dict[str, Any]]] = {}
    for idx, g in enumerate(genomes):
        groups.setdefault(find(idx), []).append(g)
    return list(groups.values())


def build_report(
    genomes: list[dict[str, Any]],
    rank_lookup: dict[str, dict[str, Any]],
    duplicate_groups: list[list[dict[str, Any]]],
    species_complex: dict[str, Any] | None,
    complex_genome_count: int,
) -> str:
    """Format a printed summary report: per-forma-specialis-group counts,
    quality metadata (provider/gene-count/assembly-level/BUSCO), duplicate-
    strain-group notes, and a species-complex note when one exists. Pure
    formatting over Tasks 1-2's data, no new querying."""
    lines = []

    # Count genomes by forma-specialis group
    group_counts: dict[str, int] = {}
    for g in genomes:
        grp = g.get("_forma_group", "unknown")
        group_counts[grp] = group_counts.get(grp, 0) + 1

    # Report per-group counts and quality metadata
    if group_counts:
        lines.append("Forma specialis groups:")
        for grp in sorted(group_counts.keys()):
            count = group_counts[grp]
            lines.append(f"  {grp}: {count}")
            # Collect quality metadata for this group
            group_genomes = [g for g in genomes if g.get("_forma_group") == grp]
            gene_counts = [g.get("protein_coding_count") for g in group_genomes if g.get("protein_coding_count") is not None]
            busco_scores = [g.get("busco_score") for g in group_genomes if g.get("busco_score") is not None]
            pipelines = set(g.get("annotation_pipeline") for g in group_genomes if g.get("annotation_pipeline"))
            assembly_levels = set(g.get("assembly_level") for g in group_genomes if g.get("assembly_level"))

            if gene_counts:
                avg_genes = sum(gene_counts) / len(gene_counts)
                lines.append(f"    Avg protein-coding genes: {avg_genes:.0f}")
            if busco_scores:
                avg_busco = sum(busco_scores) / len(busco_scores)
                lines.append(f"    Avg BUSCO: {avg_busco:.1f}")
            if pipelines:
                lines.append(f"    Pipelines: {', '.join(sorted(p for p in pipelines if p))}")
            if assembly_levels:
                lines.append(f"    Assembly levels: {', '.join(sorted(assembly_levels))}")

    # Report duplicate groups (those with more than one genome)
    multi_genome_groups = [g for g in duplicate_groups if len(g) > 1]
    if multi_genome_groups:
        if lines:
            lines.append("")
        lines.append("Duplicate/merged strains:")
        for group in multi_genome_groups:
            accessions = ", ".join(g.get("accession", "?") for g in group)
            lines.append(f"  Merged: {accessions}")

    # Report species complex if present
    if species_complex:
        if lines:
            lines.append("")
        lines.append("Species complex:")
        name = species_complex.get("name", "Unknown")
        lines.append(f"  {name} ({complex_genome_count} genomes)")

    return "\n".join(lines)


DISCOVER_HEADER = [
    "Short", "Species", "Strain", "Group", "TaxonGroup",
    "Protein_Source", "Protein_Accession", "Taxon_ID",
    "Genome_Source", "Genome_Accession", "GFF3_Source", "GFF3_Accession",
]


def discover_species(
    species: str,
    study_dir: Path,
    *,
    ingroup_groups: list[str] | None,
    outgroup_groups: list[str] | None,
    auto: bool,
    include_species_complex: bool,
) -> None:
    """Full orchestration for `bin/ni discover`: resolve taxon, query
    genomes (+ species-complex genomes if `include_species_complex`), fetch
    taxonomy ranks (batched, twice per Task 1's design), group, detect
    duplicates, collapse GCA/GCF pairs, assign Shorts, apply the selection
    mode, write species.csv, print build_report's output.

    Selection modes:
    - default (no ingroup_groups/outgroup_groups, auto=False): every row
      written with Group="" -- the human decides IN/OUT by hand afterward.
    - explicit ingroup_groups/outgroup_groups: only genomes whose forma-
      specialis group is named in one of the two lists are written, with
      Group set to "IN" or "OUT" accordingly. An unknown label is a hard
      error (sys.exit) naming the real available groups.
    - --auto: REPORT-ONLY. Prints the same report plus an extra proposal
      note, but writes every row with Group="" exactly like the default
      mode -- --auto must NEVER write IN/OUT into species.csv.

    --auto and explicit ingroup_groups/outgroup_groups are mutually
    exclusive (sys.exit if both given). Refuses to run if species.csv
    already exists (discover only seeds a brand-new study).
    """
    if (ingroup_groups or outgroup_groups) and auto:
        sys.exit("ERROR: --ingroup-groups/--outgroup-groups and --auto are mutually exclusive")

    species_csv = study_dir / "species.csv"
    if species_csv.exists():
        sys.exit(f"ERROR: {species_csv} already exists -- discover only seeds a new study")

    taxa = resolve_taxon_id(species)
    if len(taxa) != 1:
        sys.exit(f"ERROR: could not resolve '{species}' to exactly one taxon ({len(taxa)} matches)")
    taxon_id = taxa[0]["taxon_id"]

    genomes = query_species_genomes(taxon_id)

    # First batched taxonomy call: ranks for every genome's own record taxon.
    record_taxon_ids = sorted({g["record_taxon_id"] for g in genomes} | {taxon_id})
    rank_lookup = fetch_taxonomy_ranks(record_taxon_ids)

    # Second batched call: ranks for every parent id not already looked up
    # (needed so find_forma_specialis_group/find_species_complex can inspect
    # ancestors' own rank field, per Task 1's design).
    parent_ids = sorted({p for node in rank_lookup.values() for p in node.get("parents", [])} - set(rank_lookup))
    if parent_ids:
        rank_lookup.update(fetch_taxonomy_ranks(parent_ids))

    species_complex = find_species_complex(taxon_id, rank_lookup)
    complex_genome_count = 0
    if species_complex:
        if include_species_complex:
            genomes = query_species_genomes(species_complex["taxon_id"])
            extra_ids = sorted({g["record_taxon_id"] for g in genomes} - set(rank_lookup))
            if extra_ids:
                rank_lookup.update(fetch_taxonomy_ranks(extra_ids))
        else:
            all_complex_genomes = query_species_genomes(species_complex["taxon_id"])
            complex_genome_count = len(all_complex_genomes) - len(genomes)

    genomes = _collapse_paired(genomes)

    for g in genomes:
        g["_forma_group"], g["_used_fallback"] = find_forma_specialis_group(
            g["record_taxon_id"], rank_lookup, g["organism_name"]
        )

    duplicate_groups = detect_duplicate_groups(genomes)
    # One representative genome per duplicate group -- the rest are merged in.
    representatives = [grp[0] for grp in duplicate_groups]

    available_group_labels = {g["_forma_group"] for g in representatives}
    for label in (ingroup_groups or []) + (outgroup_groups or []):
        if label not in available_group_labels:
            sys.exit(f"ERROR: group label {label!r} not found -- available groups: {sorted(available_group_labels)}")

    report = build_report(representatives, rank_lookup, duplicate_groups, species_complex, complex_genome_count)
    print(report)
    if auto:
        print("\n--auto: this is a PROPOSAL to review, not a trusted grouping. "
              "Group is left blank for every row -- assign IN/OUT by hand.")

    rows = []
    used_shorts: set[str] = set()
    species_abbrev = "".join(w[:2] for w in species.split()[:2]) or "Sp"
    for g in representatives:
        group_value = ""
        if ingroup_groups and g["_forma_group"] in ingroup_groups:
            group_value = "IN"
        elif outgroup_groups and g["_forma_group"] in outgroup_groups:
            group_value = "OUT"
        elif ingroup_groups or outgroup_groups:
            continue  # explicit mode: skip genomes not in any named group
        short = pick_short(g, species_abbrev, used_shorts)
        used_shorts.add(short)
        rows.append({
            "Short": short,
            "Species": derive_species_name(g, rank_lookup),
            "Strain": g.get("strain") or g.get("isolate") or "",
            "Group": group_value,
            "TaxonGroup": g["_forma_group"],
            "Protein_Source": "ncbi",
            "Protein_Accession": g["accession"],
            "Taxon_ID": g["record_taxon_id"],
            "Genome_Source": "ncbi",
            "Genome_Accession": g["accession"],
            "GFF3_Source": "",
            "GFF3_Accession": "",
        })

    study_dir.mkdir(parents=True, exist_ok=True)
    with open(species_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=DISCOVER_HEADER)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {species_csv} ({len(rows)} species)")
