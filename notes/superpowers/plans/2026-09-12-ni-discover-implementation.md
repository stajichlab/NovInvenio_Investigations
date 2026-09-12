# `bin/ni discover` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `bin/ni discover --species "<name>" --study-dir <dir>` — given
a species name, enumerate every NCBI-annotated genome, group by
FORMA_SPECIALIS taxonomy rank, detect duplicate strain registrations and
sibling-species-complex genomes, and write `species.csv` rows (`Group` blank
by default; explicit `--ingroup-groups`/`--outgroup-groups` or a report-only
`--auto` for the other two selection modes).

**Architecture:** A new `lib/ni_discover.py` module, reusing
`lib/ni_resolve.py`'s already-shipped `_collapse_paired` (unchanged import,
no new pair-preference rule — the design spec's re-review found none is
needed) and `resolve_taxon_id` (for resolving `--species` to a taxon ID).
Everything else is new: a richer NCBI genome-record parser (more fields than
`resolve`'s `_parse_assembly_record` needs), a batched NCBI Taxonomy
rank/lineage lookup, forma-specialis grouping, duplicate-strain detection,
and the group/count report table. `bin/ni` gets a third subcommand.

**Tech Stack:** Python 3.12, `pytest`, stdlib `subprocess`/`json`/`re`, the
`datasets` CLI (already a pixi dependency). No new dependencies.

**Spec:** `notes/superpowers/specs/2026-09-11-ni-discover-design.md`

## Global Constraints

- **Never write an automatic `IN`/`OUT` split.** `--auto` prints an extended
  proposal report but writes `Group` blank, identically to the default mode
  — this is the spec's central, twice-reviewed correction. No task in this
  plan may reintroduce automatic group-assignment logic.
- Query NCBI by resolved taxon ID, reusing `resolve_taxon_id` from
  `lib/ni_resolve.py` (do not re-implement taxon resolution).
- Restricting to annotated genomes (`--annotated` on the `datasets` CLI call)
  is always on, not a flag.
- Duplicate-strain detection is `assembly_info.biosample.accession` match
  **OR** normalized strain/isolate string match, unioned — never
  "BioSample-first, string-match only as a fallback" (a fallback-only rule
  never fires on real data, since every live-checked record carries a
  BioSample).
- FORMA_SPECIALIS grouping uses NCBI Taxonomy's own rank field via a batched
  lookup (not string-parsing `organism_name` as the primary mechanism) — the
  regex is a fallback only, for a taxon with no FORMA_SPECIALIS ancestor.
- GCA/GCF pair-collapse reuses `_collapse_paired` unmodified — no new
  provider-based preference rule (verified unnecessary: real pairs share the
  same `annotation_info.provider`, RefSeq propagates the submitter's
  annotation rather than producing an independent one).
- `discover` refuses to run if `studies/<domain>/<set>/species.csv` already
  exists (`sys.exit`, clear message) — it seeds a *new* study only.
- `--ingroup-groups`/`--outgroup-groups` and `--auto` are mutually exclusive.
- Run tests via `pixi run pytest tests/test_ni_discover.py -v`.

---

### Task 1: NCBI genome query + batched taxonomy rank lookup

**Files:**
- Create: `lib/ni_discover.py`
- Test: `tests/test_ni_discover.py` (new)

**Interfaces:**
- Produces: `query_species_genomes(taxon_id: str) -> list[dict[str, Any]]` —
  richer per-genome dict than `ni_resolve`'s `_parse_assembly_record`: at
  least `accession`, `paired_accession`, `organism_name`, `record_taxon_id`,
  `strain`, `isolate`, `biosample_accession`, `assembly_level`,
  `has_annotation`, `annotation_pipeline`, `protein_coding_count`,
  `busco_score` (float or None).
- Produces: `fetch_taxonomy_ranks(taxon_ids: list[str]) -> dict[str, dict]` —
  one batched `datasets summary taxonomy taxon <id1> <id2> ...` call,
  returns `{taxon_id: {"rank": str, "name": str, "parents": list[str],
  "species_name": str}}` per id (`species_name` from that id's own
  `classification.species.name` field).
- Produces: `find_forma_specialis_group(record_taxon_id: str, rank_lookup:
  dict[str, dict], organism_name: str) -> tuple[str, bool]` — returns
  `(group_label, used_fallback)`; walks `rank_lookup[record_taxon_id]`'s
  `parents` looking for one whose own rank is `FORMA_SPECIALIS` (a **second**
  batched `fetch_taxonomy_ranks` call on the parent ids, since the first
  call's `parents` field carries no rank info of its own — confirmed live
  during this spec's review, not a guess); falls back to a `f\. sp\. (\S+)`
  regex on `organism_name` if no such ancestor exists, returning
  `used_fallback=True` in that case. Group label is the bare forma name
  (species prefix stripped), or `"no-fsp-in-name"` if neither the rank walk
  nor the regex finds one.
- Produces: `find_species_complex(taxon_id: str, rank_lookup: dict[str,
  dict]) -> dict | None` — walks the resolved species' own parents for a
  `SPECIES_GROUP`-rank ancestor; returns `{"taxon_id": ..., "name": ...}` or
  `None`.
- Later tasks import all of the above.

**Before writing code — live verification is mandatory, not optional (this
spec's own review process caught wrong-guessed field shapes twice already):**
```bash
pixi run datasets summary genome taxon "Fusarium oxysporum" --annotated --as-json-lines | head -c 3000 | python3 -m json.tool 2>&1 | head -80
pixi run datasets summary taxonomy taxon 1229664 | python3 -m json.tool
pixi run datasets summary taxonomy taxon 1229664 61366 5507 171631 | python3 -m json.tool
```
Confirm: (1) the genome-record fields this task needs — especially
`annotation_info`'s pipeline/provider/gene-count/busco field names, which
have NOT been pinned exactly by prior live checks in this session (unlike
`paired_accession`, `refseq_category`, etc., which Task 1/2 of the earlier
`resolve` plan already confirmed); (2) that `datasets summary taxonomy taxon`
accepts multiple ids in one call and what the response shape is for `rank`,
`parents`, and `classification.species.name`; (3) that a `SPECIES_GROUP`-rank
ancestor is reachable the same way. Adjust the sketch below to match reality,
not the other way around.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_ni_discover.py`. Use `tests/fixtures/fox_annotated.jsonl`
(a frozen live capture — see Step 1a) rather than hand-written JSON for the
grouping/parsing tests, per the spec's explicit requirement.

- [ ] **Step 1a: Capture the live fixtures**

```bash
mkdir -p tests/fixtures
pixi run datasets summary genome taxon "Fusarium oxysporum" --annotated --as-json-lines > tests/fixtures/fox_annotated.jsonl
pixi run datasets summary genome taxon 171631 --annotated --as-json-lines > tests/fixtures/fox_complex_annotated.jsonl
```
Verify: `wc -l tests/fixtures/fox_annotated.jsonl` should be 70 (per this
spec's own live-verified count as of 2026-09-11 — if NCBI's data has moved
since, that's fine, just note the new count in your report; the fixture is
frozen at whatever it actually is today, not forced to match an old number).
`fox_complex_annotated.jsonl` should be a superset including the
`Fusarium odoratissimum` records (74 as of the spec's writing).

Write these tests in `tests/test_ni_discover.py`:

```python
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "lib"))
sys.path.insert(0, str(REPO / "bin"))
import ni_discover as nd  # noqa: E402

FIXTURES = REPO / "tests" / "fixtures"


def _load_fixture_lines(name: str) -> str:
    return (FIXTURES / name).read_text()


def test_query_species_genomes_parses_live_fixture():
    fake_run = MagicMock(stdout=_load_fixture_lines("fox_annotated.jsonl"), returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        genomes = nd.query_species_genomes("5507")
    assert len(genomes) > 0
    # Every record must carry the fields later tasks depend on.
    for g in genomes:
        assert "accession" in g
        assert "organism_name" in g
        assert "record_taxon_id" in g
        assert "has_annotation" in g
    # Spot-check a known record from the live fixture (adjust if the live
    # capture no longer contains this exact accession -- the point is that
    # SOME record's fields are checked against known-real values, not that
    # this specific accession is eternal).
    fo47 = [g for g in genomes if g["strain"] == "Fo47"]
    assert fo47, "expected at least one Fo47 record in the live fixture"


def test_fetch_taxonomy_ranks_batched_call_shape():
    # Mirrors the real 'datasets summary taxonomy taxon <id1> <id2> ...' shape
    # you verified live in Task 1's pre-code step -- fill in with the REAL
    # response shape you found, not a guess.
    raise NotImplementedError(
        "Write this test against the real `datasets summary taxonomy taxon` "
        "JSON shape you verified live in the 'before writing code' step -- "
        "do not invent a shape. Mock subprocess.run to return that real "
        "shape (trimmed to 1-2 records), assert fetch_taxonomy_ranks parses "
        "rank/name/parents/species_name correctly for each id."
    )


def test_find_forma_specialis_group_walks_rank_lookup():
    # cubense race 1 (taxid 1229664) -> parent 61366 "Fusarium oxysporum f.
    # sp. cubense" (rank FORMA_SPECIALIS) -- live-verified in the design
    # spec's review. Build rank_lookup with exactly the real shape you
    # confirmed in fetch_taxonomy_ranks's own test above.
    rank_lookup = {
        "1229664": {"rank": "STRAIN", "name": "Fusarium oxysporum f. sp. cubense race 1", "parents": ["61366"], "species_name": "Fusarium oxysporum"},
        "61366": {"rank": "FORMA_SPECIALIS", "name": "Fusarium oxysporum f. sp. cubense", "parents": ["5507"], "species_name": "Fusarium oxysporum"},
    }
    label, used_fallback = nd.find_forma_specialis_group("1229664", rank_lookup, "Fusarium oxysporum f. sp. cubense race 1")
    assert label == "cubense"
    assert used_fallback is False


def test_find_forma_specialis_group_falls_back_to_regex():
    # A taxid with no FORMA_SPECIALIS ancestor in rank_lookup at all.
    rank_lookup = {
        "999999": {"rank": "STRAIN", "name": "Some organism f. sp. madeup", "parents": ["1"], "species_name": "Some organism"},
        "1": {"rank": "SPECIES", "name": "Some organism", "parents": [], "species_name": "Some organism"},
    }
    label, used_fallback = nd.find_forma_specialis_group("999999", rank_lookup, "Some organism f. sp. madeup")
    assert label == "madeup"
    assert used_fallback is True


def test_find_forma_specialis_group_no_fsp_anywhere():
    rank_lookup = {
        "660027": {"rank": "STRAIN", "name": "Fusarium oxysporum Fo47", "parents": ["5507"], "species_name": "Fusarium oxysporum"},
        "5507": {"rank": "SPECIES", "name": "Fusarium oxysporum", "parents": ["171631"], "species_name": "Fusarium oxysporum"},
    }
    label, used_fallback = nd.find_forma_specialis_group("660027", rank_lookup, "Fusarium oxysporum Fo47")
    assert label == "no-fsp-in-name"


def test_find_species_complex_detects_species_group_ancestor():
    rank_lookup = {
        "5507": {"rank": "SPECIES", "name": "Fusarium oxysporum", "parents": ["171631"], "species_name": "Fusarium oxysporum"},
        "171631": {"rank": "SPECIES_GROUP", "name": "Fusarium oxysporum species complex", "parents": [], "species_name": ""},
    }
    complex_info = nd.find_species_complex("5507", rank_lookup)
    assert complex_info == {"taxon_id": "171631", "name": "Fusarium oxysporum species complex"}


def test_find_species_complex_none_when_no_species_group_ancestor():
    rank_lookup = {
        "5507": {"rank": "SPECIES", "name": "Fusarium oxysporum", "parents": ["4751"], "species_name": "Fusarium oxysporum"},
        "4751": {"rank": "KINGDOM", "name": "Fungi", "parents": [], "species_name": ""},
    }
    assert nd.find_species_complex("5507", rank_lookup) is None
```

Delete the `NotImplementedError` placeholder test and replace it with a real
one once you've done the live verification — it exists only to force you to
look before writing, per this task's own instructions; a plan step that ships
with an actual `NotImplementedError` in the final code is not acceptable, it
must not survive past Step 1.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run pytest tests/test_ni_discover.py -v`
Expected: `ModuleNotFoundError: No module named 'ni_discover'`.

- [ ] **Step 3: Write `lib/ni_discover.py`**

Write the real implementation based on what Step 1's live verification found
— do not treat this outline as final field paths, only as the functions'
shapes and responsibilities:

```python
"""NCBI species-population enumeration for bin/ni discover.

See notes/superpowers/specs/2026-09-11-ni-discover-design.md for the full
design, and its two independent bioinformatics reviews (both live-verified
against real NCBI data) for why every rule here exists -- several are
corrections of an initially-plausible-but-wrong first draft.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).parent))
from ni_resolve import _collapse_paired, resolve_taxon_id  # noqa: E402

_FSP_REGEX = re.compile(r"f\.\s*sp\.\s*(\S+)", re.IGNORECASE)


def query_species_genomes(taxon_id: str) -> list[dict[str, Any]]:
    """Every NCBI-annotated genome for a taxon. Fill in the exact field paths
    from your live verification -- this is a starting shape, not a final one."""
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
            # TODO(implementer, per live verification): confirm the real
            # gene-count and BUSCO field paths -- these are placeholders
            # from the design spec's prose description, not a verified path.
            "protein_coding_count": ann.get("stats", {}).get("gene_counts", {}).get("protein_coding"),
            "busco_score": ann.get("busco", {}).get("complete"),
        })
    return genomes


def fetch_taxonomy_ranks(taxon_ids: list[str]) -> dict[str, dict[str, Any]]:
    """One batched `datasets summary taxonomy taxon <id1> <id2> ...` call.
    Fill in the real response-parsing logic from your live verification."""
    if not taxon_ids:
        return {}
    cmd = ["datasets", "summary", "taxonomy", "taxon", *taxon_ids]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    # TODO(implementer): parse result.stdout per the REAL shape you verified
    # live -- this is not necessarily JSON-lines like the genome command;
    # confirm the actual output format (could be a single JSON object with a
    # "reports" list, or something else entirely) before writing this parser.
    raise NotImplementedError("Parse per the live-verified response shape from Task 1's Step 0")


def find_forma_specialis_group(
    record_taxon_id: str, rank_lookup: dict[str, dict[str, Any]], organism_name: str
) -> tuple[str, bool]:
    seen = set()
    current = record_taxon_id
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
        current = parents[0] if parents else None

    m = _FSP_REGEX.search(organism_name)
    if m:
        return m.group(1), True
    return "no-fsp-in-name", False


def find_species_complex(taxon_id: str, rank_lookup: dict[str, dict[str, Any]]) -> dict[str, Any] | None:
    seen = set()
    current = taxon_id
    while current and current not in seen:
        seen.add(current)
        node = rank_lookup.get(current)
        if node is None:
            return None
        if node["rank"] == "SPECIES_GROUP":
            return {"taxon_id": current, "name": node["name"]}
        parents = node.get("parents", [])
        current = parents[0] if parents else None
    return None
```

Note the two `NotImplementedError`/`TODO` markers in `fetch_taxonomy_ranks`
and the gene-count/BUSCO field paths in `query_species_genomes` — these are
explicitly not-yet-verified spots this task's own Step 0 requires you to
resolve with a real query before this step is done. Do not leave either in
the committed code.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run pytest tests/test_ni_discover.py -v`
Expected: all tests PASS (at least 6 — adjust count for whatever real tests
Step 1's `fetch_taxonomy_ranks` verification produced).

- [ ] **Step 5: Commit**

```bash
git add lib/ni_discover.py tests/test_ni_discover.py tests/fixtures/fox_annotated.jsonl tests/fixtures/fox_complex_annotated.jsonl
git commit -m "Add NCBI genome query + batched taxonomy rank lookup for bin/ni discover

query_species_genomes, fetch_taxonomy_ranks, find_forma_specialis_group
(rank-based grouping with organism_name regex fallback), and
find_species_complex, per notes/superpowers/specs/2026-09-11-ni-discover-design.md.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 2: Species/Short derivation + duplicate-strain detection

**Files:**
- Modify: `lib/ni_discover.py`
- Modify: `tests/test_ni_discover.py`

**Interfaces:**
- Consumes: `query_species_genomes`, `fetch_taxonomy_ranks` (Task 1).
- Produces: `derive_species_name(genome: dict, rank_lookup: dict[str, dict])
  -> str` — from `rank_lookup[genome["record_taxon_id"]]["species_name"]`
  (Task 1 already captures this per taxon; do not re-derive from
  `organism_name` string surgery, per the spec's explicit correction).
- Produces: `pick_short(genome: dict, species_abbrev: str, used_shorts:
  set[str]) -> str` — species-prefixed, sanitized, de-duplicated identifier.
  Prefers `isolate` over `strain` when `strain` looks like a bare
  race/pathotype token (a short regex/heuristic: a strain value matching
  `^(race\s*\d+|TR\d+|R\d+)$`, case-insensitive, with a real `isolate`
  available) — otherwise uses `strain`. Sanitizes to `[A-Za-z0-9_]`, adds a
  numeric suffix on collision with an already-used Short.
- Produces: `detect_duplicate_groups(genomes: list[dict]) -> list[list[dict]]`
  — partitions genomes into duplicate-groups via a union-find (or simpler:
  build an undirected graph where two genomes are linked if
  `biosample_accession` matches (both non-None) **or** normalized
  strain/isolate strings match, then take connected components). A genome
  with no match to any other is its own singleton group. Normalization:
  lowercase, strip a leading `fo`/species-abbreviation prefix, strip
  non-alphanumerics.
- Later tasks import all of the above.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ni_discover.py`:

```python
def test_derive_species_name_uses_taxonomy_not_string_surgery():
    genome = {"record_taxon_id": "660027", "organism_name": "Fusarium oxysporum Fo47"}
    rank_lookup = {"660027": {"rank": "STRAIN", "name": "Fusarium oxysporum Fo47", "parents": ["5507"], "species_name": "Fusarium oxysporum"}}
    assert nd.derive_species_name(genome, rank_lookup) == "Fusarium oxysporum"


def test_pick_short_prefers_isolate_over_race_token_strain():
    genome = {"strain": "TR4", "isolate": "UK0001"}
    short = nd.pick_short(genome, "Foxy", set())
    assert "TR4" not in short
    assert "UK0001" in short


def test_pick_short_uses_strain_when_not_a_race_token():
    genome = {"strain": "Fo47", "isolate": None}
    short = nd.pick_short(genome, "Foxy", set())
    assert "Fo47" in short


def test_pick_short_deduplicates_on_collision():
    used = set()
    first = nd.pick_short({"strain": "X1", "isolate": None}, "Foxy", used)
    used.add(first)
    second = nd.pick_short({"strain": "X1", "isolate": None}, "Foxy", used)
    assert first != second


def test_pick_short_never_bare_numeric():
    short = nd.pick_short({"strain": "9", "isolate": None}, "Foxy", set())
    assert not short.isdigit()


def test_detect_duplicate_groups_unions_biosample_and_strain_match():
    # Fo5176 registered 3 times with 3 DIFFERENT biosample accessions
    # (live-verified in the design spec's review) -- must still collapse via
    # the strain-string match arm of the union, not the biosample arm.
    genomes = [
        {"accession": "GCA_1", "biosample_accession": "SAMN_A", "strain": "Fo5176", "isolate": None},
        {"accession": "GCA_2", "biosample_accession": "SAMN_B", "strain": "Fo5176", "isolate": None},
        {"accession": "GCA_3", "biosample_accession": "SAMN_C", "strain": "Fo5176", "isolate": None},
        {"accession": "GCA_4", "biosample_accession": "SAMN_D", "strain": "SomethingElse", "isolate": None},
    ]
    groups = nd.detect_duplicate_groups(genomes)
    sizes = sorted(len(g) for g in groups)
    assert sizes == [1, 3]


def test_detect_duplicate_groups_unions_via_biosample_when_strain_differs():
    genomes = [
        {"accession": "GCA_1", "biosample_accession": "SAMN_SHARED", "strain": "NameA", "isolate": None},
        {"accession": "GCA_2", "biosample_accession": "SAMN_SHARED", "strain": "NameB", "isolate": None},
    ]
    groups = nd.detect_duplicate_groups(genomes)
    assert len(groups) == 1
    assert len(groups[0]) == 2
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run pytest tests/test_ni_discover.py -v`
Expected: `AttributeError` on the 3 not-yet-defined functions.

- [ ] **Step 3: Add to `lib/ni_discover.py`**

```python
_RACE_TOKEN_RE = re.compile(r"^(race\s*\d+|TR\d+|R\d+)$", re.IGNORECASE)


def derive_species_name(genome: dict[str, Any], rank_lookup: dict[str, dict[str, Any]]) -> str:
    node = rank_lookup.get(genome["record_taxon_id"])
    if node and node.get("species_name"):
        return node["species_name"]
    return genome.get("organism_name", "")


def _sanitize(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "", s.replace(" ", "_"))


def pick_short(genome: dict[str, Any], species_abbrev: str, used_shorts: set[str]) -> str:
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run pytest tests/test_ni_discover.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add lib/ni_discover.py tests/test_ni_discover.py
git commit -m "Add Species/Short derivation + duplicate-strain detection to bin/ni discover

derive_species_name reads NCBI Taxonomy's own classification.species.name
(not organism_name string surgery); pick_short prefers isolate over a
bare race/pathotype-token strain and is never a bare numeric; 
detect_duplicate_groups unions BioSample-match and normalized-strain-
match as independent criteria, per the design spec's live-verified
correction (a fallback-only rule never fires on real data).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 3: Group/count report table

**Files:**
- Modify: `lib/ni_discover.py`
- Modify: `tests/test_ni_discover.py`

**Interfaces:**
- Consumes: everything from Tasks 1-2.
- Produces: `build_report(
    genomes: list[dict], rank_lookup: dict[str, dict], duplicate_groups:
    list[list[dict]], species_complex: dict | None, complex_genome_count: int
  ) -> str` — the full printed report: per-forma-specialis-group counts,
  provider/gene-count/assembly-level/N50 spread per group (or overall if
  per-group is too sparse — your call on the exact layout, but the data must
  be present), duplicate-group notes (which accessions were merged),
  species-complex note (if one exists: name, and how many additional genomes
  live there).
- No new querying in this task — pure formatting over data the earlier two
  tasks already produce.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ni_discover.py`. Use small, hand-built genome lists
(this task is formatting logic, not API-shape-sensitive, so hand-built
fixtures are fine here unlike Tasks 1-2):

```python
def test_build_report_includes_group_counts():
    genomes = [
        {"accession": "GCA_1", "record_taxon_id": "1", "organism_name": "X f. sp. a", "strain": "s1", "isolate": None, "assembly_level": "Scaffold", "annotation_pipeline": "p1", "protein_coding_count": 100, "busco_score": 90.0},
        {"accession": "GCA_2", "record_taxon_id": "2", "organism_name": "X f. sp. a", "strain": "s2", "isolate": None, "assembly_level": "Scaffold", "annotation_pipeline": "p1", "protein_coding_count": 200, "busco_score": 95.0},
        {"accession": "GCA_3", "record_taxon_id": "3", "organism_name": "X", "strain": "s3", "isolate": None, "assembly_level": "Chromosome", "annotation_pipeline": "p2", "protein_coding_count": 150, "busco_score": None},
    ]
    rank_lookup = {}  # forma_group precomputed and stashed on each genome dict for this test
    for g, grp in zip(genomes, ["a", "a", "no-fsp-in-name"]):
        g["_forma_group"] = grp
    report = nd.build_report(genomes, rank_lookup, duplicate_groups=[[g] for g in genomes], species_complex=None, complex_genome_count=0)
    assert "a" in report
    assert "2" in report  # count for group "a"
    assert "no-fsp-in-name" in report


def test_build_report_notes_duplicate_groups():
    genomes = [
        {"accession": "GCA_1", "strain": "Fo5176", "isolate": None, "_forma_group": "no-fsp-in-name"},
        {"accession": "GCA_2", "strain": "Fo5176", "isolate": None, "_forma_group": "no-fsp-in-name"},
    ]
    report = nd.build_report(genomes, {}, duplicate_groups=[genomes], species_complex=None, complex_genome_count=0)
    assert "GCA_1" in report and "GCA_2" in report
    assert "duplicate" in report.lower() or "merged" in report.lower()


def test_build_report_notes_species_complex():
    report = nd.build_report(
        [], {}, duplicate_groups=[],
        species_complex={"taxon_id": "171631", "name": "Fusarium oxysporum species complex"},
        complex_genome_count=4,
    )
    assert "Fusarium oxysporum species complex" in report
    assert "4" in report
```

Note: these tests assume `build_report` reads a precomputed `"_forma_group"`
key on each genome dict rather than recomputing it — decide during
implementation whether `build_report` takes already-grouped genomes (cleaner
separation, Task 4 computes groups once and passes them in) or recomputes
grouping itself; if you change this interface, update these tests to match
and note the change in your report (this is a legitimate implementation
judgment call the brief doesn't need to over-specify).

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run pytest tests/test_ni_discover.py -v`
Expected: `AttributeError: ... 'build_report'`.

- [ ] **Step 3: Implement `build_report` in `lib/ni_discover.py`**

Write this from scratch based on the test expectations above — the exact
formatting (column widths, ordering) is your call; the requirements are: every
forma-specialis group and its count appears, duplicate-group merges are
named (which accessions), and a species-complex note appears when one was
found, naming it and the additional genome count.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run pytest tests/test_ni_discover.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add lib/ni_discover.py tests/test_ni_discover.py
git commit -m "Add group/count report table for bin/ni discover

build_report() formats the printed summary: per-forma-specialis-group
counts and quality metadata (provider/gene-count/assembly-level/BUSCO),
duplicate-strain-group notes, and a species-complex note when one
exists -- pure formatting over Tasks 1-2's data, no new querying.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 4: Row assembly, selection modes, and `bin/ni discover` CLI

**Files:**
- Modify: `lib/ni_discover.py`
- Modify: `bin/ni`
- Modify: `tests/test_ni_discover.py`

**Interfaces:**
- Consumes: everything from Tasks 1-3.
- Produces: `discover_species(
    species: str, study_dir: Path, *, ingroup_groups: list[str] | None,
    outgroup_groups: list[str] | None, auto: bool, include_species_complex: bool
  ) -> None` — the full orchestration: resolve taxon, query genomes
  (+ species-complex genomes if `--include-species-complex`), fetch
  taxonomy ranks (batched, twice per Task 1's design), group, detect
  duplicates, collapse GCA/GCF pairs (`_collapse_paired`), assign Shorts,
  apply the selection mode, write `species.csv`, print `build_report`'s
  output (plus, for `--auto`, the extra proposal content).
- Produces: `bin/ni discover --species "<name>" --study-dir <dir>
  [--ingroup-groups G1,G2] [--outgroup-groups G3] [--auto]
  [--include-species-complex]`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ni_discover.py`. Monkeypatch `query_species_genomes`
and `fetch_taxonomy_ranks` (Tasks 1) so this task's tests don't need live
calls or the frozen fixtures again:

```python
import csv

DISCOVER_HEADER = [
    "Short", "Species", "Strain", "Group", "TaxonGroup",
    "Protein_Source", "Protein_Accession", "Taxon_ID",
    "Genome_Source", "Genome_Accession", "GFF3_Source", "GFF3_Accession",
]


def _fake_genome(accession, taxid, strain, group_after_lookup="no-fsp-in-name"):
    return {
        "accession": accession, "paired_accession": None, "organism_name": "Test species",
        "record_taxon_id": taxid, "strain": strain, "isolate": None,
        "biosample_accession": f"SAMN_{accession}", "assembly_level": "Chromosome",
        "has_annotation": True, "annotation_pipeline": "p", "protein_coding_count": 100,
        "busco_score": 95.0,
    }


def test_discover_species_default_mode_writes_blank_group(tmp_path, monkeypatch):
    study_dir = tmp_path / "studies" / "fungi" / "toy_discover"
    genomes = [_fake_genome("GCA_1", "1", "StrainA"), _fake_genome("GCA_2", "2", "StrainB")]
    monkeypatch.setattr(nd, "resolve_taxon_id", lambda s: [{"taxon_id": "999", "scientific_name": s}])
    monkeypatch.setattr(nd, "query_species_genomes", lambda tid: genomes)
    monkeypatch.setattr(nd, "fetch_taxonomy_ranks", lambda ids: {
        "1": {"rank": "STRAIN", "name": "Test species StrainA", "parents": ["999"], "species_name": "Test species"},
        "2": {"rank": "STRAIN", "name": "Test species StrainB", "parents": ["999"], "species_name": "Test species"},
        "999": {"rank": "SPECIES", "name": "Test species", "parents": [], "species_name": "Test species"},
    })

    nd.discover_species("Test species", study_dir, ingroup_groups=None, outgroup_groups=None, auto=False, include_species_complex=False)

    with open(study_dir / "species.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2
    assert all(r["Group"] == "" for r in rows)
    assert all(r["Protein_Source"] == "ncbi" for r in rows)
    assert all(r["Genome_Source"] == "ncbi" for r in rows)


def test_discover_species_explicit_groups_filters_and_assigns(tmp_path, monkeypatch):
    study_dir = tmp_path / "studies" / "fungi" / "toy_discover2"
    genomes = [
        _fake_genome("GCA_1", "1", "StrainA"),
        _fake_genome("GCA_2", "2", "StrainB"),
    ]
    monkeypatch.setattr(nd, "resolve_taxon_id", lambda s: [{"taxon_id": "999", "scientific_name": s}])
    monkeypatch.setattr(nd, "query_species_genomes", lambda tid: genomes)
    monkeypatch.setattr(nd, "fetch_taxonomy_ranks", lambda ids: {
        "1": {"rank": "FORMA_SPECIALIS", "name": "Test species f. sp. alpha", "parents": ["999"], "species_name": "Test species"},
        "2": {"rank": "STRAIN", "name": "Test species StrainB", "parents": ["999"], "species_name": "Test species"},
        "999": {"rank": "SPECIES", "name": "Test species", "parents": [], "species_name": "Test species"},
    })

    nd.discover_species("Test species", study_dir, ingroup_groups=["alpha"], outgroup_groups=["no-fsp-in-name"], auto=False, include_species_complex=False)

    with open(study_dir / "species.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2
    groups = {r["Group"] for r in rows}
    assert groups == {"IN", "OUT"}


def test_discover_species_unknown_group_label_errors(tmp_path, monkeypatch):
    study_dir = tmp_path / "studies" / "fungi" / "toy_discover3"
    genomes = [_fake_genome("GCA_1", "1", "StrainA")]
    monkeypatch.setattr(nd, "resolve_taxon_id", lambda s: [{"taxon_id": "999", "scientific_name": s}])
    monkeypatch.setattr(nd, "query_species_genomes", lambda tid: genomes)
    monkeypatch.setattr(nd, "fetch_taxonomy_ranks", lambda ids: {
        "1": {"rank": "STRAIN", "name": "Test species StrainA", "parents": ["999"], "species_name": "Test species"},
        "999": {"rank": "SPECIES", "name": "Test species", "parents": [], "species_name": "Test species"},
    })
    try:
        nd.discover_species("Test species", study_dir, ingroup_groups=["nonexistent_group"], outgroup_groups=None, auto=False, include_species_complex=False)
        assert False, "expected SystemExit"
    except SystemExit as e:
        assert "nonexistent_group" in str(e)


def test_discover_species_auto_never_writes_group(tmp_path, monkeypatch, capsys):
    study_dir = tmp_path / "studies" / "fungi" / "toy_discover4"
    genomes = [_fake_genome("GCA_1", "1", "StrainA"), _fake_genome("GCA_2", "2", "StrainB")]
    monkeypatch.setattr(nd, "resolve_taxon_id", lambda s: [{"taxon_id": "999", "scientific_name": s}])
    monkeypatch.setattr(nd, "query_species_genomes", lambda tid: genomes)
    monkeypatch.setattr(nd, "fetch_taxonomy_ranks", lambda ids: {
        "1": {"rank": "FORMA_SPECIALIS", "name": "Test species f. sp. alpha", "parents": ["999"], "species_name": "Test species"},
        "2": {"rank": "STRAIN", "name": "Test species StrainB", "parents": ["999"], "species_name": "Test species"},
        "999": {"rank": "SPECIES", "name": "Test species", "parents": [], "species_name": "Test species"},
    })

    nd.discover_species("Test species", study_dir, ingroup_groups=None, outgroup_groups=None, auto=True, include_species_complex=False)

    with open(study_dir / "species.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert all(r["Group"] == "" for r in rows)  # --auto NEVER writes Group


def test_discover_species_refuses_if_species_csv_exists(tmp_path, monkeypatch):
    study_dir = tmp_path / "studies" / "fungi" / "toy_discover5"
    study_dir.mkdir(parents=True)
    (study_dir / "species.csv").write_text("Short,Species\n")
    try:
        nd.discover_species("Test species", study_dir, ingroup_groups=None, outgroup_groups=None, auto=False, include_species_complex=False)
        assert False, "expected SystemExit"
    except SystemExit as e:
        assert "species.csv" in str(e)


def test_discover_species_both_explicit_and_auto_errors(tmp_path, monkeypatch):
    study_dir = tmp_path / "studies" / "fungi" / "toy_discover6"
    try:
        nd.discover_species("Test species", study_dir, ingroup_groups=["x"], outgroup_groups=None, auto=True, include_species_complex=False)
        assert False, "expected SystemExit"
    except SystemExit:
        pass
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run pytest tests/test_ni_discover.py -v`
Expected: `AttributeError: ... 'discover_species'`.

- [ ] **Step 3: Add `discover_species` to `lib/ni_discover.py`**

```python
def discover_species(
    species: str,
    study_dir: Path,
    *,
    ingroup_groups: list[str] | None,
    outgroup_groups: list[str] | None,
    auto: bool,
    include_species_complex: bool,
) -> None:
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
        writer = csv.DictWriter(fh, fieldnames=[
            "Short", "Species", "Strain", "Group", "TaxonGroup",
            "Protein_Source", "Protein_Accession", "Taxon_ID",
            "Genome_Source", "Genome_Accession", "GFF3_Source", "GFF3_Accession",
        ])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {species_csv} ({len(rows)} species)")
```

Add `import csv` and `from pathlib import Path` (if not already present) to
`lib/ni_discover.py`'s top-of-file imports.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run pytest tests/test_ni_discover.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Wire `discover` into `bin/ni`**

Read the current `bin/ni` (has `resolve` and `fetch` subcommands already).
Add a third subparser:

```python
discover_ap = sub.add_parser("discover", help="Enumerate a species' NCBI genomes into a new study's species.csv")
discover_ap.add_argument("--species", required=True)
discover_ap.add_argument("--study-dir", required=True)
discover_ap.add_argument("--ingroup-groups", default=None, help="Comma-separated group labels")
discover_ap.add_argument("--outgroup-groups", default=None, help="Comma-separated group labels")
discover_ap.add_argument("--auto", action="store_true")
discover_ap.add_argument("--include-species-complex", action="store_true")
```

And in `main()`'s dispatch:

```python
if args.command == "discover":
    from ni_discover import discover_species
    discover_species(
        args.species,
        Path(args.study_dir),
        ingroup_groups=args.ingroup_groups.split(",") if args.ingroup_groups else None,
        outgroup_groups=args.outgroup_groups.split(",") if args.outgroup_groups else None,
        auto=args.auto,
        include_species_complex=args.include_species_complex,
    )
    return 0
```

- [ ] **Step 6: Run the full test suite**

Run: `pixi run pytest tests/ -v`
Expected: all tests pass (Tasks 1-4's `test_ni_discover.py` tests plus the
existing `test_ni_resolve.py`/`test_build_study_config.py`/etc. suites,
unaffected by this addition).

- [ ] **Step 7: Commit**

```bash
git add lib/ni_discover.py bin/ni tests/test_ni_discover.py
git commit -m "Add bin/ni discover CLI: species-population enumeration + selection modes

discover_species() orchestrates the full pipeline (resolve taxon, query
genomes, batched taxonomy rank lookup, forma-specialis grouping,
duplicate collapsing, GCA/GCF pair-collapse, Short assignment) and
implements all three selection modes -- default (Group blank), explicit
--ingroup-groups/--outgroup-groups, and --auto (report-only, never
writes Group), per notes/superpowers/specs/2026-09-11-ni-discover-design.md.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 5: Documentation

**Files:**
- Modify: `README.md`
- Modify: `.claude/skills/new-study/SKILL.md`

**Interfaces:**
- Consumes: `bin/ni discover` (Task 4) — pure documentation.

- [ ] **Step 1: Update `README.md`**

Read the current `README.md` (has a `bin/ni` section with `resolve`/`fetch`
usage and the two pangenome/species-comparison examples already). Add a
`discover` subsection describing when to use it (starting a pangenome study
from scratch by enumerating a species' population, rather than typing out
known strains by hand) with a short worked example:

```bash
pixi run python bin/ni discover --species "Fusarium oxysporum" --study-dir studies/fungi/my_pangenome_study
```
followed by 1-2 sentences on reviewing the printed group table and assigning
`Group` by hand (or using `--ingroup-groups`/`--outgroup-groups` once you
know which forma-specialis groups you want), and a one-line caveat that
forma specialis is a pathotype label, not a phylogenetic split — link to the
design spec for the full reasoning rather than repeating it.

- [ ] **Step 2: Update the onboarding skill**

Read `.claude/skills/new-study/SKILL.md`'s current content (it documents
`bin/ni resolve` already). Add a short note that `bin/ni discover` is an
alternative starting point for pangenome/multi-strain studies when the
species list isn't already known — one paragraph, pointing at the design
spec, not a full re-explanation.

- [ ] **Step 3: Run the full test suite once more**

Run: `pixi run pytest tests/ -v`
Expected: unaffected (doc-only task), all tests still pass.

- [ ] **Step 4: Commit**

```bash
git add README.md .claude/skills/new-study/SKILL.md
git commit -m "Document bin/ni discover in the README and onboarding skill

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Self-Review Notes (for whoever executes this plan)

- Task 1's live-verification step is the highest-risk part of this whole
  plan — the `annotation_info` gene-count/BUSCO field paths and the
  `datasets summary taxonomy taxon` response shape have NOT been pinned by
  any prior live check in this session (unlike the genome-record fields
  `resolve`'s own tasks already nailed down). Do not skip Task 1's Step 0.
- If `datasets summary taxonomy taxon` turns out not to accept multiple ids
  in one invocation (contrary to what this spec's review implied), fall back
  to one call per unique id — still correct, just not literally "2 batched
  calls total"; note the deviation in your Task 1 report rather than
  silently changing the plan's stated architecture.
- Task 4's `discover_species` is long — if it grows unwieldy during
  implementation, splitting the "resolve taxon + query + rank lookups" setup
  half from the "group/filter/write rows" half into two functions is a
  reasonable, plan-consistent refactor; note it in your report rather than
  silently doing it.
