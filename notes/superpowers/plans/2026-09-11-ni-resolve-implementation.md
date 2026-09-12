# `bin/ni resolve` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `bin/ni resolve` — given a `species.csv` with `Short`/`Species`/
`Strain`/`Group`/`TaxonGroup` filled in and source columns blank, it queries
UniProt and NCBI to fill in `Protein_Source`/`Genome_Source`/`*_Accession`, or
leaves a row blank with a reported reason when resolution is genuinely
ambiguous. Also `bin/ni fetch`, a thin pass-through to the existing
`bin/build_study_config.py`.

**Architecture:** All resolution logic lives in `lib/ni_resolve.py` (pure
functions over parsed API responses, unit-testable with canned fixtures — no
live network calls in tests). `bin/ni` is a thin argparse dispatcher with two
subcommands. HTTP via stdlib `urllib.request` (matches
`bin/fetch_uniprot_proteome.py`'s existing convention — no new dependency).
NCBI queries via the `datasets` CLI through `subprocess.run` (matches
`bin/fetch_genome_assembly.py`'s existing convention).

**Tech Stack:** Python 3.12, `pytest` (`pixi run pytest tests/`), stdlib
`urllib.request`/`json`/`subprocess`, no new pixi dependencies.

**Spec:** `notes/superpowers/specs/2026-09-11-ni-dataset-resolution-design.md`

## Global Constraints

- Never guess: any ambiguity (multiple candidates a tiebreak can't settle,
  zero candidates, an unresolvable taxon name) leaves that row's source
  columns blank and reports why — never picks "first" or "most recent" as a
  silent default.
- Query UniProt/NCBI by resolved **taxon ID**, never by free-text organism
  name — a mycovirus named after its fungal host can be a `Reference proteome`
  in its own right (confirmed live during this spec's review) and a
  name-based query returns it as a false candidate.
- Collapse NCBI's GCA/GCF paired assemblies into one candidate before any
  tie-detection — both members of a pair carry
  `assembly_info.refseq_category: reference genome` and are the same
  underlying assembly.
- When an NCBI assembly is about to feed `Protein_Source=ncbi` (no UniProt
  protein for this row), it must have `annotation_info` present — an
  unannotated assembly can't produce a protein FASTA
  (`fetch_genome_assembly.py --include-protein` only warns, doesn't error, on
  a missing protein file).
- Before accepting `Protein_Source=uniprot`, verify the proteome's FTP
  directory actually exists (`HEAD` request) — the live REST API can list a
  proteome not yet in the FTP snapshot `fetch_uniprot_proteome.py` downloads
  from.
- Record `Taxon_ID` for both the `uniprot` and `ncbi` protein paths (not just
  `uniprot`).
- A row with either `Protein_Source` or `Genome_Source` already filled is
  never touched by `resolve` — safe to re-run after partial manual edits.
- Every accepted row's report line states which path produced it (default
  reference match vs. strain-targeted widened match) and flags a strain-string
  mismatch when the record's own strain doesn't contain `species.csv`'s
  `Strain` value.
- Run tests via `pixi run pytest tests/test_ni_resolve.py -v`.

---

### Task 1: UniProt resolution (Step 0 taxon lookup, Step 1 reference search, FTP check)

**Files:**
- Create: `lib/ni_resolve.py`
- Test: `tests/test_ni_resolve.py` (new)

**Interfaces:**
- Produces: `resolve_taxon_id(species: str) -> list[dict]` — each dict
  `{"taxon_id": str, "scientific_name": str}`. Empty list = unresolvable name.
- Produces: `UNIPROT_REST = "https://rest.uniprot.org"`,
  `UNIPROT_FTP_BASE` (same value as `bin/fetch_uniprot_proteome.py`'s
  constant — import it from there rather than redefining, to guarantee they
  never drift apart: `from fetch_uniprot_proteome import UNIPROT_FTP_BASE,
  DOMAIN_FOLDERS`).
- Produces: `search_uniprot_proteomes(taxon_id: str, *, reference_only: bool)
  -> list[dict]` — each dict has at least `proteome_id`, `strain`,
  `busco_score` (float or None), `genome_assembly_id` (str or None),
  `superkingdom` (str). Paginates internally (follows UniProt's cursor
  header) so callers never see a truncated result.
- Produces: `check_ftp_available(proteome_id: str, superkingdom: str) ->
  bool`.
- Produces: `strain_matches(candidate_strain: str | None, wanted_strain: str)
  -> bool` — case-insensitive substring match; `False` if `candidate_strain`
  is `None`/empty.
- Later tasks (2-4) import all of the above from `lib/ni_resolve.py`.

**Before writing code:** run a real query to confirm the exact response shape
you're parsing, since UniProt's JSON is nested and this spec's own review
caught wrong-guessed field names once already:
```bash
curl -s 'https://rest.uniprot.org/proteomes/search?query=taxonomy_id:5334+AND+proteome_type:REFERENCE&format=json' | python3 -m json.tool | head -60
curl -s 'https://rest.uniprot.org/taxonomy/search?query=Neurospora+crassa&fields=taxon_id,scientific_name&format=json' | python3 -m json.tool
```
Confirm the exact JSON paths for: proteome ID, strain (likely under
`taxonomy.strain` or similar — UniProt's schema nests strain info under the
`taxonomy` block in most proteome records), `proteomeCompletenessReport.
buscoReport.score`, `genomeAssembly.assemblyId`, `superkingdom`, and the
response's pagination — UniProt returns a `Link` header with `rel="next"`
containing the next page's URL when more results exist (`resp.headers.get
("Link")`, parse the `<url>; rel="next"` pattern) — if what you find differs
from this task's sketch below, adjust the parsing to match reality, not the
sketch.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_ni_resolve.py`:

```python
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "lib"))
sys.path.insert(0, str(REPO / "bin"))
import ni_resolve as nr  # noqa: E402


def _mock_response(payload, link_header=None):
    """A context-manager mock matching urllib.request.urlopen's return value."""
    resp = MagicMock()
    resp.__enter__.return_value = resp
    resp.read.return_value = json.dumps(payload).encode("utf-8")
    resp.headers = {"Link": link_header} if link_header else {}
    resp.headers.get = (lambda k, default=None: (link_header if k == "Link" else default))
    return resp


def test_resolve_taxon_id_single_match():
    payload = {"results": [{"taxonId": "5334", "scientificName": "Neurospora crassa"}]}
    with patch("urllib.request.urlopen", return_value=_mock_response(payload)):
        hits = nr.resolve_taxon_id("Neurospora crassa")
    assert hits == [{"taxon_id": "5334", "scientific_name": "Neurospora crassa"}]


def test_resolve_taxon_id_no_match():
    payload = {"results": []}
    with patch("urllib.request.urlopen", return_value=_mock_response(payload)):
        hits = nr.resolve_taxon_id("Xylaria sp.")
    assert hits == []


def test_search_uniprot_proteomes_excludes_mycovirus_via_taxon_scoping():
    # Taxon-ID-scoped query should never see the mycovirus in the first place --
    # this test asserts the QUERY URL built uses taxonomy_id, not organism_name,
    # which is the actual fix (the mycovirus contamination happens at the API
    # level for organism_name queries; scoping by taxon_id is what avoids it).
    payload = {
        "results": [
            {
                "id": "UP000001805",
                "taxonomy": {"strain": "ATCC 24698 / 74-OR23-1A / CBS 708.71 / DSM 1257 / FGSC 987"},
                "superkingdom": "Eukaryota",
                "genomeAssembly": {"assemblyId": "GCA_000182925.2"},
                "proteomeCompletenessReport": {"buscoReport": {"score": 98.5}},
            }
        ]
    }
    captured_url = {}

    def fake_urlopen(req, *a, **kw):
        url = req.full_url if hasattr(req, "full_url") else req
        captured_url["url"] = url
        return _mock_response(payload)

    with patch("urllib.request.urlopen", side_effect=fake_urlopen):
        candidates = nr.search_uniprot_proteomes("5334", reference_only=True)

    assert "taxonomy_id:5334" in captured_url["url"]
    assert "organism_name" not in captured_url["url"]
    assert len(candidates) == 1
    assert candidates[0]["proteome_id"] == "UP000001805"
    assert candidates[0]["genome_assembly_id"] == "GCA_000182925.2"
    assert candidates[0]["busco_score"] == 98.5


def test_search_uniprot_proteomes_paginates():
    page1 = {"results": [{"id": "UP1", "taxonomy": {}, "superkingdom": "Eukaryota"}]}
    page2 = {"results": [{"id": "UP2", "taxonomy": {}, "superkingdom": "Eukaryota"}]}
    responses = [
        _mock_response(page1, link_header='<https://rest.uniprot.org/proteomes/search?cursor=abc>; rel="next"'),
        _mock_response(page2),
    ]

    def fake_urlopen(req, *a, **kw):
        return responses.pop(0)

    with patch("urllib.request.urlopen", side_effect=fake_urlopen):
        candidates = nr.search_uniprot_proteomes("999", reference_only=False)

    assert [c["proteome_id"] for c in candidates] == ["UP1", "UP2"]


def test_check_ftp_available_true_on_200():
    resp = MagicMock()
    resp.__enter__.return_value = resp
    resp.status = 200
    with patch("urllib.request.urlopen", return_value=resp):
        assert nr.check_ftp_available("UP000001805", "Eukaryota") is True


def test_check_ftp_available_false_on_404():
    import urllib.error
    with patch("urllib.request.urlopen", side_effect=urllib.error.HTTPError(
        "url", 404, "not found", {}, None
    )):
        assert nr.check_ftp_available("UP001497681", "Eukaryota") is False


def test_strain_matches_substring_case_insensitive():
    assert nr.strain_matches("H4-8 / ATCC 38548", "h4-8") is True
    assert nr.strain_matches("Tattone D", "H4-8") is False
    assert nr.strain_matches(None, "H4-8") is False
    assert nr.strain_matches("", "H4-8") is False
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run pytest tests/test_ni_resolve.py -v`
Expected: `ModuleNotFoundError: No module named 'ni_resolve'`.

- [ ] **Step 3: Write `lib/ni_resolve.py`**

Start the file with this structure (fill in the exact JSON field paths after
running the live queries above — the sketch below is this plan's best-known
shape, not guaranteed byte-for-byte against the live schema):

```python
"""UniProt/NCBI dataset resolution for bin/ni resolve.

See notes/superpowers/specs/2026-09-11-ni-dataset-resolution-design.md for the
full design and the bioinformatics rationale behind each rule enforced here.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from fetch_uniprot_proteome import DOMAIN_FOLDERS, UNIPROT_FTP_BASE  # noqa: E402

UNIPROT_REST = "https://rest.uniprot.org"


def resolve_taxon_id(species: str) -> list[dict[str, str]]:
    """Resolve a species name to NCBI taxon ID(s) via UniProt's taxonomy search.
    Empty list means no resolvable taxon -- caller reports and leaves the row blank."""
    url = (
        f"{UNIPROT_REST}/taxonomy/search"
        f"?query={urllib.parse.quote(species)}&fields=taxon_id,scientific_name&format=json"
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
    """Search UniProt proteomes for a taxon ID. reference_only=True restricts to
    proteome_type:REFERENCE (Step 1); False searches all proteomes for this taxon
    (Step 1b, strain-targeted widening)."""
    query = f"taxonomy_id:{taxon_id}"
    if reference_only:
        query += " AND proteome_type:REFERENCE"
    url = f"{UNIPROT_REST}/proteomes/search?query={urllib.parse.quote(query)}&format=json&size=500"

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
                "strain": r.get("taxonomy", {}).get("strain"),
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
```

Add `import urllib.parse` to the import block (used by the quote calls above).

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run pytest tests/test_ni_resolve.py -v`
Expected: all 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add lib/ni_resolve.py tests/test_ni_resolve.py
git commit -m "Add UniProt taxon/proteome resolution for bin/ni resolve

Step 0 (taxon lookup) and Step 1 (reference-proteome search, taxon-ID
scoped to avoid mycovirus contamination) and the FTP-availability check,
per notes/superpowers/specs/2026-09-11-ni-dataset-resolution-design.md.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 2: NCBI resolution (Step 2 query, GCA/GCF pair collapse, ranking)

**Files:**
- Modify: `lib/ni_resolve.py`
- Modify: `tests/test_ni_resolve.py`

**Interfaces:**
- Consumes: nothing from Task 1 directly (independent query path), but shares
  the module and `strain_matches()`.
- Produces: `query_ncbi_assemblies(taxon_id: str, *, require_reference: bool)
  -> list[dict]` — each dict has `accession` (GCA), `paired_accession` (GCF or
  None), `refseq_category` (str or None), `assembly_level` (str), `strain`
  (str or None), `has_annotation` (bool), `assembly_status` (str), `submitter`
  (str or None), `submission_date` (str or None).
- Produces: `rank_ncbi_candidates(candidates: list[dict], *,
  require_annotation: bool) -> list[dict]` — returns the candidates in the
  top surviving tier (could be 0, 1, or several if truly tied); caller decides
  accept/tie/none from the length of what's returned.
- Later tasks import both.

**Before writing code:** run a real query to confirm the exact JSON field
paths (these are this plan's best guess, pin them against reality):
```bash
datasets summary genome taxon 5334 --as-json-lines | python3 -m json.tool | head -80
```
Confirm the paths for: assembly accession, `assembly_info.refseq_category`,
`assembly_info.assembly_level`, `assembly_info.paired_accession` (or
equivalent pairing field — it may be named differently; find it),
`assembly_info.assembly_status`, `organism.infraspecific_names.strain`,
`annotation_info` (presence = annotated), `assembly_info.submitter`,
`assembly_info.submission_date`. Adjust the sketch below to match.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ni_resolve.py`:

```python
NCRASSA_PAIR_JSONL = "\n".join([
    json.dumps({
        "accession": "GCA_000182925.2",
        "assembly_info": {
            "refseq_category": "reference genome",
            "assembly_level": "Chromosome",
            "paired_accession": "GCF_000182925.2",
            "assembly_status": "current",
            "submitter": "Broad Institute",
            "submission_date": "2014-05-01",
        },
        "organism": {"infraspecific_names": {"strain": "OR74A"}},
        "annotation_info": {"name": "GCF_000182925.2-RS_2014_05"},
    }),
    json.dumps({
        "accession": "GCF_000182925.2",
        "assembly_info": {
            "refseq_category": "reference genome",
            "assembly_level": "Chromosome",
            "paired_accession": "GCA_000182925.2",
            "assembly_status": "current",
            "submitter": "Broad Institute",
            "submission_date": "2014-05-01",
        },
        "organism": {"infraspecific_names": {"strain": "OR74A"}},
        "annotation_info": {"name": "GCF_000182925.2-RS_2014_05"},
    }),
])


def test_query_ncbi_assemblies_parses_jsonl():
    fake_run = MagicMock(stdout=NCRASSA_PAIR_JSONL, returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        candidates = nr.query_ncbi_assemblies("5334", require_reference=False)
    assert len(candidates) == 2
    accessions = {c["accession"] for c in candidates}
    assert accessions == {"GCA_000182925.2", "GCF_000182925.2"}


def test_rank_ncbi_candidates_collapses_gca_gcf_pair_to_one():
    fake_run = MagicMock(stdout=NCRASSA_PAIR_JSONL, returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        candidates = nr.query_ncbi_assemblies("5334", require_reference=False)
    top = nr.rank_ncbi_candidates(candidates, require_annotation=False)
    assert len(top) == 1  # NOT 2 -- this is the pair-collapse fix
    assert top[0]["accession"] == "GCA_000182925.2"  # GCA preferred for recording


def test_rank_ncbi_candidates_requires_annotation_when_asked():
    unannotated = json.dumps({
        "accession": "GCA_999999999.1",
        "assembly_info": {
            "refseq_category": None,
            "assembly_level": "Complete Genome",
            "paired_accession": None,
            "assembly_status": "current",
            "submitter": "Some Lab",
            "submission_date": "2025-01-01",
        },
        "organism": {"infraspecific_names": {}},
    })  # no annotation_info key at all
    fake_run = MagicMock(stdout=unannotated, returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        candidates = nr.query_ncbi_assemblies("1", require_reference=False)
    assert nr.rank_ncbi_candidates(candidates, require_annotation=True) == []
    assert len(nr.rank_ncbi_candidates(candidates, require_annotation=False)) == 1


def test_rank_ncbi_candidates_ties_when_no_reference_designation():
    two_unranked = "\n".join([
        json.dumps({
            "accession": f"GCA_{i}00000000.1",
            "assembly_info": {
                "refseq_category": None, "assembly_level": "Chromosome",
                "paired_accession": None, "assembly_status": "current",
                "submitter": f"Lab {i}", "submission_date": "2025-01-01",
            },
            "organism": {"infraspecific_names": {"strain": f"strain{i}"}},
            "annotation_info": {"name": "x"},
        }) for i in (1, 2)
    ])
    fake_run = MagicMock(stdout=two_unranked, returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        candidates = nr.query_ncbi_assemblies("1", require_reference=False)
    tied = nr.rank_ncbi_candidates(candidates, require_annotation=False)
    assert len(tied) == 2  # genuine tie, not collapsed -- different accessions, no ranking signal
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run pytest tests/test_ni_resolve.py -v`
Expected: 4 new failures (`AttributeError`/`ImportError` on the not-yet-defined functions).

- [ ] **Step 3: Add to `lib/ni_resolve.py`**

```python
def query_ncbi_assemblies(taxon_id: str, *, require_reference: bool) -> list[dict[str, Any]]:
    """Query NCBI Datasets for every assembly of this taxon. require_reference=True
    passes --reference (Step 2's default path); False queries all assemblies
    (Step 2b, strain-targeted widening)."""
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
            "paired_accession": info.get("paired_accession"),
            "refseq_category": info.get("refseq_category"),
            "assembly_level": info.get("assembly_level"),
            "assembly_status": info.get("assembly_status"),
            "strain": r.get("organism", {}).get("infraspecific_names", {}).get("strain"),
            "has_annotation": "annotation_info" in r,
            "submitter": info.get("submitter"),
            "submission_date": info.get("submission_date"),
        })
    return candidates


_ASSEMBLY_LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 0, "Scaffold": 1, "Contig": 2}


def _collapse_paired(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One entry per unique underlying assembly (a GCA/GCF pair collapses to
    one, keyed by whichever of the pair is the GCA -- INSDC accessions start
    with 'GCA_', matching UniProt's own genomeAssembly.assemblyId convention)."""
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run pytest tests/test_ni_resolve.py -v`
Expected: all tests (Task 1's 6 + Task 2's 4) PASS, 10 total.

- [ ] **Step 5: Commit**

```bash
git add lib/ni_resolve.py tests/test_ni_resolve.py
git commit -m "Add NCBI Datasets assembly resolution for bin/ni resolve

query_ncbi_assemblies + rank_ncbi_candidates: collapses GCA/GCF pairs
before tie detection (both carry refseq_category=reference genome for
the same underlying assembly), requires annotation_info when feeding
the Protein_Source=ncbi path, ranks reference > representative > best
complete/chromosome-level.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 3: Strain-targeted widening + per-row orchestration

**Files:**
- Modify: `lib/ni_resolve.py`
- Modify: `tests/test_ni_resolve.py`

**Interfaces:**
- Consumes: `resolve_taxon_id`, `search_uniprot_proteomes`, `check_ftp_available`,
  `strain_matches` (Task 1); `query_ncbi_assemblies`, `rank_ncbi_candidates`
  (Task 2).
- Produces: a `ResolvedRow` dataclass:
  ```python
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
  ```
- Produces: `resolve_row(short: str, species: str, strain: str) -> ResolvedRow`
  — orchestrates Steps 0, 1, 1b, 2, 2b for one row. This is the function
  Task 4's CSV-driving code calls per blank row.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ni_resolve.py`:

```python
def test_resolve_row_uniprot_reference_match(monkeypatch):
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [{"taxon_id": "5334", "scientific_name": s}])
    monkeypatch.setattr(nr, "search_uniprot_proteomes", lambda tid, reference_only: (
        [{"proteome_id": "UP000001805", "strain": "OR74A / FGSC 987", "busco_score": 98.5,
          "genome_assembly_id": "GCA_000182925.2", "superkingdom": "Eukaryota"}]
        if reference_only else []
    ))
    monkeypatch.setattr(nr, "check_ftp_available", lambda pid, sk: True)
    monkeypatch.setattr(nr, "query_ncbi_assemblies", lambda tid, require_reference: [])

    row = nr.resolve_row("Ncra", "Neurospora crassa", "OR74A")

    assert row.resolved is True
    assert row.protein_source == "uniprot"
    assert row.protein_accession == "UP000001805"
    assert row.taxon_id == "5334"
    assert row.genome_source == "ncbi"
    assert row.genome_accession == "GCA_000182925.2"


def test_resolve_row_uniprot_ftp_not_yet_available_falls_back_to_ncbi(monkeypatch):
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [{"taxon_id": "1", "scientific_name": s}])
    monkeypatch.setattr(nr, "search_uniprot_proteomes", lambda tid, reference_only: (
        [{"proteome_id": "UP001497681", "strain": "Tattone D", "busco_score": 97.0,
          "genome_assembly_id": "GCA_023508785.1", "superkingdom": "Eukaryota"}]
        if reference_only else []
    ))
    monkeypatch.setattr(nr, "check_ftp_available", lambda pid, sk: False)  # not in FTP snapshot yet
    monkeypatch.setattr(nr, "query_ncbi_assemblies", lambda tid, require_reference: [
        {"accession": "GCA_023508785.1", "paired_accession": "GCF_023508785.1",
         "refseq_category": "reference genome", "assembly_level": "Chromosome",
         "assembly_status": "current", "strain": "H4-8", "has_annotation": True,
         "submitter": "x", "submission_date": "2022-01-01"},
        {"accession": "GCF_023508785.1", "paired_accession": "GCA_023508785.1",
         "refseq_category": "reference genome", "assembly_level": "Chromosome",
         "assembly_status": "current", "strain": "H4-8", "has_annotation": True,
         "submitter": "x", "submission_date": "2022-01-01"},
    ] if not require_reference else [])

    row = nr.resolve_row("Scom", "Schizophyllum commune", "H4-8")

    assert row.resolved is True
    assert row.protein_source == "ncbi"  # NOT uniprot -- FTP wasn't available
    assert row.protein_accession == "GCA_023508785.1"
    assert row.genome_source == "ncbi"
    assert row.genome_accession == "GCA_023508785.1"
    assert any("FTP" in line for line in row.report_lines)


def test_resolve_row_no_match_anywhere_leaves_blank(monkeypatch):
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [{"taxon_id": "1", "scientific_name": s}])
    monkeypatch.setattr(nr, "search_uniprot_proteomes", lambda tid, reference_only: [])
    monkeypatch.setattr(nr, "query_ncbi_assemblies", lambda tid, require_reference: [])

    row = nr.resolve_row("Xsp", "Xylaria sp.", "")

    assert row.resolved is False
    assert row.protein_source == ""
    assert row.genome_source == ""


def test_resolve_row_unresolvable_taxon_name(monkeypatch):
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [])

    row = nr.resolve_row("Xsp", "Xylaria sp.", "")

    assert row.resolved is False
    assert any("taxon" in line.lower() for line in row.report_lines)


def test_resolve_row_strain_targeted_widening(monkeypatch):
    # Default reference-only UniProt search finds nothing; strain is given;
    # widened (non-reference) search matches by strain name.
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [{"taxon_id": "1", "scientific_name": s}])

    def fake_search(tid, reference_only):
        if reference_only:
            return []
        return [
            {"proteome_id": "UP1", "strain": "Fo47", "busco_score": 95.0,
             "genome_assembly_id": None, "superkingdom": "Eukaryota"},
            {"proteome_id": "UP2", "strain": "race 4", "busco_score": 93.0,
             "genome_assembly_id": None, "superkingdom": "Eukaryota"},
        ]

    monkeypatch.setattr(nr, "search_uniprot_proteomes", fake_search)
    monkeypatch.setattr(nr, "check_ftp_available", lambda pid, sk: True)
    monkeypatch.setattr(nr, "query_ncbi_assemblies", lambda tid, require_reference: [])

    row = nr.resolve_row("Foxy1", "Fusarium oxysporum", "Fo47")

    assert row.resolved is True
    assert row.protein_source == "uniprot"
    assert row.protein_accession == "UP1"
    assert any("strain-targeted" in line.lower() for line in row.report_lines)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run pytest tests/test_ni_resolve.py -v`
Expected: 5 new failures (`AttributeError: module 'ni_resolve' has no attribute 'resolve_row'` etc.)

- [ ] **Step 3: Add to `lib/ni_resolve.py`**

```python
from dataclasses import dataclass, field


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
            widened = query_ncbi_assemblies(taxon_id, require_reference=False)
            matches = [c for c in widened if strain_matches(c["strain"], strain)]
            if len(matches) == 1:
                top = matches
                report.append(f"strain-targeted NCBI match: {strain!r}")

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

    resolved = bool(protein_source and genome_source)
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


def _pick_uniprot_candidate(candidates: list[dict[str, Any]], strain: str, report: list[str]) -> dict[str, Any] | None:
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        matches = [c for c in candidates if strain_matches(c["strain"], strain)]
        if len(matches) == 1:
            return matches[0]
        report.append(f"{len(candidates)} ambiguous UniProt reference proteomes, strain did not disambiguate")
    return None
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run pytest tests/test_ni_resolve.py -v`
Expected: all tests (15 total) PASS.

- [ ] **Step 5: Commit**

```bash
git add lib/ni_resolve.py tests/test_ni_resolve.py
git commit -m "Add per-row resolution orchestration + strain-targeted widening

resolve_row() drives Steps 0/1/1b/2/2b for one species.csv row: UniProt
reference match with FTP-availability check, strain-targeted widening to
non-reference proteomes/assemblies when the default search finds nothing
and Strain is given (the pangenome/multi-strain use case), and NCBI
fallback with the annotation-required rule. Never guesses -- ambiguity
leaves the row unresolved with a reported reason.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 4: `bin/ni` CLI + `species.csv` I/O + report formatting

**Files:**
- Create: `bin/ni`
- Modify: `lib/ni_resolve.py` (add `resolve_species_csv`)
- Modify: `tests/test_ni_resolve.py`

**Interfaces:**
- Consumes: `resolve_row` (Task 3).
- Produces: `resolve_species_csv(study_dir: Path) -> None` — reads
  `study_dir/species.csv`, calls `resolve_row` for every row with both source
  columns blank, rewrites the file in place, prints the summary to stdout.
- Produces: `bin/ni`'s CLI surface: `bin/ni resolve --study-dir <dir>`,
  `bin/ni fetch --study-dir <dir> [-- extra build_study_config.py args]`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ni_resolve.py`:

```python
import csv


HEADER = [
    "Short", "Species", "Strain", "Group", "TaxonGroup",
    "Protein_Source", "Protein_Accession", "Taxon_ID",
    "Genome_Source", "Genome_Accession", "GFF3_Source", "GFF3_Accession",
]


def _write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER)
        w.writeheader()
        w.writerows(rows)


def test_resolve_species_csv_fills_blank_rows_only(tmp_path, monkeypatch, capsys):
    study_dir = tmp_path / "studies" / "fungi" / "toy"
    study_dir.mkdir(parents=True)
    _write_csv(study_dir / "species.csv", [
        {  # already filled -- must NOT be touched
            "Short": "Already", "Species": "X", "Strain": "", "Group": "IN", "TaxonGroup": "X",
            "Protein_Source": "local_faa", "Protein_Accession": "/some/path.faa", "Taxon_ID": "",
            "Genome_Source": "local_genome", "Genome_Accession": "/some/path.fna",
            "GFF3_Source": "", "GFF3_Accession": "",
        },
        {  # blank -- should resolve
            "Short": "Ncra", "Species": "Neurospora crassa", "Strain": "OR74A", "Group": "OUT",
            "TaxonGroup": "Pezizomycotina",
            "Protein_Source": "", "Protein_Accession": "", "Taxon_ID": "",
            "Genome_Source": "", "Genome_Accession": "", "GFF3_Source": "", "GFF3_Accession": "",
        },
    ])

    def fake_resolve_row(short, species, strain):
        assert short == "Ncra"
        return nr.ResolvedRow(
            short=short, resolved=True, protein_source="uniprot",
            protein_accession="UP000001805", taxon_id="5334",
            genome_source="ncbi", genome_accession="GCA_000182925.2",
            report_lines=[],
        )

    monkeypatch.setattr(nr, "resolve_row", fake_resolve_row)
    nr.resolve_species_csv(study_dir)

    with open(study_dir / "species.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))

    already = next(r for r in rows if r["Short"] == "Already")
    assert already["Protein_Source"] == "local_faa"  # untouched

    ncra = next(r for r in rows if r["Short"] == "Ncra")
    assert ncra["Protein_Source"] == "uniprot"
    assert ncra["Protein_Accession"] == "UP000001805"
    assert ncra["Genome_Source"] == "ncbi"
    assert ncra["Genome_Accession"] == "GCA_000182925.2"

    out = capsys.readouterr().out
    assert "Resolved 1 of 1" in out


def test_resolve_species_csv_leaves_unresolved_blank_and_reports(tmp_path, monkeypatch, capsys):
    study_dir = tmp_path / "studies" / "fungi" / "toy2"
    study_dir.mkdir(parents=True)
    _write_csv(study_dir / "species.csv", [
        {"Short": "Xsp", "Species": "Xylaria sp.", "Strain": "", "Group": "OUT", "TaxonGroup": "X",
         "Protein_Source": "", "Protein_Accession": "", "Taxon_ID": "",
         "Genome_Source": "", "Genome_Accession": "", "GFF3_Source": "", "GFF3_Accession": ""},
    ])

    monkeypatch.setattr(nr, "resolve_row", lambda short, species, strain: nr.ResolvedRow(
        short=short, resolved=False, report_lines=["not a resolvable NCBI taxon name"]
    ))
    nr.resolve_species_csv(study_dir)

    with open(study_dir / "species.csv", newline="") as fh:
        row = next(csv.DictReader(fh))
    assert row["Protein_Source"] == ""

    out = capsys.readouterr().out
    assert "need a decision" in out
    assert "not a resolvable NCBI taxon name" in out
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run pytest tests/test_ni_resolve.py -v`
Expected: 2 new failures (`AttributeError: ... 'resolve_species_csv'`).

- [ ] **Step 3: Add `resolve_species_csv` to `lib/ni_resolve.py`**

```python
import csv as _csv


def resolve_species_csv(study_dir: Path) -> None:
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
```

All of `resolve_species_csv`'s output goes to stdout (plain `print()`, no
`file=sys.stderr`) — this is a report a human reads, not a progress log
alongside piped data, so it doesn't need `bin/build_study_config.py`'s
stdout/stderr split.

- [ ] **Step 4: Create `bin/ni`**

```python
#!/usr/bin/env python3
"""ni: dataset resolution and config assembly for a new NII study.

See notes/superpowers/specs/2026-09-11-ni-dataset-resolution-design.md.

Usage:
  bin/ni resolve --study-dir studies/<domain>/<set>
  bin/ni fetch    --study-dir studies/<domain>/<set> [-- extra build_study_config.py args]
"""
import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from ni_resolve import resolve_species_csv  # noqa: E402

BIN = Path(__file__).parent


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    resolve_ap = sub.add_parser("resolve", help="Fill in species.csv's source columns")
    resolve_ap.add_argument("--study-dir", required=True)

    fetch_ap = sub.add_parser("fetch", help="Pass through to bin/build_study_config.py")
    fetch_ap.add_argument("--study-dir", required=True)
    fetch_ap.add_argument("extra_args", nargs=argparse.REMAINDER)

    args = ap.parse_args()

    if args.command == "resolve":
        resolve_species_csv(Path(args.study_dir))
        return 0

    if args.command == "fetch":
        cmd = [sys.executable, str(BIN / "build_study_config.py"), "--study-dir", args.study_dir]
        cmd.extend(a for a in args.extra_args if a != "--")
        result = subprocess.run(cmd)
        return result.returncode

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
```

```bash
chmod +x bin/ni
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `pixi run pytest tests/test_ni_resolve.py -v`
Expected: all tests PASS (17 total).

- [ ] **Step 6: Manual smoke test of `bin/ni fetch`'s pass-through**

```bash
pixi run python bin/ni fetch --study-dir studies/bacteria/UHM_lachnoNovelclade -- --skip-fetch
```
Expected: identical output to running `bin/build_study_config.py --study-dir studies/bacteria/UHM_lachnoNovelclade --skip-fetch` directly (this study's data is already fetched, so `--skip-fetch` makes this a fast, side-effect-free check that the pass-through wiring is correct).

- [ ] **Step 7: Commit**

```bash
git add bin/ni lib/ni_resolve.py tests/test_ni_resolve.py
git commit -m "Add bin/ni CLI: resolve (species.csv fill-in) and fetch (pass-through)

bin/ni resolve drives resolve_species_csv() over an existing species.csv,
touching only blank rows and printing a resolved/unresolved summary.
bin/ni fetch execs the existing bin/build_study_config.py unchanged --
no new download logic, per the spec's explicit resolve-vs-fetch separation.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 5: Documentation — wire `bin/ni` into the onboarding skill

**Files:**
- Modify: `.claude/skills/new-study/SKILL.md`
- Modify: `DESIGN.md` (add Phase 2 to "Explicitly deferred" section)

**Interfaces:**
- Consumes: `bin/ni resolve`/`bin/ni fetch` (Task 4) — pure documentation, no
  new code.

- [ ] **Step 1: Update the onboarding skill**

Read `.claude/skills/new-study/SKILL.md`'s current step 2 (source
classification) and step 5 (`bin/build_study_config.py` invocation). Add a
new step between them (renumber accordingly): once `species.csv` has
`Short`/`Species`/`Strain`/`Group`/`TaxonGroup` filled in, run `bin/ni resolve
--study-dir studies/<domain>/<set_name>` to auto-fill the source columns for
any species classified as needing a UniProt/NCBI lookup (as opposed to
`local_faa`/`local_genome`, which the skill still classifies by hand per its
existing step 2 — `bin/ni resolve` only touches blank rows, so this is
additive, not a replacement for the manual classification of local-file
species). Note explicitly: rows `bin/ni resolve` reports as needing a
decision must be resolved by hand (following the spec's own guidance — check
the reported candidates, pick one, or determine the species genuinely has no
usable public data yet) before moving on to `bin/ni fetch`.

- [ ] **Step 2: Update `DESIGN.md`**

Find the "9. Explicitly deferred / open items" section. Add a bullet noting
Phase 2 (open-ended clade/taxon discovery — "find good Ascomycota outgroups")
is deferred, with a pointer to
`notes/superpowers/specs/2026-09-11-ni-dataset-resolution-design.md`'s own
"Scope: Phase 1 only" section for the reasoning.

- [ ] **Step 3: Run the full test suite once more**

Run: `pixi run pytest tests/ -v`
Expected: all tests pass (this task changes no code, but confirms nothing
about the doc edits broke anything, and closes out the plan with a clean
suite).

- [ ] **Step 4: Commit**

```bash
git add .claude/skills/new-study/SKILL.md DESIGN.md
git commit -m "Wire bin/ni resolve into the new-study onboarding skill

Adds bin/ni resolve as the step between manual species classification
and bin/build_study_config.py, for species needing a UniProt/NCBI lookup
rather than an already-known local file. Notes Phase 2 (open-ended clade
discovery) as deferred in DESIGN.md, per the ni spec's own scope section.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Self-Review Notes (for whoever executes this plan)

- Task 1/2's "before writing code" live-query steps are load-bearing, not
  optional — this spec's own review process already caught one wrong-guessed
  UniProt query shape (`organism_name` vs `taxonomy_id`) and one wrong
  assumption about NCBI's `refseq_category` (GCA/GCF pairs both carrying
  `reference genome`). Do not skip verifying the sketch's field paths against
  a real live call before treating them as correct.
- No task in this plan performs a live end-to-end resolution against real
  UniProt/NCBI data — every test mocks the network/subprocess boundary. A
  reasonable follow-up (not blocking this plan) is a one-time manual smoke
  test of `bin/ni resolve` against a real, small `species.csv` (e.g. a copy
  of a few `UHM_lachnoNovelclade` rows with source columns blanked) to confirm
  the live API shapes still match what Tasks 1-2's fixtures assume, before
  relying on `bin/ni resolve` for a real new study.
