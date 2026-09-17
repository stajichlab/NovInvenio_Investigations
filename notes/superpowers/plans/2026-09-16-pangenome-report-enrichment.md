# Pangenome Report Enrichment Bundle Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add 4 independently-gated report-enrichment features to `nf_NovInvenio`'s
pangenome island/Pfam-enrichment step: InterPro hotlinks, Pfam2GO GO-term
annotation, island genomic-locus assignment, and a per-strain statistical
outlier flag.

**Architecture:** Two in-place augmentations of existing script outputs
(InterPro links, outlier flag) plus one new script + new gated Nextflow process
(Pfam2GO) plus one new function folded into the existing `REPORT_TABLES`
process (island locus). One final wiring task threads a new
`enrichment_for_report` channel through `REPORT_TABLES`/`REPORT_RENDER` and adds
two new inputs (`cluster_tsv`, `gene_positions`) to `REPORT_TABLES`.

**Tech Stack:** Python 3 (argparse, csv, statistics/numpy), Nextflow DSL2,
pytest, pixi.

**Spec:** `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/notes/superpowers/specs/2026-09-16-pangenome-report-enrichment-design.md`

## Global Constraints

- A Nextflow process can only be invoked ONCE per workflow — multiplicity via
  channel fan-out (`Channel.fromList`), never a Groovy loop.
- New gating param `pangenome_pfam2go` must not collide with any existing param.
- No hand-written narrative/prose in the generated report — every new line in
  `render_report_markdown()` is templated from real data (a literal listing
  driven by a TSV, not invented text).
- Real 10-strain execution smoke test required for the final wiring task
  (Task 5) — `--help` alone is not sufficient. No committed 10-strain fixture
  exists in the repo; Task 5 must build one (recipe given in Task 5).
- `pfam_accession` stays versioned everywhere it already exists (do not break
  `tests/test_pangenome_domain_enrichment.py:96`'s `"PF13577.9"` assertion); a
  new shared `bare_pfam_accession()` helper handles stripping wherever a bare
  accession is actually needed (URL building, pfam2go lookup).
- Repo sentinel convention for "value could not be resolved" is `-`, not `""`
  or a crash.
- `pangenome_id_sep` (existing param, default `'|'`) must be threaded through
  anywhere a family/member ID gets split on a separator — never hardcode `"|"`
  (this is the F6 fix already applied elsewhere in this codebase; do not
  reintroduce the bug it fixed).
- `lib/compressed_io.py` is untouched by this plan.

---

### Task 1: InterPro hotlinks for Pfam domains

**Files:**
- Modify: `bin/pangenome_domain_enrichment.py`
- Modify: `bin/pangenome_report_render.py`
- Test: `tests/test_pangenome_domain_enrichment.py`
- Test: `tests/test_pangenome_report_render.py`

**Interfaces:**
- Produces: `bare_pfam_accession(accession: str) -> str` in
  `pangenome_domain_enrichment.py` — importable by Task 2's
  `pangenome_pfam2go.py`. A new `pfam_url` column in
  `island_pfam_enrichment.tsv`, immediately after the existing `pfam_accession`
  column.
- Consumes: nothing new (uses the existing `pfam_accession` column already
  produced by `parse_domtblout_accessions()`).

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_pangenome_domain_enrichment.py`:

```python
from pangenome_domain_enrichment import bare_pfam_accession


def test_bare_pfam_accession_strips_version_suffix():
    assert bare_pfam_accession("PF13577.9") == "PF13577"
    assert bare_pfam_accession("PF00001.1") == "PF00001"


def test_bare_pfam_accession_passes_through_sentinel():
    assert bare_pfam_accession("-") == "-"
```

Update the existing header-assertion test (currently asserting the 8-column
header ending in `fisher_p`, `fdr_q`) to expect `pfam_url` immediately after
`pfam_accession`:

```python
    assert header == [
        "domain", "pfam_accession", "pfam_url", "n_with_domain_in_islands",
        "n_with_domain_in_background", "n_island_families",
        "n_background_families", "fisher_p", "fdr_q",
    ]
    data_row = lines[1].split("\t")
    assert data_row[0] == "SnoaL_2"
    assert data_row[1] == "PF13577.9"
    assert data_row[2] == "https://www.ebi.ac.uk/interpro/entry/pfam/PF13577/"
```

Add to `tests/test_pangenome_report_render.py`:

```python
def test_domain_table_renders_pfam_url_as_link():
    md = render_report_markdown(
        {}, {}, {}, [{"domain": "SnoaL_2", "pfam_accession": "PF13577.9",
                      "pfam_url": "https://www.ebi.ac.uk/interpro/entry/pfam/PF13577/",
                      "fisher_p": "1e-5", "fdr_q": "2e-5"}],
        0, None, None, [], None,
    )
    assert "[SnoaL_2](https://www.ebi.ac.uk/interpro/entry/pfam/PF13577/)" in md


def test_domain_table_falls_back_to_bare_name_without_pfam_url_key():
    # Existing-style row with no pfam_url key at all -- must not KeyError.
    md = render_report_markdown(
        {}, {}, {}, [{"domain": "SnoaL_2", "fisher_p": "1e-5", "fdr_q": "2e-5"}],
        0, None, None, [], None,
    )
    assert "| SnoaL_2 |" in md


def test_domain_table_falls_back_to_bare_name_when_pfam_url_is_sentinel():
    md = render_report_markdown(
        {}, {}, {}, [{"domain": "SnoaL_2", "pfam_url": "-",
                      "fisher_p": "1e-5", "fdr_q": "2e-5"}],
        0, None, None, [], None,
    )
    assert "| SnoaL_2 |" in md
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest tests/test_pangenome_domain_enrichment.py tests/test_pangenome_report_render.py -v`
Expected: FAIL (`bare_pfam_accession` not defined, header mismatch, no `pfam_url` rendering).

- [ ] **Step 3: Implement `bare_pfam_accession()` and the `pfam_url` column**

In `bin/pangenome_domain_enrichment.py`, add near the top (after imports):

```python
def bare_pfam_accession(accession: str) -> str:
    """Strip a Pfam accession's version suffix (PF13577.9 -> PF13577).
    Passes '-' (the unresolved-accession sentinel) through unchanged."""
    if accession == "-" or "." not in accession:
        return accession
    return accession.split(".", 1)[0]
```

In `main()`, change the output-writing block to add `pfam_url`:

```python
    with open(args.output, "w") as out:
        out.write(
            "domain\tpfam_accession\tpfam_url\tn_with_domain_in_islands\t"
            "n_with_domain_in_background\tn_island_families\tn_background_families\t"
            "fisher_p\tfdr_q\n"
        )
        for r in enrichment:
            accession = domain_accessions.get(r["domain"], "-")
            if accession == "-":
                pfam_url = "-"
            else:
                pfam_url = f"https://www.ebi.ac.uk/interpro/entry/pfam/{bare_pfam_accession(accession)}/"
            out.write(
                f"{r['domain']}\t{accession}\t{pfam_url}\t"
                f"{r['n_with_domain_in_islands']}\t{r['n_with_domain_in_background']}\t"
                f"{r['n_island_families']}\t{r['n_background_families']}\t"
                f"{r['fisher_p']:.3e}\t{r['fdr_q']:.3e}\n"
            )
```

- [ ] **Step 4: Implement the render-side link formatting**

In `bin/pangenome_report_render.py`'s `render_report_markdown()`, replace the
"Pfam domain enrichment" table's row-building loop (currently
`f"| {row['domain']} | {float(row['fisher_p']):.2e} | {float(row['fdr_q']):.2e} |"`)
with:

```python
        for row in top_domains:
            pfam_url = row.get("pfam_url")
            domain_cell = f"[{row['domain']}]({pfam_url})" if pfam_url and pfam_url != "-" else row["domain"]
            lines.append(f"| {domain_cell} | {float(row['fisher_p']):.2e} | {float(row['fdr_q']):.2e} |")
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pixi run pytest tests/test_pangenome_domain_enrichment.py tests/test_pangenome_report_render.py -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add bin/pangenome_domain_enrichment.py bin/pangenome_report_render.py \
        tests/test_pangenome_domain_enrichment.py tests/test_pangenome_report_render.py
git commit -m "pangenome: add InterPro hotlinks for Pfam domains in the report"
```

---

### Task 2: Pfam2GO annotation script

**Files:**
- Create: `bin/pangenome_pfam2go.py`
- Test: `tests/test_pangenome_pfam2go.py`
- Modify: `bin/pangenome_report_render.py` (GO-column rendering, unit-testable
  independent of any Nextflow wiring — Task 5 wires the actual pipeline)
- Test: `tests/test_pangenome_report_render.py`

**Interfaces:**
- Consumes: `bare_pfam_accession()` from Task 1's `pangenome_domain_enrichment.py`.
- Produces: CLI `pangenome_pfam2go.py --island_pfam_enrichment PATH --pfam2go PATH --output PATH`,
  writing the input table plus two new columns `n_go_terms`, `go_terms`. Task 5
  wires this into the workflow as a new `PFAM2GO` process.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_pangenome_pfam2go.py`:

```python
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
import pangenome_pfam2go


PFAM2GO_FIXTURE = """\
!date: today
Pfam:PF13577 SnoaL_2 > GO:oxidoreductase activity ; GO:0016491
Pfam:PF13577 SnoaL_2 > GO:metabolic process ; GO:0008152
Pfam:PF00001 7tm_1 > GO:G protein-coupled receptor activity ; GO:0004930
"""


def test_parse_pfam2go_maps_bare_accession_to_go_ids(tmp_path):
    p = tmp_path / "pfam2go"
    p.write_text(PFAM2GO_FIXTURE)
    mapping = pangenome_pfam2go.parse_pfam2go(str(p))
    assert mapping["PF13577"] == ["GO:0016491", "GO:0008152"]
    assert mapping["PF00001"] == ["GO:0004930"]


def test_main_strips_version_suffix_before_lookup(tmp_path):
    pfam2go_path = tmp_path / "pfam2go"
    pfam2go_path.write_text(PFAM2GO_FIXTURE)
    enrichment_path = tmp_path / "island_pfam_enrichment.tsv"
    enrichment_path.write_text(
        "domain\tpfam_accession\tpfam_url\tfisher_p\tfdr_q\n"
        "SnoaL_2\tPF13577.9\thttp://example/PF13577/\t1e-5\t2e-5\n"
        "unknown_domain\t-\t-\t1e-2\t2e-2\n"
    )
    output_path = tmp_path / "output.tsv"
    argv = [
        "pangenome_pfam2go.py",
        "--island_pfam_enrichment", str(enrichment_path),
        "--pfam2go", str(pfam2go_path),
        "--output", str(output_path),
    ]
    import sys as _sys
    old_argv = _sys.argv
    _sys.argv = argv
    try:
        pangenome_pfam2go.main()
    finally:
        _sys.argv = old_argv

    with open(output_path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["domain"] == "SnoaL_2"
    assert rows[0]["n_go_terms"] == "2"
    assert rows[0]["go_terms"] == "GO:0016491;GO:0008152"
    assert rows[1]["domain"] == "unknown_domain"
    assert rows[1]["n_go_terms"] == "0"
    assert rows[1]["go_terms"] == "-"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest tests/test_pangenome_pfam2go.py -v`
Expected: FAIL (`pangenome_pfam2go` module not found).

- [ ] **Step 3: Implement `bin/pangenome_pfam2go.py`**

```python
#!/usr/bin/env python3
"""GO-term annotation of Pfam domains already found by this pipeline's own
hmmscan (island_pfam_enrichment.tsv) via a standard pfam2go mapping file
(http://current.geneontology.org/ontology/external2go/pfam2go). Ported from
studies/fungi/Afumigatus_pangenome/bin/map_pfam_to_go.py -- annotates only
domains already discovered upstream, does not run a fresh InterProScan.

Usage:
  pangenome_pfam2go.py --island_pfam_enrichment island_pfam_enrichment.tsv \\
      --pfam2go pfam2go --output island_pfam_enrichment.go.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pangenome_domain_enrichment import bare_pfam_accession  # noqa: E402

PFAM2GO_LINE_RE = re.compile(r"^Pfam:(PF\d+)\s.*?;\s*(GO:\d+)\s*$")


def parse_pfam2go(path: str) -> dict[str, list[str]]:
    """{bare_pfam_accession: [go_id, ...]} from a standard pfam2go file.
    One line per Pfam-accession/GO-term pair -- a domain with N GO terms
    appears on N separate lines."""
    mapping: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            m = PFAM2GO_LINE_RE.match(line.strip())
            if not m:
                continue
            accession, go_id = m.group(1), m.group(2)
            mapping.setdefault(accession, []).append(go_id)
    return mapping


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--island_pfam_enrichment", required=True)
    ap.add_argument("--pfam2go", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    pfam2go = parse_pfam2go(args.pfam2go)
    print(f"pangenome_pfam2go: {len(pfam2go)} Pfam accessions with >=1 GO term in {args.pfam2go}",
          file=sys.stderr)

    with open(args.island_pfam_enrichment, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or []) + ["n_go_terms", "go_terms"]

    n_annotated = 0
    for row in rows:
        go_ids = pfam2go.get(bare_pfam_accession(row.get("pfam_accession", "-")), [])
        row["n_go_terms"] = str(len(go_ids))
        row["go_terms"] = ";".join(go_ids) if go_ids else "-"
        if go_ids:
            n_annotated += 1

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"pangenome_pfam2go: {n_annotated}/{len(rows)} domains got >=1 GO term",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

Make it executable: `chmod +x bin/pangenome_pfam2go.py`.

- [ ] **Step 4: Add the render-side GO-column test and implementation**

Add to `tests/test_pangenome_report_render.py`:

```python
def test_domain_table_adds_go_columns_when_present():
    md = render_report_markdown(
        {}, {}, {}, [{"domain": "SnoaL_2", "fisher_p": "1e-5", "fdr_q": "2e-5",
                      "n_go_terms": "2", "go_terms": "GO:0016491;GO:0008152"}],
        0, None, None, [], None,
    )
    assert "GO terms" in md
    assert "GO:0016491;GO:0008152" in md


def test_domain_table_omits_go_columns_when_absent():
    md = render_report_markdown(
        {}, {}, {}, [{"domain": "SnoaL_2", "fisher_p": "1e-5", "fdr_q": "2e-5"}],
        0, None, None, [], None,
    )
    assert "GO terms" not in md
```

In `bin/pangenome_report_render.py`'s `render_report_markdown()`, extend the
"Pfam domain enrichment" table to add GO columns when present on ANY row:

```python
    lines += ["## Pfam domain enrichment", ""]
    if not top_domains:
        lines += ["No significantly enriched Pfam domains found.", ""]
    else:
        lines += ["![Top enriched domains](figures/island_domain_enrichment.png)", ""]
        has_go = any(row.get("go_terms") for row in top_domains)
        if has_go:
            lines += ["| Domain | Fisher p | FDR q | GO terms |", "|---|---|---|---|"]
        else:
            lines += ["| Domain | Fisher p | FDR q |", "|---|---|---|"]
        for row in top_domains:
            pfam_url = row.get("pfam_url")
            domain_cell = f"[{row['domain']}]({pfam_url})" if pfam_url and pfam_url != "-" else row["domain"]
            cells = f"| {domain_cell} | {float(row['fisher_p']):.2e} | {float(row['fdr_q']):.2e} |"
            if has_go:
                cells += f" {row.get('go_terms', '-')} |"
            lines.append(cells)
        lines.append("")
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pixi run pytest tests/test_pangenome_pfam2go.py tests/test_pangenome_report_render.py -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add bin/pangenome_pfam2go.py bin/pangenome_report_render.py \
        tests/test_pangenome_pfam2go.py tests/test_pangenome_report_render.py
git commit -m "pangenome: add Pfam2GO GO-term annotation script"
```

---

### Task 3: Island genomic-locus assignment (folded into REPORT_TABLES)

**Files:**
- Modify: `bin/pangenome_report_tables.py`
- Modify: `bin/pangenome_report_render.py` (top-islands table)
- Test: `tests/test_pangenome_report_tables.py`
- Test: `tests/test_pangenome_report_render.py`

**Interfaces:**
- Consumes: `lib/pangenome_matrix.read_cluster_tsv` (existing).
- Produces: `add_island_locus(islands_rows, member_to_rep, gene_positions, id_sep="|")`
  in `pangenome_report_tables.py` — adds 6 columns to each row:
  `locus_id, locus_contig, locus_start, locus_end, n_members_with_coordinates,
  n_contigs_in_locus`. New CLI args on `pangenome_report_tables.py`:
  `--cluster_tsv`, `--gene_positions`, `--id_sep` (default `"|"`). Task 5's
  wiring consumes these.

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_pangenome_report_tables.py`:

```python
def test_add_island_locus_computes_span_from_gene_positions():
    islands_rows = [{
        "n_strains": "1", "example_strain": "S1", "island_size": "2",
        "member_families": "famA|famB", "n_supporting_pairs": "1",
        "classifications": "starship_explained",
    }]
    # famA's rep is a different strain's protein; S1's own copy is protS1_a.
    member_to_rep = {"protS1_a": "famA_rep", "protS1_b": "famB_rep",
                      "famA_rep": "famA_rep", "famB_rep": "famB_rep"}
    gene_positions = {
        ("S1", "protS1_a"): {"contig": "contig1", "start": 100, "end": 200},
        ("S1", "protS1_b"): {"contig": "contig1", "start": 300, "end": 400},
    }
    rows = add_island_locus(islands_rows, member_to_rep, gene_positions, id_sep="|")
    assert rows[0]["locus_id"] == "S1:contig1:100-400"
    assert rows[0]["locus_contig"] == "contig1"
    assert rows[0]["n_members_with_coordinates"] == 2
    assert rows[0]["n_contigs_in_locus"] == 1


def test_add_island_locus_flags_multi_contig():
    islands_rows = [{
        "n_strains": "1", "example_strain": "S1", "island_size": "2",
        "member_families": "famA|famB", "n_supporting_pairs": "1",
        "classifications": "starship_explained",
    }]
    member_to_rep = {"protS1_a": "famA_rep", "protS1_b": "famB_rep"}
    gene_positions = {
        ("S1", "protS1_a"): {"contig": "contig1", "start": 100, "end": 200},
        ("S1", "protS1_b"): {"contig": "contig2", "start": 10, "end": 50},
    }
    rows = add_island_locus(islands_rows, member_to_rep, gene_positions, id_sep="|")
    assert rows[0]["n_contigs_in_locus"] == 2


def test_add_island_locus_sentinel_on_total_failure():
    islands_rows = [{
        "n_strains": "1", "example_strain": "S1", "island_size": "1",
        "member_families": "famA", "n_supporting_pairs": "0",
        "classifications": "unexplained_physical",
    }]
    rows = add_island_locus(islands_rows, {}, {}, id_sep="|")
    assert rows[0]["locus_id"] == "-"
    assert rows[0]["n_members_with_coordinates"] == 0


def test_add_island_locus_respects_custom_id_sep():
    islands_rows = [{
        "n_strains": "1", "example_strain": "S1", "island_size": "1",
        "member_families": "famA", "n_supporting_pairs": "0",
        "classifications": "unexplained_physical",
    }]
    member_to_rep = {"protS1_a": "famA_rep"}
    gene_positions = {("S1", "protS1_a"): {"contig": "contig1", "start": 5, "end": 50}}
    rows = add_island_locus(islands_rows, member_to_rep, gene_positions, id_sep="_")
    assert rows[0]["locus_id"] == "S1:contig1:5-50"
```

Add to `tests/test_pangenome_report_render.py`:

```python
def test_accessory_islands_section_adds_top_islands_table_when_locus_present():
    md = render_report_markdown(
        {}, {}, {}, [], 1, None, None, [], None,
        islands_with_domains_rows=[{
            "locus_id": "S1:contig1:100-400", "island_size": "5",
            "n_strains": "2", "pfam_domains": "SnoaL_2",
        }],
    )
    assert "S1:contig1:100-400" in md


def test_accessory_islands_section_omits_top_islands_table_without_locus():
    md = render_report_markdown(
        {}, {}, {}, [], 1, None, None, [], None,
        islands_with_domains_rows=[{"locus_id": "-", "island_size": "5",
                                     "n_strains": "2", "pfam_domains": "-"}],
    )
    assert "Top islands" not in md
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest tests/test_pangenome_report_tables.py tests/test_pangenome_report_render.py -v`
Expected: FAIL (`add_island_locus` not defined; `render_report_markdown` doesn't
accept `islands_with_domains_rows`).

- [ ] **Step 3: Implement `add_island_locus()` in `bin/pangenome_report_tables.py`**

```python
def add_island_locus(
    islands_rows: list[dict],
    member_to_rep: dict[str, str],
    gene_positions: dict[tuple[str, str], dict],
    id_sep: str = "|",
) -> list[dict]:
    """Adds locus_id/locus_contig/locus_start/locus_end/
    n_members_with_coordinates/n_contigs_in_locus to each island row, by
    resolving member_families (rep-protein IDs) down to the island's own
    example_strain's actual proteins via member_to_rep (rep -> member
    inverted, restricted to that strain), then spanning gene_positions.
    Members with no resolvable coordinate (e.g. a rescue-pass genome_only
    call with no annotated protein_id) are excluded from the span and
    counted, not treated as an error. n_contigs_in_locus > 1 is reported,
    not silently collapsed (probable paralog-copy pull-in)."""
    # Invert member_to_rep (member -> rep) into rep -> [members on this strain].
    rep_to_members: dict[str, list[str]] = {}
    for member, rep in member_to_rep.items():
        rep_to_members.setdefault(rep, []).append(member)

    out = []
    for row in islands_rows:
        strain = row["example_strain"]
        starts, ends, contigs = [], [], set()
        n_resolved = 0
        for family in row["member_families"].split(id_sep):
            resolved_member = None
            for member in rep_to_members.get(family, [family]):
                if (strain, member) in gene_positions:
                    resolved_member = member
                    break
            if resolved_member is None:
                continue
            pos = gene_positions[(strain, resolved_member)]
            starts.append(pos["start"])
            ends.append(pos["end"])
            contigs.add(pos["contig"])
            n_resolved += 1

        new_row = dict(row)
        if n_resolved == 0:
            new_row.update({
                "locus_id": "-", "locus_contig": "-", "locus_start": "-",
                "locus_end": "-", "n_members_with_coordinates": 0,
                "n_contigs_in_locus": 0,
            })
        else:
            contig = sorted(contigs)[0]
            new_row.update({
                "locus_id": f"{strain}:{contig}:{min(starts)}-{max(ends)}",
                "locus_contig": contig, "locus_start": min(starts),
                "locus_end": max(ends), "n_members_with_coordinates": n_resolved,
                "n_contigs_in_locus": len(contigs),
            })
        out.append(new_row)
    return out
```

- [ ] **Step 4: Wire `add_island_locus()` into `main()`, add CLI args, fix the
zero-island fallback header**

In `bin/pangenome_report_tables.py`'s `main()`, add args:

```python
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--gene_positions", required=True)
    ap.add_argument("--id_sep", default="|")
```

Load `cluster_tsv` (rep -> members, restricted per strain via
`gene_positions`'s own `Short` column) and `gene_positions` into the shapes
`add_island_locus()` expects, then call it right after
`annotate_islands_with_domains()` and before writing `islands_with_domains.tsv`:

```python
    from lib.pangenome_matrix import read_cluster_tsv  # or existing import style in this file
    member_to_rep = read_cluster_tsv(args.cluster_tsv)

    gene_positions: dict[tuple[str, str], dict] = {}
    with open(args.gene_positions, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene_positions[(row["Short"], row["protein_id"])] = {
                "contig": row["contig"], "start": int(row["start"]), "end": int(row["end"]),
            }

    annotated = annotate_islands_with_domains(islands_rows, family_domains)
    annotated = add_island_locus(annotated, member_to_rep, gene_positions, id_sep=args.id_sep)
```

Update the zero-island fallback header (currently a hardcoded list) to include
the 6 new columns:

```python
        fieldnames = list(annotated[0].keys()) if annotated else [
            "n_strains", "example_strain", "island_size", "member_families",
            "n_supporting_pairs", "classifications", "pfam_domains",
            "locus_id", "locus_contig", "locus_start", "locus_end",
            "n_members_with_coordinates", "n_contigs_in_locus",
        ]
```

(Check the actual import path/style already used in this file for
`lib.pangenome_matrix` — match the existing `sys.path.insert` + bare-import
convention rather than a package-style `from lib.pangenome_matrix import`.)

- [ ] **Step 5: Implement the top-islands table in `render_report_markdown()`**

Add a new optional parameter `islands_with_domains_rows: list[dict] | None = None`
to `render_report_markdown()`'s signature, and in the "## Accessory islands"
section (after the existing count/figure lines), add:

```python
    top_islands = [r for r in (islands_with_domains_rows or []) if r.get("locus_id", "-") != "-"]
    if top_islands:
        top_islands.sort(key=lambda r: -int(r["island_size"]))
        lines += ["", "**Top islands (by size):**", "",
                  "| Locus | Size | Strains | Pfam domains |", "|---|---|---|---|"]
        for row in top_islands[:20]:
            lines.append(f"| {row['locus_id']} | {row['island_size']} | "
                          f"{row['n_strains']} | {row.get('pfam_domains', '-')} |")
        lines.append("")
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `pixi run pytest tests/test_pangenome_report_tables.py tests/test_pangenome_report_render.py -v`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add bin/pangenome_report_tables.py bin/pangenome_report_render.py \
        tests/test_pangenome_report_tables.py tests/test_pangenome_report_render.py
git commit -m "pangenome: add island genomic-locus assignment to REPORT_TABLES"
```

---

### Task 4: Per-strain outlier flag (median/MAD modified z-score)

**Files:**
- Modify: `bin/pangenome_report_tables.py`
- Modify: `bin/pangenome_report_render.py`
- Test: `tests/test_pangenome_report_tables.py`
- Test: `tests/test_pangenome_report_render.py`

**Interfaces:**
- Produces: `per_strain_summary()` gains two new columns per row:
  `singleton_z`, `is_outlier`.

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_pangenome_report_tables.py`:

```python
def test_per_strain_summary_flags_clear_outlier():
    # 5 strains, singleton counts 10,11,9,10,150 -- strain E is the outlier.
    totals = {
        "A": {"Short": "A", "n_families": 100, "core": 80, "soft_core": 5, "shell": 3, "cloud": 2, "singleton": 10},
        "B": {"Short": "B", "n_families": 101, "core": 80, "soft_core": 5, "shell": 3, "cloud": 2, "singleton": 11},
        "C": {"Short": "C", "n_families": 99, "core": 80, "soft_core": 5, "shell": 3, "cloud": 2, "singleton": 9},
        "D": {"Short": "D", "n_families": 100, "core": 80, "soft_core": 5, "shell": 3, "cloud": 2, "singleton": 10},
        "E": {"Short": "E", "n_families": 240, "core": 80, "soft_core": 5, "shell": 3, "cloud": 2, "singleton": 150},
    }
    rows = add_outlier_flags(list(totals.values()))
    by_short = {r["Short"]: r for r in rows}
    assert by_short["E"]["is_outlier"] == "Y"
    assert by_short["A"]["is_outlier"] == "N"


def test_per_strain_summary_no_outliers_all_n():
    totals = [
        {"Short": s, "singleton": v} for s, v in
        [("A", 10), ("B", 11), ("C", 9), ("D", 10), ("E", 12)]
    ]
    rows = add_outlier_flags(totals)
    assert all(r["is_outlier"] == "N" for r in rows)


def test_per_strain_summary_below_min_n_emits_sentinel():
    totals = [{"Short": "A", "singleton": 10}, {"Short": "B", "singleton": 500}]
    rows = add_outlier_flags(totals)
    assert all(r["singleton_z"] == "-" and r["is_outlier"] == "-" for r in rows)


def test_per_strain_summary_zero_mad_emits_sentinel():
    totals = [{"Short": s, "singleton": 10} for s in ("A", "B", "C", "D")]
    rows = add_outlier_flags(totals)
    assert all(r["singleton_z"] == "-" and r["is_outlier"] == "-" for r in rows)
```

Add to `tests/test_pangenome_report_render.py`:

```python
def test_flagged_outlier_strains_line_present():
    md = render_report_markdown(
        {}, {}, {}, [], 0, None, None, [100, 101, 99, 100, 240],
        per_strain_rows=[
            {"Short": "A", "is_outlier": "N"}, {"Short": "E", "is_outlier": "Y"},
        ],
    )
    assert "Outlier strains" in md
    assert "E" in md


def test_flagged_outlier_strains_line_absent_when_none_flagged():
    md = render_report_markdown(
        {}, {}, {}, [], 0, None, None, [100, 101],
        per_strain_rows=[{"Short": "A", "is_outlier": "N"}, {"Short": "B", "is_outlier": "N"}],
    )
    assert "Outlier strains" not in md
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest tests/test_pangenome_report_tables.py tests/test_pangenome_report_render.py -v`
Expected: FAIL

- [ ] **Step 3: Implement `add_outlier_flags()` in `bin/pangenome_report_tables.py`**

```python
def add_outlier_flags(totals: list[dict], mad_multiplier: float = 0.6745, threshold: float = 3.5) -> list[dict]:
    """Modified z-score (Iglewicz-Hoaglin) on each strain's `singleton` count
    across the whole cohort -- NOT a mean/stdev z-score, which is bounded
    (max |z| = sqrt(n-1)) and cannot exceed ~3.0 at n=10, making a fixed >3
    cutoff unreachable for small cohorts. singleton_z/is_outlier are `-`
    (not a fabricated number, not a crash) when n<3 or MAD==0 (most strains
    share the same singleton count -- statistic is undefined)."""
    values = [t["singleton"] for t in totals]
    out = [dict(t) for t in totals]
    if len(values) < 3:
        for row in out:
            row["singleton_z"] = "-"
            row["is_outlier"] = "-"
        return out

    sorted_values = sorted(values)
    n = len(sorted_values)
    median = sorted_values[n // 2] if n % 2 else (sorted_values[n // 2 - 1] + sorted_values[n // 2]) / 2
    deviations = sorted([abs(v - median) for v in values])
    mad = deviations[n // 2] if n % 2 else (deviations[n // 2 - 1] + deviations[n // 2]) / 2

    if mad == 0:
        for row in out:
            row["singleton_z"] = "-"
            row["is_outlier"] = "-"
        return out

    for row in out:
        z = mad_multiplier * (row["singleton"] - median) / mad
        row["singleton_z"] = f"{z:.2f}"
        row["is_outlier"] = "Y" if abs(z) > threshold else "N"
    return out
```

Wire it into `per_strain_summary()`'s existing return (call
`add_outlier_flags()` on the list of per-strain dicts before returning), and
add `singleton_z`/`is_outlier` to the `fieldnames` list in `main()`'s
`per_strain_summary.tsv` writer.

- [ ] **Step 4: Implement the flagged-strains line in `render_report_markdown()`**

Add a new optional parameter `per_strain_rows: list[dict] | None = None`, and in
the "## Per-strain summary" section, after the existing
min/median/max/`per_strain_summary.tsv` line:

```python
    flagged = [r["Short"] for r in (per_strain_rows or []) if r.get("is_outlier") == "Y"]
    if flagged:
        lines += [f"**Outlier strains (singleton-count modified z-score beyond threshold):** "
                  f"{', '.join(flagged)}", ""]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pixi run pytest tests/test_pangenome_report_tables.py tests/test_pangenome_report_render.py -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add bin/pangenome_report_tables.py bin/pangenome_report_render.py \
        tests/test_pangenome_report_tables.py tests/test_pangenome_report_render.py
git commit -m "pangenome: add per-strain outlier flag (median/MAD modified z-score)"
```

---

### Task 5: Workflow wiring, Pfam2GO module, docs, and real-execution smoke test

**Files:**
- Create: `modules/pangenome/pfam2go.nf`
- Modify: `modules/pangenome/report.nf`
- Modify: `workflows/pangenome_profile.nf`
- Modify: `nextflow.config`
- Modify: `pangenome.nf` (help text)
- Modify: `README.md`

**Interfaces:**
- Consumes: Task 1-4's script changes; existing `CLUSTER_TIER1.out.cluster_tsv`,
  `GENE_POSITIONS.out.positions` emits.
- Produces: new `pangenome_pfam2go = null` param; a new `enrichment_for_report`
  channel variable feeding both `REPORT_TABLES` and `REPORT_RENDER`.

- [ ] **Step 1: Create `modules/pangenome/pfam2go.nf`**

```groovy
nextflow.enable.dsl=2

// PFAM2GO -- optional GO-term annotation of island_pfam_enrichment.tsv via a
// standard pfam2go mapping file. Only invoked when params.pangenome_pfam2go
// is set (an operator-supplied local file -- see pangenome.nf's help text for
// where to obtain it; this pipeline does not auto-fetch it).
process PFAM2GO {
    label 'low_cpu'
    tag "pfam2go"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(island_pfam_enrichment)
    path(pfam2go)

    output:
    path("island_pfam_enrichment.go.tsv"), emit: annotated

    script:
    """
    pangenome_pfam2go.py \
        --island_pfam_enrichment ${island_pfam_enrichment} \
        --pfam2go ${pfam2go} \
        --output island_pfam_enrichment.go.tsv
    """
}
```

- [ ] **Step 2: Add `--cluster_tsv`/`--gene_positions`/`--id_sep` to
`REPORT_TABLES` in `modules/pangenome/report.nf`**

Add `path(cluster_tsv)` and `path(gene_positions)` to `REPORT_TABLES`'s
`input:` block (after `domtblout`), and add to its script block:

```groovy
        --cluster_tsv ${cluster_tsv} \
        --gene_positions ${gene_positions} \
        --id_sep '${params.pangenome_id_sep}' \
```

- [ ] **Step 3: Wire `enrichment_for_report`, the new `REPORT_TABLES` inputs,
and the `PFAM2GO` include into `workflows/pangenome_profile.nf`**

Add to the `include` block near the other pangenome module includes:

```groovy
include { PFAM2GO } from '../modules/pangenome/pfam2go'
```

In the `if (params.pangenome_island_pfam_hmm) { ... }` block, after
`DOMAIN_ENRICHMENT(...)` and before `REPORT_TABLES(...)`, add:

```groovy
        if (params.pangenome_pfam2go) {
            PFAM2GO(DOMAIN_ENRICHMENT.out.enrichment, file(params.pangenome_pfam2go))
            enrichment_for_report = PFAM2GO.out.annotated
        } else {
            enrichment_for_report = DOMAIN_ENRICHMENT.out.enrichment
        }
```

Update the `REPORT_TABLES(...)` call to use `enrichment_for_report` (in place
of the direct `DOMAIN_ENRICHMENT.out.enrichment` reference) and add the two new
trailing args:

```groovy
        REPORT_TABLES(
            BUILD_ISLANDS.out.islands,
            enrichment_for_report,
            PAIR_CLASSIFICATION.out.classification,
            rescued_matrix,
            FREQUENCY_BINS.out.table,
            FAMILY_PFAM_SCAN.out.domtblout,
            CLUSTER_TIER1.out.cluster_tsv,
            GENE_POSITIONS.out.positions,
        )
```

Update the `REPORT_RENDER(...)` call's `island_pfam_enrichment` argument
(currently `DOMAIN_ENRICHMENT.out.enrichment`) to `enrichment_for_report`:

```groovy
        REPORT_RENDER(
            FREQUENCY_BINS.out.table,
            rescued_matrix,
            REPORT_TABLES.out.islands_with_domains,
            REPORT_TABLES.out.size_distribution,
            REPORT_TABLES.out.classification_counts,
            enrichment_for_report,
            REPORT_TABLES.out.marker_summary,
            REPORT_TABLES.out.per_strain_summary,
        )
```

- [ ] **Step 4: Add the new param to `nextflow.config`**

In the "Accessory-island construction + Pfam functional enrichment" params
block (near `pangenome_island_pfam_hmm`), add:

```groovy
    pangenome_pfam2go = null               // path to a local pfam2go mapping file (see
                                            // pangenome.nf --help for where to obtain it);
                                            // GO-term annotation skipped entirely when unset
```

- [ ] **Step 5: Update `pangenome.nf`'s help text and `README.md`**

In `pangenome.nf`'s help text, near the existing
`--pangenome_island_pfam_hmm`/`--pangenome_marker_*` entries, add:

```
      --pangenome_pfam2go              Path to a local pfam2go mapping file -- enables GO-term
                                       annotation of enriched Pfam domains (off by default).
                                       Download from:
                                       http://current.geneontology.org/ontology/external2go/pfam2go
                                       This pipeline does NOT auto-fetch it; obtain it once and
                                       pass its local path.
```

Add an equivalent short note to `README.md` wherever the existing pangenome
params are documented (match whatever section format Task 9 of the prior plan
used for `--pangenome_island_pfam_hmm`, if `README.md` documents it there —
check first; if `README.md` has no pangenome section at all, as the prior
plan's Task 9 found, skip this file and note so in the report).

- [ ] **Step 6: Build the 10-strain smoke-test fixture**

No committed fixture exists (confirmed absent from the repo by grep). Build one
from the already-onboarded Coccidioides study's own config, matching the prior
plan's Task 8 approach:

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
SMOKE_DIR="${SCRATCH:?}/report_enrichment_smoketest"
mkdir -p "$SMOKE_DIR"
head -1 /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/config.csv > "$SMOKE_DIR/config_smoketest.csv"
tail -n +2 /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/config.csv | head -10 >> "$SMOKE_DIR/config_smoketest.csv"
```

For the Pfam2GO path, do NOT download the real (multi-MB) pfam2go file for
this smoke test. Instead, run once WITHOUT `--pangenome_pfam2go` to see which
real Pfam accessions appear in this 10-strain run's `island_pfam_enrichment.tsv`
(if any — a 10-strain run may find zero significant islands, matching the
prior plan's own real result), and build a tiny synthetic pfam2go-format fixture
covering at least one of those real accessions (or, if zero domains are
enriched at this scale, covering the accession from Task 2's own test fixture)
so the with-pfam2go run has something deterministic to annotate.

- [ ] **Step 7: Run both smoke-test paths**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
rm -rf work .nextflow.log* .nextflow

# Without --pangenome_pfam2go
nextflow run pangenome.nf \
    --pangenome_samplesheet "$SMOKE_DIR/config_smoketest.csv" \
    --pangenome_data_dir /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/data_dir \
    --pangenome_island_pfam_hmm /bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm \
    --pangenome_project report_enrichment_smoketest_nopfam2go \
    --outdir "$SMOKE_DIR/out_nopfam2go" \
    -profile local -with-trace "$SMOKE_DIR/out_nopfam2go/trace.txt"

# With --pangenome_pfam2go
nextflow run pangenome.nf \
    --pangenome_samplesheet "$SMOKE_DIR/config_smoketest.csv" \
    --pangenome_data_dir /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/data_dir \
    --pangenome_island_pfam_hmm /bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm \
    --pangenome_pfam2go "$SMOKE_DIR/pfam2go_test_fixture" \
    --pangenome_project report_enrichment_smoketest_withpfam2go \
    --outdir "$SMOKE_DIR/out_withpfam2go" \
    -profile local -with-trace "$SMOKE_DIR/out_withpfam2go/trace.txt"
```

Verify for BOTH runs: `completed == failed_count + completed_count` with
`failed == 0` (check `trace.txt` or Nextflow's own summary), `report.md` exists,
`islands_with_domains.tsv` has the 6 new locus columns in its header, and
`per_strain_summary.tsv` has `singleton_z`/`is_outlier` columns. Verify
SPECIFICALLY for the with-pfam2go run: `report.md` contains a "GO terms" column
header (proving Component 2's channel wiring actually reaches the render step —
this is exactly the class of bug B1 in the design review caught, so do not skip
this check). Verify for the no-pfam2go run: `report.md` does NOT contain a "GO
terms" column header.

If either run has zero significant islands (a real, expected possibility at
this scale, per the prior plan's own documented result), confirm the
`islands_with_domains.tsv` zero-row header still contains all 6 locus columns
(this is exactly what Task 3's zero-island fallback-header fix targets) and
that the report gracefully omits the top-islands table and the outlier-strains
line rather than erroring.

- [ ] **Step 8: Run the full test suite**

Run: `pixi run pytest tests/ -v`
Expected: PASS, 0 failed.

- [ ] **Step 9: Commit**

```bash
git add modules/pangenome/pfam2go.nf modules/pangenome/report.nf \
        workflows/pangenome_profile.nf nextflow.config pangenome.nf README.md
git commit -m "pangenome: wire Pfam2GO + island-locus + outlier-flag report enrichment into the workflow"
```

---

### Task 6: Full regression check

**Files:** none (verification only)

- [ ] **Step 1: Run the full test suite one final time**

Run: `pixi run pytest tests/ -v`
Expected: PASS, 0 failed, matching or exceeding Task 5's count.

- [ ] **Step 2: Report final counts**

No commit — this task is a final confirmation gate before the whole-branch
review, per `superpowers:subagent-driven-development`'s process.
