# Pangenome Island + Pfam Enrichment Pipeline Step Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generalize the Afumigatus study's manual accessory-island + Pfam
functional-enrichment scripts into a permanent, flag-gated step in
`nf_NovInvenio`'s `pangenome_profile.nf` subworkflow.

**Architecture:** Seven new modules (`BUILD_ISLANDS`, a new `MARKER_HMMSEARCH`
process for 0+ named marker searches, `SELECT_BACKGROUND_REPS`,
`FAMILY_PFAM_SCAN`, `DOMAIN_ENRICHMENT`, `REPORT_TABLES`, `REPORT_RENDER`), all in
`nf_NovInvenio`'s shared `bin/`/`modules/pangenome/` (pipeline code, not
study-specific), gated behind a new, distinct `--pangenome_island_pfam_hmm`
param (not the pre-existing `--pangenome_pfam_hmm`, which stays load-bearing
for the captain-by-name branch only). All new scripts port proven, real
algorithms from `studies/fungi/Afumigatus_pangenome/bin/` verbatim where the
logic is already correct, with fixes from two rounds of Fable review and two
rounds of Opus review folded in: one Pfam scan (not two, split at enrichment
time), named marker columns (not one collapsing `has_marker_gene` column), a
dedicated `MARKER_HMMSEARCH` process (not an illegal repeated invocation of
the existing captain-gene modules), and all five of the spec's required
report figures restored (an earlier draft silently dropped three of them).

**Tech Stack:** Python 3 (stdlib `csv`/`argparse`, `numpy`, `scipy.stats` for
Fisher exact + BH correction, `matplotlib` `Agg` backend for figures), Nextflow
DSL2, `pytest` (this repo's existing `tests/` convention, `pixi run pytest`).

**Spec:** `notes/superpowers/specs/2026-09-16-pangenome-island-pfam-enrichment-design.md`

## Global Constraints

- All new files live in `nf_NovInvenio` (`/bigdata/stajichlab/jstajich/projects/NovInvenio`,
  branch `pangenome-profiling-module`), not the NII study repo — this is pipeline
  code, reusable by any future pangenome study.
- Every new script gets a `tests/test_<script>.py` following this repo's existing
  pattern (`sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))`, plain
  `pytest`, no fixtures framework beyond `tmp_path`).
- Run `pixi run pytest tests/ -v` after each task (ambient `python3` is 3.9 and
  can't even collect this repo's existing tests — always use `pixi run`).
- No hand-written narrative in any generated report — `REPORT_RENDER`'s Markdown
  must be fully templated so it works for any future study automatically.
- The 10%-missing-from-`is_core` warning in `accessory_islands()` is a real
  data-integrity safety check (catches GFF3-ID-vs-presence-matrix-ID mismatches)
  and must be ported unchanged, not dropped as "just a warning."
- Background for Pfam enrichment is always shell+cloud families only (never
  core/soft_core/singleton) — matches `pangenome_cooccurrence.py`'s own
  co-occurrence-eligibility selection; singleton island members are excluded from
  enrichment testing, not treated as an error.

---

### Task 1: `bin/pangenome_build_islands.py` — island construction + significance gate

**Files:**
- Create: `bin/pangenome_build_islands.py`
- Test: `tests/test_pangenome_build_islands.py`

**Interfaces:**
- Produces: `accessory_islands(gene_order: list[tuple], is_core: dict[str, bool]) -> list[list[tuple]]`,
  `build_pair_index(significant_pairs: dict[frozenset, str]) -> dict[str, list[tuple[frozenset, str]]]`,
  `find_significant_islands(strain_gene_orders: dict[str, list[tuple]], is_core: dict[str, bool], significant_pairs: dict[frozenset, str]) -> list[dict]`,
  `load_hit_families(tblout_path: str, member_to_rep: dict[str, str], id_sep: str = "|") -> set[str]`
  — Task 3's Nextflow module wiring and any later task calling this script's CLI
  rely on these exact names/signatures and on the CLI flags defined in Step 7
  below (`--marker_tblout NAME=PATH`, repeatable).

- [ ] **Step 1: Write the failing tests for `accessory_islands()`**

```python
# tests/test_pangenome_build_islands.py
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_build_islands import accessory_islands, build_pair_index, find_significant_islands


def test_accessory_islands_merges_consecutive_noncore_same_contig():
    gene_order = [
        ("g1", "chr1", 0, 0), ("g2", "chr1", 1, 1), ("g3", "chr1", 2, 2),
        ("g4", "chr1", 3, 3), ("g5", "chr1", 4, 4),
    ]
    is_core = {"g1": True, "g2": False, "g3": False, "g4": False, "g5": True}
    islands = accessory_islands(gene_order, is_core)
    assert len(islands) == 1
    assert [g[0] for g in islands[0]] == ["g2", "g3", "g4"]


def test_accessory_islands_splits_at_contig_boundary():
    gene_order = [
        ("g1", "chr1", 0, 0), ("g2", "chr1", 1, 1),
        ("g3", "chr2", 0, 0), ("g4", "chr2", 1, 1),
    ]
    is_core = {"g1": False, "g2": False, "g3": False, "g4": False}
    islands = accessory_islands(gene_order, is_core)
    assert len(islands) == 2
    assert [g[0] for g in islands[0]] == ["g1", "g2"]
    assert [g[0] for g in islands[1]] == ["g3", "g4"]


def test_accessory_islands_defaults_missing_gene_to_core():
    gene_order = [("g1", "chr1", 0, 0), ("g2", "chr1", 1, 1), ("g3", "chr1", 2, 2)]
    is_core = {"g1": False, "g3": False}  # g2 missing -> defaults to core, splits the run
    islands = accessory_islands(gene_order, is_core)
    assert len(islands) == 2
    assert [g[0] for g in islands[0]] == ["g1"]
    assert [g[0] for g in islands[1]] == ["g3"]


def test_accessory_islands_warns_above_10pct_missing():
    gene_order = [(f"g{i}", "chr1", i, i) for i in range(10)]
    is_core = {"g0": False}  # 9/10 missing -> 90% > 10% threshold
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        accessory_islands(gene_order, is_core)
    assert any("absent from is_core" in str(w.message) for w in caught)


def test_build_pair_index_indexes_by_both_families():
    pairs = {frozenset({"famA", "famB"}): "unexplained_physical"}
    index = build_pair_index(pairs)
    assert index["famA"] == [(frozenset({"famA", "famB"}), "unexplained_physical")]
    assert index["famB"] == [(frozenset({"famA", "famB"}), "unexplained_physical")]


def test_find_significant_islands_keeps_only_islands_with_a_supporting_pair():
    strain_gene_orders = {
        "s1": [
            ("famA", "chr1", 0, 0), ("famB", "chr1", 1, 1),  # supported island
            ("core1", "chr1", 2, 2),
            ("famC", "chr1", 3, 3), ("famD", "chr1", 4, 4),  # unsupported island
        ],
    }
    is_core = {"famA": False, "famB": False, "core1": True, "famC": False, "famD": False}
    significant_pairs = {frozenset({"famA", "famB"}): "unexplained_physical"}
    islands = find_significant_islands(strain_gene_orders, is_core, significant_pairs)
    assert len(islands) == 1
    assert set(islands[0]["members"]) == {"famA", "famB"}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
pixi run pytest tests/test_pangenome_build_islands.py -v
```

Expected: FAIL (`ModuleNotFoundError: No module named 'pangenome_build_islands'`).

- [ ] **Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""Accessory-island construction + statistical-significance gate, generalized
from studies/fungi/Afumigatus_pangenome/bin/synteny_windows.py +
find_accessory_islands.py (NovInvenio_Investigations repo) into a permanent
nf_NovInvenio pipeline step -- see
notes/superpowers/specs/2026-09-16-pangenome-island-pfam-enrichment-design.md.

Consumes family_positions.tsv's already-computed per-strain `rank` column
directly (no GFF3 re-parsing needed -- accessory_islands() only needs
(gene_id, contig) pairs in traversal order, and rank already IS that order).
"Significant" island membership is gated entirely by pair_classification.tsv's
`classification` column (starship_explained/unexplained_physical/
ambiguous_linkage) -- that column already reflects the upstream FDR-filtered,
exact-permutation-tested co-occurrence result; this script does not re-derive
or recompute any significance test.

Usage:
  pangenome_build_islands.py --family_positions family_positions.tsv \\
      --frequency_table frequency_table.tsv \\
      --pair_classification pair_classification.tsv \\
      --cluster_tsv tier1_cluster.tsv \\
      [--marker_tblout captain=captain_vs_study.tblout] \\
      [--marker_tblout sm_backbone=SM_backbone_vs_study.tblout] \\
      --output significant_islands.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import read_cluster_tsv  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402


PHYSICAL_CLASSIFICATIONS = frozenset(
    {"starship_explained", "unexplained_physical", "ambiguous_linkage"}
)


def accessory_islands(
    gene_order: list[tuple], is_core: dict[str, bool]
) -> list[list[tuple]]:
    """Maximal runs of consecutive non-core genes, never crossing a contig
    boundary. A gene_id absent from `is_core` defaults to core (conservative:
    an unrecognized gene more likely reflects a GFF3-ID-vs-presence-matrix-ID
    mismatch than a genuine accessory gene) -- warns if >10% of genes seen are
    missing from `is_core`, since that's the real failure signature of such a
    mismatch (every island silently collapsing to nothing, no exception
    raised)."""
    islands: list[list[tuple]] = []
    current: list[tuple] = []
    prev_contig = None
    missing = 0
    for gene in gene_order:
        gene_id, contig = gene[0], gene[1]
        if gene_id not in is_core:
            missing += 1
        core = is_core.get(gene_id, True)
        if contig != prev_contig and current:
            islands.append(current)
            current = []
        if not core:
            current.append(gene)
        elif current:
            islands.append(current)
            current = []
        prev_contig = contig
    if current:
        islands.append(current)
    if gene_order and missing / len(gene_order) > 0.1:
        warnings.warn(
            f"accessory_islands: {missing}/{len(gene_order)} genes "
            f"({missing / len(gene_order):.0%}) were absent from is_core and "
            "defaulted to core -- check for a GFF3 ID vs presence-matrix "
            "protein ID format mismatch before trusting these islands.",
            stacklevel=2,
        )
    return islands


def load_is_core(frequency_table_path: str) -> dict[str, bool]:
    """{family: True if core/soft_core, False if shell/cloud/singleton}."""
    is_core: dict[str, bool] = {}
    with open(frequency_table_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            is_core[row["family"]] = row["bin"] in ("core", "soft_core")
    return is_core


def load_strain_gene_orders(family_positions_path: str) -> dict[str, list[tuple]]:
    """{Short: [(family, contig, rank, rank), ...]}, sorted by rank."""
    by_strain: dict[str, list[tuple]] = {}
    with open(family_positions_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            short, family, contig, rank = row["Short"], row["family"], row["contig"], int(row["rank"])
            by_strain.setdefault(short, []).append((family, contig, rank, rank))
    for rows in by_strain.values():
        rows.sort(key=lambda r: r[2])
    return by_strain


def load_significant_physical_pairs(pair_classification_path: str) -> dict[frozenset, str]:
    """{frozenset({family_a, family_b}): classification} for every physically-linked pair."""
    pairs: dict[frozenset, str] = {}
    with open(pair_classification_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["classification"] in PHYSICAL_CLASSIFICATIONS:
                pairs[frozenset({row["family_a"], row["family_b"]})] = row["classification"]
    return pairs


def load_hit_families(tblout_path: str, member_to_rep: dict[str, str], id_sep: str = "|") -> set[str]:
    """Families with >=1 hmmsearch/hmmscan hit anywhere in the cohort --
    generic across any tblout (captain, SM-backbone, or any future named
    marker search)."""
    families: set[str] = set()
    with open_maybe_compressed(tblout_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            target = line.split()[0]
            if id_sep not in target:
                continue
            short, protein_id = target.split(id_sep, 1)
            family = member_to_rep.get(f"{short}{id_sep}{protein_id}")
            if family is not None:
                families.add(family)
    return families


def build_pair_index(significant_pairs: dict[frozenset, str]) -> dict[str, list[tuple[frozenset, str]]]:
    """{family: [(pair, classification), ...]} -- avoids an O(islands x
    total_pairs) full scan; indexed lookup is O(islands x island_size x
    avg_pairs_per_family), the difference between finishing and not (a naive
    full-scan version was killed after 12+ minutes with no output on the
    Afumigatus real-data run)."""
    index: dict[str, list[tuple[frozenset, str]]] = {}
    for pair, classification in significant_pairs.items():
        for family in pair:
            index.setdefault(family, []).append((pair, classification))
    return index


def find_significant_islands(
    strain_gene_orders: dict[str, list[tuple]],
    is_core: dict[str, bool],
    significant_pairs: dict[frozenset, str],
) -> list[dict]:
    """Per strain: build accessory islands, keep only islands containing >=1
    significant physically-linked pair (both members present in that SAME
    island in that strain)."""
    pairs_by_family = build_pair_index(significant_pairs)
    results: list[dict] = []
    for short, gene_order in strain_gene_orders.items():
        islands = accessory_islands(gene_order, is_core)
        for island in islands:
            if len(island) < 2:
                continue
            members = [g[0] for g in island]
            member_set = set(members)
            candidates: dict[frozenset, str] = {}
            for family in member_set:
                for pair, classification in pairs_by_family.get(family, ()):
                    if pair <= member_set:
                        candidates[pair] = classification
            if not candidates:
                continue
            results.append({
                "strain": short,
                "contig": island[0][1],
                "size": len(island),
                "members": members,
                "supporting_pairs": list(candidates.items()),
            })
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family_positions", required=True)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--pair_classification", required=True)
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--min_island_size", type=int, default=2)
    ap.add_argument(
        "--marker_tblout", action="append", default=[], metavar="NAME=PATH",
        help="Repeatable. e.g. --marker_tblout captain=captain.tblout "
             "--marker_tblout sm_backbone=sm_backbone.tblout",
    )
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    markers: dict[str, str] = {}
    for entry in args.marker_tblout:
        if "=" not in entry:
            sys.exit(f"ERROR: --marker_tblout must be NAME=PATH, got: {entry!r}")
        name, path = entry.split("=", 1)
        markers[name] = path

    is_core = load_is_core(args.frequency_table)
    strain_gene_orders = load_strain_gene_orders(args.family_positions)
    significant_pairs = load_significant_physical_pairs(args.pair_classification)
    print(f"pangenome_build_islands: {len(significant_pairs)} significant physically-linked "
          f"pairs, {len(strain_gene_orders)} strains", file=sys.stderr)

    islands = find_significant_islands(strain_gene_orders, is_core, significant_pairs)
    islands = [isl for isl in islands if isl["size"] >= args.min_island_size]
    print(f"pangenome_build_islands: {len(islands)} significant islands found across all strains",
          file=sys.stderr)

    marker_families: dict[str, set[str]] = {name: set() for name in markers}
    if markers:
        member_to_rep = read_cluster_tsv(args.cluster_tsv)
        for name, path in markers.items():
            marker_families[name] = load_hit_families(path, member_to_rep)

    marker_names = sorted(markers)
    with open(args.output, "w") as out:
        header = (
            "n_strains\texample_strain\tisland_size\tmember_families\t"
            "n_supporting_pairs\tclassifications"
        )
        for name in marker_names:
            header += f"\thas_{name}"
        out.write(header + "\n")

        by_key: dict[tuple, dict] = {}
        for island in islands:
            member_set = frozenset(island["members"])
            classifications = frozenset(c for _, c in island["supporting_pairs"])
            key = (member_set, classifications)
            entry = by_key.setdefault(key, {
                "strains": [], "size": island["size"], "members": island["members"],
                "classifications": classifications,
                "n_supporting_pairs": len(island["supporting_pairs"]),
            })
            entry["strains"].append(island["strain"])

        for (member_set, classifications), entry in sorted(
            by_key.items(), key=lambda kv: -len(kv[1]["strains"])
        ):
            row = (
                f"{len(entry['strains'])}\t{entry['strains'][0]}\t{entry['size']}\t"
                f"{','.join(entry['members'])}\t{entry['n_supporting_pairs']}\t"
                f"{','.join(sorted(classifications))}"
            )
            for name in marker_names:
                has_marker = bool(member_set & marker_families[name])
                row += f"\t{'Y' if has_marker else 'N'}"
            out.write(row + "\n")

    print(f"pangenome_build_islands: {len(by_key)} distinct significant islands "
          f"(deduped across strains) written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pixi run pytest tests/test_pangenome_build_islands.py -v
```

Expected: PASS (6 tests).

- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_build_islands.py tests/test_pangenome_build_islands.py
git commit -m "$(cat <<'EOF'
pangenome: add island construction + significance gate

Generalizes Afumigatus study's synteny_windows.py/find_accessory_islands.py
into a permanent pipeline script. Named marker tblouts (--marker_tblout
NAME=PATH, repeatable) replace the study's hardcoded captain/sm_backbone
flags, producing one has_<name> column per marker instead of a single
column that would collapse independent signals.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 2: `bin/pangenome_select_background_reps.py` — Pfam-eligible family rep selection

**Files:**
- Create: `bin/pangenome_select_background_reps.py`
- Test: `tests/test_pangenome_select_background_reps.py`

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces: `select_background_families(frequency_table_path: str) -> set[str]`
  (family IDs in `shell`/`cloud` bins only) and `write_background_fasta(rep_fasta_path: str, background_families: set[str], out_path: str) -> int`
  (returns count written) — Task 4's Nextflow wiring calls this script's CLI.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_pangenome_select_background_reps.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_select_background_reps import select_background_families, write_background_fasta


def test_select_background_families_keeps_only_shell_and_cloud(tmp_path):
    freq = tmp_path / "frequency_table.tsv"
    freq.write_text(
        "family\tfrequency\tstrain_count\tbin\n"
        "famA\t0.9\t90\tcore\n"
        "famB\t0.3\t30\tshell\n"
        "famC\t0.05\t5\tcloud\n"
        "famD\t0.01\t1\tsingleton\n"
    )
    result = select_background_families(str(freq))
    assert result == {"famB", "famC"}


def test_write_background_fasta_filters_by_family_set(tmp_path):
    rep_fasta = tmp_path / "tier1_rep_seq.fasta"
    rep_fasta.write_text(
        ">famA some description\nMSEQA\n"
        ">famB\nMSEQB\n"
        ">famC\nMSEQC\n"
    )
    out_path = tmp_path / "background_reps.fa"
    n = write_background_fasta(str(rep_fasta), {"famB", "famC"}, str(out_path))
    assert n == 2
    content = out_path.read_text()
    assert ">famB" in content and ">famC" in content and ">famA" not in content
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pixi run pytest tests/test_pangenome_select_background_reps.py -v
```

Expected: FAIL (`ModuleNotFoundError`).

- [ ] **Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""Filter tier1_rep_seq.fasta down to the families actually eligible for
Pfam-domain enrichment testing: shell+cloud bins only, matching
pangenome_cooccurrence.py's own co-occurrence-eligibility selection (never
core/soft_core/singleton -- core was never eligible for co-occurrence
testing in the first place, and a whole-genome background would spuriously
enrich for "accessory-genome-typical" domains regardless of which specific
island is tested).

Usage:
  pangenome_select_background_reps.py --rep_fasta tier1_rep_seq.fasta \\
      --frequency_table frequency_table.tsv --output background_reps.fa
"""
from __future__ import annotations

import argparse
import csv
import sys


def select_background_families(frequency_table_path: str) -> set[str]:
    """All families in the shell or cloud bin."""
    background = set()
    with open(frequency_table_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["bin"] in ("shell", "cloud"):
                background.add(row["family"])
    return background


def write_background_fasta(rep_fasta_path: str, background_families: set[str], out_path: str) -> int:
    """Writes only the FASTA records whose header ID is in
    `background_families`. Returns the number of records written."""
    n_written = 0
    writing = False
    with open(rep_fasta_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                family_id = line[1:].split(None, 1)[0].rstrip()
                writing = family_id in background_families
                if writing:
                    n_written += 1
            if writing:
                fout.write(line)
    return n_written


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rep_fasta", required=True)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    background = select_background_families(args.frequency_table)
    n = write_background_fasta(args.rep_fasta, background, args.output)
    print(f"pangenome_select_background_reps: {n} shell+cloud family reps written to {args.output}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pixi run pytest tests/test_pangenome_select_background_reps.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_select_background_reps.py tests/test_pangenome_select_background_reps.py
git commit -m "$(cat <<'EOF'
pangenome: add Pfam-eligible background rep selection

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 3: `bin/pangenome_domain_enrichment.py` — Fisher exact + BH domain enrichment

**Files:**
- Create: `bin/pangenome_domain_enrichment.py`
- Test: `tests/test_pangenome_domain_enrichment.py`

**Interfaces:**
- Consumes: Task 1's `significant_islands.tsv` format (columns
  `member_families` comma-list), Task 2's background-family concept.
- Produces: `parse_domtblout(paths: list[str], max_ievalue: float = 1e-3) -> dict[str, set[str]]`,
  `domain_enrichment(island_member_families: set[str], background_families: set[str], family_domains: dict[str, set[str]]) -> list[dict]`
  (each dict has keys `domain, n_with_domain_in_islands, n_with_domain_in_background, n_island_families, n_background_families, fisher_p, fdr_q`)
  — Task 5's `island_pfam_enrichment.tsv` schema and Task 7's report-tables join
  rely on these exact column names.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_pangenome_domain_enrichment.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_domain_enrichment import parse_domtblout, domain_enrichment


def test_parse_domtblout_filters_by_ievalue(tmp_path):
    dtbl = tmp_path / "test.domtblout"
    # hmmscan domtblout columns (0-indexed): 0 target_name, 1 target_acc,
    # 2 tlen, 3 query_name, 4 query_acc, 5 qlen, 6 full_evalue, 7 full_score,
    # 8 full_bias, 9 dom#, 10 dom_of, 11 c-Evalue, 12 i-Evalue, 13 dom_score,
    # 14 dom_bias, 15-20 hmm/ali/env coords, 21 acc. parse_domtblout reads
    # parts[12] -- the i-Evalue must actually sit there, or this test can't
    # distinguish pass/fail (caught during Opus review: an earlier draft of
    # this fixture put the discriminating value at index 9/10 instead).
    dtbl.write_text(
        "# comment line\n"
        "PF00001 - 100 famA - 50 1.0e-10 100.0 20.0 1 1 1.0e-10 1.0e-10 100.0 20.0 10 20 10 20 10 20 0.99\n"
        "PF00002 - 100 famA - 50 5.0e-02 100.0 20.0 1 1 5.0e-02 5.0e-02 100.0 20.0 10 20 10 20 10 20 0.99\n"
    )
    hits = parse_domtblout([str(dtbl)], max_ievalue=1e-3)
    assert hits == {"famA": {"PF00001"}}


def test_domain_enrichment_fisher_and_bh(tmp_path):
    # famA (island) has domainX; 3 background families total (famA, famB, famC),
    # only famA has domainX -> strong enrichment signal.
    island_members = {"famA"}
    background = {"famA", "famB", "famC"}
    family_domains = {"famA": {"domainX"}}
    result = domain_enrichment(island_members, background, family_domains)
    assert len(result) == 1
    row = result[0]
    assert row["domain"] == "domainX"
    assert row["n_with_domain_in_islands"] == 1
    assert row["n_with_domain_in_background"] == 1
    assert "fdr_q" in row


def test_domain_enrichment_excludes_island_members_outside_background(capsys):
    # famZ is an island member but a singleton (not in shell+cloud background)
    island_members = {"famA", "famZ"}
    background = {"famA", "famB"}
    family_domains = {"famA": {"domainX"}}
    result = domain_enrichment(island_members, background, family_domains)
    captured = capsys.readouterr()
    assert "outside the eligible" in captured.err
    # famZ excluded -> n_island_families should reflect only famA
    assert all(row["n_island_families"] == 1 for row in result)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pixi run pytest tests/test_pangenome_domain_enrichment.py -v
```

Expected: FAIL (`ModuleNotFoundError`).

- [ ] **Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""Pfam-domain enrichment test for accessory-island member families vs. the
CORRECT background: all families actually eligible for co-occurrence testing
(shell+cloud bins, matching pangenome_cooccurrence.py's own selection) --
never the whole genome, since core genes were never eligible for testing in
the first place and a whole-genome background would spuriously enrich for
"accessory-genome-typical" domains regardless of which specific island is
being tested.

Generalized from studies/fungi/Afumigatus_pangenome/bin/
summarize_island_functions.py's domain_enrichment() -- ported verbatim
(same Fisher-exact + BH-FDR statistics), split out of that study's combined
island-annotation-plus-enrichment script into a standalone enrichment-only
step (island annotation now lives in pangenome_report_tables.py, Task 7).

Usage:
  pangenome_domain_enrichment.py --significant_islands significant_islands.tsv \\
      --domtblout background_reps_vs_pfam.domtblout \\
      --frequency_table frequency_table.tsv \\
      --output island_pfam_enrichment.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys

import numpy as np
from scipy.stats import fisher_exact, false_discovery_control


def parse_domtblout(paths: list[str], max_ievalue: float = 1e-3) -> dict[str, set[str]]:
    """{family_id: {pfam_domain_name, ...}} from one or more hmmscan
    --domtblout files. Domain-level i-Evalue re-checked explicitly (can be
    looser than the sequence-level cutoff hmmscan's own -E already applied)."""
    hits: dict[str, set[str]] = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 13:
                    continue
                target_name, query_name = parts[0], parts[3]
                try:
                    i_evalue = float(parts[12])
                except ValueError:
                    continue
                if i_evalue > max_ievalue:
                    continue
                hits.setdefault(query_name, set()).add(target_name)
    return hits


def domain_enrichment(
    island_member_families: set[str],
    background_families: set[str],
    family_domains: dict[str, set[str]],
) -> list[dict]:
    """One-sided Fisher's exact test per Pfam domain, BH-FDR corrected across
    every domain tested. `island_member_families` is intersected with
    `background_families` first: singleton island members were never eligible
    for co-occurrence testing (shell/cloud only), so they're excluded from
    enrichment testing entirely -- a real, expected occurrence, not an
    error."""
    excluded = island_member_families - background_families
    island_member_families = island_member_families & background_families
    if excluded:
        print(
            f"domain_enrichment: {len(excluded)} island-member families are outside "
            "the eligible (shell+cloud) background -- most likely singleton families "
            "an island happened to include; excluded from the enrichment test "
            f"({len(island_member_families)} remain)", file=sys.stderr,
        )

    all_domains: set[str] = set()
    for family in island_member_families:
        all_domains |= family_domains.get(family, set())

    n_background = len(background_families)
    n_island = len(island_member_families)
    rows = []
    pvalues = []
    for domain in sorted(all_domains):
        families_with_domain = {f for f in background_families if domain in family_domains.get(f, ())}
        both = len(families_with_domain & island_member_families)
        island_only = n_island - both
        domain_only = len(families_with_domain) - both
        neither = n_background - n_island - domain_only
        _, p = fisher_exact([[both, island_only], [domain_only, neither]], alternative="greater")
        rows.append({
            "domain": domain,
            "n_with_domain_in_islands": both,
            "n_with_domain_in_background": len(families_with_domain),
            "n_island_families": n_island,
            "n_background_families": n_background,
            "fisher_p": p,
        })
        pvalues.append(p)

    if pvalues:
        qvalues = false_discovery_control(np.asarray(pvalues), method="bh")
        for row, q in zip(rows, qvalues):
            row["fdr_q"] = float(q)
    rows.sort(key=lambda r: r["fisher_p"])
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--significant_islands", required=True)
    ap.add_argument("--domtblout", required=True, action="append")
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    from pangenome_select_background_reps import select_background_families

    family_domains = parse_domtblout(args.domtblout)
    print(f"pangenome_domain_enrichment: {len(family_domains)} families with >=1 Pfam domain hit",
          file=sys.stderr)

    background = select_background_families(args.frequency_table)
    print(f"pangenome_domain_enrichment: {len(background)} eligible (shell+cloud) background families",
          file=sys.stderr)

    with open(args.significant_islands, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    island_member_families: set[str] = set()
    for row in rows:
        island_member_families.update(row["member_families"].split(","))

    enrichment = domain_enrichment(island_member_families, background, family_domains)
    with open(args.output, "w") as out:
        out.write(
            "domain\tn_with_domain_in_islands\tn_with_domain_in_background\t"
            "n_island_families\tn_background_families\tfisher_p\tfdr_q\n"
        )
        for r in enrichment:
            out.write(
                f"{r['domain']}\t{r['n_with_domain_in_islands']}\t{r['n_with_domain_in_background']}\t"
                f"{r['n_island_families']}\t{r['n_background_families']}\t"
                f"{r['fisher_p']:.3e}\t{r['fdr_q']:.3e}\n"
            )

    print(f"pangenome_domain_enrichment: {len(enrichment)} domains tested, "
          f"{sum(1 for r in enrichment if r['fdr_q'] < 0.05)} significant at FDR<0.05",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pixi run pytest tests/test_pangenome_domain_enrichment.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_domain_enrichment.py tests/test_pangenome_domain_enrichment.py
git commit -m "$(cat <<'EOF'
pangenome: add Pfam domain enrichment (Fisher exact + BH)

Ported from studies/fungi/Afumigatus_pangenome/bin/
summarize_island_functions.py's domain_enrichment(), split out of that
study's combined island-annotation-plus-enrichment script into a
standalone step (island annotation moved to pangenome_report_tables.py).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 4: `modules/pangenome/islands.nf` and `modules/pangenome/pfam_enrichment.nf` — Nextflow module wiring

**Files:**
- Create: `modules/pangenome/islands.nf`
- Create: `modules/pangenome/pfam_enrichment.nf`

**Interfaces:**
- Consumes: Task 1's `pangenome_build_islands.py` CLI (`--family_positions
  --frequency_table --pair_classification --cluster_tsv --min_island_size
  --marker_tblout NAME=PATH (repeatable) --output`), Task 2's
  `pangenome_select_background_reps.py` CLI, Task 3's
  `pangenome_domain_enrichment.py` CLI.
- Produces: `BUILD_ISLANDS`, `MARKER_HMMSEARCH`, `SELECT_BACKGROUND_REPS`,
  `FAMILY_PFAM_SCAN`, `DOMAIN_ENRICHMENT` Nextflow processes — Task 8's
  workflow wiring calls these by exact name.

**Correction from an earlier draft (caught by Opus review before implementation
started): the original design's marker-search plan — aliasing the existing
`HMMFETCH_CAPTAIN`/`CAPTAIN_HMMSEARCH` (`modules/pangenome/captain.nf`) and
invoking the alias once per `--pangenome_marker_names` entry inside a Groovy
`.eachWithIndex` loop — is illegal Nextflow DSL2** (a process, or its alias,
can only be invoked once per workflow; DSL2 raises "Process 'X' has been
already used" on a second call). It also can't preserve the marker's name
alongside its output (`CAPTAIN_HMMSEARCH`'s output filename is hardcoded to
`captain_vs_study.tblout`, not parameterized). Fixed by writing one new,
small process (`MARKER_HMMSEARCH`, below) with a `tuple` input/output so
Nextflow's normal channel-based fan-out runs it once per marker via a single
invocation — this is idiomatic Nextflow multiplicity, not a loop — and its
output filename is name-templated, which also incidentally fixes the
`publishDir` collision the earlier draft's caveat flagged. This process does
**not** touch `modules/pangenome/captain.nf` (untouched, unchanged, still used
by the existing captain-gene branch on its own).

- [ ] **Step 1: Write `modules/pangenome/islands.nf`**

```groovy
// BUILD_ISLANDS -- accessory-island construction + statistical significance
// gate (pangenome_build_islands.py). Named marker hit tables (0+, e.g.
// captain/sm_backbone) are passed as two PARALLEL lists -- marker_names
// (val) and marker_tblout_files (path, so Nextflow actually stages them
// into this task's work directory) -- rather than one combined list of
// tuples, to sidestep Nextflow's default flattening behavior on
// channel.collect() over tuples (a real bug an earlier draft of this task
// had: [[n1,p1],[n2,p2]].collect() flattens to [n1,p1,n2,p2] unless you
// explicitly ask for a nested list). See
// notes/superpowers/specs/2026-09-16-pangenome-island-pfam-enrichment-design.md.
process BUILD_ISLANDS {
    label 'low_cpu'
    tag "build_islands"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(family_positions)
    path(frequency_table)
    path(pair_classification)
    path(cluster_tsv)
    val(marker_names)             // parallel list of marker names, may be empty
    path(marker_tblout_files)     // parallel list of tblout files, same length, may be empty

    output:
    path("significant_islands.tsv"), emit: islands

    script:
    def marker_pairs = [marker_names, marker_tblout_files].transpose()
    def marker_args = marker_pairs.collect { name, path -> "--marker_tblout ${name}=${path}" }.join(' ')
    """
    pangenome_build_islands.py \
        --family_positions ${family_positions} \
        --frequency_table ${frequency_table} \
        --pair_classification ${pair_classification} \
        --cluster_tsv ${cluster_tsv} \
        --min_island_size ${params.pangenome_island_min_size} \
        ${marker_args} \
        --output significant_islands.tsv
    """
}

// MARKER_HMMSEARCH -- new, small, dedicated process for named marker
// searches (e.g. captain/DUF3435, sm_backbone/PKS-NRPS, or any future
// named marker per notes/superpowers/specs/
// 2026-09-16-pangenome-island-pfam-enrichment-design.md). Deliberately
// NOT an alias of modules/pangenome/captain.nf's CAPTAIN_HMMSEARCH (that
// process's output filename is hardcoded and it can't be invoked more
// than once per workflow) -- this is genuinely new, tiny code: same
// hmmsearch command, `tuple` input/output so Nextflow's channel-based
// multiplicity runs it once per marker via ONE process invocation (Task 8
// wires this with Channel.fromList, not a loop), and a name-templated
// output filename.
process MARKER_HMMSEARCH {
    label 'med_cpu'
    tag { marker_name }
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome/markers" }, mode: 'copy'

    input:
    tuple val(marker_name), path(marker_hmm)
    path(all_strains_fa)

    output:
    tuple val(marker_name), path("${marker_name}_vs_study.tblout"), emit: result

    script:
    """
    hmmsearch --tblout ${marker_name}_vs_study.tblout \
        -E ${params.pangenome_marker_evalue} \
        --cpu ${task.cpus} ${marker_hmm} ${all_strains_fa} > /dev/null
    """
}
```

- [ ] **Step 2: Write `modules/pangenome/pfam_enrichment.nf`**

```groovy
// SELECT_BACKGROUND_REPS / FAMILY_PFAM_SCAN / DOMAIN_ENRICHMENT -- Pfam
// functional-enrichment testing for accessory-island member families.
// ONE hmmscan over the full shell+cloud-eligible background (not a
// separate island-vs-background scan pair) -- island/background
// partitioning happens inside DOMAIN_ENRICHMENT itself, matching
// summarize_island_functions.py's actual logic (background must be a
// superset of island members for the Fisher test to be valid).
process SELECT_BACKGROUND_REPS {
    label 'low_cpu'
    tag "select_background_reps"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"

    input:
    path(rep_fasta)
    path(frequency_table)

    output:
    path("background_reps.fa"), emit: fasta

    script:
    """
    pangenome_select_background_reps.py \
        --rep_fasta ${rep_fasta} \
        --frequency_table ${frequency_table} \
        --output background_reps.fa
    """
}

process FAMILY_PFAM_SCAN {
    label 'med_cpu'
    tag "family_pfam_scan"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(background_reps_fasta)
    path(pfam_hmm)

    output:
    path("pfam.domtblout"), emit: domtblout

    script:
    """
    hmmscan --domtblout pfam.domtblout \
        -E ${params.pangenome_pfam_domain_evalue} \
        --cpu ${task.cpus} \
        ${pfam_hmm} ${background_reps_fasta} > /dev/null
    """
}

process DOMAIN_ENRICHMENT {
    label 'low_cpu'
    tag "domain_enrichment"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(significant_islands)
    path(domtblout)
    path(frequency_table)

    output:
    path("island_pfam_enrichment.tsv"), emit: enrichment

    script:
    """
    pangenome_domain_enrichment.py \
        --significant_islands ${significant_islands} \
        --domtblout ${domtblout} \
        --frequency_table ${frequency_table} \
        --output island_pfam_enrichment.tsv
    """
}
```

- [ ] **Step 4: Verify Nextflow parses these files with no syntax errors**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
nextflow run pangenome.nf --help > /tmp/nf_syntax_check.out 2>&1; echo "exit: $?"
```

Expected: `exit: 0` (these two new files aren't `include`d into `pangenome_profile.nf`
yet — Task 8 does that — so this only checks that `pangenome.nf` itself still
parses; a real syntax check of the two new files happens once Task 8 includes
them).

- [ ] **Step 5: Commit**

```bash
git add modules/pangenome/islands.nf modules/pangenome/pfam_enrichment.nf
git commit -m "$(cat <<'EOF'
pangenome: add BUILD_ISLANDS/MARKER_HMMSEARCH/SELECT_BACKGROUND_REPS/FAMILY_PFAM_SCAN/DOMAIN_ENRICHMENT modules

MARKER_HMMSEARCH is new, tuple-based code (not an alias of the existing
captain.nf modules) -- a repeated alias invocation is illegal Nextflow
DSL2, caught by Opus review before implementation started. Not yet wired
into pangenome_profile.nf (Task 8).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 5: `bin/pangenome_report_tables.py` — tidy report tables (no matplotlib)

**Files:**
- Create: `bin/pangenome_report_tables.py`
- Test: `tests/test_pangenome_report_tables.py`

**Interfaces:**
- Consumes: Task 1's `significant_islands.tsv` schema, Task 3's
  `island_pfam_enrichment.tsv` schema (columns `domain, n_with_domain_in_islands,
  n_with_domain_in_background, n_island_families, n_background_families,
  fisher_p, fdr_q`), and `presence_matrix.tsv` (the rescued matrix, wide
  format: header `family` + one column per strain; each row a family with
  `present`/`genome_only`/`absent` per strain).
- Produces: `annotate_islands_with_domains(islands_rows: list[dict], family_domains: dict[str, set[str]]) -> list[dict]`
  (adds a `pfam_domains` key, comma-joined sorted domain names or `"-"`),
  `island_size_distribution(islands_rows: list[dict]) -> dict[int, int]`
  (`{size: count}`), `classification_counts(pair_classification_path: str) -> dict[str, int]`,
  `per_strain_summary(presence_matrix_path: str, family_bin: dict[str, str]) -> list[dict]`
  (one dict per strain: `Short, n_genes, core, soft_core, shell, cloud, singleton`)
  — Task 7's `REPORT_TABLES` module and Task 6's `pangenome_report_render.py`
  rely on these names/shapes.

**Addition from Opus review (report-completeness pass): `per_strain_summary`.**
A genuinely comprehensive pangenome report surfaces per-strain outliers (this
is literally how the Afumigatus study found its one real accessory-count
outlier strain, by hand — `REPORT.md`'s QC section). Streams
`presence_matrix.tsv` row-by-row (one row per family) rather than loading a
full family×strain array into memory, since only per-strain totals are
needed here — Task 6's presence/absence raster figure is the one place that
genuinely needs the full 2D array, and it loads its own copy via
`lib/pangenome_matrix.PresenceMatrix`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_pangenome_report_tables.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_report_tables import (
    annotate_islands_with_domains,
    island_size_distribution,
    classification_counts,
    per_strain_summary,
)


def test_annotate_islands_with_domains_unions_member_domains():
    islands_rows = [{"member_families": "famA,famB", "island_size": "2"}]
    family_domains = {"famA": {"PF00001"}, "famB": {"PF00002"}}
    result = annotate_islands_with_domains(islands_rows, family_domains)
    assert result[0]["pfam_domains"] == "PF00001,PF00002"


def test_annotate_islands_with_domains_no_hits_is_dash():
    islands_rows = [{"member_families": "famA", "island_size": "1"}]
    result = annotate_islands_with_domains(islands_rows, {})
    assert result[0]["pfam_domains"] == "-"


def test_island_size_distribution_counts_by_size():
    islands_rows = [{"island_size": "2"}, {"island_size": "2"}, {"island_size": "5"}]
    assert island_size_distribution(islands_rows) == {2: 2, 5: 1}


def test_classification_counts_tallies_column(tmp_path):
    pc = tmp_path / "pair_classification.tsv"
    pc.write_text(
        "family_a\tfamily_b\tclassification\n"
        "f1\tf2\ttrans\n"
        "f3\tf4\ttrans\n"
        "f5\tf6\tambiguous_linkage\n"
    )
    assert classification_counts(str(pc)) == {"trans": 2, "ambiguous_linkage": 1}


def test_per_strain_summary_counts_genes_and_bins(tmp_path):
    pm = tmp_path / "presence_matrix.tsv"
    pm.write_text(
        "family\ts1\ts2\n"
        "famA\tpresent\tpresent\n"    # core
        "famB\tpresent\tabsent\n"     # shell
        "famC\tgenome_only\tabsent\n" # cloud -- genome_only counts as present
    )
    family_bin = {"famA": "core", "famB": "shell", "famC": "cloud"}
    result = per_strain_summary(str(pm), family_bin)
    by_strain = {r["Short"]: r for r in result}
    assert by_strain["s1"] == {"Short": "s1", "n_genes": 3, "core": 1, "soft_core": 0, "shell": 1, "cloud": 1, "singleton": 0}
    assert by_strain["s2"] == {"Short": "s2", "n_genes": 1, "core": 1, "soft_core": 0, "shell": 0, "cloud": 0, "singleton": 0}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pixi run pytest tests/test_pangenome_report_tables.py -v
```

Expected: FAIL.

- [ ] **Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""Tidy, matplotlib-free aggregation tables for the pangenome island+Pfam
report -- joins significant_islands.tsv to island_pfam_enrichment.tsv,
computes size distribution, classification counts, and per-strain gene/bin
summaries. Deliberately kept free of any plotting dependency so these
tables are independently useful (e.g. to this repo's docs/ publishing
pipeline) and independently testable. Generalizes
studies/fungi/Afumigatus_pangenome/bin/annotate_islands_with_enrichment.py's
join logic, plus a new per_strain_summary (Opus review: this is how that
study's own accessory-count outlier strain was found, by hand -- generating
it automatically here).

Usage:
  pangenome_report_tables.py --significant_islands significant_islands.tsv \\
      --island_pfam_enrichment island_pfam_enrichment.tsv \\
      --pair_classification pair_classification.tsv \\
      --presence_matrix presence_matrix.tsv \\
      --frequency_table frequency_table.tsv \\
      --domtblout background_reps_vs_pfam.domtblout \\
      --out_dir report_tables/
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pangenome_domain_enrichment import parse_domtblout  # noqa: E402


def annotate_islands_with_domains(islands_rows: list[dict], family_domains: dict[str, set[str]]) -> list[dict]:
    """Adds a `pfam_domains` key (comma-joined sorted domain names, or `-`
    if none) to each island row, from the union of its member families'
    domains."""
    out = []
    for row in islands_rows:
        row = dict(row)
        members = row["member_families"].split(",")
        domains: set[str] = set()
        for m in members:
            domains |= family_domains.get(m, set())
        row["pfam_domains"] = ",".join(sorted(domains)) if domains else "-"
        out.append(row)
    return out


def island_size_distribution(islands_rows: list[dict]) -> dict[int, int]:
    """{island_size: count}."""
    dist: dict[int, int] = {}
    for row in islands_rows:
        size = int(row["island_size"])
        dist[size] = dist.get(size, 0) + 1
    return dist


def classification_counts(pair_classification_path: str) -> dict[str, int]:
    """{classification: count} across all rows in pair_classification.tsv."""
    counts: dict[str, int] = {}
    with open(pair_classification_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            c = row["classification"]
            counts[c] = counts.get(c, 0) + 1
    return counts


def per_strain_summary(presence_matrix_path: str, family_bin: dict[str, str]) -> list[dict]:
    """One row per strain: total genes present, plus counts broken down by
    frequency bin (core/soft_core/shell/cloud/singleton). `genome_only`
    counts as present (matches lib/pangenome_matrix.PresenceMatrix.is_present's
    semantics) -- rescued genome-only calls are real presence evidence, not
    a weaker state."""
    with open(presence_matrix_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        strains = header[1:]
        totals = {
            s: {"Short": s, "n_genes": 0, "core": 0, "soft_core": 0, "shell": 0, "cloud": 0, "singleton": 0}
            for s in strains
        }
        for row in reader:
            family = row[0]
            b = family_bin.get(family)
            for strain, call in zip(strains, row[1:]):
                if call != "absent":
                    totals[strain]["n_genes"] += 1
                    if b in totals[strain]:
                        totals[strain][b] += 1
    return list(totals.values())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--significant_islands", required=True)
    ap.add_argument("--island_pfam_enrichment", required=True)
    ap.add_argument("--pair_classification", required=True)
    ap.add_argument("--presence_matrix", required=True)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--domtblout", required=True, action="append")
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.significant_islands, newline="") as fh:
        islands_rows = list(csv.DictReader(fh, delimiter="\t"))

    family_domains = parse_domtblout(args.domtblout)
    annotated = annotate_islands_with_domains(islands_rows, family_domains)
    with open(out_dir / "islands_with_domains.tsv", "w", newline="") as out:
        fieldnames = list(annotated[0].keys()) if annotated else []
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in sorted(annotated, key=lambda r: -int(r["island_size"])):
            writer.writerow(row)

    dist = island_size_distribution(islands_rows)
    with open(out_dir / "island_size_distribution.tsv", "w") as out:
        out.write("island_size\tcount\n")
        for size, count in sorted(dist.items()):
            out.write(f"{size}\t{count}\n")

    counts = classification_counts(args.pair_classification)
    with open(out_dir / "classification_counts.tsv", "w") as out:
        out.write("classification\tcount\n")
        for classification, count in sorted(counts.items(), key=lambda kv: -kv[1]):
            out.write(f"{classification}\t{count}\n")

    family_bin: dict[str, str] = {}
    with open(args.frequency_table, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            family_bin[row["family"]] = row["bin"]
    strain_rows = per_strain_summary(args.presence_matrix, family_bin)
    with open(out_dir / "per_strain_summary.tsv", "w", newline="") as out:
        fieldnames = ["Short", "n_genes", "core", "soft_core", "shell", "cloud", "singleton"]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in sorted(strain_rows, key=lambda r: r["n_genes"]):
            writer.writerow(row)

    print(f"pangenome_report_tables: wrote islands_with_domains.tsv "
          f"({len(annotated)} islands), island_size_distribution.tsv, "
          f"classification_counts.tsv, per_strain_summary.tsv "
          f"({len(strain_rows)} strains) to {out_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pixi run pytest tests/test_pangenome_report_tables.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_report_tables.py tests/test_pangenome_report_tables.py
git commit -m "$(cat <<'EOF'
pangenome: add tidy report-tables aggregation + per-strain summary

per_strain_summary added per Opus report-completeness review: automates
the kind of by-hand outlier discovery the Afumigatus study's own QC
section relied on.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 6: `bin/pangenome_report_render.py` — figures + templated Markdown report

**Files:**
- Create: `bin/pangenome_report_render.py`
- Test: `tests/test_pangenome_report_render.py`

**Correction from an earlier draft (caught by Opus review, report-completeness
pass): the earlier draft of this task only implemented 3 figures (composition
pie, island-size distribution, classification bar) and silently dropped 3 of
the 5 figures the spec itself required** (design spec lines ~246-254:
"ports `plot_pangenome_summary.py`'s matplotlib figure generation --
frequency histogram, core/shell/cloud pie, presence/absence raster,
rarefaction curve, classification bar chart, **plus** new figures for
island-size distribution and top enriched domains"). This was an
undocumented oversight, not a deliberate scope decision — this revision
restores the frequency histogram, presence/absence raster, and
accumulation/rarefaction curve (ported from
`studies/fungi/Afumigatus_pangenome/bin/plot_pangenome_summary.py`), adds the
spec-required top-enriched-domains figure (which Afumigatus itself never had
a checked-in generator for — `results/figures/island_domain_enrichment.png`
was hand-made and is not reproducible today; this makes it reproducible for
the first time), and adds a Heaps'-law openness fit + core-genome asymptote
fit as a genuine scientific upgrade over the Afumigatus precedent (which only
asserted "pangenome is open" from the curve's visual shape, with no fitted
statistic).

**Files consumed now include `presence_matrix.tsv`** (new — needed for the
raster and accumulation-curve figures) **and `per_strain_summary.tsv`** (Task 5's
new output).

**Interfaces:**
- Consumes: Task 5's four output TSVs (`islands_with_domains.tsv`,
  `island_size_distribution.tsv`, `classification_counts.tsv`,
  `per_strain_summary.tsv`) plus `frequency_table.tsv` and `presence_matrix.tsv`.
- Produces: `render_report_markdown(counts: dict, size_dist: dict,
  classification_counts_dict: dict, top_domains: list[dict], n_islands: int,
  heaps_fit: dict, core_decay: dict, strain_gene_counts: list[int]) -> str`
  (the full Markdown text) — no later task depends on further functions from
  this file, but keep this one testable-without-matplotlib by having `main()`
  call it separately from the figure-generation calls.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_pangenome_report_render.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_report_render import render_report_markdown, fit_heaps_law, fit_core_decay
import numpy as np


def test_render_report_markdown_includes_key_sections():
    md = render_report_markdown(
        counts={"core": 100, "soft_core": 10, "shell": 50, "cloud": 200, "singleton": 300},
        size_dist={2: 50, 3: 20, 10: 5},
        classification_counts_dict={"trans": 40, "unexplained_physical": 10},
        top_domains=[{"domain": "PF00001", "fisher_p": 1e-5, "fdr_q": 1e-4}],
        n_islands=75,
        heaps_fit={"kappa": 500.0, "gamma": 0.4, "r_squared": 0.95, "is_open": True},
        core_decay={"core_inf": 90.0, "tau": 20.0, "fit_ok": True},
        strain_gene_counts=[8000, 8100, 8050],
    )
    assert "# Pangenome Island + Pfam Enrichment Report" in md
    assert "## Pangenome composition" in md
    assert "## Accessory islands" in md
    assert "## Pfam domain enrichment" in md
    assert "## Pangenome openness" in md
    assert "PF00001" in md
    assert "75" in md  # n_islands appears somewhere
    assert "open" in md.lower()


def test_render_report_markdown_handles_zero_enriched_domains():
    md = render_report_markdown(
        counts={"core": 1, "soft_core": 0, "shell": 0, "cloud": 0, "singleton": 0},
        size_dist={},
        classification_counts_dict={},
        top_domains=[],
        n_islands=0,
        heaps_fit=None,
        core_decay=None,
        strain_gene_counts=[],
    )
    assert "No significantly enriched" in md


def test_fit_heaps_law_detects_open_pangenome():
    # A slowly-still-growing pangenome curve (gamma < 1 -> open)
    n = np.arange(1, 21)
    pan_mean = 100 * n**0.3  # gamma=0.3, kappa=100 by construction
    result = fit_heaps_law(pan_mean)
    assert result["is_open"] is True
    assert abs(result["gamma"] - 0.3) < 0.05
    assert result["r_squared"] > 0.99


def test_fit_core_decay_returns_asymptote():
    n = np.arange(1, 21).astype(float)
    core_mean = 50 + 200 * np.exp(-n / 5)  # decays toward asymptote 50
    result = fit_core_decay(core_mean)
    assert result["fit_ok"] is True
    assert abs(result["core_inf"] - 50) < 5
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pixi run pytest tests/test_pangenome_report_render.py -v
```

Expected: FAIL.

- [ ] **Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""Figures (matplotlib, PNG+PDF per figure) + templated Markdown report for
the pangenome island+Pfam enrichment step. Fully automated -- no
hand-written narrative -- so it works for any future pangenome study, not
just Afumigatus (whose REPORT.md was hand-assembled prose; this is not
that). Figure functions generalize
studies/fungi/Afumigatus_pangenome/bin/plot_pangenome_summary.py's five
figures (frequency histogram, composition pie, presence/absence raster,
accumulation/rarefaction curve, classification bar), plus two new figures
(top enriched Pfam domains, which Afumigatus embedded by hand with no
checked-in generator) and a Heaps'-law openness fit + core-genome
asymptote fit that Afumigatus's own report never computed (it only
asserted openness from the curve's visual shape).

Usage:
  pangenome_report_render.py --frequency_table frequency_table.tsv \\
      --presence_matrix presence_matrix.tsv \\
      --islands_with_domains islands_with_domains.tsv \\
      --island_size_distribution island_size_distribution.tsv \\
      --classification_counts classification_counts.tsv \\
      --island_pfam_enrichment island_pfam_enrichment.tsv \\
      --per_strain_summary per_strain_summary.tsv \\
      --n_permutations 20 --seed 0 \\
      --out_dir report/
"""
from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix  # noqa: E402

BAND_ORDER = ["core", "soft_core", "shell", "cloud", "singleton"]
BAND_COLORS = {
    "core": "#2c7fb8", "soft_core": "#7fcdbb", "shell": "#fed976",
    "cloud": "#fd8d3c", "singleton": "#bdbdbd",
}


def _savefig_both(fig, out_dir: Path, name: str) -> None:
    """Writes both <name>.png (into out_dir/figures/) and <name>.pdf
    (into out_dir/figures_pdf/) -- matches Afumigatus's convention of a
    PNG for the Markdown embed and a parallel PDF vector copy for print."""
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)
    (out_dir / "figures_pdf").mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "figures" / f"{name}.png", dpi=150)
    fig.savefig(out_dir / "figures_pdf" / f"{name}.pdf")


def plot_frequency_distribution(frequency_table_rows: list[dict], out_dir: Path) -> None:
    """Ported from plot_pangenome_summary.py's plot_frequency_distribution --
    family-frequency histogram, one series per band."""
    fig, ax = plt.subplots(figsize=(8, 5))
    for band in BAND_ORDER:
        freqs = [float(r["frequency"]) for r in frequency_table_rows if r["bin"] == band]
        if freqs:
            ax.hist(freqs, bins=40, range=(0, 1), color=BAND_COLORS[band],
                     label=f"{band} (n={len(freqs)})", alpha=0.85)
    ax.set_xlabel("Fraction of strains carrying the family")
    ax.set_ylabel("Number of gene families")
    ax.set_title("Family-frequency distribution")
    ax.legend()
    fig.tight_layout()
    _savefig_both(fig, out_dir, "frequency_distribution")
    plt.close(fig)


def plot_frequency_bins(counts: dict[str, int], out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    labels = [k for k in BAND_ORDER if counts.get(k, 0) > 0]
    values = [counts[k] for k in labels]
    ax.pie(values, labels=labels, autopct="%1.1f%%", colors=[BAND_COLORS[b] for b in labels])
    ax.set_title("Pangenome composition")
    fig.tight_layout()
    _savefig_both(fig, out_dir, "core_shell_cloud_pie")
    plt.close(fig)


def plot_presence_absence_matrix(matrix: PresenceMatrix, frequency_table_rows: list[dict], out_dir: Path) -> None:
    """Ported from plot_pangenome_summary.py's plot_presence_absence_matrix
    -- families (rows, frequency-sorted) x strains raster via imshow, kept
    fast at tens-of-thousands-of-families scale by rasterizing rather than
    drawing one element per cell."""
    freq_by_family = {r["family"]: float(r["frequency"]) for r in frequency_table_rows}
    family_order = sorted((f for f in matrix.families if f in freq_by_family), key=lambda f: -freq_by_family[f])
    strain_order = list(matrix.strains)
    arr = np.zeros((len(family_order), len(strain_order)), dtype=bool)
    for i, fam in enumerate(family_order):
        for j, strain in enumerate(strain_order):
            arr[i, j] = matrix.is_present(fam, strain)

    fig, ax = plt.subplots(figsize=(8, 10))
    ax.imshow(arr, aspect="auto", cmap="Greys", interpolation="nearest")
    ax.set_xlabel(f"Strains (n={arr.shape[1]})")
    ax.set_ylabel(f"Gene families (n={arr.shape[0]}), sorted by frequency")
    ax.set_title("Presence/absence matrix")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    _savefig_both(fig, out_dir, "presence_absence_matrix")
    plt.close(fig)


def accumulation_curve(matrix: PresenceMatrix, n_permutations: int = 20, seed: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Ported from plot_pangenome_summary.py's accumulation_curve --
    (pan_mean, pan_std, core_mean, core_std) averaged over `n_permutations`
    random strain orderings."""
    n_strains = len(matrix.strains)
    n_families = len(matrix.families)
    arr = np.zeros((n_families, n_strains), dtype=bool)
    for i, fam in enumerate(matrix.families):
        for j, strain in enumerate(matrix.strains):
            arr[i, j] = matrix.is_present(fam, strain)

    rng = random.Random(seed)
    pan_runs = np.zeros((n_permutations, n_strains), dtype=int)
    core_runs = np.zeros((n_permutations, n_strains), dtype=int)
    strain_indices = list(range(n_strains))
    for p in range(n_permutations):
        order = strain_indices[:]
        rng.shuffle(order)
        cum_any = np.zeros(n_families, dtype=bool)
        cum_all = np.ones(n_families, dtype=bool)
        for k, s in enumerate(order):
            cum_any |= arr[:, s]
            cum_all &= arr[:, s]
            pan_runs[p, k] = cum_any.sum()
            core_runs[p, k] = cum_all.sum()
    return (
        pan_runs.mean(axis=0), pan_runs.std(axis=0),
        core_runs.mean(axis=0), core_runs.std(axis=0),
    )


def fit_heaps_law(pan_mean: np.ndarray) -> dict:
    """Fit Heaps' law n = kappa * N^gamma to the pangenome accumulation
    curve via log-log linear regression. gamma < 1 conventionally indicates
    an OPEN pangenome (new families keep appearing as more genomes are
    added); gamma >= 1 indicates a closed one. This is the fitted openness
    statistic Afumigatus's own report never computed -- it only asserted
    openness from the curve's visual shape."""
    n_strains = np.arange(1, len(pan_mean) + 1)
    log_n = np.log(n_strains)
    log_pan = np.log(pan_mean)
    gamma, log_kappa = np.polyfit(log_n, log_pan, 1)
    kappa = np.exp(log_kappa)
    predicted = log_kappa + gamma * log_n
    ss_res = np.sum((log_pan - predicted) ** 2)
    ss_tot = np.sum((log_pan - log_pan.mean()) ** 2)
    r_squared = float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    return {"kappa": float(kappa), "gamma": float(gamma), "r_squared": r_squared, "is_open": bool(gamma < 1)}


def fit_core_decay(core_mean: np.ndarray) -> dict:
    """Fit an exponential-decay model core(n) = core_inf + amplitude *
    exp(-n / tau) to the core-genome accumulation curve, reporting the
    extrapolated asymptotic core-genome size (core_inf) -- the standard
    companion statistic to Heaps' law openness."""
    from scipy.optimize import curve_fit

    n_strains = np.arange(1, len(core_mean) + 1).astype(float)

    def model(n, core_inf, amplitude, tau):
        return core_inf + amplitude * np.exp(-n / tau)

    try:
        popt, _ = curve_fit(
            model, n_strains, core_mean,
            p0=[core_mean[-1], core_mean[0] - core_mean[-1], max(len(core_mean) / 4, 1)],
            maxfev=10000,
        )
        core_inf, _amplitude, tau = popt
        return {"core_inf": float(core_inf), "tau": float(tau), "fit_ok": True}
    except RuntimeError:
        return {"core_inf": float(core_mean[-1]), "tau": float("nan"), "fit_ok": False}


def plot_accumulation_curve(pan_mean, pan_std, core_mean, core_std, heaps_fit: dict, out_dir: Path) -> None:
    n = len(pan_mean)
    x = np.arange(1, n + 1)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, pan_mean, color="#2c7fb8", label="Pangenome (union)")
    ax.fill_between(x, pan_mean - pan_std, pan_mean + pan_std, color="#2c7fb8", alpha=0.2)
    ax.plot(x, core_mean, color="#d95f02", label="Core (intersection)")
    ax.fill_between(x, core_mean - core_std, core_mean + core_std, color="#d95f02", alpha=0.2)
    ax.set_xlabel("Number of strains sampled")
    ax.set_ylabel("Number of gene families")
    openness = "open" if heaps_fit["is_open"] else "closed"
    ax.set_title(f"Accumulation curve (Heaps' γ={heaps_fit['gamma']:.2f}, {openness})")
    ax.legend()
    fig.tight_layout()
    _savefig_both(fig, out_dir, "accumulation_curve")
    plt.close(fig)


def plot_classification_counts(classification_counts_dict: dict[str, int], out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = list(classification_counts_dict.keys())
    values = [classification_counts_dict[k] for k in labels]
    ax.bar(labels, values)
    ax.set_ylabel("Number of pairs")
    ax.set_title("Pair classification breakdown")
    fig.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    _savefig_both(fig, out_dir, "pair_classification_summary")
    plt.close(fig)


def plot_island_size_distribution(size_dist: dict[int, int], out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    sizes = sorted(size_dist)
    counts = [size_dist[s] for s in sizes]
    ax.bar([str(s) for s in sizes], counts)
    ax.set_xlabel("Island size (genes)")
    ax.set_ylabel("Number of islands")
    ax.set_title("Accessory island size distribution")
    fig.tight_layout()
    _savefig_both(fig, out_dir, "island_size_distribution")
    plt.close(fig)


def plot_domain_enrichment(top_domains: list[dict], out_dir: Path, top_n: int = 20) -> None:
    """Horizontal bar chart of the top-N (by fdr_q) significantly enriched
    Pfam domains among island-member families -- the figure
    (`island_domain_enrichment.png`) Afumigatus's own REPORT.md embedded by
    hand with NO checked-in generator; this makes it reproducible."""
    domains = top_domains[:top_n]
    if not domains:
        return
    fig, ax = plt.subplots(figsize=(8, max(4, 0.3 * len(domains))))
    labels = [d["domain"] for d in reversed(domains)]
    neg_log_q = [-np.log10(float(d["fdr_q"])) if float(d["fdr_q"]) > 0 else 300.0 for d in reversed(domains)]
    ax.barh(labels, neg_log_q, color="#2c7fb8")
    ax.set_xlabel("-log10(FDR q)")
    ax.set_title(f"Top {len(domains)} enriched Pfam domains in accessory islands")
    fig.tight_layout()
    _savefig_both(fig, out_dir, "island_domain_enrichment")
    plt.close(fig)


def render_report_markdown(
    counts: dict[str, int],
    size_dist: dict[int, int],
    classification_counts_dict: dict[str, int],
    top_domains: list[dict],
    n_islands: int,
    heaps_fit: dict | None,
    core_decay: dict | None,
    strain_gene_counts: list[int],
) -> str:
    total_families = sum(counts.values())
    lines = ["# Pangenome Island + Pfam Enrichment Report", ""]
    lines += ["## Pangenome composition", ""]
    lines += [f"Total families: {total_families}", ""]
    for label in BAND_ORDER:
        if counts.get(label, 0):
            pct = 100 * counts[label] / total_families if total_families else 0
            lines.append(f"- **{label}**: {counts[label]} ({pct:.1f}%)")
    lines += ["", "![Composition](figures/core_shell_cloud_pie.png)", ""]
    lines += ["![Frequency distribution](figures/frequency_distribution.png)", ""]

    if heaps_fit is not None:
        lines += ["## Pangenome openness", ""]
        openness = "open" if heaps_fit["is_open"] else "closed"
        lines += [f"Heaps' law fit: κ={heaps_fit['kappa']:.1f}, γ={heaps_fit['gamma']:.3f} "
                  f"(R²={heaps_fit['r_squared']:.3f}) -- pangenome is **{openness}** "
                  f"(γ {'<' if heaps_fit['is_open'] else '>='} 1).", ""]
        if core_decay is not None and core_decay["fit_ok"]:
            lines += [f"Extrapolated asymptotic core-genome size: {core_decay['core_inf']:.0f} families.", ""]
        lines += ["![Accumulation curve](figures/accumulation_curve.png)", ""]
        lines += ["![Presence/absence matrix](figures/presence_absence_matrix.png)", ""]

    if strain_gene_counts:
        lines += ["## Per-strain summary", ""]
        lines += [f"Genes per strain: min={min(strain_gene_counts)}, "
                  f"median={sorted(strain_gene_counts)[len(strain_gene_counts)//2]}, "
                  f"max={max(strain_gene_counts)} (n={len(strain_gene_counts)} strains). "
                  "See per_strain_summary.tsv for outliers.", ""]

    lines += ["## Accessory islands", ""]
    lines += [f"{n_islands} statistically significant accessory islands found "
              "(built from adjacency of non-core genes, gated by containing "
              "at least one FDR-significant physically-linked pair).", ""]
    if size_dist:
        lines += ["![Island sizes](figures/island_size_distribution.png)", ""]

    lines += ["## Pair classification breakdown", ""]
    for classification, count in sorted(classification_counts_dict.items(), key=lambda kv: -kv[1]):
        lines.append(f"- **{classification}**: {count}")
    lines += ["", "![Classification breakdown](figures/pair_classification_summary.png)", ""]

    lines += ["## Pfam domain enrichment", ""]
    if not top_domains:
        lines += ["No significantly enriched Pfam domains found.", ""]
    else:
        lines += ["![Top enriched domains](figures/island_domain_enrichment.png)", ""]
        lines += ["| Domain | Fisher p | FDR q |", "|---|---|---|"]
        for row in top_domains:
            lines.append(f"| {row['domain']} | {float(row['fisher_p']):.2e} | {float(row['fdr_q']):.2e} |")
        lines.append("")

    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--presence_matrix", required=True)
    ap.add_argument("--islands_with_domains", required=True)
    ap.add_argument("--island_size_distribution", required=True)
    ap.add_argument("--classification_counts", required=True)
    ap.add_argument("--island_pfam_enrichment", required=True)
    ap.add_argument("--per_strain_summary", required=True)
    ap.add_argument("--fdr_threshold", type=float, default=0.05)
    ap.add_argument("--n_permutations", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.frequency_table, newline="") as fh:
        frequency_table_rows = list(csv.DictReader(fh, delimiter="\t"))
    counts: dict[str, int] = {}
    for row in frequency_table_rows:
        counts[row["bin"]] = counts.get(row["bin"], 0) + 1

    size_dist: dict[int, int] = {}
    with open(args.island_size_distribution, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            size_dist[int(row["island_size"])] = int(row["count"])

    classification_counts_dict: dict[str, int] = {}
    with open(args.classification_counts, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            classification_counts_dict[row["classification"]] = int(row["count"])

    with open(args.islands_with_domains, newline="") as fh:
        n_islands = sum(1 for _ in csv.DictReader(fh, delimiter="\t"))

    top_domains: list[dict] = []
    with open(args.island_pfam_enrichment, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if float(row["fdr_q"]) < args.fdr_threshold:
                top_domains.append(row)
    top_domains.sort(key=lambda r: float(r["fisher_p"]))

    strain_gene_counts: list[int] = []
    with open(args.per_strain_summary, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            strain_gene_counts.append(int(row["n_genes"]))

    plot_frequency_distribution(frequency_table_rows, out_dir)
    plot_frequency_bins(counts, out_dir)
    if size_dist:
        plot_island_size_distribution(size_dist, out_dir)
    if classification_counts_dict:
        plot_classification_counts(classification_counts_dict, out_dir)
    plot_domain_enrichment(top_domains, out_dir)

    matrix = PresenceMatrix.from_tsv(args.presence_matrix)
    plot_presence_absence_matrix(matrix, frequency_table_rows, out_dir)
    pan_mean, pan_std, core_mean, core_std = accumulation_curve(
        matrix, n_permutations=args.n_permutations, seed=args.seed,
    )
    heaps_fit = fit_heaps_law(pan_mean)
    core_decay = fit_core_decay(core_mean)
    plot_accumulation_curve(pan_mean, pan_std, core_mean, core_std, heaps_fit, out_dir)

    with open(out_dir / "pangenome_openness.tsv", "w") as out:
        out.write("kappa\tgamma\tr_squared\tis_open\tcore_inf\ttau\n")
        out.write(f"{heaps_fit['kappa']:.4f}\t{heaps_fit['gamma']:.4f}\t{heaps_fit['r_squared']:.4f}\t"
                  f"{heaps_fit['is_open']}\t{core_decay['core_inf']:.4f}\t{core_decay['tau']:.4f}\n")

    markdown = render_report_markdown(
        counts, size_dist, classification_counts_dict, top_domains, n_islands,
        heaps_fit, core_decay, strain_gene_counts,
    )
    (out_dir / "report.md").write_text(markdown)

    print(f"pangenome_report_render: wrote {out_dir}/report.md, pangenome_openness.tsv, "
          f"and figures/, figures_pdf/", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pixi run pytest tests/test_pangenome_report_render.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add bin/pangenome_report_render.py tests/test_pangenome_report_render.py
git commit -m "$(cat <<'EOF'
pangenome: restore dropped figures, add Heaps' law openness fit

An earlier draft of this task silently dropped 3 of the 5 spec-required
figures (frequency histogram, presence/absence raster, accumulation
curve) and never added the spec-required top-enriched-domains figure --
caught by Opus review before implementation started. This restores all
five, adds a domain-enrichment figure Afumigatus itself only ever made by
hand (no checked-in generator), and adds a Heaps' law + core-genome-decay
fit as a genuine statistical upgrade over that study's own report, which
only asserted openness from the curve's visual shape.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 7: `modules/pangenome/report.nf` — REPORT_TABLES / REPORT_RENDER modules

**Files:**
- Create: `modules/pangenome/report.nf`

**Interfaces:**
- Consumes: Task 5's `pangenome_report_tables.py` CLI, Task 6's
  `pangenome_report_render.py` CLI.
- Produces: `REPORT_TABLES`, `REPORT_RENDER` Nextflow processes — Task 8's
  workflow wiring calls these by exact name.

- [ ] **Step 1: Write `modules/pangenome/report.nf`**

**Correction from an earlier draft**: `REPORT_TABLES` and `REPORT_RENDER` now
also take `presence_matrix`/`frequency_table` (for `per_strain_summary.tsv`
and the restored figures) — see Task 5/6's corrections.

```groovy
// REPORT_TABLES / REPORT_RENDER -- tidy aggregation tables, then
// figures+Markdown, for the pangenome island+Pfam enrichment step. Split
// into two processes (not one) so table aggregation stays testable without
// a matplotlib dependency and independently reusable (e.g. by this repo's
// docs/ publishing pipeline).
process REPORT_TABLES {
    label 'low_cpu'
    tag "report_tables"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome/report_tables" }, mode: 'copy'

    input:
    path(significant_islands)
    path(island_pfam_enrichment)
    path(pair_classification)
    path(presence_matrix)
    path(frequency_table)
    path(domtblout)

    output:
    path("islands_with_domains.tsv"), emit: islands_with_domains
    path("island_size_distribution.tsv"), emit: size_distribution
    path("classification_counts.tsv"), emit: classification_counts
    path("per_strain_summary.tsv"), emit: per_strain_summary

    script:
    """
    pangenome_report_tables.py \
        --significant_islands ${significant_islands} \
        --island_pfam_enrichment ${island_pfam_enrichment} \
        --pair_classification ${pair_classification} \
        --presence_matrix ${presence_matrix} \
        --frequency_table ${frequency_table} \
        --domtblout ${domtblout} \
        --out_dir .
    """
}

process REPORT_RENDER {
    label 'low_cpu'
    tag "report_render"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/pangenome" }, mode: 'copy'

    input:
    path(frequency_table)
    path(presence_matrix)
    path(islands_with_domains)
    path(island_size_distribution)
    path(classification_counts)
    path(island_pfam_enrichment)
    path(per_strain_summary)

    output:
    path("report/report.md"), emit: report
    path("report/pangenome_openness.tsv"), emit: openness
    path("report/figures/*"), emit: figures_png
    path("report/figures_pdf/*"), emit: figures_pdf

    script:
    """
    pangenome_report_render.py \
        --frequency_table ${frequency_table} \
        --presence_matrix ${presence_matrix} \
        --islands_with_domains ${islands_with_domains} \
        --island_size_distribution ${island_size_distribution} \
        --classification_counts ${classification_counts} \
        --island_pfam_enrichment ${island_pfam_enrichment} \
        --per_strain_summary ${per_strain_summary} \
        --n_permutations ${params.pangenome_accumulation_permutations} \
        --seed ${params.pangenome_accumulation_seed} \
        --out_dir report
    """
}
```

- [ ] **Step 2: Commit**

```bash
git add modules/pangenome/report.nf
git commit -m "$(cat <<'EOF'
pangenome: add REPORT_TABLES/REPORT_RENDER modules

Wires presence_matrix/frequency_table through for the restored figures
(Task 5/6 corrections) alongside the original island/enrichment inputs.
Not yet wired into pangenome_profile.nf (Task 8).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 8: Wire everything into `workflows/pangenome_profile.nf` + new params

**Files:**
- Modify: `workflows/pangenome_profile.nf`
- Modify: `nextflow.config` (add the new params from the spec's Params table)
- Modify: `conf/ucr_hpcc_slurm.config` (add resource labels for the new
  processes, following this file's existing `withName:` pattern)

**Interfaces:**
- Consumes: `BUILD_ISLANDS`, `MARKER_HMMSEARCH`, `SELECT_BACKGROUND_REPS`,
  `FAMILY_PFAM_SCAN`, `DOMAIN_ENRICHMENT` (Task 4), `REPORT_TABLES`,
  `REPORT_RENDER` (Task 7). Does **not** touch or import from
  `modules/pangenome/captain.nf` — see Task 4's correction.

**Three corrections from an earlier draft, all caught by Opus review before
implementation started (verified directly against the real code, not just
asserted):**

1. **Param name collision.** `params.pangenome_pfam_hmm` **already exists**
   in `nextflow.config` (confirmed: `nextflow.config:251`,
   `pangenome_pfam_hmm = null // ...this larger Pfam-style HMM database`)
   and is genuinely load-bearing for the *existing* captain-by-name branch
   (`workflows/pangenome_profile.nf`'s `captain_requested` check:
   `params.pangenome_captain_hmm_name && params.pangenome_pfam_hmm`). An
   earlier draft of this task reused that same name as the gate for this
   new step — that would both re-declare an existing config key and
   silently arm the whole new island+Pfam branch for every existing
   captain-by-name study. Fixed: a **new, distinct** param,
   `pangenome_island_pfam_hmm`, gates this step; the pre-existing
   `pangenome_pfam_hmm` is untouched and unrelated.
2. **Illegal repeated process invocation.** See Task 4's correction — the
   marker search is now `MARKER_HMMSEARCH`, invoked ONCE via a
   `Channel.fromList` of `[name, hmm]` tuples (normal Nextflow multiplicity),
   not a Groovy loop calling a process repeatedly (illegal DSL2).
3. **`collect()` flattens tuples by default.** An earlier draft's
   `marker_tblouts_ch.collect()` over a channel of `[name, tblout]` pairs
   would silently flatten to `[name1, tblout1, name2, tblout2, ...]`, not
   `[[name1, tblout1], [name2, tblout2]]`. Fixed by using `.toList()`
   (preserves nesting) and splitting into two **parallel** lists
   (`marker_names`, `marker_tblout_files`) before passing to `BUILD_ISLANDS`
   — also sidesteps a related staging concern (a `path`-qualified input
   is what actually gets Nextflow to stage the files into the task's work
   directory; a `val` containing embedded file objects is a subtler,
   harder-to-verify pattern this avoids entirely).

- [ ] **Step 1: Add new params to `nextflow.config`**, near the other
  `pangenome_*` params (find them with `grep -n "pangenome_" nextflow.config`
  first, and add these alongside in the same block — do **not** touch the
  existing `pangenome_pfam_hmm` line):

```groovy
    pangenome_island_pfam_hmm = null       // enables the whole island+Pfam step when set (distinct from the existing pangenome_pfam_hmm, used only by the captain-by-name branch)
    pangenome_island_min_size = 2
    pangenome_pfam_domain_evalue = 1e-3
    pangenome_marker_names = ''            // comma list, e.g. 'captain,sm_backbone'
    pangenome_marker_hmm_paths = ''        // parallel comma list of HMM paths
    pangenome_marker_evalue = 1e-5
    pangenome_accumulation_permutations = 20
    pangenome_accumulation_seed = 0
```

- [ ] **Step 2: Read `workflows/pangenome_profile.nf`'s existing captain-branch
  code** (the `captain_requested` block) and confirm the `rescued_matrix`
  local variable name (set earlier in the same workflow body, in the rescue-pass
  section) — this task's new branch reuses it directly for
  `presence_matrix`, matching what `FREQUENCY_BINS`/`COOCCURRENCE` already
  consume, rather than the pre-rescue `PRESENCE_MATRIX.out.matrix`.

- [ ] **Step 3: Add the `include` lines to the top of `workflows/pangenome_profile.nf`**
  (alongside the existing `include` block — Nextflow `include` statements are
  script-level only and cannot appear inside the `workflow { }` body):

```groovy
include { BUILD_ISLANDS; MARKER_HMMSEARCH }                                from '../modules/pangenome/islands'
include { SELECT_BACKGROUND_REPS; FAMILY_PFAM_SCAN; DOMAIN_ENRICHMENT }     from '../modules/pangenome/pfam_enrichment'
include { REPORT_TABLES; REPORT_RENDER }                                   from '../modules/pangenome/report'
```

- [ ] **Step 4: Add the new island+Pfam branch inside the `workflow { }` body**,
  after the existing `GENE_POSITIONS`/`FAMILY_POSITIONS` section (check the
  file's own numbered-comment convention, e.g. `// --- 6. ... ---`, and
  continue that numbering):

```groovy
// --- N. Accessory islands + Pfam functional enrichment (optional) --------
// Named marker searches (0+): ONE MARKER_HMMSEARCH invocation over a
// Channel.fromList of [name, hmm_path] tuples runs one task per marker via
// Nextflow's normal channel-based fan-out -- not a loop calling the process
// repeatedly (DSL2 forbids invoking a process more than once per workflow).
// CONCAT_PROTEOMES.out.fasta is a single-emission channel; `.first()` turns
// it into a value reused for every marker task rather than being consumed
// by only the first one.
if (params.pangenome_island_pfam_hmm) {
    def marker_names_list = params.pangenome_marker_names ? params.pangenome_marker_names.split(',') as List : []
    def marker_hmm_paths_list = params.pangenome_marker_hmm_paths ? params.pangenome_marker_hmm_paths.split(',') as List : []
    if (marker_names_list.size() != marker_hmm_paths_list.size()) {
        error "ERROR: --pangenome_marker_names and --pangenome_marker_hmm_paths must have " +
              "the same number of comma-separated entries (got ${marker_names_list.size()} names, " +
              "${marker_hmm_paths_list.size()} paths)"
    }

    if (marker_names_list) {
        marker_input_ch = Channel.fromList(
            [marker_names_list, marker_hmm_paths_list.collect { file(it) }].transpose()
        )
        MARKER_HMMSEARCH(marker_input_ch, CONCAT_PROTEOMES.out.fasta.first())
        // .toList() (NOT .collect(), which flattens [[n,p],[n,p]] to
        // [n,p,n,p] by default) keeps each [name, tblout] pair intact as one
        // list-of-pairs value; .collect{it[0]}/.collect{it[1]} then split
        // that into the two PARALLEL lists BUILD_ISLANDS's val+path inputs
        // expect (see Task 4).
        marker_pairs_ch = MARKER_HMMSEARCH.out.result.toList()
        marker_names_ch = marker_pairs_ch.map { pairs -> pairs.collect { it[0] } }
        marker_tblout_files_ch = marker_pairs_ch.map { pairs -> pairs.collect { it[1] } }
    } else {
        marker_names_ch = Channel.value([])
        marker_tblout_files_ch = Channel.value([])
    }

    BUILD_ISLANDS(
        FAMILY_POSITIONS.out.positions,
        FREQUENCY_BINS.out.table,
        PAIR_CLASSIFICATION.out.classification,
        CLUSTER_TIER1.out.cluster_tsv,
        marker_names_ch,
        marker_tblout_files_ch,
    )

    SELECT_BACKGROUND_REPS(CLUSTER_TIER1.out.rep_fasta, FREQUENCY_BINS.out.table)
    FAMILY_PFAM_SCAN(SELECT_BACKGROUND_REPS.out.fasta, file(params.pangenome_island_pfam_hmm))
    DOMAIN_ENRICHMENT(BUILD_ISLANDS.out.islands, FAMILY_PFAM_SCAN.out.domtblout, FREQUENCY_BINS.out.table)

    REPORT_TABLES(
        BUILD_ISLANDS.out.islands,
        DOMAIN_ENRICHMENT.out.enrichment,
        PAIR_CLASSIFICATION.out.classification,
        rescued_matrix,
        FREQUENCY_BINS.out.table,
        FAMILY_PFAM_SCAN.out.domtblout,
    )
    REPORT_RENDER(
        FREQUENCY_BINS.out.table,
        rescued_matrix,
        REPORT_TABLES.out.islands_with_domains,
        REPORT_TABLES.out.size_distribution,
        REPORT_TABLES.out.classification_counts,
        DOMAIN_ENRICHMENT.out.enrichment,
        REPORT_TABLES.out.per_strain_summary,
    )
}
```

Note: the exact upstream channel names above (`FAMILY_POSITIONS.out.positions`,
`CLUSTER_TIER1.out.cluster_tsv`, `CLUSTER_TIER1.out.rep_fasta`,
`FREQUENCY_BINS.out.table`, `PAIR_CLASSIFICATION.out.classification`,
`CONCAT_PROTEOMES.out.fasta`, `rescued_matrix`) were verified directly against
each module's real `output:`/`emit:` block (`grep -n "emit:" modules/pangenome/*.nf`)
and the workflow's own local variable assignment before writing this task —
all confirmed correct as written above. If a future edit to any of those
modules renames an `emit:`, this step will fail loudly with a Nextflow
compile error naming the missing channel, not silently.

- [ ] **Step 5: Add resource-label overrides for the new heavy processes**
  (`FAMILY_PFAM_SCAN` is a full Pfam-A.hmm scan over potentially thousands of
  family reps; `REPORT_RENDER`'s accumulation-curve permutation loop over the
  full presence matrix is also non-trivial at real scale) to
  `conf/ucr_hpcc_slurm.config`, in the same `// --- Pangenome cluster-profiling
  subworkflow (pangenome.nf) ---` section as the other pangenome overrides:

```groovy
    withName: '.*FAMILY_PFAM_SCAN' {
        cpus   = 8
        memory = '16.GB'
        time   = { task.attempt > 1 ? '8.h' : '2.h' }
    }
    withName: '.*REPORT_RENDER' {
        memory = '8.GB'
        time   = { task.attempt > 1 ? '4.h' : '1.h' }
    }
```

- [ ] **Step 6: Verify the whole pipeline still parses**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
nextflow run pangenome.nf --help > /tmp/nf_syntax_check2.out 2>&1; echo "exit: $?"
cat /tmp/nf_syntax_check2.out | tail -20
```

Expected: `exit: 0`, help text prints. If a channel-name mismatch error
appears, fix the specific `.out.<name>` reference and re-run.

- [ ] **Step 7: Real execution smoke test (not just `--help`) — required, not
  optional.** The marker-search channel wiring in Step 4 (tuple-based
  fan-out, `.toList()`, splitting into parallel lists) is genuinely subtle
  Nextflow channel semantics that a `--help` parse check cannot validate —
  it only catches syntax errors, not runtime channel-shape mismatches. Reuse
  the small samplesheet/data_dir already built and verified working for the
  Coccidioides pangenome study's own small-subset validation
  (`studies/fungi/coccidioides_pangenome/data_dir`,
  `studies/fungi/coccidioides_pangenome/config.csv` — build a fresh 10-strain
  subset the same way that study's Task 6 did: first 5 rows per species from
  `config.csv`, written to a temp file) and a real Pfam-A.hmm
  (`/bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm`,
  confirmed present on disk). Run **twice** — once with zero markers
  configured (the empty-list path), once with one marker configured (the
  `Channel.fromList`/fan-out path) — since these are the two structurally
  different branches in Step 4's code:

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio

# Zero-marker run:
nextflow run pangenome.nf \
    --pangenome_samplesheet /path/to/config_smoketest.csv \
    --pangenome_data_dir /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/data_dir \
    --pangenome_island_pfam_hmm /bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm \
    --pangenome_project island_pfam_smoketest_nomarker \
    --outdir "$SCRATCH/island_pfam_smoketest_nomarker" \
    -profile local

# One-marker run (any real single-profile HMM works for this check; the
# content doesn't matter, only that the fan-out wiring executes correctly):
nextflow run pangenome.nf \
    --pangenome_samplesheet /path/to/config_smoketest.csv \
    --pangenome_data_dir /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/data_dir \
    --pangenome_island_pfam_hmm /bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm \
    --pangenome_marker_names test_marker \
    --pangenome_marker_hmm_paths /bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm \
    --pangenome_project island_pfam_smoketest_onemarker \
    --outdir "$SCRATCH/island_pfam_smoketest_onemarker" \
    -profile local
```

Expected for both: completes without error; `significant_islands.tsv` has a
`has_test_marker` column in the second run only (zero columns of that shape
in the first); `report/report.md` and `report/figures/*.png` (7 figures:
frequency_distribution, core_shell_cloud_pie, presence_absence_matrix,
accumulation_curve, pair_classification_summary, island_size_distribution,
island_domain_enrichment — the last two only if any islands/enriched domains
exist at this small scale, which may legitimately be zero, matching the
same "small-N gives no significant results" pattern the original onboarding
plan's Task 6 already established for `PAIR_CLASSIFICATION`) exist.

- [ ] **Step 8: Commit**

```bash
git add workflows/pangenome_profile.nf nextflow.config conf/ucr_hpcc_slurm.config
git commit -m "$(cat <<'EOF'
pangenome: wire island+Pfam enrichment step into pangenome_profile.nf

Flag-gated behind the new --pangenome_island_pfam_hmm (distinct from the
pre-existing --pangenome_pfam_hmm, which stays load-bearing for the
captain-by-name branch only). Named marker searches
(--pangenome_marker_names/--pangenome_marker_hmm_paths) run via the new
MARKER_HMMSEARCH process, invoked once over a Channel.fromList fan-out --
not a loop, which DSL2 forbids. Verified via a real 10-strain execution
(not just --help) in both the zero-marker and one-marker code paths.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 9: Update `pangenome.nf`'s help text and README/DESIGN docs

**Files:**
- Modify: `pangenome.nf` (the `print_help()` function)
- Modify: `README.md` (if it documents `pangenome.nf`'s params — check first)

**Interfaces:** None — documentation only.

- [ ] **Step 1: Check what `pangenome.nf`'s `print_help()` currently documents**

```bash
grep -n "print_help\|Optional arguments" pangenome.nf
```

- [ ] **Step 2: Add the new params to the help text**, in the same
  `Optional arguments:` block as `--pangenome_captain_hmm`:

```groovy
  --pangenome_island_pfam_hmm      Path to Pfam-A.hmm -- enables the accessory-island
                                    + Pfam functional-enrichment step (off by default).
                                    Distinct from --pangenome_pfam_hmm, which is only
                                    used by the captain-by-name branch above.
  --pangenome_island_min_size      Minimum island size to report (default: 2).
  --pangenome_pfam_domain_evalue   hmmscan domain-level E-value cutoff (default: 1e-3).
  --pangenome_marker_names         Comma list of named marker searches (e.g.
                                    'captain,sm_backbone'), run via MARKER_HMMSEARCH.
  --pangenome_marker_hmm_paths     Parallel comma list of HMM paths for each named marker.
  --pangenome_marker_evalue        hmmsearch E-value cutoff for marker searches (default: 1e-5).
  --pangenome_accumulation_permutations  Random strain-order permutations for the
                                    rarefaction/accumulation curve (default: 20).
  --pangenome_accumulation_seed    RNG seed for the accumulation curve (default: 0).
```

- [ ] **Step 3: Check and update `README.md` if it lists `pangenome.nf` params**

```bash
grep -n "pangenome_captain_hmm\|pangenome_dereplicate" README.md
```

If found, add the same params documented in Step 2 nearby, matching that
section's existing format.

- [ ] **Step 4: Commit**

```bash
git add pangenome.nf README.md
git commit -m "$(cat <<'EOF'
pangenome: document the new island+Pfam enrichment params

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01JocL68YDR4Rpfb5bqnrJ87
EOF
)"
```

---

### Task 10: Full test suite regression check

**Files:** None new — verification only.

- [ ] **Step 1: Run the complete pytest suite**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
pixi run pytest tests/ -v 2>&1 | tail -60
```

Expected: all tests pass (the pre-existing suite plus this plan's ~20 new tests
across Tasks 1, 2, 3, 5, 6), zero failures, zero new errors.

- [ ] **Step 2: If any test fails, diagnose before fixing** — per this repo's
  own convention (see `.claude` guidance elsewhere in this session), check the
  actual failure output/traceback first; do not guess at a fix.

- [ ] **Step 3: No commit needed for this task** (verification only) — but if
  Step 2 required a fix, commit that fix with its own message before moving on.

---

## Explicitly out of scope for this plan

- Running this new step against the real Coccidioides pangenome data (a
  separate execution step, once this plan's tasks are all committed and the
  test suite is green — needs an actual `Pfam-A.hmm` path available and a
  decision on whether Coccidioides has any named marker HMMs worth searching
  for, per the spec's "out of scope" section).
- Publishing a Claude Artifact from any report output (an agent-side action
  per-study, not a pipeline capability — see the spec).
- The MULE/DDE-transposase marker idea (`nf_NovInvenio/todo/
  mule-dde-transposase-island-marker.md`) — no specific Pfam accessions are
  chosen in this plan.
