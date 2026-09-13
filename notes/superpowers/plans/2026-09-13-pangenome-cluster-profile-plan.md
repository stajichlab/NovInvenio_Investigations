# Pangenome Cluster-Profile Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone pangenome analysis for `studies/fungi/Afumigatus_pangenome` (293 *A. fumigatus* strains) that bins gene families into core/soft-core/shell/cloud/singleton by presence frequency, detects statistically supported co-occurring gain/loss between accessory families, tests whether co-varying families are physically clustered in each strain's own genome (Starship/giant-transposon mobility signature), and validates all of this against a battery of known-Starship positive/negative controls from Gluck-Thaler et al. 2025 (mBio) before running the same method on unpublished candidate clusters — including a first concrete target, the HAC (hrmA-Associated Cluster) gene family and the independent `hacA` locus.

**Architecture:** A chain of small, single-purpose Python scripts under `studies/fungi/Afumigatus_pangenome/bin/`, sharing one data structure (`PresenceMatrix`, in a new study-local `lib/pangenome_matrix.py`) and NovInvenio's existing `config.csv`/GFF3 conventions. Each script is runnable standalone on the real 293-strain dataset and unit-tested on small synthetic fixtures (no NovInvenio-style end-to-end pytest gate — the benchmark suite itself, run against real published Starship data, is the biological validation, per the spec's "Testing / validation plan" section).

**Tech Stack:** Python 3.12 (NII's existing pixi environment) + pandas, biopython, scipy (already NII dependencies) + mmseqs2, diamond, blast+ (tblastn), mash (new pixi dependencies for this study — see Task 1; all invoked as bare subprocess commands, not pixi-specific paths, so they keep working unchanged if this later moves into the existing `ghcr.io/stajichlab/novinvenio` container, which already ships mmseqs2/diamond/blast, or a future dedicated container — only `mash` would need adding to that image).

**Spec:** `notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md` (see especially its "Review disposition" section — every design choice below follows directly from it, including the bioinformatics-review revisions).

## Global Constraints

- Scripts live under `studies/fungi/Afumigatus_pangenome/bin/` (study-specific per NII's own convention), a shared module under `studies/fungi/Afumigatus_pangenome/lib/`, tests under `studies/fungi/Afumigatus_pangenome/tests/`.
- Tier-1 clustering (the unit used for presence/frequency/co-occurrence/synteny): mmseqs2 `--min-seq-id 0.9 -c 0.8 --cov-mode 0 --cluster-reassign`; diamond equivalent `diamond cluster --approx-id 90 --member-cover 80`.
- Tier-2 clustering (superfamily label only, never used for presence/frequency): re-cluster tier-1 representatives at 30-50% identity.
- Presence calls are three-state: `present` (protein-model hit), `genome_only` (tblastn rescue hit, no protein model), `absent`. Both `present` and `genome_only` count as "this strain carries the family" for frequency/co-occurrence/synteny purposes (`PresenceMatrix.is_present()`); only `absent` means truly missing.
- Every external-tool subprocess call uses the bare command name (`mmseqs`, `diamond`, `tblastn`, `mash`) and never a hardcoded environment-specific path.
- Co-occurrence pairs require a minimum frequency floor (>=5 dereplicated strains carrying each family) before testing, use Fisher's-exact p-values with Benjamini-Hochberg FDR correction (not a raw similarity threshold), and are polarized to gain/loss using the two outgroup strains (`Aslen_ref`, `Neofi_ref`) rather than labeled "co-loss"/"co-gain" from correlation alone.
- All new Python modules use `from __future__ import annotations` and `sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))` for local imports, matching NovInvenio's existing `bin/` convention.

---

### Task 1: Environment setup — add clustering/rescue/dereplication tools

**Files:**
- Modify: `pixi.toml`

**Interfaces:**
- Consumes: nothing (foundational task)
- Produces: `mmseqs`, `diamond`, `tblastn`/`blastp` (from the `blast` package), `mash` all resolve on `PATH` inside `pixi run`

- [ ] **Step 1: Add dependencies to `pixi.toml`**

Add to the `[dependencies]` table (matching NovInvenio's own pinned versions where this study's dependency overlaps, so results are comparable if code is later ported):

```toml
mmseqs2 = ">=18.8cc5c,<19"
diamond = ">=2.2.0,<3"
blast = ">=2.17.0,<3"
mash = ">=2.3,<3"
```

- [ ] **Step 2: Install and verify**

Run:
```bash
pixi install
pixi run mmseqs version
pixi run diamond --version
pixi run tblastn -version
pixi run mash --version
```
Expected: each prints a version string with no error.

- [ ] **Step 3: Commit**

```bash
git add pixi.toml pixi.lock
git commit -m "Add mmseqs2/diamond/blast/mash for the pangenome cluster-profile analysis"
```

---

### Task 2: Shared `PresenceMatrix` data structure

**Files:**
- Create: `studies/fungi/Afumigatus_pangenome/lib/__init__.py` (empty)
- Create: `studies/fungi/Afumigatus_pangenome/lib/pangenome_matrix.py`
- Test: `studies/fungi/Afumigatus_pangenome/tests/test_pangenome_matrix.py`

**Interfaces:**
- Consumes: nothing
- Produces: `PRESENT`, `GENOME_ONLY`, `ABSENT` constants; `read_cluster_tsv(path) -> dict[str,str]`; `build_families(member_to_rep) -> dict[str,list[str]]`; `PresenceMatrix` class with `families: list[str]`, `strains: list[str]`, `set_call(family, strain, state, copies=0)`, `call(family, strain) -> str`, `is_present(family, strain) -> bool`, `presence_vector(family) -> list[bool]`, `frequency(family) -> float`, `strain_count(family) -> int`, `to_tsv(path)`, `PresenceMatrix.from_tsv(path)` classmethod. Every later task depends on this exact interface.

- [ ] **Step 1: Write the failing tests**

```python
# studies/fungi/Afumigatus_pangenome/tests/test_pangenome_matrix.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import (
    PRESENT, GENOME_ONLY, ABSENT,
    read_cluster_tsv, build_families, PresenceMatrix,
)


def test_read_cluster_tsv_parses_rep_member_pairs(tmp_path):
    p = tmp_path / "clusters.tsv"
    p.write_text("repA\trepA\nrepA\tmemberA2\nrepB\trepB\n")
    result = read_cluster_tsv(p)
    assert result == {"repA": "repA", "memberA2": "repA", "repB": "repB"}


def test_build_families_groups_by_rep_including_singletons():
    member_to_rep = {"repA": "repA", "memberA2": "repA", "repB": "repB"}
    families = build_families(member_to_rep)
    assert families == {"repA": ["memberA2", "repA"], "repB": ["repB"]}


def test_presence_matrix_set_and_query():
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2", "s3"])
    pm.set_call("famA", "s1", PRESENT)
    pm.set_call("famA", "s2", GENOME_ONLY)
    # famA absent in s3 (never set) -- default ABSENT
    assert pm.call("famA", "s1") == PRESENT
    assert pm.call("famA", "s3") == ABSENT
    assert pm.is_present("famA", "s1") is True
    assert pm.is_present("famA", "s2") is True   # genome_only still counts as present
    assert pm.is_present("famA", "s3") is False
    assert pm.presence_vector("famA") == [True, True, False]
    assert pm.strain_count("famA") == 2
    assert pm.frequency("famA") == 2 / 3


def test_presence_matrix_tsv_roundtrip(tmp_path):
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", PRESENT)
    pm.set_call("famB", "s2", GENOME_ONLY)
    out = tmp_path / "matrix.tsv"
    pm.to_tsv(out)
    loaded = PresenceMatrix.from_tsv(out)
    assert loaded.strains == ["s1", "s2"]
    assert loaded.families == ["famA", "famB"]
    assert loaded.call("famA", "s1") == PRESENT
    assert loaded.call("famB", "s2") == GENOME_ONLY
    assert loaded.call("famA", "s2") == ABSENT


def test_set_call_rejects_invalid_state():
    pm = PresenceMatrix(families=["famA"], strains=["s1"])
    try:
        pm.set_call("famA", "s1", "not_a_real_state")
        assert False, "expected ValueError"
    except ValueError:
        pass
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_pangenome_matrix.py -v`
Expected: FAIL/ERROR — `pangenome_matrix` module does not exist yet.

- [ ] **Step 3: Write the implementation**

```python
# studies/fungi/Afumigatus_pangenome/lib/__init__.py
```
(empty file)

```python
# studies/fungi/Afumigatus_pangenome/lib/pangenome_matrix.py
"""Shared data structures for the pangenome cluster-profile analysis
(notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md).

PresenceMatrix holds one row per gene family (a tier-1 mmseqs/diamond
cluster representative) and one column per strain, with a three-state
call per cell -- see the Global Constraints in this plan's header for
what each state means.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

PRESENT = "present"
GENOME_ONLY = "genome_only"
ABSENT = "absent"
STATES = (PRESENT, GENOME_ONLY, ABSENT)


def read_cluster_tsv(path: str | Path) -> dict[str, str]:
    """Parse an mmseqs/diamond cluster TSV (rep\\tmember per line) into
    {member_id: rep_id}."""
    member_to_rep: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            member_to_rep[member] = rep
    return member_to_rep


def build_families(member_to_rep: dict[str, str]) -> dict[str, list[str]]:
    """Return {rep_id: sorted [member_ids]} for every cluster, including
    singletons -- unlike NovInvenio's lib/clusters.py::build_families,
    which drops singleton clusters (irrelevant there; here a singleton
    family is a real, reportable frequency bin)."""
    families: dict[str, list[str]] = {}
    for member, rep in member_to_rep.items():
        families.setdefault(rep, []).append(member)
    return {rep: sorted(members) for rep, members in families.items()}


@dataclass
class PresenceMatrix:
    families: list[str]
    strains: list[str]
    calls: dict[tuple[str, str], str] = field(default_factory=dict)
    copy_number: dict[tuple[str, str], int] = field(default_factory=dict)

    def set_call(self, family: str, strain: str, state: str, copies: int = 0) -> None:
        if state not in STATES:
            raise ValueError(f"invalid state {state!r}, must be one of {STATES}")
        self.calls[(family, strain)] = state
        if copies:
            self.copy_number[(family, strain)] = copies

    def call(self, family: str, strain: str) -> str:
        return self.calls.get((family, strain), ABSENT)

    def is_present(self, family: str, strain: str) -> bool:
        return self.call(family, strain) in (PRESENT, GENOME_ONLY)

    def presence_vector(self, family: str) -> list[bool]:
        return [self.is_present(family, s) for s in self.strains]

    def strain_count(self, family: str) -> int:
        return sum(self.presence_vector(family))

    def frequency(self, family: str) -> float:
        if not self.strains:
            return 0.0
        return self.strain_count(family) / len(self.strains)

    def to_tsv(self, path: str | Path) -> None:
        with open(path, "w") as fh:
            fh.write("family\t" + "\t".join(self.strains) + "\n")
            for fam in self.families:
                row = [self.call(fam, s) for s in self.strains]
                fh.write(fam + "\t" + "\t".join(row) + "\n")

    @classmethod
    def from_tsv(cls, path: str | Path) -> "PresenceMatrix":
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            strains = header[1:]
            families: list[str] = []
            calls: dict[tuple[str, str], str] = {}
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                fam = parts[0]
                families.append(fam)
                for strain, state in zip(strains, parts[1:]):
                    calls[(fam, strain)] = state
        pm = cls(families=families, strains=strains)
        pm.calls = calls
        return pm
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_pangenome_matrix.py -v`
Expected: PASS (5 tests)

- [ ] **Step 5: Commit**

```bash
git add studies/fungi/Afumigatus_pangenome/lib/__init__.py \
        studies/fungi/Afumigatus_pangenome/lib/pangenome_matrix.py \
        studies/fungi/Afumigatus_pangenome/tests/test_pangenome_matrix.py
git commit -m "Add shared PresenceMatrix data structure for pangenome cluster-profile analysis"
```

---

### Task 3: Strain inventory, assembly-quality proxy, and dereplication

**Files:**
- Create: `studies/fungi/Afumigatus_pangenome/bin/dereplicate_strains.py`
- Test: `studies/fungi/Afumigatus_pangenome/tests/test_dereplicate_strains.py`

**Interfaces:**
- Consumes: nothing new (reads `config.csv`, strain DNA FASTAs from `data_dir/dna/`)
- Produces: `compute_assembly_stats(fasta_path) -> dict` (keys `n_contigs`, `n50`, `total_length`); `parse_mash_dist(lines: list[str], threshold: float) -> list[set[str]]`; `choose_representatives(dedup_groups: list[set[str]], assembly_stats: dict[str, dict]) -> dict[str, str]` (member Short -> representative Short, picking the highest-N50 member per group); output file `strain_inventory.tsv` (columns: `Short, n_contigs, n50, total_length, dedup_group, is_representative`) consumed by every later task that needs "which strains actually count."

This directly addresses two open controls from the review: assembly completeness/fragmentation (via N50/contig count, a directly measurable proxy — the review's "draft vs. long-read" concern is really about fragmentation risk, which N50 measures without needing to guess sequencing technology from metadata) and strain dereplication (via mash).

- [ ] **Step 1: Write the failing tests**

```python
# studies/fungi/Afumigatus_pangenome/tests/test_dereplicate_strains.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from dereplicate_strains import (
    compute_assembly_stats, parse_mash_dist, choose_representatives,
)


def test_compute_assembly_stats_n50_and_contig_count(tmp_path):
    fa = tmp_path / "strain.dna.fa"
    # two contigs: 100bp and 300bp -> total 400, N50 is the 300bp contig
    fa.write_text(">contig1\n" + "A" * 100 + "\n>contig2\n" + "C" * 300 + "\n")
    stats = compute_assembly_stats(fa)
    assert stats["n_contigs"] == 2
    assert stats["total_length"] == 400
    assert stats["n50"] == 300


def test_parse_mash_dist_groups_below_threshold():
    # mash dist -t output: query, ref1_dist, ref2_dist, ... one row per query
    # here: s1~s2 near-identical (0.0005), s3 distinct (0.05)
    lines = [
        "#query\ts1\ts2\ts3",
        "s1\t0.0000\t0.0005\t0.0500",
        "s2\t0.0005\t0.0000\t0.0520",
        "s3\t0.0500\t0.0520\t0.0000",
    ]
    groups = parse_mash_dist(lines, threshold=0.001)
    assert {"s1", "s2"} in groups
    assert {"s3"} in groups


def test_choose_representatives_picks_highest_n50():
    groups = [{"s1", "s2"}, {"s3"}]
    assembly_stats = {
        "s1": {"n50": 500_000, "n_contigs": 20, "total_length": 29_000_000},
        "s2": {"n50": 2_000_000, "n_contigs": 8, "total_length": 29_100_000},
        "s3": {"n50": 1_000_000, "n_contigs": 15, "total_length": 28_900_000},
    }
    reps = choose_representatives(groups, assembly_stats)
    assert reps == {"s1": "s2", "s2": "s2", "s3": "s3"}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_dereplicate_strains.py -v`
Expected: FAIL — `dereplicate_strains` module does not exist yet.

- [ ] **Step 3: Write the implementation**

```python
# studies/fungi/Afumigatus_pangenome/bin/dereplicate_strains.py
#!/usr/bin/env python3
"""Strain inventory: per-strain assembly-quality proxy (N50/contig count)
and mash-based dereplication, feeding every later step of the pangenome
cluster-profile analysis (see notes/superpowers/specs/
2026-09-13-pangenome-cluster-profile-design.md, component 1b).

Usage:
  dereplicate_strains.py --config config.csv --data_dir data_dir \\
      --output strain_inventory.tsv
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
NOVINVENIO_LIB = Path(__file__).resolve().parents[4] / "NovInvenio" / "lib"
if NOVINVENIO_LIB.exists():
    sys.path.insert(0, str(NOVINVENIO_LIB))
from config_parser import parse_config  # noqa: E402


def compute_assembly_stats(fasta_path: str | Path) -> dict:
    """Return {n_contigs, n50, total_length} for a DNA FASTA, without
    external tools -- a fast, dependency-free fragmentation proxy."""
    lengths = []
    current = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
        if current:
            lengths.append(current)
    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2
    running = 0
    n50 = 0
    for length in lengths:
        running += length
        if running >= half:
            n50 = length
            break
    return {"n_contigs": len(lengths), "n50": n50, "total_length": total}


def parse_mash_dist(lines: list[str], threshold: float) -> list[set[str]]:
    """Parse `mash dist -t` output (a query x reference distance matrix,
    tab-separated, first row is '#query\\tref1\\tref2\\t...') into
    dereplication groups: strains whose pairwise mash distance is below
    `threshold` are grouped together (union-find over the threshold graph)."""
    header = lines[0].lstrip("#").split("\t")
    names = header[1:]
    parent = {name: name for name in names}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for row in lines[1:]:
        parts = row.split("\t")
        query = parts[0]
        for ref, dist_str in zip(names, parts[1:]):
            if query == ref:
                continue
            if float(dist_str) < threshold:
                union(query, ref)

    groups: dict[str, set[str]] = {}
    for name in names:
        root = find(name)
        groups.setdefault(root, set()).add(name)
    return list(groups.values())


def choose_representatives(
    dedup_groups: list[set[str]], assembly_stats: dict[str, dict]
) -> dict[str, str]:
    """For every strain in every group, map it to the group's chosen
    representative (the member with the highest N50 -- the least
    fragmented assembly, per the review's fragmentation-risk concern)."""
    result: dict[str, str] = {}
    for group in dedup_groups:
        rep = max(group, key=lambda s: assembly_stats[s]["n50"])
        for member in group:
            result[member] = rep
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--mash_threshold", type=float, default=0.001)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    samples = parse_config(args.config)
    data_dir = Path(args.data_dir)
    dna_paths = {s.short: data_dir / "dna" / s.dna for s in samples if s.dna}

    stats = {short: compute_assembly_stats(p) for short, p in dna_paths.items()}

    sketch_prefix = "strain_sketches"
    subprocess.run(
        ["mash", "sketch", "-o", sketch_prefix] + [str(p) for p in dna_paths.values()],
        check=True,
    )
    dist_out = subprocess.run(
        ["mash", "dist", "-t", f"{sketch_prefix}.msh", f"{sketch_prefix}.msh"],
        check=True, capture_output=True, text=True,
    ).stdout
    groups = parse_mash_dist(dist_out.splitlines(), args.mash_threshold)
    reps = choose_representatives(groups, stats)
    dedup_group_id = {}
    for i, group in enumerate(groups):
        for member in group:
            dedup_group_id[member] = i

    with open(args.output, "w") as fh:
        fh.write("Short\tn_contigs\tn50\ttotal_length\tdedup_group\tis_representative\n")
        for short in sorted(dna_paths):
            s = stats[short]
            fh.write(
                f"{short}\t{s['n_contigs']}\t{s['n50']}\t{s['total_length']}\t"
                f"{dedup_group_id[short]}\t{int(reps[short] == short)}\n"
            )


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_dereplicate_strains.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Run against the real 293-strain dataset**

```bash
cd studies/fungi/Afumigatus_pangenome
chmod +x bin/dereplicate_strains.py
pixi run bin/dereplicate_strains.py --config config.csv --data_dir data_dir \
    --output strain_inventory.tsv
```
Inspect `strain_inventory.tsv`: confirm every one of the 293 ingroup strains appears exactly once, and manually spot-check a couple of dedup groups (if any) against strain names in `config.csv` for plausibility (e.g. two rows with obviously similar strain-name spellings).

- [ ] **Step 6: Commit**

```bash
git add studies/fungi/Afumigatus_pangenome/bin/dereplicate_strains.py \
        studies/fungi/Afumigatus_pangenome/tests/test_dereplicate_strains.py
git commit -m "Add strain inventory: assembly-quality proxy + mash dereplication"
```

---

### Task 4: Two-tier clustering backend (mmseqs2 + diamond)

**Files:**
- Create: `studies/fungi/Afumigatus_pangenome/bin/cluster_backend.py`
- Test: `studies/fungi/Afumigatus_pangenome/tests/test_cluster_backend.py`

**Interfaces:**
- Consumes: `read_cluster_tsv`, `build_families` from Task 2's `lib/pangenome_matrix.py`
- Produces: `two_tier_families(tier1_cluster_tsv, tier2_cluster_tsv) -> dict[str, dict]` mapping `tier1_rep -> {"members": [...], "superfamily": tier2_rep_of_tier1_rep}`; CLI subcommands `mmseqs-tier1`, `mmseqs-tier2`, `diamond-tier1`, `diamond-tier2` that shell out to the real tools. Later tasks (5-8, 10, 11) consume the tier-1 family membership; the benchmark scorecard (Task 10) consumes `superfamily` to test whether known-related paralogs (e.g. the HAC family) land under one superfamily label without collapsing their individual tier-1 presence calls.

- [ ] **Step 1: Write the failing test**

```python
# studies/fungi/Afumigatus_pangenome/tests/test_cluster_backend.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from cluster_backend import two_tier_families


def test_two_tier_families_attaches_superfamily_label(tmp_path):
    tier1 = tmp_path / "tier1_cluster.tsv"
    # two tier-1 families: repA (with member memberA2) and repB (singleton)
    tier1.write_text("repA\trepA\nrepA\tmemberA2\nrepB\trepB\n")
    tier2 = tmp_path / "tier2_cluster.tsv"
    # tier-2 re-clusters the tier-1 reps: repA and repB fall under superfamily repA
    tier2.write_text("repA\trepA\nrepA\trepB\n")

    result = two_tier_families(tier1, tier2)

    assert result["repA"]["members"] == ["memberA2", "repA"]
    assert result["repA"]["superfamily"] == "repA"
    assert result["repB"]["members"] == ["repB"]
    assert result["repB"]["superfamily"] == "repA"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_cluster_backend.py -v`
Expected: FAIL — `cluster_backend` module does not exist yet.

- [ ] **Step 3: Write the implementation**

```python
# studies/fungi/Afumigatus_pangenome/bin/cluster_backend.py
#!/usr/bin/env python3
"""Two-tier, swappable (mmseqs2 | diamond) clustering backend for the
pangenome cluster-profile analysis (notes/superpowers/specs/
2026-09-13-pangenome-cluster-profile-design.md, component 1). Tier 1 is
the allele/ortholog unit used by every downstream step (presence,
frequency, co-occurrence, synteny); tier 2 re-clusters tier-1
representatives loosely, purely as a superfamily annotation label -- it
is never used for presence/frequency calls.

Usage:
  cluster_backend.py mmseqs-tier1 --fasta all_ingroup.fa --out_prefix tier1
  cluster_backend.py mmseqs-tier2 --fasta tier1_rep_seq.fasta --out_prefix tier2
  cluster_backend.py diamond-tier1 --fasta all_ingroup.fa --out_prefix tier1
  cluster_backend.py diamond-tier2 --fasta tier1_rep_seq.fasta --out_prefix tier2
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import build_families, read_cluster_tsv  # noqa: E402


def run_mmseqs_cluster(
    fasta: Path, out_prefix: str, min_seq_id: float, cov: float, cluster_reassign: bool
) -> Path:
    cmd = [
        "mmseqs", "easy-cluster", str(fasta), out_prefix, "tmp_mmseqs",
        "--min-seq-id", str(min_seq_id), "-c", str(cov), "--cov-mode", "0",
    ]
    if cluster_reassign:
        cmd.append("--cluster-reassign")
    subprocess.run(cmd, check=True)
    return Path(f"{out_prefix}_cluster.tsv")


def run_diamond_cluster(fasta: Path, out_prefix: str, approx_id: float, member_cover: float) -> Path:
    cmd = [
        "diamond", "cluster", "-d", str(fasta), "-o", f"{out_prefix}_cluster.tsv",
        "--approx-id", str(approx_id), "--member-cover", str(member_cover),
    ]
    subprocess.run(cmd, check=True)
    return Path(f"{out_prefix}_cluster.tsv")


def two_tier_families(tier1_cluster_tsv: str | Path, tier2_cluster_tsv: str | Path) -> dict:
    """Combine tier-1 family membership with the tier-2 superfamily label
    for each tier-1 representative. Returns
    {tier1_rep: {"members": [...], "superfamily": tier2_rep_of_tier1_rep}}."""
    tier1_families = build_families(read_cluster_tsv(tier1_cluster_tsv))
    tier1_rep_to_tier2_rep = read_cluster_tsv(tier2_cluster_tsv)
    return {
        rep: {
            "members": members,
            "superfamily": tier1_rep_to_tier2_rep.get(rep, rep),
        }
        for rep, members in tier1_families.items()
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)
    for name in ("mmseqs-tier1", "mmseqs-tier2", "diamond-tier1", "diamond-tier2"):
        p = sub.add_parser(name)
        p.add_argument("--fasta", required=True, type=Path)
        p.add_argument("--out_prefix", required=True)
    args = ap.parse_args()

    if args.command == "mmseqs-tier1":
        run_mmseqs_cluster(args.fasta, args.out_prefix, min_seq_id=0.9, cov=0.8, cluster_reassign=True)
    elif args.command == "mmseqs-tier2":
        run_mmseqs_cluster(args.fasta, args.out_prefix, min_seq_id=0.4, cov=0.8, cluster_reassign=False)
    elif args.command == "diamond-tier1":
        run_diamond_cluster(args.fasta, args.out_prefix, approx_id=90, member_cover=80)
    elif args.command == "diamond-tier2":
        run_diamond_cluster(args.fasta, args.out_prefix, approx_id=40, member_cover=80)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_cluster_backend.py -v`
Expected: PASS

- [ ] **Step 5: Run against the real dataset, both backends**

```bash
cd studies/fungi/Afumigatus_pangenome
chmod +x bin/cluster_backend.py
# build one FASTA of all ingroup proteins first (reuse NovInvenio's extract_candidates.py
# pattern, or a simple cat of data_dir/pep/*.pep.fa with Short-prefixed headers)
pixi run bin/cluster_backend.py mmseqs-tier1 --fasta all_ingroup.fa --out_prefix mmseqs_tier1
pixi run bin/cluster_backend.py mmseqs-tier2 --fasta mmseqs_tier1_rep_seq.fasta --out_prefix mmseqs_tier2
pixi run bin/cluster_backend.py diamond-tier1 --fasta all_ingroup.fa --out_prefix diamond_tier1
pixi run bin/cluster_backend.py diamond-tier2 --fasta diamond_tier1_rep_seq.fasta --out_prefix diamond_tier2
```
Confirm both produce non-empty `*_cluster.tsv` files, and that `two_tier_families()` loads each pair without error.

- [ ] **Step 6: Commit**

```bash
git add studies/fungi/Afumigatus_pangenome/bin/cluster_backend.py \
        studies/fungi/Afumigatus_pangenome/tests/test_cluster_backend.py
git commit -m "Add two-tier mmseqs2/diamond clustering backend"
```

---

### Task 5: Genome-level rescue pass (tblastn) + isoform collapse wiring

**Files:**
- Create: `studies/fungi/Afumigatus_pangenome/bin/rescue_pass.py`
- Test: `studies/fungi/Afumigatus_pangenome/tests/test_rescue_pass.py`

**Interfaces:**
- Consumes: `PresenceMatrix`, `PRESENT`, `GENOME_ONLY`, `ABSENT` from Task 2
- Produces: `parse_tblastn_hits(tblastn_tsv_lines, min_pident=90.0, min_qcov=80.0) -> set[tuple[str,str]]` (family_rep, strain) pairs with a qualifying genomic hit; `apply_rescue(matrix, rescue_hits) -> None` (mutates matrix in place, only upgrading `ABSENT` calls — never downgrades an existing `PRESENT`). Task 6 (frequency binning) and everything downstream consumes the matrix after this rescue pass is applied, not before.

- [ ] **Step 1: Write the failing tests**

```python
# studies/fungi/Afumigatus_pangenome/tests/test_rescue_pass.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY, ABSENT
from rescue_pass import parse_tblastn_hits, apply_rescue


def test_parse_tblastn_hits_filters_by_identity_and_coverage():
    # outfmt6 with qcovs appended: qseqid sseqid pident length ... qcovs
    # family "famA" query, subject header encodes strain as "s2|contig1"
    lines = [
        "famA\ts2|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",
        "famA\ts3|contig1\t60.0\t100\t0\t0\t1\t100\t500\t600\t1e-10\t80\t95",  # too low identity
        "famA\ts4|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t50",  # too low coverage
    ]
    hits = parse_tblastn_hits(lines, min_pident=90.0, min_qcov=80.0)
    assert hits == {("famA", "s2")}


def test_apply_rescue_only_upgrades_absent_calls():
    pm = PresenceMatrix(families=["famA"], strains=["s1", "s2", "s3"])
    pm.set_call("famA", "s1", PRESENT)
    # s2, s3 default to ABSENT
    apply_rescue(pm, {("famA", "s1"), ("famA", "s2")})
    assert pm.call("famA", "s1") == PRESENT       # unchanged, was already PRESENT
    assert pm.call("famA", "s2") == GENOME_ONLY   # upgraded from ABSENT
    assert pm.call("famA", "s3") == ABSENT         # no rescue hit, stays ABSENT
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_rescue_pass.py -v`
Expected: FAIL — `rescue_pass` module does not exist yet.

- [ ] **Step 3: Write the implementation**

```python
# studies/fungi/Afumigatus_pangenome/bin/rescue_pass.py
#!/usr/bin/env python3
"""Genome-level rescue pass: tblastn every tier-1 family representative
against each strain's own genome to catch coverage-failed protein-model
absences (fragmented/split gene models, draft-assembly artifacts) --
notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md,
component 1b, Fable review finding 2.

Usage:
  rescue_pass.py --matrix presence_matrix.tsv --tblastn_tsv all_vs_all.tblastn.tsv \\
      --output presence_matrix.rescued.tsv

Upstream isoform collapse: before building the protein FASTA fed into
Task 4's clustering, run NovInvenio's existing bin/collapse_isoforms.py
per strain (see NovInvenio's CLAUDE.md "collapse_isoforms.py" entry) so
per-strain copy counts reflect genes, not alternative transcripts:

  NovInvenio/bin/collapse_isoforms.py \\
      --protein-fasta data_dir/pep/<Short>.pep.fa \\
      --feature-table <Short>_feature_table.txt.gz \\
      --output data_dir/pep_collapsed/<Short>.pep.fa
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import ABSENT, GENOME_ONLY, PresenceMatrix  # noqa: E402


def parse_tblastn_hits(
    lines: list[str], min_pident: float = 90.0, min_qcov: float = 80.0
) -> set[tuple[str, str]]:
    """Parse tblastn outfmt6 lines (with a trailing qcovs column) into the
    set of (family_rep, strain) pairs with a qualifying genomic hit.
    Subject IDs are expected as '<Short>|<contig>' (set via -subject_besthit
    or a renamed genome DB, so the strain is recoverable from the hit)."""
    hits: set[tuple[str, str]] = set()
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        family, subject, pident, qcovs = parts[0], parts[1], float(parts[2]), float(parts[-1])
        if pident >= min_pident and qcovs >= min_qcov:
            strain = subject.split("|", 1)[0]
            hits.add((family, strain))
    return hits


def apply_rescue(matrix: PresenceMatrix, rescue_hits: set[tuple[str, str]]) -> None:
    """Upgrade ABSENT calls to GENOME_ONLY wherever a qualifying rescue
    hit exists. Never downgrades an existing PRESENT call."""
    for family, strain in rescue_hits:
        if matrix.call(family, strain) == ABSENT:
            matrix.set_call(family, strain, GENOME_ONLY)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--tblastn_tsv", required=True)
    ap.add_argument("--min_pident", type=float, default=90.0)
    ap.add_argument("--min_qcov", type=float, default=80.0)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)
    with open(args.tblastn_tsv) as fh:
        hits = parse_tblastn_hits(fh.readlines(), args.min_pident, args.min_qcov)
    apply_rescue(matrix, hits)
    matrix.to_tsv(args.output)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_rescue_pass.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Run against the real dataset**

Build one BLAST-formatted genome DB per strain with subject headers prefixed `<Short>|`, tblastn every tier-1 representative against every strain's genome DB with `-outfmt "6 std qcovs"`, concatenate results, then:

```bash
pixi run bin/rescue_pass.py --matrix presence_matrix.tsv \
    --tblastn_tsv all_vs_all.tblastn.tsv --output presence_matrix.rescued.tsv
```
Spot-check: pick 2-3 families that were `ABSENT` in a strain before rescue and confirm they're now `GENOME_ONLY` only where a real qualifying hit exists (inspect the raw tblastn line for that pair).

- [ ] **Step 6: Commit**

```bash
git add studies/fungi/Afumigatus_pangenome/bin/rescue_pass.py \
        studies/fungi/Afumigatus_pangenome/tests/test_rescue_pass.py
git commit -m "Add genome-level tblastn rescue pass for coverage-failed presence calls"
```

---

### Task 6: Frequency binning

**Files:**
- Create: `studies/fungi/Afumigatus_pangenome/bin/frequency_bins.py`
- Test: `studies/fungi/Afumigatus_pangenome/tests/test_frequency_bins.py`

**Interfaces:**
- Consumes: `PresenceMatrix` from Task 2
- Produces: `assign_bin(freq, strain_count, core_cutoff=0.95, softcore_cutoff=0.90, shell_cutoff=0.15) -> str` (one of `"core"`, `"soft_core"`, `"shell"`, `"cloud"`, `"singleton"`); `compute_frequency_table(matrix) -> list[dict]` (one row per family: `family, frequency, strain_count, bin`). Components 4/5 (co-occurrence, synteny) consume this to select the `shell`+`cloud` subset and exclude `singleton`.

- [ ] **Step 1: Write the failing tests**

```python
# studies/fungi/Afumigatus_pangenome/tests/test_frequency_bins.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT
from frequency_bins import assign_bin, compute_frequency_table


def test_assign_bin_boundaries():
    assert assign_bin(freq=1.0, strain_count=100) == "core"
    assert assign_bin(freq=0.95, strain_count=95) == "core"
    assert assign_bin(freq=0.92, strain_count=92) == "soft_core"
    assert assign_bin(freq=0.50, strain_count=50) == "shell"
    assert assign_bin(freq=0.05, strain_count=5) == "cloud"
    assert assign_bin(freq=0.01, strain_count=1) == "singleton"


def test_compute_frequency_table():
    pm = PresenceMatrix(families=["famCore", "famSingleton"], strains=["s1", "s2", "s3", "s4"])
    for s in pm.strains:
        pm.set_call("famCore", s, PRESENT)
    pm.set_call("famSingleton", "s1", PRESENT)

    table = compute_frequency_table(pm)
    by_family = {row["family"]: row for row in table}
    assert by_family["famCore"]["bin"] == "core"
    assert by_family["famCore"]["strain_count"] == 4
    assert by_family["famSingleton"]["bin"] == "singleton"
    assert by_family["famSingleton"]["strain_count"] == 1
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_frequency_bins.py -v`
Expected: FAIL — `frequency_bins` module does not exist yet.

- [ ] **Step 3: Write the implementation**

```python
# studies/fungi/Afumigatus_pangenome/bin/frequency_bins.py
#!/usr/bin/env python3
"""Core/soft-core/shell/cloud/singleton frequency binning -- notes/
superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md,
component 3.

Usage:
  frequency_bins.py --matrix presence_matrix.rescued.tsv --output frequency_table.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix  # noqa: E402


def assign_bin(
    freq: float,
    strain_count: int,
    core_cutoff: float = 0.95,
    softcore_cutoff: float = 0.90,
    shell_cutoff: float = 0.15,
) -> str:
    if strain_count <= 1:
        return "singleton"
    if freq >= core_cutoff:
        return "core"
    if freq >= softcore_cutoff:
        return "soft_core"
    if freq >= shell_cutoff:
        return "shell"
    return "cloud"


def compute_frequency_table(
    matrix: PresenceMatrix,
    core_cutoff: float = 0.95,
    softcore_cutoff: float = 0.90,
    shell_cutoff: float = 0.15,
) -> list[dict]:
    rows = []
    for fam in matrix.families:
        freq = matrix.frequency(fam)
        count = matrix.strain_count(fam)
        rows.append({
            "family": fam,
            "frequency": freq,
            "strain_count": count,
            "bin": assign_bin(freq, count, core_cutoff, softcore_cutoff, shell_cutoff),
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--core_cutoff", type=float, default=0.95)
    ap.add_argument("--softcore_cutoff", type=float, default=0.90)
    ap.add_argument("--shell_cutoff", type=float, default=0.15)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)
    table = compute_frequency_table(
        matrix, args.core_cutoff, args.softcore_cutoff, args.shell_cutoff
    )
    with open(args.output, "w") as fh:
        fh.write("family\tfrequency\tstrain_count\tbin\n")
        for row in table:
            fh.write(f"{row['family']}\t{row['frequency']:.4f}\t{row['strain_count']}\t{row['bin']}\n")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_frequency_bins.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Run against the real dataset and inspect the histogram**

```bash
pixi run bin/frequency_bins.py --matrix presence_matrix.rescued.tsv --output frequency_table.tsv
pixi run python3 -c "
import pandas as pd
df = pd.read_csv('frequency_table.tsv', sep='\t')
print(df['bin'].value_counts())
print(df['frequency'].describe())
"
```
Per the spec, use this real histogram to sanity-check (and adjust if warranted) the default `--core_cutoff`/`--softcore_cutoff`/`--shell_cutoff` values before treating them as final — record any change and why in `PANGENOME_CLUSTER_PROFILE_NOTES.md`.

- [ ] **Step 6: Commit**

```bash
git add studies/fungi/Afumigatus_pangenome/bin/frequency_bins.py \
        studies/fungi/Afumigatus_pangenome/tests/test_frequency_bins.py
git commit -m "Add core/soft-core/shell/cloud/singleton frequency binning"
```

---

### Task 7: Co-occurrence (co-loss/co-gain) statistics

**Files:**
- Create: `studies/fungi/Afumigatus_pangenome/bin/cooccurrence.py`
- Test: `studies/fungi/Afumigatus_pangenome/tests/test_cooccurrence.py`

**Interfaces:**
- Consumes: `PresenceMatrix` from Task 2, `compute_frequency_table`/`assign_bin` from Task 6 (to restrict to `shell`+`cloud` families)
- Produces: `jaccard(a, b) -> float`; `fisher_pvalue(a, b) -> float`; `benjamini_hochberg(pvalues) -> list[float]`; `polarize_direction(present_in_outgroup, freq_in_ingroup) -> str` (one of `"gain"`, `"loss"`, `"ambiguous"`); `clade_composition(strains_present, clade_of_strain) -> dict[str,int]`; `permutation_null_pvalue(a, b, clade_of_strain, n_perms, rng) -> float`; `find_cooccurring_pairs(matrix, frequency_table, clade_of_strain, outgroup_presence, min_strain_count=5, fdr_alpha=0.05, n_perms=1000, seed=0) -> list[dict]` — the integration function producing the final reportable pairs. Task 10 (benchmark scorecard) calls `find_cooccurring_pairs` to check recovery of/false positives against known pairs.

- [ ] **Step 1: Write the failing tests**

```python
# studies/fungi/Afumigatus_pangenome/tests/test_cooccurrence.py
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT
from cooccurrence import (
    jaccard, fisher_pvalue, benjamini_hochberg, polarize_direction,
    clade_composition, permutation_null_pvalue, find_cooccurring_pairs,
)


def test_jaccard_identical_and_disjoint():
    assert jaccard([True, True, False], [True, True, False]) == 1.0
    assert jaccard([True, False], [False, True]) == 0.0


def test_fisher_pvalue_significant_for_perfect_cooccurrence():
    a = [True, True, True, True, False, False, False, False]
    b = [True, True, True, True, False, False, False, False]
    p = fisher_pvalue(a, b)
    assert p < 0.05


def test_benjamini_hochberg_orders_correctly():
    pvals = [0.01, 0.02, 0.03, 0.9]
    qvals = benjamini_hochberg(pvals)
    assert len(qvals) == 4
    assert qvals[0] <= qvals[1] <= qvals[2]
    assert qvals[3] >= qvals[2]


def test_polarize_direction():
    assert polarize_direction(present_in_outgroup=True, freq_in_ingroup=0.6) == "loss"
    assert polarize_direction(present_in_outgroup=False, freq_in_ingroup=0.6) == "gain"


def test_clade_composition_counts_per_clade():
    clade_of_strain = {"s1": "Clade_1", "s2": "Clade_1", "s3": "Clade_2"}
    counts = clade_composition(["s1", "s2", "s3"], clade_of_strain)
    assert counts == {"Clade_1": 2, "Clade_2": 1}


def test_permutation_null_pvalue_high_when_confound_is_purely_clade():
    # a and b are both simply "is this strain in Clade_1" -- co-occurrence
    # is entirely explained by clade membership, so the within-clade
    # permutation null should NOT find this significant.
    clade_of_strain = {f"s{i}": ("Clade_1" if i < 4 else "Clade_2") for i in range(8)}
    a = [clade_of_strain[f"s{i}"] == "Clade_1" for i in range(8)]
    b = list(a)
    rng = random.Random(0)
    p = permutation_null_pvalue(a, b, [clade_of_strain[f"s{i}"] for i in range(8)], n_perms=200, rng=rng)
    assert p > 0.05


def test_find_cooccurring_pairs_applies_frequency_floor_and_fdr():
    strains = [f"s{i}" for i in range(10)]
    pm = PresenceMatrix(families=["famA", "famB", "famRare1", "famRare2"], strains=strains)
    # famA/famB co-occur in strains 0-5 (6 strains, clears the floor)
    for s in strains[:6]:
        pm.set_call("famA", s, PRESENT)
        pm.set_call("famB", s, PRESENT)
    # famRare1/famRare2 co-occur only in strains 0-1 (below the floor of 5)
    pm.set_call("famRare1", "s0", PRESENT)
    pm.set_call("famRare1", "s1", PRESENT)
    pm.set_call("famRare2", "s0", PRESENT)
    pm.set_call("famRare2", "s1", PRESENT)

    frequency_table = [
        {"family": "famA", "bin": "shell"}, {"family": "famB", "bin": "shell"},
        {"family": "famRare1", "bin": "cloud"}, {"family": "famRare2", "bin": "cloud"},
    ]
    clade_of_strain = {s: "Clade_1" for s in strains}
    outgroup_presence = {"famA": False, "famB": False, "famRare1": False, "famRare2": False}

    pairs = find_cooccurring_pairs(
        pm, frequency_table, clade_of_strain, outgroup_presence,
        min_strain_count=5, fdr_alpha=0.05, n_perms=100, seed=0,
    )
    reported = {(p["family_a"], p["family_b"]) for p in pairs}
    assert ("famA", "famB") in reported or ("famB", "famA") in reported
    assert ("famRare1", "famRare2") not in reported and ("famRare2", "famRare1") not in reported
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_cooccurrence.py -v`
Expected: FAIL — `cooccurrence` module does not exist yet.

- [ ] **Step 3: Write the implementation**

```python
# studies/fungi/Afumigatus_pangenome/bin/cooccurrence.py
#!/usr/bin/env python3
"""Co-occurrence (co-loss/co-gain) statistics between shell+cloud gene
families -- notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-
design.md, component 4 (revised per the Fable bioinformatics review:
frequency floor + Fisher/BH-FDR instead of a raw Jaccard threshold,
outgroup-based gain/loss polarization, within-clade permutation null as a
cheap stand-in for a full phylogenetic correction).

Usage:
  cooccurrence.py --matrix presence_matrix.rescued.tsv \\
      --frequency_table frequency_table.tsv --config config.csv \\
      --output cooccurring_pairs.tsv
"""
from __future__ import annotations

import argparse
import itertools
import random
import sys
from pathlib import Path

from scipy.stats import fisher_exact, false_discovery_control

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix  # noqa: E402


def jaccard(a: list[bool], b: list[bool]) -> float:
    intersection = sum(1 for x, y in zip(a, b) if x and y)
    union = sum(1 for x, y in zip(a, b) if x or y)
    return intersection / union if union else 0.0


def fisher_pvalue(a: list[bool], b: list[bool]) -> float:
    both = sum(1 for x, y in zip(a, b) if x and y)
    a_only = sum(1 for x, y in zip(a, b) if x and not y)
    b_only = sum(1 for x, y in zip(a, b) if y and not x)
    neither = sum(1 for x, y in zip(a, b) if not x and not y)
    _, p = fisher_exact([[both, a_only], [b_only, neither]], alternative="greater")
    return p


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    if not pvalues:
        return []
    return list(false_discovery_control(pvalues, method="bh"))


def polarize_direction(present_in_outgroup: bool, freq_in_ingroup: float) -> str:
    """A family present in the outgroup that most ingroup strains also
    carry implies the strains lacking it lost it; a family absent from
    the outgroup implies the strains carrying it gained it."""
    if present_in_outgroup:
        return "loss"
    return "gain"


def clade_composition(strains_present: list[str], clade_of_strain: dict[str, str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for s in strains_present:
        clade = clade_of_strain.get(s, "unknown")
        counts[clade] = counts.get(clade, 0) + 1
    return counts


def permutation_null_pvalue(
    a: list[bool], b: list[bool], clades: list[str], n_perms: int, rng: random.Random
) -> float:
    """Shuffle `b` WITHIN each clade group (not globally) n_perms times,
    recompute Fisher's p-value each time, and return the fraction of
    permutations at least as extreme as the observed statistic -- a
    pair significant only because both families mark the same clade
    will look unremarkable under this null."""
    observed = fisher_pvalue(a, b)
    indices_by_clade: dict[str, list[int]] = {}
    for i, clade in enumerate(clades):
        indices_by_clade.setdefault(clade, []).append(i)

    at_least_as_extreme = 0
    for _ in range(n_perms):
        b_perm = list(b)
        for indices in indices_by_clade.values():
            values = [b_perm[i] for i in indices]
            rng.shuffle(values)
            for i, v in zip(indices, values):
                b_perm[i] = v
        if fisher_pvalue(a, b_perm) <= observed:
            at_least_as_extreme += 1
    return at_least_as_extreme / n_perms


def find_cooccurring_pairs(
    matrix: PresenceMatrix,
    frequency_table: list[dict],
    clade_of_strain: dict[str, str],
    outgroup_presence: dict[str, bool],
    min_strain_count: int = 5,
    fdr_alpha: float = 0.05,
    n_perms: int = 1000,
    seed: int = 0,
) -> list[dict]:
    eligible = [
        row["family"] for row in frequency_table
        if row["bin"] in ("shell", "cloud") and matrix.strain_count(row["family"]) >= min_strain_count
    ]
    clades = [clade_of_strain.get(s, "unknown") for s in matrix.strains]

    candidates = []
    for fam_a, fam_b in itertools.combinations(eligible, 2):
        vec_a = matrix.presence_vector(fam_a)
        vec_b = matrix.presence_vector(fam_b)
        p = fisher_pvalue(vec_a, vec_b)
        candidates.append((fam_a, fam_b, vec_a, vec_b, p))

    if not candidates:
        return []
    qvalues = benjamini_hochberg([c[4] for c in candidates])

    rng = random.Random(seed)
    results = []
    for (fam_a, fam_b, vec_a, vec_b, p), q in zip(candidates, qvalues):
        if q >= fdr_alpha:
            continue
        strains_present_a = [s for s, present in zip(matrix.strains, vec_a) if present]
        results.append({
            "family_a": fam_a,
            "family_b": fam_b,
            "jaccard": jaccard(vec_a, vec_b),
            "fisher_p": p,
            "fdr_q": q,
            "permutation_p": permutation_null_pvalue(vec_a, vec_b, clades, n_perms, rng),
            "direction_a": polarize_direction(outgroup_presence.get(fam_a, False), sum(vec_a) / len(vec_a)),
            "clade_composition": clade_composition(strains_present_a, clade_of_strain),
        })
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--min_strain_count", type=int, default=5)
    ap.add_argument("--fdr_alpha", type=float, default=0.05)
    ap.add_argument("--n_perms", type=int, default=1000)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    NOVINVENIO_LIB = Path(__file__).resolve().parents[4] / "NovInvenio" / "lib"
    sys.path.insert(0, str(NOVINVENIO_LIB))
    from config_parser import parse_config  # noqa: E402

    samples = parse_config(args.config)
    clade_of_strain = {s.short: s.taxon_group for s in samples}
    outgroup_shorts = {s.short for s in samples if s.group == "OUT"}

    matrix = PresenceMatrix.from_tsv(args.matrix)
    with open(args.frequency_table) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        frequency_table = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]

    outgroup_presence = {
        fam: any(matrix.is_present(fam, s) for s in outgroup_shorts) for fam in matrix.families
    }

    pairs = find_cooccurring_pairs(
        matrix, frequency_table, clade_of_strain, outgroup_presence,
        args.min_strain_count, args.fdr_alpha, args.n_perms,
    )
    with open(args.output, "w") as fh:
        fh.write("family_a\tfamily_b\tjaccard\tfisher_p\tfdr_q\tpermutation_p\tdirection_a\tclade_composition\n")
        for row in pairs:
            fh.write(
                f"{row['family_a']}\t{row['family_b']}\t{row['jaccard']:.4f}\t"
                f"{row['fisher_p']:.2e}\t{row['fdr_q']:.2e}\t{row['permutation_p']:.4f}\t"
                f"{row['direction_a']}\t{row['clade_composition']}\n"
            )


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_cooccurrence.py -v`
Expected: PASS (7 tests)

- [ ] **Step 5: Run against the real dataset**

```bash
pixi run bin/cooccurrence.py --matrix presence_matrix.rescued.tsv \
    --frequency_table frequency_table.tsv --config config.csv \
    --output cooccurring_pairs.tsv
```
Inspect the top hits by `fdr_q`: for each, check `clade_composition` — a pair whose presence is concentrated in one `TaxonGroup` clade is a candidate clade-confound even if it clears the permutation null (report both numbers, don't silently drop it).

- [ ] **Step 6: Commit**

```bash
git add studies/fungi/Afumigatus_pangenome/bin/cooccurrence.py \
        studies/fungi/Afumigatus_pangenome/tests/test_cooccurrence.py
git commit -m "Add co-occurrence statistics: Fisher/BH-FDR, outgroup polarization, clade-stratified null"
```

---

### Task 8: Synteny / physical-clustering (accessory-island) analysis

**Files:**
- Create: `studies/fungi/Afumigatus_pangenome/bin/synteny_windows.py`
- Test: `studies/fungi/Afumigatus_pangenome/tests/test_synteny_windows.py`

**Interfaces:**
- Consumes: `PresenceMatrix` from Task 2
- Produces: `parse_gff3_gene_order(gff3_path) -> list[tuple[str,str,int,int]]` (`gene_id, contig_id, start, end`, sorted by contig then start); `sliding_windows(gene_order, n) -> list[list[str]]` (excludes cross-contig windows); `accessory_islands(gene_order, is_core) -> list[list[str]]` (maximal runs of consecutive non-core genes per contig); `linkage_fraction(family_a_genes, family_b_genes, gene_position, k=10) -> float` (fraction of co-carrying strains where the two families' genes are within k genes on the same contig). Task 10 (benchmark scorecard) uses these against Table S14-16's known region boundaries.

- [ ] **Step 1: Write the failing tests**

```python
# studies/fungi/Afumigatus_pangenome/tests/test_synteny_windows.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from synteny_windows import parse_gff3_gene_order, sliding_windows, accessory_islands, linkage_fraction


def test_parse_gff3_gene_order_sorts_by_contig_and_start(tmp_path):
    gff3 = tmp_path / "strain.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "contig1\tsrc\tgene\t500\t600\t.\t+\t.\tID=geneB\n"
        "contig1\tsrc\tgene\t100\t200\t.\t+\t.\tID=geneA\n"
        "contig2\tsrc\tgene\t1\t50\t.\t+\t.\tID=geneC\n"
        "contig1\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=geneA.t1;Parent=geneA\n"
    )
    order = parse_gff3_gene_order(gff3)
    assert order == [
        ("geneA", "contig1", 100, 200),
        ("geneB", "contig1", 500, 600),
        ("geneC", "contig2", 1, 50),
    ]


def test_sliding_windows_excludes_cross_contig():
    gene_order = [
        ("g1", "c1", 1, 10), ("g2", "c1", 20, 30), ("g3", "c1", 40, 50),
        ("g4", "c2", 1, 10), ("g5", "c2", 20, 30),
    ]
    windows = sliding_windows(gene_order, n=3)
    gene_id_windows = [[g[0] for g in w] for w in windows]
    assert ["g1", "g2", "g3"] in gene_id_windows
    # no window of size 3 crosses from c1 into c2
    assert not any("g4" in w and "g3" in w for w in gene_id_windows)


def test_accessory_islands_finds_maximal_noncore_runs():
    gene_order = [
        ("g1", "c1", 1, 10), ("g2", "c1", 20, 30), ("g3", "c1", 40, 50),
        ("g4", "c1", 60, 70), ("g5", "c1", 80, 90),
    ]
    is_core = {"g1": True, "g2": False, "g3": False, "g4": False, "g5": True}
    islands = accessory_islands(gene_order, is_core)
    assert [g[0] for g in islands[0]] == ["g2", "g3", "g4"]


def test_linkage_fraction_measures_physical_proximity():
    gene_position = {
        "s1": {"famA": ("c1", 5), "famB": ("c1", 7)},   # 2 genes apart, within k
        "s2": {"famA": ("c1", 5), "famB": ("c1", 500)},  # far apart, same contig
        "s3": {"famA": ("c1", 5)},                       # famB absent in s3
    }
    frac = linkage_fraction("famA", "famB", gene_position, k=10)
    assert frac == 0.5  # only s1 of {s1, s2} (both-present strains) is within k
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_synteny_windows.py -v`
Expected: FAIL — `synteny_windows` module does not exist yet.

- [ ] **Step 3: Write the implementation**

```python
# studies/fungi/Afumigatus_pangenome/bin/synteny_windows.py
#!/usr/bin/env python3
"""Per-strain gene-order synteny scanning for Starship-driven physical
clustering -- notes/superpowers/specs/2026-09-13-pangenome-cluster-
profile-design.md, component 5 (refined per the Fable bioinformatics
review: per-strain window enumeration, contig-edge exclusion, an
accessory-island test alongside fixed-size windows, and a direct
linkage-fraction statistic per correlated family pair).

Usage:
  synteny_windows.py --gff3_dir data_dir/gff3 --config config.csv \\
      --matrix presence_matrix.rescued.tsv --window_sizes 3,10 \\
      --output synteny_windows.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))


_ID_RE = re.compile(r"ID=([^;\n]+)")


def parse_gff3_gene_order(gff3_path: str | Path) -> list[tuple[str, str, int, int]]:
    """Return [(gene_id, contig_id, start, end), ...] for every 'gene'
    feature, sorted by contig then start position."""
    genes = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            contig, start, end, attrs = fields[0], int(fields[3]), int(fields[4]), fields[8]
            m = _ID_RE.search(attrs)
            if not m:
                continue
            genes.append((m.group(1), contig, start, end))
    genes.sort(key=lambda g: (g[1], g[2]))
    return genes


def sliding_windows(gene_order: list[tuple[str, str, int, int]], n: int) -> list[list[tuple]]:
    windows = []
    for i in range(len(gene_order) - n + 1):
        window = gene_order[i:i + n]
        contigs = {g[1] for g in window}
        if len(contigs) == 1:
            windows.append(window)
    return windows


def accessory_islands(
    gene_order: list[tuple[str, str, int, int]], is_core: dict[str, bool]
) -> list[list[tuple]]:
    """Maximal runs of consecutive non-core genes, never crossing a
    contig boundary."""
    islands: list[list[tuple]] = []
    current: list[tuple] = []
    prev_contig = None
    for gene in gene_order:
        gene_id, contig = gene[0], gene[1]
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
    return islands


def linkage_fraction(
    family_a: str,
    family_b: str,
    gene_position: dict[str, dict[str, tuple[str, int]]],
    k: int = 10,
) -> float:
    """Fraction of strains carrying BOTH families where their gene-order
    positions are within k genes of each other on the same contig."""
    both_present = [
        s for s, positions in gene_position.items()
        if family_a in positions and family_b in positions
    ]
    if not both_present:
        return 0.0
    linked = 0
    for s in both_present:
        contig_a, pos_a = gene_position[s][family_a]
        contig_b, pos_b = gene_position[s][family_b]
        if contig_a == contig_b and abs(pos_a - pos_b) <= k:
            linked += 1
    return linked / len(both_present)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gff3_dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--window_sizes", default="3,10")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    # Full per-strain window enumeration + correlated-window detection
    # across all 293 strains' GFF3s is run interactively (Step 5 below)
    # rather than duplicated here -- this CLI's job is exposing the
    # tested building blocks above for that run.
    print(
        "Use parse_gff3_gene_order/sliding_windows/accessory_islands/"
        "linkage_fraction directly for the real run; see plan Task 8 Step 5.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_synteny_windows.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Run against the real dataset**

Write a short interactive/notebook-style script (not committed as a formal CLI, since this step's exact shape depends on what Task 7's co-occurring pairs actually look like) that: for each strain, calls `parse_gff3_gene_order` on its GFF3, builds `gene_position` (gene_id -> family via the tier-1 cluster membership from Task 4, then family -> (contig, index-in-order)), and for every co-occurring pair from Task 7's `cooccurring_pairs.tsv`, computes `linkage_fraction` at k=3 and k=10. Report which pairs are both statistically co-occurring (Task 7) AND physically linked (this task) — those are the strongest Starship-cargo candidates.

- [ ] **Step 6: Commit**

```bash
git add studies/fungi/Afumigatus_pangenome/bin/synteny_windows.py \
        studies/fungi/Afumigatus_pangenome/tests/test_synteny_windows.py
git commit -m "Add per-strain synteny windows, accessory islands, and linkage-fraction statistic"
```

---

### Task 9: Starship-supplement ground truth extraction + ID crosswalk

**Files:**
- Create: `studies/fungi/Afumigatus_pangenome/bin/parse_starship_supplement.py`
- Create: `studies/fungi/Afumigatus_pangenome/bin/id_crosswalk.py`
- Test: `studies/fungi/Afumigatus_pangenome/tests/test_parse_starship_supplement.py`
- Test: `studies/fungi/Afumigatus_pangenome/tests/test_id_crosswalk.py`

**Interfaces:**
- Consumes: nothing new (reads `mbio.01092-25-s0002.xlsx`)
- Produces: `high_confidence_starships(table_s6_df, table_s21_df) -> dict[str, dict]` (`name_id -> {"population_freq": float, "presence": {isolate_id: bool}}`); `cargo_gene_sets(table_s12_or_s13_df) -> dict[str, set[str]]` (`starship_id -> {gene_id, ...}`); `parse_diamond_blastp_besthits(lines) -> dict[str,str]` (query_id -> best-hit subject_id, one row per query, highest bitscore kept). Task 10 consumes all three to build ground truth and crosswalk this study's own IDs to the paper's.

- [ ] **Step 1: Write the failing tests**

```python
# studies/fungi/Afumigatus_pangenome/tests/test_parse_starship_supplement.py
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from parse_starship_supplement import high_confidence_starships, cargo_gene_sets


def test_high_confidence_starships_combines_frequency_and_presence():
    table_s6 = pd.DataFrame({"nameID": ["Gnosis-h1", "Osiris-h4"], "freq": [0.599, 0.20]})
    table_s21 = pd.DataFrame({
        "isolateID": ["iso1", "iso2", "iso1", "iso2"],
        "nameID": ["Gnosis-h1", "Gnosis-h1", "Osiris-h4", "Osiris-h4"],
        "presence/absence": [1, 0, 0, 1],
    })
    result = high_confidence_starships(table_s6, table_s21)
    assert result["Gnosis-h1"]["population_freq"] == 0.599
    assert result["Gnosis-h1"]["presence"] == {"iso1": True, "iso2": False}
    assert result["Osiris-h4"]["presence"] == {"iso1": False, "iso2": True}


def test_cargo_gene_sets_groups_by_starship_id():
    table_s13 = pd.DataFrame({
        "starshipID": ["47-10_s00011", "47-10_s00011", "1F1SW-F4_e00006"],
        "geneID": ["47-10_000766", "47-10_000767", "1F1SW-F4_g001"],
    })
    result = cargo_gene_sets(table_s13)
    assert result["47-10_s00011"] == {"47-10_000766", "47-10_000767"}
    assert result["1F1SW-F4_e00006"] == {"1F1SW-F4_g001"}
```

```python
# studies/fungi/Afumigatus_pangenome/tests/test_id_crosswalk.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from id_crosswalk import parse_diamond_blastp_besthits


def test_parse_diamond_blastp_besthits_keeps_highest_bitscore():
    # outfmt6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    lines = [
        "XP_748727.1\tstudy_protein_A\t99.0\t500\t0\t0\t1\t500\t1\t500\t0.0\t950",
        "XP_748727.1\tstudy_protein_B\t60.0\t300\t0\t0\t1\t300\t1\t300\t1e-20\t150",
        "AFUB_079030\tstudy_protein_C\t95.0\t400\t0\t0\t1\t400\t1\t400\t0.0\t800",
    ]
    result = parse_diamond_blastp_besthits(lines)
    assert result == {"XP_748727.1": "study_protein_A", "AFUB_079030": "study_protein_C"}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_parse_starship_supplement.py studies/fungi/Afumigatus_pangenome/tests/test_id_crosswalk.py -v`
Expected: FAIL — neither module exists yet.

- [ ] **Step 3: Write the implementation**

```python
# studies/fungi/Afumigatus_pangenome/bin/parse_starship_supplement.py
#!/usr/bin/env python3
"""Extract structured ground truth from mbio.01092-25-s0002.xlsx (Gluck-
Thaler et al. 2025, doi:10.1128/mbio.01092-25) for the benchmark
scorecard -- notes/superpowers/specs/2026-09-13-pangenome-cluster-
profile-design.md, component 6.

Usage:
  parse_starship_supplement.py --xlsx mbio.01092-25-s0002.xlsx \\
      --output_prefix starship_ground_truth
"""
from __future__ import annotations

import argparse

import pandas as pd


def high_confidence_starships(table_s6: pd.DataFrame, table_s21: pd.DataFrame) -> dict:
    """Combine Table S6 (nameID, freq) with Table S21 (isolateID, nameID,
    presence/absence) into {name_id: {population_freq, presence: {isolate: bool}}}."""
    result = {}
    freq_by_name = dict(zip(table_s6["nameID"], table_s6["freq"]))
    for name_id, freq in freq_by_name.items():
        subset = table_s21[table_s21["nameID"] == name_id]
        presence = dict(zip(subset["isolateID"], subset["presence/absence"].astype(bool)))
        result[name_id] = {"population_freq": freq, "presence": presence}
    return result


def cargo_gene_sets(table_s12_or_s13: pd.DataFrame) -> dict[str, set[str]]:
    """Group Table S12/S13's (starshipID, geneID) rows into
    {starship_id: {gene_id, ...}}."""
    result: dict[str, set[str]] = {}
    for starship_id, gene_id in zip(table_s12_or_s13["starshipID"], table_s12_or_s13["geneID"]):
        result.setdefault(starship_id, set()).add(gene_id)
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--output_prefix", required=True)
    args = ap.parse_args()

    xl = pd.ExcelFile(args.xlsx)
    table_s6 = xl.parse("Table S6", header=1)
    table_s21 = xl.parse("Table S21", header=1)
    table_s13 = xl.parse("Table S13", header=1)

    starships = high_confidence_starships(table_s6, table_s21)
    cargo = cargo_gene_sets(table_s13)

    with open(f"{args.output_prefix}_starships.tsv", "w") as fh:
        fh.write("nameID\tpopulation_freq\tpresence_json\n")
        import json
        for name_id, info in starships.items():
            fh.write(f"{name_id}\t{info['population_freq']}\t{json.dumps(info['presence'])}\n")

    with open(f"{args.output_prefix}_cargo.tsv", "w") as fh:
        fh.write("starshipID\tgeneID\n")
        for starship_id, genes in cargo.items():
            for gene in sorted(genes):
                fh.write(f"{starship_id}\t{gene}\n")


if __name__ == "__main__":
    main()
```

```python
# studies/fungi/Afumigatus_pangenome/bin/id_crosswalk.py
#!/usr/bin/env python3
"""Sequence-based crosswalk between the paper's gene/protein IDs (AFUB_*,
Afu*g*, per-strain g#, XP_* accessions) and this study's own protein IDs
-- notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-
design.md, component 6's required ID-crosswalk control. String matching
does not work (the two ID schemes are unrelated); this runs the paper's
cited protein sequences (fetched separately, e.g. via NCBI efetch for the
accessions in Table S19/S3-S16, into a FASTA passed as --paper_fasta)
through diamond blastp against this study's own proteomes and keeps the
best hit per paper ID.

Usage:
  diamond makedb --in all_ingroup.fa -d all_ingroup
  diamond blastp -q paper_reference_proteins.fa -d all_ingroup \\
      -o paper_vs_study.tsv --outfmt 6
  id_crosswalk.py --diamond_tsv paper_vs_study.tsv --output crosswalk.tsv
"""
from __future__ import annotations

import argparse


def parse_diamond_blastp_besthits(lines: list[str]) -> dict[str, str]:
    """Parse diamond blastp outfmt6 lines and keep, for each query
    (paper gene/protein ID), the subject (this study's protein ID) with
    the highest bitscore."""
    best_bitscore: dict[str, float] = {}
    best_subject: dict[str, str] = {}
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        query, subject, bitscore = parts[0], parts[1], float(parts[11])
        if bitscore > best_bitscore.get(query, -1.0):
            best_bitscore[query] = bitscore
            best_subject[query] = subject
    return best_subject


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--diamond_tsv", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open(args.diamond_tsv) as fh:
        crosswalk = parse_diamond_blastp_besthits(fh.readlines())
    with open(args.output, "w") as fh:
        fh.write("paper_id\tstudy_protein_id\n")
        for paper_id, study_id in sorted(crosswalk.items()):
            fh.write(f"{paper_id}\t{study_id}\n")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_parse_starship_supplement.py studies/fungi/Afumigatus_pangenome/tests/test_id_crosswalk.py -v`
Expected: PASS (4 tests total)

- [ ] **Step 5: Run against the real dataset**

```bash
chmod +x bin/parse_starship_supplement.py bin/id_crosswalk.py
pixi run bin/parse_starship_supplement.py --xlsx mbio.01092-25-s0002.xlsx \
    --output_prefix starship_ground_truth
```
For the crosswalk: fetch the paper's cited protein accessions (e.g. `XP_748727.1` for hacA, plus every distinct ID referenced in `starship_ground_truth_cargo.tsv`) into a FASTA via NCBI efetch, then run the `diamond makedb`/`diamond blastp`/`id_crosswalk.py` pipeline documented in the script's docstring. Also build the strain-overlap table here: cross-reference `config.csv`'s `Short`/`Strain` against Table S2/S1's isolate names (manual/semi-manual matching — record the mapping in `PANGENOME_CLUSTER_PROFILE_NOTES.md`, not committed as guessed code).

- [ ] **Step 6: Commit**

```bash
git add studies/fungi/Afumigatus_pangenome/bin/parse_starship_supplement.py \
        studies/fungi/Afumigatus_pangenome/bin/id_crosswalk.py \
        studies/fungi/Afumigatus_pangenome/tests/test_parse_starship_supplement.py \
        studies/fungi/Afumigatus_pangenome/tests/test_id_crosswalk.py
git commit -m "Add Starship-supplement ground-truth extraction and sequence-based ID crosswalk"
```

---

### Task 10: Benchmark scorecard (positive + negative controls, mmseqs vs. diamond)

**Files:**
- Create: `studies/fungi/Afumigatus_pangenome/bin/benchmark_scorecard.py`
- Test: `studies/fungi/Afumigatus_pangenome/tests/test_benchmark_scorecard.py`

**Interfaces:**
- Consumes: `PresenceMatrix` (Task 2); `two_tier_families` (Task 4); ground truth from Task 9 (`high_confidence_starships`, `cargo_gene_sets`, crosswalk TSV)
- Produces: `score_presence_recovery(predicted, truth) -> dict` (`tp, fp, tn, fn, accuracy, jaccard`); `score_cargo_grouping(predicted_family_of_gene, truth_cargo_sets) -> dict[str, dict]` (per-Starship `purity`, `completeness`); `run_scorecard(...) -> list[dict]` (one row per control x backend). This is the empirical answer to the mmseqs-vs-diamond question and the gate before Task 11's HAC screen is trusted.

- [ ] **Step 1: Write the failing tests**

```python
# studies/fungi/Afumigatus_pangenome/tests/test_benchmark_scorecard.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from benchmark_scorecard import score_presence_recovery, score_cargo_grouping


def test_score_presence_recovery_computes_confusion_counts():
    predicted = {"iso1": True, "iso2": False, "iso3": True, "iso4": False}
    truth =     {"iso1": True, "iso2": True,  "iso3": True, "iso4": False}
    result = score_presence_recovery(predicted, truth)
    assert result["tp"] == 2
    assert result["fn"] == 1
    assert result["fp"] == 0
    assert result["tn"] == 1
    assert result["accuracy"] == 0.75


def test_score_cargo_grouping_purity_and_completeness():
    # our clustering put gene1, gene2 in family "famX", and gene3 alone in "famY"
    predicted_family_of_gene = {"gene1": "famX", "gene2": "famX", "gene3": "famY"}
    # ground truth: gene1, gene2, gene3 all belong to the same Starship's cargo
    truth_cargo_sets = {"starship1": {"gene1", "gene2", "gene3"}}
    result = score_cargo_grouping(predicted_family_of_gene, truth_cargo_sets)
    # famX correctly grouped 2 of the 3 true cargo genes together (some completeness),
    # and every gene it contains is a true cargo member (perfect purity)
    assert result["starship1"]["purity"] == 1.0
    assert round(result["starship1"]["completeness"], 4) == round(2 / 3, 4)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_benchmark_scorecard.py -v`
Expected: FAIL — `benchmark_scorecard` module does not exist yet.

- [ ] **Step 3: Write the implementation**

```python
# studies/fungi/Afumigatus_pangenome/bin/benchmark_scorecard.py
#!/usr/bin/env python3
"""Score each clustering backend (mmseqs2 vs. diamond) against known
Starship positive controls (Tables S6/S21/S12/S13/S7/S14-16) plus
negative controls (conserved non-mobile SM clusters, random family
pairs) -- notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-
design.md, component 6.

Usage:
  benchmark_scorecard.py --crosswalk crosswalk.tsv \\
      --ground_truth_starships starship_ground_truth_starships.tsv \\
      --ground_truth_cargo starship_ground_truth_cargo.tsv \\
      --matrix_mmseqs presence_matrix.mmseqs.rescued.tsv \\
      --matrix_diamond presence_matrix.diamond.rescued.tsv \\
      --output benchmark_scorecard.tsv
"""
from __future__ import annotations

import argparse


def score_presence_recovery(predicted: dict[str, bool], truth: dict[str, bool]) -> dict:
    tp = fp = tn = fn = 0
    for key, truth_val in truth.items():
        pred_val = predicted.get(key, False)
        if truth_val and pred_val:
            tp += 1
        elif truth_val and not pred_val:
            fn += 1
        elif not truth_val and pred_val:
            fp += 1
        else:
            tn += 1
    total = tp + fp + tn + fn
    accuracy = (tp + tn) / total if total else 0.0
    union = tp + fp + fn
    jaccard = tp / union if union else 0.0
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn, "accuracy": accuracy, "jaccard": jaccard}


def score_cargo_grouping(
    predicted_family_of_gene: dict[str, str], truth_cargo_sets: dict[str, set[str]]
) -> dict[str, dict]:
    """For each true cargo set, find the predicted family containing the
    largest overlap with it, and report purity (fraction of that
    family's members that are true cargo members) and completeness
    (fraction of true cargo members recovered in that family)."""
    result = {}
    for starship_id, true_genes in truth_cargo_sets.items():
        family_hits: dict[str, set[str]] = {}
        for gene in true_genes:
            fam = predicted_family_of_gene.get(gene)
            if fam is not None:
                family_hits.setdefault(fam, set()).add(gene)
        if not family_hits:
            result[starship_id] = {"purity": 0.0, "completeness": 0.0}
            continue
        best_family = max(family_hits, key=lambda f: len(family_hits[f]))
        family_members = {
            g for g, f in predicted_family_of_gene.items() if f == best_family
        }
        overlap = family_hits[best_family]
        purity = len(overlap) / len(family_members) if family_members else 0.0
        completeness = len(overlap) / len(true_genes) if true_genes else 0.0
        result[starship_id] = {"purity": purity, "completeness": completeness}
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--crosswalk", required=True)
    ap.add_argument("--ground_truth_starships", required=True)
    ap.add_argument("--ground_truth_cargo", required=True)
    ap.add_argument("--matrix_mmseqs", required=True)
    ap.add_argument("--matrix_diamond", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    # The real run wires PresenceMatrix.from_tsv() for each backend, the
    # crosswalk TSV, and Task 9's ground-truth TSVs into calls to
    # score_presence_recovery/score_cargo_grouping per control Starship;
    # see Task 10 Step 5 in the plan for the exact sequence, since it
    # depends on the crosswalk actually being populated first (Task 9).
    print("See plan Task 10 Step 5 for the real-data scorecard run.", )


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_benchmark_scorecard.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Run the real benchmark**

Using the crosswalk (Task 9) to translate each control Starship's cargo genes (Table S12/S13) and per-strain presence (Table S21) into this study's own IDs (where the strain-overlap table says a given control strain is actually in this study's 293), for each backend (mmseqs, diamond):
1. `score_presence_recovery` — predicted presence (from that backend's rescued `PresenceMatrix`) vs. Table S21's truth, for every overlapping strain x control Starship.
2. `score_cargo_grouping` — predicted tier-1 family membership vs. Table S12/S13's cargo sets.
3. Negative controls: pick 2-3 conserved secondary-metabolite gene clusters (present in ~all strains — cross-reference published *A. fumigatus* SM cluster gene lists, e.g. fumagillin/pseurotin, gliotoxin) and confirm `score_presence_recovery`/synteny (Task 8) do NOT flag them as variably mobile; and score a handful of random shell-family pairs (matched for frequency) through Task 7's `find_cooccurring_pairs` to estimate its real false-positive rate.

Write results to `benchmark_scorecard.tsv` (one row per control x backend x metric) and record the mmseqs-vs-diamond conclusion, plus any parameter adjustment this suggests (e.g. tier-1 identity), in `PANGENOME_CLUSTER_PROFILE_NOTES.md`.

- [ ] **Step 6: Commit**

```bash
git add studies/fungi/Afumigatus_pangenome/bin/benchmark_scorecard.py \
        studies/fungi/Afumigatus_pangenome/tests/test_benchmark_scorecard.py
git commit -m "Add benchmark scorecard: mmseqs vs. diamond against known-Starship positive/negative controls"
```

---

### Task 11: HAC (hrmA-Associated Cluster) + hacA targeted screen

**Files:**
- Create: `studies/fungi/Afumigatus_pangenome/bin/hac_screen.py`
- Test: `studies/fungi/Afumigatus_pangenome/tests/test_hac_screen.py`

**Interfaces:**
- Consumes: `PresenceMatrix` (Task 2); crosswalk TSV (Task 9); benchmark scorecard result (Task 10) as the gate before trusting this screen's output
- Produces: `screen_family(matrix, family_id, strain_to_starship) -> list[dict]` (one row per strain: `strain, state, copy_number, starship`, `starship` = `"unknown"` when no paper Starship-assignment data exists for that strain); a `hac_screen_report.tsv` with one section for the HAC family (PF11001/IPR047092 paralogs) and one for the independent `hacA` locus.

- [ ] **Step 1: Write the failing test**

```python
# studies/fungi/Afumigatus_pangenome/tests/test_hac_screen.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY
from hac_screen import screen_family


def test_screen_family_reports_state_copy_number_and_starship():
    pm = PresenceMatrix(families=["hacFamily"], strains=["s1", "s2", "s3"])
    pm.set_call("hacFamily", "s1", PRESENT, copies=2)
    pm.set_call("hacFamily", "s2", GENOME_ONLY)
    # s3 left ABSENT (default)
    strain_to_starship = {"s1": "Nebuchadnezzar-h1"}

    rows = screen_family(pm, "hacFamily", strain_to_starship)
    by_strain = {r["strain"]: r for r in rows}

    assert by_strain["s1"]["state"] == PRESENT
    assert by_strain["s1"]["copy_number"] == 2
    assert by_strain["s1"]["starship"] == "Nebuchadnezzar-h1"
    assert by_strain["s2"]["state"] == GENOME_ONLY
    assert by_strain["s2"]["starship"] == "unknown"
    assert by_strain["s3"]["copy_number"] == 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_hac_screen.py -v`
Expected: FAIL — `hac_screen` module does not exist yet.

- [ ] **Step 3: Write the implementation**

```python
# studies/fungi/Afumigatus_pangenome/bin/hac_screen.py
#!/usr/bin/env python3
"""Targeted screen for the HAC (hrmA-Associated Cluster, PF11001/
IPR047092 paralog family) and the independent hacA (Afu3g04070) locus
across all 293 strains -- notes/superpowers/specs/2026-09-13-pangenome-
cluster-profile-design.md, component 7. Run only after Task 10's
benchmark scorecard has validated the clustering backend being used
here (per the spec's explicit ask to test recovery of known Starships
before trusting the method on this target).

Usage:
  hac_screen.py --matrix presence_matrix.rescued.tsv \\
      --hac_family_id <tier1_rep_for_the_PF11001_family> \\
      --haca_family_id <tier1_rep_for_hacA> \\
      --strain_starship_map strain_to_starship.tsv \\
      --output hac_screen_report.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix  # noqa: E402


def screen_family(
    matrix: PresenceMatrix, family_id: str, strain_to_starship: dict[str, str]
) -> list[dict]:
    rows = []
    for strain in matrix.strains:
        state = matrix.call(family_id, strain)
        rows.append({
            "strain": strain,
            "state": state,
            "copy_number": matrix.copy_number.get((family_id, strain), 0),
            "starship": strain_to_starship.get(strain, "unknown"),
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--hac_family_id", required=True)
    ap.add_argument("--haca_family_id", required=True)
    ap.add_argument("--strain_starship_map", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)
    strain_to_starship: dict[str, str] = {}
    with open(args.strain_starship_map) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                strain_to_starship[parts[0]] = parts[1]

    hac_rows = screen_family(matrix, args.hac_family_id, strain_to_starship)
    haca_rows = screen_family(matrix, args.haca_family_id, strain_to_starship)

    with open(args.output, "w") as fh:
        fh.write("locus\tstrain\tstate\tcopy_number\tstarship\n")
        for row in hac_rows:
            fh.write(f"HAC\t{row['strain']}\t{row['state']}\t{row['copy_number']}\t{row['starship']}\n")
        for row in haca_rows:
            fh.write(f"hacA\t{row['strain']}\t{row['state']}\t{row['copy_number']}\t{row['starship']}\n")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run pytest studies/fungi/Afumigatus_pangenome/tests/test_hac_screen.py -v`
Expected: PASS

- [ ] **Step 5: Run against the real dataset**

Identify `--hac_family_id`/`--haca_family_id` (the tier-1 representative IDs whose members include the crosswalked PF11001-family genes / `hacA` respectively, from Task 9's crosswalk), build `strain_to_starship.tsv` for whichever strains overlap the paper's Starship-annotated set, then:

```bash
chmod +x bin/hac_screen.py
pixi run bin/hac_screen.py --matrix presence_matrix.rescued.tsv \
    --hac_family_id <id> --haca_family_id <id> \
    --strain_starship_map strain_to_starship.tsv \
    --output hac_screen_report.tsv
```
Summarize: how many of the 293 strains carry the HAC family at all, how many carry more than one copy (candidate paralog expansion), and — for strains with known Starship assignment — which Starship each currently rides in, cross-checked against the paper's documented instances (`Nebuchadnezzar-h1`, `Osiris-h3`, `Logos-h1`, `Gnosis-h2`, `Logos-h2`, `navis10-var35`).

- [ ] **Step 6: Commit**

```bash
git add studies/fungi/Afumigatus_pangenome/bin/hac_screen.py \
        studies/fungi/Afumigatus_pangenome/tests/test_hac_screen.py
git commit -m "Add HAC/hacA targeted screen"
```

---

## Self-Review

**Spec coverage:** Component 1 (two-tier clustering) -> Task 4; component 1b (rescue + dereplication) -> Tasks 3, 5; component 2 (presence matrix) -> Task 2; component 3 (frequency binning) -> Task 6; component 4 (co-occurrence) -> Task 7; component 5 (synteny) -> Task 8; component 6 (benchmark suite, positive + negative controls, ID crosswalk, strain overlap) -> Tasks 9, 10; component 7 (HAC/hacA screen) -> Task 11. Open controls (draft/long-read fraction, dereplication, core-cutoff sensitivity) are addressed in Tasks 3 and 6's "run against real data" steps rather than as separate tasks, since they're properties of the real dataset a script can't unit-test in isolation.

**Placeholder scan:** No TBD/TODO markers. Two CLI `main()` bodies (Tasks 8, 10) intentionally print a pointer to the plan's "Step 5" real-data instructions rather than duplicating a long interactive analysis inline — this is a documented design choice (the real run's exact shape depends on upstream outputs not yet produced when the script is written), not an unfinished implementation; the tested, reusable functions in each module are complete.

**Type/interface consistency:** `PresenceMatrix`, `PRESENT`/`GENOME_ONLY`/`ABSENT`, `read_cluster_tsv`, `build_families` (Task 2) are imported with the same names/signatures in Tasks 3-11. `two_tier_families()`'s return shape (`{"members": [...], "superfamily": ...}`, Task 4) matches how Task 10 describes consuming it. `compute_frequency_table()`'s row shape (`family, frequency, strain_count, bin`, Task 6) matches Task 7's `frequency_table` parameter usage.

**Scope:** Eleven tasks is a lot for one plan, but they chain linearly (each depends on the previous task's real output, not just its interface) and the spec explicitly rejected splitting this into a multi-phase NovInvenio workflow in favor of one standalone analysis — decomposing further would just be renaming the same dependency chain.
