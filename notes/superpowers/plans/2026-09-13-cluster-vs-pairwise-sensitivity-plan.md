# Cluster-vs-Pairwise Sensitivity Investigation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the standalone scripts and run the real comparisons needed to
measure, for `pezizo_set1` (and a second clade, Agaricales/agaricomycetes),
how much novelty/loss-detection recall the mmseqs cluster pathway loses
relative to the pairwise pathway, whether raw cluster membership alone (no
HMM) does better or worse, and whether a targeted within-family paralog
-split refinement (Tier R) recovers the difference without inflating false
novelty calls — producing one report with a clear recommendation.

**Architecture:** Four presence-calling tiers (P/C+H/C/R) scored through the
same two code paths: `nf_NovInvenio/bin/score_controls.py` (extended with a
`--presence-mode` switch, Tiers C/R reuse it unchanged by pointing at
different `--cluster-tsv` input) for the 16/2-control recall/FP comparison,
and three new `NovInvenio_Investigations/bin/` scripts for genome-wide
concordance, compute-cost accounting, and report assembly. `refine_ambiguous_
families.py` (new, `nf_NovInvenio/bin/`) produces Tier R's refined
`families_cluster.tsv`-shaped input by running one small within-family
diamond all-vs-all (not a reuse of Tier P's existing search) and cutting
edges the pipeline's own already-published per-genome `self_hits/*.paralog_
cutoffs.tsv` marks as a registered paralog pair. Everything here is
standalone scripts against already-published results — no `nf_NovInvenio`
Nextflow/DSL2 changes (that is this spec's explicit follow-on, not this
plan).

**Tech Stack:** Python 3 (pandas), pytest, diamond (already a pipeline
dependency, invoked via `subprocess`), existing `nf_NovInvenio/lib/` modules
(`config_parser.py`, `hits.py`, `busco.py`).

**Spec:** `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/notes/superpowers/specs/2026-09-13-cluster-vs-pairwise-sensitivity-design.md`

## Global Constraints

- No changes to `nf_NovInvenio`'s Nextflow workflows/modules (`workflows/*.nf`,
  `modules/*.nf`) — everything here is `bin/` scripts, run by hand against
  already-published `results/` outputs.
- `score_controls.py`'s existing `hmm`-mode behavior must be byte-identical
  after this plan (default `--presence-mode hmm`, no argument changes to
  existing call sites) — every existing caller of this script keeps working
  unchanged.
- Tier R must not reuse Tier P's own diamond search results (`search_cache/
  *.diamond.tsv.gz`) even though they already exist for `pezizo_set1` — it
  runs its own small within-family diamond search, because the point of Tier
  R is to measure whether it's viable in a deployment where Tier P never ran.
- Tier R's ambiguous-family definition (verified against real data, not the
  earlier wrong assumption in the spec's first draft): a family is ambiguous
  iff (a) it is not in `oversized_families.tsv`, (b) at least one ingroup
  species contributes ≥2 members (species-duplication), and (c) its ingroup
  presence fraction is `>= --ingroup-min-frac` **and** its outgroup presence
  fraction is `> --other-max-frac` in the *existing* (Tier C+H) `presence_
  matrix.tsv` — i.e. it is already a near-miss novelty candidate. Do not
  reintroduce the "spans ingroup and outgroup at the raw-membership level"
  definition — verified false (mmseqs clusters the ingroup only).
- All new scripts follow this repo's existing test convention: import the
  `bin/` script as a module via `importlib.util.spec_from_file_location`
  (nf_NovInvenio side) or `sys.path.insert` (NovInvenio_Investigations side,
  matching `tests/test_build_study_config.py`), unit-test pure functions
  directly, plus one end-to-end subprocess test per script against small
  fixture files (no fixture ever depends on the real 11-genome data).

---

## File Structure

```
nf_NovInvenio/
  bin/score_controls.py              MODIFY — add --presence-mode
  tests/test_score_controls.py       MODIFY — cluster_membership-mode tests
  bin/refine_ambiguous_families.py   NEW — Tier R family-splitting script
  tests/test_refine_ambiguous_families.py  NEW

NovInvenio_Investigations/
  bin/generate_busco_map.py          NEW — busco_id\tprotein_id from a full_table.tsv
  tests/test_generate_busco_map.py   NEW
  bin/compare_cluster_tiers.py       NEW — Analysis 1 runner (4 tiers x N control sets)
  tests/test_compare_cluster_tiers.py NEW
  bin/genome_wide_concordance.py     NEW — Analysis 2 (concordance + inflation check)
  tests/test_genome_wide_concordance.py NEW
  bin/trace_cost_report.py           NEW — Analysis 3 (nextflow trace -> CPU-hours)
  tests/test_trace_cost_report.py    NEW
  studies/fungi/pezizo_set1_cluster/tier_comparison/   OUTPUT — real run results + report.md
```

---

### Task 1: `score_controls.py` — add `--presence-mode cluster_membership`

**Files:**
- Modify: `/bigdata/stajichlab/jstajich/projects/NovInvenio/bin/score_controls.py:216-286` (`score_controls()`, `main()`)
- Test: `/bigdata/stajichlab/jstajich/projects/NovInvenio/tests/test_score_controls.py`

**Interfaces:**
- Produces: `build_cluster_membership_presence(rep_to_members: dict[str, list[str]], protein_to_proteome: dict[str, str], proteome_cols: list[str]) -> dict[str, dict[str, int]]`
  — `rep -> {proteome_short: 0/1}`, used by later tasks (Tier C and Tier R
  scoring both call `score_controls()` with `presence_mode='cluster_membership'`).
- `score_controls(..., presence_mode='hmm', protein_to_proteome=None)` — new
  keyword args, default preserves current behavior exactly.

- [ ] **Step 1: Write the failing tests**

```python
# Append to nf_NovInvenio/tests/test_score_controls.py

def test_build_cluster_membership_presence_from_membership_alone():
    # fam rep pA1 has members pA1(In1), pA2(In2) -- present in In1+In2, absent Out1.
    # fam rep pB1 has members pB1(In1), pB1x(Out1) -- present in In1+Out1.
    rep_to_members = {'pA1': ['pA1', 'pA2'], 'pB1': ['pB1', 'pB1x']}
    protein_to_proteome = {'pA1': 'In1', 'pA2': 'In2', 'pB1': 'In1', 'pB1x': 'Out1'}
    cols = ['In1', 'In2', 'Out1']
    presence = sc.build_cluster_membership_presence(rep_to_members, protein_to_proteome, cols)
    assert presence == {
        'pA1': {'In1': 1, 'In2': 1, 'Out1': 0},
        'pB1': {'In1': 1, 'In2': 0, 'Out1': 1},
    }


def test_score_controls_cluster_membership_mode_ignores_matrix_presence_columns(tmp_path):
    # Matrix says pA1/pA2 are present everywhere (as if HMM over-called) but raw
    # cluster membership (pA-family has only In1/In2 members) must win in this mode.
    matrix = pd.read_csv(pd.io.common.StringIO(
        "protein_id\tsource_proteome\tIn1\tIn2\tOut1\n"
        "pA1\tIn1\t1\t1\t1\n"
        "pA2\tIn2\t1\t1\t1\n"
    ), sep='\t')
    member_to_rep = {'pA1': 'pA1', 'pA2': 'pA1'}
    rep_to_members = {'pA1': ['pA1', 'pA2']}
    from config_parser import Sample
    samples = [
        Sample(group='IN', species='s1', strain='', protein='', dna='', short='In1', taxon_group=''),
        Sample(group='IN', species='s2', strain='', protein='', dna='', short='In2', taxon_group=''),
        Sample(group='OUT', species='s3', strain='', protein='', dna='', short='Out1', taxon_group=''),
    ]
    controls = [{'control_id': 'POS1', 'class': 'positive', 'expected_call': 'novel',
                 'anchor_type': 'protein_id', 'anchor': 'pA1'}]
    protein_to_proteome = {'pA1': 'In1', 'pA2': 'In2'}
    results = sc.score_controls(
        controls, matrix, member_to_rep, rep_to_members, samples,
        ingroup_min_frac=0.75, other_max_frac=0.0, busco_map={}, controls_dir=tmp_path,
        profiles_hmm=None, cpus=1, presence_mode='cluster_membership',
        protein_to_proteome=protein_to_proteome,
    )
    assert results[0]['actual_call'] == 'novel'
    assert results[0]['outcome'] == 'hit'
```

- [ ] **Step 2: Run tests to verify they fail**

```
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
pixi run pytest tests/test_score_controls.py -k cluster_membership -v
```
Expected: FAIL — `AttributeError: module 'score_controls' has no attribute
'build_cluster_membership_presence'` (first test) and `TypeError:
score_controls() got an unexpected keyword argument 'presence_mode'` (second).

- [ ] **Step 3: Implement**

Replace `score_controls()` and the surrounding block (lines 216-341 —
`resolve_anchor` stays unchanged; only `score_controls()` and `main()` change)
in `bin/score_controls.py`:

```python
def build_cluster_membership_presence(rep_to_members, protein_to_proteome, proteome_cols):
    """rep -> {proteome_short: 0/1} from raw cluster membership alone (no HMM/matrix).

    Tier C (raw mmseqs membership) and Tier R (refined membership) both call
    score_controls() with presence_mode='cluster_membership' and this as the
    precomputed presence source -- the only difference between the two tiers is
    which --cluster-tsv/--families files were loaded into rep_to_members upstream.
    """
    presence = {}
    for rep, members in rep_to_members.items():
        proteomes = {protein_to_proteome[m] for m in members if m in protein_to_proteome}
        presence[rep] = {p: int(p in proteomes) for p in proteome_cols}
    return presence


def score_controls(controls, matrix, member_to_rep, rep_to_members, samples,
                   ingroup_min_frac, other_max_frac, busco_map, controls_dir,
                   profiles_hmm, cpus, presence_mode='hmm', protein_to_proteome=None):
    proteome_cols = [c for c in matrix.columns if c not in META_COLS]
    ingroup_ids = [s.short for s in samples if s.group in INGROUP_ROLES and s.short in proteome_cols]
    outgroup_ids = [s.short for s in samples if s.group in OUTGROUP_ROLES and s.short in proteome_cols]

    cluster_presence = None
    if presence_mode == 'cluster_membership':
        cluster_presence = build_cluster_membership_presence(
            rep_to_members, protein_to_proteome or {}, proteome_cols)

    results = []
    for row in controls:
        cid = (row.get('control_id') or '').strip()
        cls = (row.get('class') or '').strip().lower()
        expected = (row.get('expected_call') or '').strip().lower()
        rep, note = resolve_anchor(row, member_to_rep, busco_map, controls_dir,
                                   profiles_hmm, cpus)

        actual = 'unresolved'
        if rep is not None:
            if presence_mode == 'cluster_membership':
                presence = cluster_presence.get(rep)
                if presence is None:
                    rep, note = None, 'family has no cluster members'
            else:
                presence = family_presence_vector(matrix, rep_to_members.get(rep, []),
                                                  proteome_cols)
                if presence is None:
                    rep, note = None, 'family has no matrix rows'
            if presence is not None:
                actual = family_call(presence, ingroup_ids, outgroup_ids,
                                     ingroup_min_frac, other_max_frac)

        if actual == 'unresolved':
            outcome = 'unresolved'
        elif cls == 'positive':
            outcome = 'hit' if actual == 'novel' else 'miss'
        elif cls == 'negative':
            outcome = 'fp' if actual == 'novel' else 'tn'
        else:
            outcome = 'n/a'

        results.append({
            'control_id': cid,
            'class': cls,
            'expected_call': expected,
            'anchor_type': (row.get('anchor_type') or '').strip(),
            'resolved_family': rep or '',
            'actual_call': actual,
            'outcome': outcome,
            'note': note,
        })
    return results
```

In `main()`, add the argument and pass the new keywords (insert after the
existing `--paralog-competition-scope`-style block, i.e. right after the
`--cpus` argument, before `--output`):

```python
    ap.add_argument('--presence-mode', choices=['hmm', 'cluster_membership'],
                    default='hmm', dest='presence_mode',
                    help="'hmm' (default): read presence from --matrix (unchanged "
                         "behavior). 'cluster_membership': ignore --matrix's presence "
                         "columns and build presence purely from --cluster-tsv "
                         "membership (Tier C/Tier R scoring) -- requires --cluster-tsv/"
                         "--families (family mode only).")
```

and, in `main()`'s body, right after `member_to_rep, rep_to_members =
identity_membership(matrix)` / the `if args.cluster_tsv:` block:

```python
    if args.presence_mode == 'cluster_membership' and not args.cluster_tsv:
        sys.exit('--presence-mode cluster_membership requires --cluster-tsv/--families '
                 '(family mode only)')
    protein_to_proteome = dict(zip(matrix['protein_id'], matrix['source_proteome']))
```

and update the `score_controls(...)` call site to:

```python
    results = score_controls(
        controls, matrix, member_to_rep, rep_to_members, samples,
        args.ingroup_min_frac, args.other_max_frac, busco_map, controls_dir,
        args.profiles, args.cpus, presence_mode=args.presence_mode,
        protein_to_proteome=protein_to_proteome,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

```
pixi run pytest tests/test_score_controls.py -v
```
Expected: all PASS, including every pre-existing test (no regressions —
`presence_mode` defaults to `'hmm'` everywhere it isn't explicitly passed).

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add bin/score_controls.py tests/test_score_controls.py
git commit -m "$(cat <<'EOF'
Add --presence-mode cluster_membership to score_controls.py

Lets a run's family presence be scored from raw mmseqs cluster membership
alone (no HMM), instead of presence_matrix.tsv -- Tier C and Tier R of the
NovInvenio_Investigations cluster-vs-pairwise sensitivity investigation both
reuse this unchanged, differing only in which --cluster-tsv they point at.
Default --presence-mode hmm preserves existing behavior exactly.
EOF
)"
```

---

### Task 2: `refine_ambiguous_families.py` — ambiguous-family detection

**Files:**
- Create: `/bigdata/stajichlab/jstajich/projects/NovInvenio/bin/refine_ambiguous_families.py`
- Test: `/bigdata/stajichlab/jstajich/projects/NovInvenio/tests/test_refine_ambiguous_families.py`

**Interfaces:**
- Consumes: `nf_NovInvenio/lib/config_parser.py`'s `parse_config`, `INGROUP_ROLES`,
  `OUTGROUP_ROLES` (unchanged, already used by `score_controls.py`).
- Produces: `species_of(protein_id: str, protein_to_proteome: dict) -> str | None`;
  `has_species_duplication(members: list[str], protein_to_proteome: dict) -> bool`;
  `family_fractions(rep: str, members: list[str], protein_to_proteome: dict, ingroup_ids: set, outgroup_ids: set) -> tuple[float, float]`
  returning `(ingroup_frac, outgroup_frac)`;
  `detect_ambiguous_families(fam_members: dict[str, list[str]], protein_to_proteome: dict, ingroup_ids: set, outgroup_ids: set, oversized_reps: set, ingroup_min_frac: float, other_max_frac: float) -> set[str]`
  — used unchanged by Task 4's `main()`.

- [ ] **Step 1: Write the failing tests**

```python
# nf_NovInvenio/tests/test_refine_ambiguous_families.py
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'refine_ambiguous_families.py'
sys.path.insert(0, str(REPO / 'lib'))
_spec = importlib.util.spec_from_file_location('refine_ambiguous_families', BIN)
raf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(raf)


def test_has_species_duplication():
    p2p = {'a': 'In1', 'b': 'In2', 'c': 'In1'}
    assert raf.has_species_duplication(['a', 'b', 'c'], p2p)  # In1 appears twice
    assert not raf.has_species_duplication(['a', 'b'], p2p)


def test_family_fractions():
    p2p = {'a': 'In1', 'b': 'In2', 'c': 'Out1'}
    ing, out = raf.family_fractions('a', ['a', 'b', 'c'], p2p, {'In1', 'In2'}, {'Out1', 'Out2'})
    assert ing == 1.0   # both In1, In2 present
    assert out == 0.5   # Out1 present, Out2 absent


def test_detect_ambiguous_families_matches_hex1_ada1_split():
    # HEX1-like family: species-duplicated (In1 x2) AND a near-miss novelty candidate
    # (ingroup_frac=1.0 >= 0.75, outgroup_frac=1.0 > 0.0) -> ambiguous.
    # ADA1-like family: no duplication -> never ambiguous regardless of presence.
    fam_members = {
        'hex1_rep': ['hex1', 'eif5a', 'orth_in2'],  # In1, In1, In2
        'ada1_rep': ['ada1_in1', 'ada1_in2'],        # In1, In2
    }
    p2p = {'hex1': 'In1', 'eif5a': 'In1', 'orth_in2': 'In2',
           'ada1_in1': 'In1', 'ada1_in2': 'In2'}
    ingroup_ids, outgroup_ids = {'In1', 'In2'}, {'Out1'}
    # Fake presence: both families present in all ingroup + all outgroup (Out1).
    # (family_fractions is computed from cluster membership + a presence lookup in the
    # real pipeline via presence_matrix.tsv; here we inject it directly via a stub.)
    def fake_fractions(rep, members, p2p, ingroup_ids, outgroup_ids):
        return (1.0, 1.0)
    raf.family_fractions = fake_fractions
    ambiguous = raf.detect_ambiguous_families(
        fam_members, p2p, ingroup_ids, outgroup_ids, oversized_reps=set(),
        ingroup_min_frac=0.75, other_max_frac=0.0)
    assert ambiguous == {'hex1_rep'}


def test_detect_ambiguous_families_excludes_oversized():
    fam_members = {'big_rep': ['a', 'b', 'c']}
    p2p = {'a': 'In1', 'b': 'In1', 'c': 'In2'}
    def fake_fractions(rep, members, p2p, ingroup_ids, outgroup_ids):
        return (1.0, 1.0)
    raf.family_fractions = fake_fractions
    ambiguous = raf.detect_ambiguous_families(
        fam_members, p2p, {'In1', 'In2'}, {'Out1'}, oversized_reps={'big_rep'},
        ingroup_min_frac=0.75, other_max_frac=0.0)
    assert ambiguous == set()
```

- [ ] **Step 2: Run tests to verify they fail**

```
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
pixi run pytest tests/test_refine_ambiguous_families.py -v
```
Expected: FAIL — `refine_ambiguous_families.py` does not exist yet
(`FileNotFoundError` / import error from `spec_from_file_location`).

- [ ] **Step 3: Implement**

```python
#!/usr/bin/env python3
"""Tier R: split ambiguous mmseqs families using a targeted within-family diamond
search, so a small, cheap refinement pass can be measured against Tier C (raw
cluster membership) and Tier P (full pairwise) -- see notes/superpowers/specs/
2026-09-13-cluster-vs-pairwise-sensitivity-design.md for the full design and the
verified definition of "ambiguous family" this module implements.

A family is ambiguous iff:
  1. it is not in --oversized-families (excluded -- refinement cost is O(k^2) per
     family and a pathologically large cluster would dominate the whole run), AND
  2. at least one ingroup species contributes >=2 members (species-duplication --
     the only signal available before any HMM exists; confirmed exact on HEX1:
     Neurospora crassa contributes both HEX1_NEUCR and IF5A_NEUCR to one family), AND
  3. it is already a near-miss novelty candidate in the existing (Tier C+H)
     presence_matrix.tsv: ingroup presence fraction >= --ingroup-min-frac AND
     outgroup presence fraction > --other-max-frac (i.e. currently rejected only
     because of outgroup presence).

Ambiguous families are refined by a single within-family diamond all-vs-all (over
the union of every ambiguous family's members, not one diamond call per family --
see main()), then split via find_paralog_edges_to_cut() + connected_components(),
using each member's own genome-local registered paralog (from the already-published
self_hits/<Short>.paralog_cutoffs.tsv, one per genome -- no new self-search) as the
forced split boundary. Non-ambiguous families pass through unmodified so the output
has the exact same rep<TAB>member / families.tsv-index contract as the pipeline's
own families_cluster.tsv + families.tsv, and score_controls.py's existing
--presence-mode cluster_membership path scores it unchanged.
"""
import argparse
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import INGROUP_ROLES, OUTGROUP_ROLES, parse_config  # noqa: E402
from hits import PARSERS, open_input  # noqa: E402

DEFAULT_DIAMOND_EVALUE = 1e-5


# --------------------------------------------------------------------- ambiguity
def species_of(protein_id, protein_to_proteome):
    return protein_to_proteome.get(protein_id)


def has_species_duplication(members, protein_to_proteome):
    species = [species_of(m, protein_to_proteome) for m in members]
    species = [s for s in species if s is not None]
    return len(species) != len(set(species))


def family_fractions(rep, members, protein_to_proteome, ingroup_ids, outgroup_ids):
    """(ingroup_frac, outgroup_frac) present, from raw cluster membership -- mirrors
    build_cluster_membership_presence()'s notion of presence (member's own genome
    counts as present), i.e. this is Tier C's own presence call for this family,
    reused here purely to find near-miss novelty candidates worth refining."""
    proteomes = {species_of(m, protein_to_proteome) for m in members}
    proteomes.discard(None)
    ing = len(proteomes & ingroup_ids) / len(ingroup_ids) if ingroup_ids else 0.0
    out = len(proteomes & outgroup_ids) / len(outgroup_ids) if outgroup_ids else 0.0
    return ing, out


def detect_ambiguous_families(fam_members, protein_to_proteome, ingroup_ids, outgroup_ids,
                              oversized_reps, ingroup_min_frac, other_max_frac):
    ambiguous = set()
    for rep, members in fam_members.items():
        if rep in oversized_reps:
            continue
        if not has_species_duplication(members, protein_to_proteome):
            continue
        ing_frac, out_frac = family_fractions(rep, members, protein_to_proteome,
                                              ingroup_ids, outgroup_ids)
        if ing_frac >= ingroup_min_frac and out_frac > other_max_frac:
            ambiguous.add(rep)
    return ambiguous


# ---------------------------------------------------------------- graph splitting
def load_paralog_map(cutoff_files):
    """protein_ID -> paralog_protein_ID, from one or more per-genome
    self_hits/<Short>.paralog_cutoffs.tsv files (already published by Tier P --
    see parse_self_hits.py; no new self-search)."""
    paralog_of = {}
    for path in cutoff_files:
        with open(path) as fh:
            header = fh.readline()
            del header
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    paralog_of[parts[0]] = parts[1]
    return paralog_of


def build_within_family_edges(hit_lines, evalue_cutoff):
    """{frozenset({a, b}), ...} for every diamond hit pair passing evalue_cutoff."""
    edges = set()
    for hit in hit_lines:
        if hit.query_id == hit.target_id:
            continue
        if hit.evalue < evalue_cutoff:
            edges.add(frozenset({hit.query_id, hit.target_id}))
    return edges


def find_paralog_edges_to_cut(members, paralog_map):
    """Edges to remove before connected components: X-P where P is X's own
    registered within-genome paralog and P is a co-member of the same family --
    mirrors build_presence_matrix.py's target-scope paralog-competition test
    (the one that already rescues HEX-1 in Tier P), applied as a forced graph
    split instead of a per-hit disqualification."""
    member_set = set(members)
    cuts = set()
    for m in members:
        p = paralog_map.get(m)
        if p and p in member_set and p != m:
            cuts.add(frozenset({m, p}))
    return cuts


def connected_components(members, edges):
    """List of member-id lists, one per connected component (a lone member with no
    surviving edge becomes its own singleton component)."""
    parent = {m: m for m in members}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for edge in edges:
        a, b = tuple(edge)
        if a in parent and b in parent:
            union(a, b)

    groups = defaultdict(list)
    for m in members:
        groups[find(m)].append(m)
    return list(groups.values())


def split_family(members, edges, paralog_map):
    """members -> list of subfamily member-lists, after removing forced-split
    (registered-paralog) edges and taking connected components on what remains."""
    cuts = find_paralog_edges_to_cut(members, paralog_map)
    remaining = edges - cuts
    return connected_components(members, remaining)


if __name__ == '__main__':
    pass  # main() wired in Task 4
```

- [ ] **Step 4: Run tests to verify they pass**

```
pixi run pytest tests/test_refine_ambiguous_families.py -v
```
Expected: PASS (all 4 tests).

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add bin/refine_ambiguous_families.py tests/test_refine_ambiguous_families.py
git commit -m "$(cat <<'EOF'
Add refine_ambiguous_families.py's ambiguous-family detection + graph split

Pure functions only in this commit (detect_ambiguous_families,
has_species_duplication, family_fractions, split_family/connected_components,
find_paralog_edges_to_cut, load_paralog_map). Implements the verified
ambiguous-family definition from the cluster-vs-pairwise sensitivity spec:
species-duplicated AND a near-miss novelty candidate in the existing
presence_matrix.tsv, excluding oversized_families.tsv. Sequence extraction,
the real diamond invocation, and CLI wiring land in the next commit.
EOF
)"
```

---

### Task 3: `refine_ambiguous_families.py` — sequence extraction, diamond, CLI, output

**Files:**
- Modify: `/bigdata/stajichlab/jstajich/projects/NovInvenio/bin/refine_ambiguous_families.py` (append `main()` + helpers)
- Test: `/bigdata/stajichlab/jstajich/projects/NovInvenio/tests/test_refine_ambiguous_families.py` (append)

**Interfaces:**
- Consumes: `PARSERS['diamond']`, `open_input` from `lib/hits.py` (already imported
  in Task 2); `detect_ambiguous_families`, `split_family`, `load_paralog_map`
  (Task 2).
- Produces: `extract_needed_sequences(pep_paths: dict[str, Path], needed_ids: set[str]) -> dict[str, str]`
  (`protein_id -> sequence`); `write_fasta(seqs: dict[str, str], path: Path) -> None`;
  `write_refined_families(fam_members: dict[str, list[str]], ambiguous_reps: set[str], subfamilies_by_rep: dict[str, list[list[str]]], out_cluster_tsv: Path, out_families_tsv: Path) -> None`
  — the two output files, in the exact same shape as `families_cluster.tsv`/
  `families.tsv`, consumed unchanged by `score_controls.py`'s `--cluster-tsv`/
  `--families`.

- [ ] **Step 1: Write the failing tests**

```python
# append to nf_NovInvenio/tests/test_refine_ambiguous_families.py

def test_extract_needed_sequences(tmp_path):
    fa1 = tmp_path / 'In1.pep.fa'
    fa1.write_text(">a desc one\nMKV\n>b desc two\nMKL\n")
    fa2 = tmp_path / 'In2.pep.fa'
    fa2.write_text(">c desc three\nMKA\n>d desc four\nMKD\n")
    seqs = raf.extract_needed_sequences({'In1': fa1, 'In2': fa2}, needed_ids={'a', 'd'})
    assert seqs == {'a': 'MKV', 'd': 'MKD'}


def test_write_fasta(tmp_path):
    out = tmp_path / 'out.fa'
    raf.write_fasta({'a': 'MKV', 'b': 'MKL'}, out)
    text = out.read_text()
    assert text == ">a\nMKV\n>b\nMKL\n"


def test_write_refined_families_passthrough_and_split(tmp_path):
    fam_members = {'rep1': ['x', 'y'], 'hex1_rep': ['hex1', 'eif5a', 'orth_in2']}
    ambiguous = {'hex1_rep'}
    subfamilies = {'hex1_rep': [['hex1', 'orth_in2'], ['eif5a']]}
    out_cluster = tmp_path / 'refined_cluster.tsv'
    out_families = tmp_path / 'refined_families.tsv'
    raf.write_refined_families(fam_members, ambiguous, subfamilies, out_cluster, out_families)

    cluster_lines = set(out_cluster.read_text().splitlines())
    # rep1 passes through unchanged (not ambiguous).
    assert 'rep1\tx' in cluster_lines
    assert 'rep1\ty' in cluster_lines
    # hex1_rep is split: each subfamily keyed by its own first member as new rep.
    assert 'hex1\thex1' in cluster_lines
    assert 'hex1\torth_in2' in cluster_lines
    assert 'eif5a\teif5a' in cluster_lines
    assert 'hex1_rep\thex1' not in cluster_lines  # old rep gone for the split family

    families_lines = out_families.read_text().splitlines()
    assert families_lines[0] == 'family_index\trepresentative_id\tn_members'
    body = {line.split('\t')[1]: line.split('\t')[2] for line in families_lines[1:]}
    assert body['rep1'] == '2'
    assert body['hex1'] == '2'
    assert body['eif5a'] == '1'
```

- [ ] **Step 2: Run tests to verify they fail**

```
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
pixi run pytest tests/test_refine_ambiguous_families.py -k "extract_needed or write_fasta or write_refined" -v
```
Expected: FAIL — `AttributeError: module has no attribute 'extract_needed_sequences'`
(and similarly for the other two).

- [ ] **Step 3: Implement**

Append to `bin/refine_ambiguous_families.py` (replacing the `if __name__ ==
'__main__': pass` stub from Task 2):

```python
# ------------------------------------------------------------------- sequences
def extract_needed_sequences(pep_paths, needed_ids):
    """protein_id -> sequence, scanning each --Short.pep.fa once, keeping only ids
    in needed_ids (the union of every ambiguous family's members) -- avoids loading
    whole proteomes when only a small fraction of their proteins are needed."""
    seqs = {}
    remaining = set(needed_ids)
    for path in pep_paths.values():
        if not remaining:
            break
        current_id, chunks = None, []
        with open(path) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if current_id in remaining:
                        seqs[current_id] = ''.join(chunks)
                        remaining.discard(current_id)
                    current_id = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)
            if current_id in remaining:
                seqs[current_id] = ''.join(chunks)
                remaining.discard(current_id)
    return seqs


def write_fasta(seqs, path):
    with open(path, 'w') as fh:
        for pid, seq in seqs.items():
            fh.write(f'>{pid}\n{seq}\n')


def run_diamond_within_family(fasta_path, cpus=1):
    """One diamond makedb + blastp all-vs-all over the combined ambiguous-family
    FASTA (not per-family -- process-startup overhead would dominate at ~1000+
    families). Returns the path to the raw outfmt-6 hits file."""
    db_path = fasta_path.with_suffix('.dmnd')
    subprocess.run(['diamond', 'makedb', '--in', str(fasta_path), '--db', str(db_path)],
                    check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    hits_path = fasta_path.with_suffix('.hits.tsv')
    subprocess.run([
        'diamond', 'blastp', '--very-sensitive', '--threads', str(cpus),
        '--query', str(fasta_path), '--db', str(db_path),
        '--outfmt', '6', 'qseqid', 'sseqid', 'evalue', 'bitscore',
        '--out', str(hits_path),
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return hits_path


# --------------------------------------------------------------------- output
def write_refined_families(fam_members, ambiguous_reps, subfamilies_by_rep,
                           out_cluster_tsv, out_families_tsv):
    """Same rep<TAB>member / families.tsv-index shape as the pipeline's own
    families_cluster.tsv + families.tsv, so score_controls.py's existing
    --cluster-tsv/--families loaders consume this file unmodified. Non-ambiguous
    families pass through verbatim; ambiguous ones are replaced by their split
    subfamilies, each keyed by its own first member as the new representative."""
    with open(out_cluster_tsv, 'w') as cfh, open(out_families_tsv, 'w') as ffh:
        ffh.write('family_index\trepresentative_id\tn_members\n')
        idx = 0
        for rep, members in fam_members.items():
            if rep not in ambiguous_reps:
                idx += 1
                ffh.write(f'fam_{idx:06d}\t{rep}\t{len(members)}\n')
                for m in members:
                    cfh.write(f'{rep}\t{m}\n')
                continue
            for subfamily in subfamilies_by_rep.get(rep, [members]):
                new_rep = subfamily[0]
                idx += 1
                ffh.write(f'fam_{idx:06d}\t{new_rep}\t{len(subfamily)}\n')
                for m in subfamily:
                    cfh.write(f'{new_rep}\t{m}\n')


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                    help="run's families_cluster.tsv (rep<TAB>member)")
    ap.add_argument('--families', required=True,
                    help="run's families.tsv (profiled families index)")
    ap.add_argument('--matrix', required=True,
                    help="run's presence_matrix.tsv (for protein_id -> source_proteome "
                         "and the near-miss-novelty check)")
    ap.add_argument('--oversized-families', default=None, dest='oversized_families',
                    help="run's oversized_families.tsv (excluded from refinement)")
    ap.add_argument('--config', required=True, help='Analysis description CSV')
    ap.add_argument('--pep', nargs='+', required=True,
                    help='Short=path.pep.fa pairs for every ingroup proteome')
    ap.add_argument('--self-hits', nargs='+', required=True, dest='self_hits',
                    help='self_hits/<Short>.paralog_cutoffs.tsv files (already '
                         'published by the pairwise run for these same genomes)')
    ap.add_argument('--ingroup-min-frac', type=float, default=0.75, dest='ingroup_min_frac')
    ap.add_argument('--other-max-frac', type=float, default=0.0, dest='other_max_frac')
    ap.add_argument('--diamond-evalue', type=float, default=DEFAULT_DIAMOND_EVALUE,
                    dest='diamond_evalue')
    ap.add_argument('--cpus', type=int, default=1)
    ap.add_argument('--tmp-dir', default=None, dest='tmp_dir')
    ap.add_argument('--output-cluster-tsv', required=True, dest='output_cluster_tsv')
    ap.add_argument('--output-families', required=True, dest='output_families')
    args = ap.parse_args()

    import pandas as pd  # local import: only main() needs it, unlike the pure fns above

    samples = parse_config(args.config)
    ingroup_ids = {s.short for s in samples if s.group in INGROUP_ROLES}
    outgroup_ids = {s.short for s in samples if s.group in OUTGROUP_ROLES}
    pep_paths = {}
    for pair in args.pep:
        short, path = pair.split('=', 1)
        pep_paths[short] = Path(path)

    matrix = pd.read_csv(args.matrix, sep='\t')
    protein_to_proteome = dict(zip(matrix['protein_id'], matrix['source_proteome']))

    fam_members = defaultdict(list)
    with open(args.cluster_tsv) as fh:
        for line in fh:
            rep, member = line.rstrip('\n').split('\t')[:2]
            fam_members[rep].append(member)
    fam_members = dict(fam_members)

    oversized_reps = set()
    if args.oversized_families:
        with open(args.oversized_families) as fh:
            header = fh.readline()
            del header
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if parts and parts[0]:
                    oversized_reps.add(parts[0])

    ambiguous = detect_ambiguous_families(
        fam_members, protein_to_proteome, ingroup_ids, outgroup_ids,
        oversized_reps, args.ingroup_min_frac, args.other_max_frac)
    print(f'{len(ambiguous)} / {len(fam_members)} families flagged ambiguous', file=sys.stderr)

    paralog_map = load_paralog_map(args.self_hits)

    tmp_dir = Path(args.tmp_dir) if args.tmp_dir else Path(tempfile.mkdtemp())
    tmp_dir.mkdir(parents=True, exist_ok=True)
    needed_ids = {m for rep in ambiguous for m in fam_members[rep]}
    seqs = extract_needed_sequences(pep_paths, needed_ids)
    combined_fasta = tmp_dir / 'ambiguous_family_members.fa'
    write_fasta(seqs, combined_fasta)

    subfamilies_by_rep = {}
    if ambiguous:
        hits_path = run_diamond_within_family(combined_fasta, cpus=args.cpus)
        with open_input(hits_path) as fh:
            all_edges = build_within_family_edges(PARSERS['diamond'](fh), args.diamond_evalue)
        for rep in ambiguous:
            members = fam_members[rep]
            member_set = set(members)
            family_edges = {e for e in all_edges if e <= member_set}
            subfamilies_by_rep[rep] = split_family(members, family_edges, paralog_map)

    write_refined_families(fam_members, ambiguous, subfamilies_by_rep,
                           args.output_cluster_tsv, args.output_families)


if __name__ == '__main__':
    main()
```

Also remove the now-superseded `if __name__ == '__main__': pass` line left
over from Task 2 (replaced by the block above).

- [ ] **Step 4: Run tests to verify they pass**

```
pixi run pytest tests/test_refine_ambiguous_families.py -v
```
Expected: all PASS (7 tests total across Tasks 2-3).

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio
git add bin/refine_ambiguous_families.py tests/test_refine_ambiguous_families.py
git commit -m "$(cat <<'EOF'
Wire refine_ambiguous_families.py's diamond step, CLI, and output writer

Completes Tier R: extract only the sequences ambiguous families actually
need, run one combined within-family diamond all-vs-all (not one process per
family), split each ambiguous family on its self_hits-registered paralog
edges, and emit a families_cluster.tsv/families.tsv-shaped refined output
that score_controls.py's --presence-mode cluster_membership path consumes
unchanged.
EOF
)"
```

---

### Task 4: `generate_busco_map.py` — fill the `--busco-map` gap

**Files:**
- Create: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/bin/generate_busco_map.py`
- Test: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/tests/test_generate_busco_map.py`

**Interfaces:**
- Consumes: `nf_NovInvenio/lib/busco.py`'s `parse_busco_full_table(path, species) -> Iterator[tuple[str, str, str, int|None]]`
  (already exists, unchanged).
- Produces: `filter_wanted(rows: list[tuple], wanted_ids: set[str]) -> dict[str, str]`
  (`busco_id -> protein_id`, first Complete row wins per id); a 2-column
  `busco_id\tprotein_id` TSV matching `score_controls.py`'s `load_busco_map()`
  format exactly.

- [ ] **Step 1: Write the failing test**

```python
# NovInvenio_Investigations/tests/test_generate_busco_map.py
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
sys.path.insert(0, '/bigdata/stajichlab/jstajich/projects/NovInvenio/lib')
import generate_busco_map as gbm  # noqa: E402


def test_filter_wanted_keeps_only_requested_ids():
    rows = [
        ('100036at4751', 'Ncra', 'protA', 500),
        ('999999at4751', 'Ncra', 'protZ', 200),  # not wanted -- must be dropped
        ('100149at4751', 'Ncra', 'protB', 300),
    ]
    wanted = {'100036at4751', '100149at4751'}
    mapping = gbm.filter_wanted(rows, wanted)
    assert mapping == {'100036at4751': 'protA', '100149at4751': 'protB'}


def test_main_writes_two_column_tsv(tmp_path):
    full_table = tmp_path / 'full_table.tsv'
    full_table.write_text(
        "# BUSCO\n"
        "100036at4751\tComplete\tprotA:1-100\t0\t500\n"
        "999999at4751\tComplete\tprotZ:1-50\t0\t200\n"
    )
    wanted_csv = tmp_path / 'controls.csv'
    wanted_csv.write_text(
        "control_id,class,expected_call,anchor_type,anchor,proteome_short,gene_name,expected_origin,source,notes\n"
        "NEG_BUSCO01,negative,core,busco,100036at4751,,x,y,z,note\n"
    )
    out = tmp_path / 'out.tsv'
    import subprocess
    subprocess.run([
        sys.executable, str(REPO / 'bin' / 'generate_busco_map.py'),
        '--full-table', str(full_table), '--species', 'Ncra',
        '--controls', str(wanted_csv), '--output', str(out),
    ], check=True)
    assert out.read_text() == '100036at4751\tprotA\n'
```

- [ ] **Step 2: Run test to verify it fails**

```
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
pixi run pytest tests/test_generate_busco_map.py -v
```
Expected: FAIL — `ModuleNotFoundError: No module named 'generate_busco_map'`.

- [ ] **Step 3: Implement**

```python
#!/usr/bin/env python3
"""Build the 2-column busco_id<TAB>protein_id map nf_NovInvenio/bin/score_controls.py's
--busco-map expects, from one species' BUSCO full_table.tsv, restricted to the busco
anchor ids an actual controls CSV references -- fixes the gap where score_controls.py
was run without --busco-map at all, silently reporting every busco-anchored negative
control as unresolved instead of scoring it (see notes/superpowers/specs/2026-09-13-
cluster-vs-pairwise-sensitivity-design.md, Analysis 1).
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, '/bigdata/stajichlab/jstajich/projects/NovInvenio/lib')
from busco import parse_busco_full_table  # noqa: E402


def filter_wanted(rows, wanted_ids):
    """First Complete row per busco_id, restricted to wanted_ids. rows: iterable of
    (busco_id, species, protein_id, length)."""
    mapping = {}
    for busco_id, _species, protein_id, _length in rows:
        if busco_id in wanted_ids and busco_id not in mapping:
            mapping[busco_id] = protein_id
    return mapping


def load_wanted_busco_ids(controls_csv):
    wanted = set()
    with open(controls_csv) as fh:
        for row in csv.DictReader(fh):
            if (row.get('anchor_type') or '').strip() == 'busco':
                wanted.add((row.get('anchor') or '').strip())
    return wanted


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--full-table', required=True, dest='full_table',
                    help="one species' BUSCO full_table.tsv")
    ap.add_argument('--species', required=True, help='label passed through to parse_busco_full_table')
    ap.add_argument('--controls', required=True, help='controls CSV to restrict ids to')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    wanted = load_wanted_busco_ids(args.controls)
    rows = list(parse_busco_full_table(args.full_table, args.species))
    mapping = filter_wanted(rows, wanted)

    missing = wanted - mapping.keys()
    if missing:
        print(f'WARNING: {len(missing)} busco id(s) not Complete in {args.species}: '
              f'{sorted(missing)}', file=sys.stderr)

    with open(args.output, 'w') as fh:
        for busco_id in sorted(mapping):
            fh.write(f'{busco_id}\t{mapping[busco_id]}\n')


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Run test to verify it passes**

```
pixi run pytest tests/test_generate_busco_map.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add bin/generate_busco_map.py tests/test_generate_busco_map.py
git commit -m "$(cat <<'EOF'
Add generate_busco_map.py to fill score_controls.py's --busco-map gap

Reuses nf_NovInvenio/lib/busco.py's parse_busco_full_table against an
existing BUSCO full_table.tsv, restricted to the busco anchor ids a given
controls CSV actually references, so all 10 BUSCO negatives can be scored
instead of 5 silently reporting unresolved for want of this file.
EOF
)"
```

---

### Task 5: `compare_cluster_tiers.py` + real run on `pezizo_set1`'s 16 controls

**Files:**
- Create: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/bin/compare_cluster_tiers.py`
- Test: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/tests/test_compare_cluster_tiers.py`
- Output (real run, not a fixture): `studies/fungi/pezizo_set1_cluster/tier_comparison/pezizo_set1.per_control.tsv`,
  `.../pezizo_set1.tier_summary.tsv`

**Interfaces:**
- Consumes: `score_controls.py` as a subprocess (all four tiers below), and
  `refine_ambiguous_families.py` as a subprocess (Tier R only).
- Produces: `run_score_controls(*, matrix, controls, config, output, cluster_tsv=None, families=None, presence_mode='hmm', busco_map=None, extra_args=()) -> None`
  (thin subprocess wrapper); `parse_summary_tsv(path: Path) -> dict[str, str]`;
  `build_tier_comparison(tier_summaries: dict[str, dict], tier_per_control: dict[str, list[dict]]) -> tuple[list[dict], list[dict]]`
  → `(per_control_rows, tier_summary_rows)`, each row tagged with a `tier`
  column, ready to write straight to TSV.

- [ ] **Step 1: Write the failing tests**

```python
# NovInvenio_Investigations/tests/test_compare_cluster_tiers.py
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
import compare_cluster_tiers as cct  # noqa: E402


def test_parse_summary_tsv(tmp_path):
    p = tmp_path / 'summary.tsv'
    p.write_text('metric\tvalue\nrecall\t0.4\nfp_rate\t0.0\n')
    assert cct.parse_summary_tsv(p) == {'recall': '0.4', 'fp_rate': '0.0'}


def test_build_tier_comparison_tags_rows_with_tier():
    tier_summaries = {
        'P': {'recall': '1.0', 'fp_rate': '0.0'},
        'C+H': {'recall': '0.4', 'fp_rate': '0.0'},
    }
    tier_per_control = {
        'P': [{'control_id': 'POS_HEX1', 'outcome': 'hit'}],
        'C+H': [{'control_id': 'POS_HEX1', 'outcome': 'miss'}],
    }
    per_control, summaries = cct.build_tier_comparison(tier_summaries, tier_per_control)
    assert {'control_id': 'POS_HEX1', 'outcome': 'hit', 'tier': 'P'} in per_control
    assert {'control_id': 'POS_HEX1', 'outcome': 'miss', 'tier': 'C+H'} in per_control
    assert {'tier': 'P', 'recall': '1.0', 'fp_rate': '0.0'} in summaries
    assert {'tier': 'C+H', 'recall': '0.4', 'fp_rate': '0.0'} in summaries
```

- [ ] **Step 2: Run tests to verify they fail**

```
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
pixi run pytest tests/test_compare_cluster_tiers.py -v
```
Expected: FAIL — module does not exist.

- [ ] **Step 3: Implement**

```python
#!/usr/bin/env python3
"""Analysis 1: run pairwise (Tier P), cluster+HMM (Tier C+H), cluster-membership-only
(Tier C), and refined-cluster (Tier R) scoring through nf_NovInvenio/bin/score_controls.py
against one controls CSV, and assemble one per-control + one per-tier-summary table --
see notes/superpowers/specs/2026-09-13-cluster-vs-pairwise-sensitivity-design.md.
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

SCORE_CONTROLS = '/bigdata/stajichlab/jstajich/projects/NovInvenio/bin/score_controls.py'


def parse_summary_tsv(path):
    with open(path) as fh:
        r = csv.reader(fh, delimiter='\t')
        next(r)  # header
        return {row[0]: row[1] for row in r if row}


def parse_per_control_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def run_score_controls(*, matrix, controls, config, output, cluster_tsv=None,
                       families=None, presence_mode='hmm', busco_map=None, cpus=1):
    cmd = [sys.executable, SCORE_CONTROLS, '--controls', controls, '--matrix', matrix,
           '--config', config, '--output', str(output), '--presence-mode', presence_mode,
           '--cpus', str(cpus)]
    if cluster_tsv:
        cmd += ['--cluster-tsv', cluster_tsv, '--families', families]
    if busco_map:
        cmd += ['--busco-map', busco_map]
    subprocess.run(cmd, check=True)


def build_tier_comparison(tier_summaries, tier_per_control):
    per_control_rows = []
    for tier, rows in tier_per_control.items():
        for row in rows:
            per_control_rows.append({**row, 'tier': tier})
    summary_rows = [{'tier': tier, **summary} for tier, summary in tier_summaries.items()]
    return per_control_rows, summary_rows


def write_tsv(rows, path):
    if not rows:
        Path(path).write_text('')
        return
    fieldnames = list(rows[0].keys())
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--controls', required=True)
    ap.add_argument('--config', required=True)
    ap.add_argument('--pairwise-matrix', required=True, dest='pairwise_matrix')
    ap.add_argument('--cluster-hmm-matrix', required=True, dest='cluster_hmm_matrix')
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv')
    ap.add_argument('--families', required=True)
    ap.add_argument('--refined-cluster-tsv', required=True, dest='refined_cluster_tsv')
    ap.add_argument('--refined-families', required=True, dest='refined_families')
    ap.add_argument('--busco-map', default=None, dest='busco_map')
    ap.add_argument('--out-dir', required=True, dest='out_dir')
    ap.add_argument('--label', required=True, help='output filename prefix, e.g. pezizo_set1')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tiers = {
        'P': dict(matrix=args.pairwise_matrix, cluster_tsv=None, families=None, presence_mode='hmm'),
        'C+H': dict(matrix=args.cluster_hmm_matrix, cluster_tsv=args.cluster_tsv,
                    families=args.families, presence_mode='hmm'),
        'C': dict(matrix=args.cluster_hmm_matrix, cluster_tsv=args.cluster_tsv,
                  families=args.families, presence_mode='cluster_membership'),
        'R': dict(matrix=args.cluster_hmm_matrix, cluster_tsv=args.refined_cluster_tsv,
                  families=args.refined_families, presence_mode='cluster_membership'),
    }

    tier_summaries, tier_per_control = {}, {}
    for tier, cfg in tiers.items():
        output = out_dir / f'{args.label}.{tier.replace("+", "")}.controls_scored.tsv'
        run_score_controls(matrix=cfg['matrix'], controls=args.controls, config=args.config,
                           output=output, cluster_tsv=cfg['cluster_tsv'],
                           families=cfg['families'], presence_mode=cfg['presence_mode'],
                           busco_map=args.busco_map)
        summary_path = output.with_suffix('').with_suffix('.summary.tsv')
        tier_summaries[tier] = parse_summary_tsv(summary_path)
        tier_per_control[tier] = parse_per_control_tsv(output)

    per_control_rows, summary_rows = build_tier_comparison(tier_summaries, tier_per_control)
    write_tsv(per_control_rows, out_dir / f'{args.label}.per_control.tsv')
    write_tsv(summary_rows, out_dir / f'{args.label}.tier_summary.tsv')
    print(f'Wrote {out_dir}/{args.label}.per_control.tsv and .tier_summary.tsv',
          file=sys.stderr)


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```
pixi run pytest tests/test_compare_cluster_tiers.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add bin/compare_cluster_tiers.py tests/test_compare_cluster_tiers.py
git commit -m "$(cat <<'EOF'
Add compare_cluster_tiers.py: run Tiers P/C+H/C/R through score_controls.py

Orchestrates the four presence-calling tiers against one controls CSV and
assembles a per-control and per-tier-summary TSV -- Analysis 1 of the
cluster-vs-pairwise sensitivity investigation.
EOF
)"
```

- [ ] **Step 6: Run for real against `pezizo_set1`'s 16 controls**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations

# 6a. Generate the busco map (Task 4), fixing the 5-of-10-unresolved gap.
python3 bin/generate_busco_map.py \
  --full-table /bigdata/stajichlab/jstajich/projects/NovInvenio/busco_pezizo5/Ncra.busco/run_fungi_odb12/full_table.tsv \
  --species Ncra \
  --controls /bigdata/stajichlab/jstajich/projects/NovInvenio/configs/controls/pezizo_set1.controls.csv \
  --output studies/fungi/pezizo_set1_cluster/tier_comparison/pezizo_set1.busco_map.tsv
# Expect stderr to report 0 missing ids (10/10 covered) -- if not, stop and
# investigate before proceeding; a partial map would silently understate FP-rate
# coverage the same way the original gap did.

# 6b. Build Tier R's refined families (gain side only for this task; loss side
# is Task 7's job, reusing the same script against loss_families/).
python3 /bigdata/stajichlab/jstajich/projects/NovInvenio/bin/refine_ambiguous_families.py \
  --cluster-tsv results/pezizo_set1_cluster/families/families_cluster.tsv \
  --families results/pezizo_set1_cluster/families/families.tsv \
  --matrix results/pezizo_set1_cluster/presence_matrix.tsv \
  --oversized-families results/pezizo_set1_cluster/families/oversized_families.tsv \
  --config studies/fungi/pezizo_set1/config.csv \
  --pep Afum=studies/fungi/pezizo_set1/data_dir/pep/Afum.pep.fa \
        Amega=studies/fungi/pezizo_set1/data_dir/pep/Amega.pep.fa \
        Ztri=studies/fungi/pezizo_set1/data_dir/pep/Ztri.pep.fa \
        Ncra=studies/fungi/pezizo_set1/data_dir/pep/Ncra.pep.fa \
        Cimm=studies/fungi/pezizo_set1/data_dir/pep/Cimm.pep.fa \
  --self-hits results/pezizo_set1/self_hits/*.paralog_cutoffs.tsv \
  --output-cluster-tsv studies/fungi/pezizo_set1_cluster/tier_comparison/gain.refined_cluster.tsv \
  --output-families studies/fungi/pezizo_set1_cluster/tier_comparison/gain.refined_families.tsv

# 6c. Run all four tiers.
python3 bin/compare_cluster_tiers.py \
  --controls /bigdata/stajichlab/jstajich/projects/NovInvenio/configs/controls/pezizo_set1.controls.csv \
  --config studies/fungi/pezizo_set1/config.csv \
  --pairwise-matrix results/pezizo_set1/presence_matrix.tsv \
  --cluster-hmm-matrix results/pezizo_set1_cluster/presence_matrix.tsv \
  --cluster-tsv results/pezizo_set1_cluster/families/families_cluster.tsv \
  --families results/pezizo_set1_cluster/families/families.tsv \
  --refined-cluster-tsv studies/fungi/pezizo_set1_cluster/tier_comparison/gain.refined_cluster.tsv \
  --refined-families studies/fungi/pezizo_set1_cluster/tier_comparison/gain.refined_families.tsv \
  --busco-map studies/fungi/pezizo_set1_cluster/tier_comparison/pezizo_set1.busco_map.tsv \
  --out-dir studies/fungi/pezizo_set1_cluster/tier_comparison \
  --label pezizo_set1
```

Inspect `studies/fungi/pezizo_set1_cluster/tier_comparison/pezizo_set1.tier_summary.tsv`
and `.per_control.tsv`. Confirm, against the spec's already-verified facts:
- Tier P: recall 1.0 (6/6), fp_rate 0.0 (10/10 now resolved, not 5).
- Tier C+H: recall 0.4 (matches the pre-existing number), fp_rate 0.0.
- Tier C: `POS_ADA1` and `POS_HAM5` now score `hit` (raw cluster membership
  was already correct for both — see spec's Analysis 1 prediction);
  `POS_HEX1` still `miss` (clustering itself merged the paralog in); `POS_LAH`
  still `unresolved` (singleton, unaffected by any presence-calling tier).
- Tier R: `POS_HEX1` should flip to `hit` if the paralog-competition split
  worked (Ncra's own `eif5a`/`hex1` registered-paralog pair is a co-member of
  that family, so the forced-split edge should separate them) — if it does
  not, **do not edit the test suite to make it pass**; inspect
  `studies/fungi/pezizo_set1_cluster/tier_comparison/pezizo_set1.R.controls_scored.tsv`'s
  `note` column for `POS_HEX1` and the intermediate `*.hits.tsv`/refined
  cluster TSV to find why (e.g. the diamond hit between HEX1 and eIF-5A
  itself didn't pass `--diamond-evalue`, or a third member's edge kept them
  connected through a different path) and report the actual finding — this
  is real data, not a fixture, and an unexpected result here is exactly the
  kind of finding the investigation exists to surface.

- [ ] **Step 7: Commit the real run's output**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add studies/fungi/pezizo_set1_cluster/tier_comparison/
git commit -m "$(cat <<'EOF'
Run Tiers P/C+H/C/R against pezizo_set1's 16 controls

First real-data run of the cluster-vs-pairwise sensitivity investigation's
Analysis 1. See tier_comparison/pezizo_set1.tier_summary.tsv and
.per_control.tsv for the recall/FP numbers and per-control outcomes.
EOF
)"
```

---

### Task 6: Extend to the Agaricales/agaricomycetes control set

**Files:**
- Modify: none (Task 5's scripts are already parametrized by controls CSV /
  config / matrix paths)
- Output: `studies/fungi/pezizo_set1_cluster/tier_comparison/agaricomycetes.per_control.tsv`,
  `.../agaricomycetes.tier_summary.tsv`, plus the Agaricales-side refined
  families files.

- [ ] **Step 1: Build Tier R's refined families for the agaricomycetes set**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
python3 /bigdata/stajichlab/jstajich/projects/NovInvenio/bin/refine_ambiguous_families.py \
  --cluster-tsv results/agaricomycetes_mmseqs/families/families_cluster.tsv \
  --families results/agaricomycetes_mmseqs/families/families.tsv \
  --matrix results/agaricomycetes_mmseqs/presence_matrix.tsv \
  --oversized-families results/agaricomycetes_mmseqs/families/oversized_families.tsv \
  --config studies/fungi/agaricomycetes_pairwise/config.csv \
  --pep Scom=studies/fungi/agaricomycetes_pairwise/data_dir/pep/Scom.pep.fa \
        Ccin=studies/fungi/agaricomycetes_pairwise/data_dir/pep/Ccin.pep.fa \
        Agbi=studies/fungi/agaricomycetes_pairwise/data_dir/pep/Agbi.pep.fa \
        Lbic=studies/fungi/agaricomycetes_pairwise/data_dir/pep/Lbic.pep.fa \
  --self-hits results/agaricomycetes_pairwise/self_hits/*.paralog_cutoffs.tsv \
  --output-cluster-tsv studies/fungi/pezizo_set1_cluster/tier_comparison/agaricomycetes.gain.refined_cluster.tsv \
  --output-families studies/fungi/pezizo_set1_cluster/tier_comparison/agaricomycetes.gain.refined_families.tsv
```
If any `--pep` path or `--self-hits` glob does not resolve (check with `ls`
first — the agaricomycetes study's `data_dir`/ingroup species list may not
exactly match `pezizo_set1`'s), correct the paths from
`studies/fungi/agaricomycetes_pairwise/config.csv`'s own `IN` rows rather
than guessing; do not skip species silently.

- [ ] **Step 2: Run all four tiers against `Agaricales.controls.csv`**

```bash
python3 bin/compare_cluster_tiers.py \
  --controls /bigdata/stajichlab/jstajich/projects/NovInvenio/configs/controls/Agaricales.controls.csv \
  --config studies/fungi/agaricomycetes_pairwise/config.csv \
  --pairwise-matrix results/agaricomycetes_pairwise/presence_matrix.tsv \
  --cluster-hmm-matrix results/agaricomycetes_mmseqs/presence_matrix.tsv \
  --cluster-tsv results/agaricomycetes_mmseqs/families/families_cluster.tsv \
  --families results/agaricomycetes_mmseqs/families/families.tsv \
  --refined-cluster-tsv studies/fungi/pezizo_set1_cluster/tier_comparison/agaricomycetes.gain.refined_cluster.tsv \
  --refined-families studies/fungi/pezizo_set1_cluster/tier_comparison/agaricomycetes.gain.refined_families.tsv \
  --out-dir studies/fungi/pezizo_set1_cluster/tier_comparison \
  --label agaricomycetes
```
Note: no `--busco-map` here — `Agaricales.controls.csv` has no `busco`-typed
rows (both controls are `fasta` anchors, spc14/spc33), so the argument is
simply omitted (matches `compare_cluster_tiers.py`'s `default=None` handling,
no code change needed).

Inspect the output. With only 2 positive controls and 0 negatives, recall is
either 0.0, 0.5, or 1.0 per tier — record whichever it is verbatim; do not
round to whatever "looks like a clean result."

- [ ] **Step 3: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add studies/fungi/pezizo_set1_cluster/tier_comparison/agaricomycetes*
git commit -m "$(cat <<'EOF'
Run Tiers P/C+H/C/R against the Agaricales spc14/spc33 controls

Second, independent ingroup/outgroup split for Analysis 1 (agaricomycetes_
pairwise vs agaricomycetes_mmseqs), so the pezizo_set1 findings aren't a
one-clade artifact.
EOF
)"
```

---

### Task 7: `genome_wide_concordance.py` — Analysis 2 (gains and losses)

**Files:**
- Create: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/bin/genome_wide_concordance.py`
- Test: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/tests/test_genome_wide_concordance.py`

**Interfaces:**
- Produces: `load_candidates(path: Path) -> set[str]` (one `source_proteome::
  protein_id` id per line, matching `candidates.txt`/`loss_candidates.txt`'s
  existing format — both `build_presence_matrix.py` and `profile_to_matrix.py`
  already emit this exact shape); `jaccard(a: set, b: set) -> float`;
  `precision_recall(predicted: set, gold: set) -> tuple[float, float]`
  (precision, recall, `gold` = Tier P's candidates); `count_inflation(baseline: set, refined: set) -> dict`
  (`{'baseline_count', 'refined_count', 'delta', 'delta_pct'}` — the Tier R
  over-splitting check from the spec's Analysis 2).

- [ ] **Step 1: Write the failing tests**

```python
# NovInvenio_Investigations/tests/test_genome_wide_concordance.py
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
import genome_wide_concordance as gwc  # noqa: E402


def test_load_candidates(tmp_path):
    p = tmp_path / 'candidates.txt'
    p.write_text('Ncra::pA1\nAfum::pB1\n')
    assert gwc.load_candidates(p) == {'Ncra::pA1', 'Afum::pB1'}


def test_load_candidates_empty_file(tmp_path):
    p = tmp_path / 'candidates.txt'
    p.write_text('')
    assert gwc.load_candidates(p) == set()


def test_jaccard():
    assert gwc.jaccard({'a', 'b'}, {'b', 'c'}) == 1 / 3
    assert gwc.jaccard(set(), set()) == 1.0  # both empty -> perfect agreement, not 0/0


def test_precision_recall():
    predicted = {'a', 'b', 'c'}
    gold = {'a', 'b', 'd'}
    precision, recall = gwc.precision_recall(predicted, gold)
    assert precision == 2 / 3
    assert recall == 2 / 3


def test_count_inflation():
    baseline = {'a', 'b'}
    refined = {'a', 'b', 'c', 'd', 'e'}
    result = gwc.count_inflation(baseline, refined)
    assert result == {'baseline_count': 2, 'refined_count': 5, 'delta': 3, 'delta_pct': 150.0}
```

- [ ] **Step 2: Run tests to verify they fail**

```
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
pixi run pytest tests/test_genome_wide_concordance.py -v
```
Expected: FAIL — module does not exist.

- [ ] **Step 3: Implement**

```python
#!/usr/bin/env python3
"""Analysis 2: genome-wide precision/recall/Jaccard of Tier C/C+H/R's novelty (or
loss) candidate lists against Tier P's, plus Tier R's candidate-count-inflation
check (the over-splitting signature the 16-control set can't see) -- see
notes/superpowers/specs/2026-09-13-cluster-vs-pairwise-sensitivity-design.md.
"""
import argparse
import csv


def load_candidates(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def jaccard(a, b):
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def precision_recall(predicted, gold):
    if not predicted:
        precision = 0.0 if gold else 1.0
    else:
        precision = len(predicted & gold) / len(predicted)
    recall = (len(predicted & gold) / len(gold)) if gold else 1.0
    return precision, recall


def count_inflation(baseline, refined):
    delta = len(refined) - len(baseline)
    delta_pct = (delta / len(baseline) * 100) if baseline else (0.0 if not refined else float('inf'))
    return {
        'baseline_count': len(baseline),
        'refined_count': len(refined),
        'delta': delta,
        'delta_pct': round(delta_pct, 1),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gold', required=True, help="Tier P's candidates.txt/loss_candidates.txt")
    ap.add_argument('--tier', action='append', nargs=2, metavar=('NAME', 'PATH'),
                    dest='tiers', required=True,
                    help='Repeatable: --tier C path/to/candidates.txt --tier "C+H" ... '
                         '--tier R path/to/refined_candidates.txt')
    ap.add_argument('--baseline-tier', default='C', dest='baseline_tier',
                    help='Which --tier name is the count_inflation baseline for Tier R '
                         '(default: C, i.e. is R inflated relative to unrefined clustering)')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    gold = load_candidates(args.gold)
    tier_candidates = {name: load_candidates(path) for name, path in args.tiers}

    rows = []
    baseline = tier_candidates.get(args.baseline_tier)
    for name, candidates in tier_candidates.items():
        precision, recall = precision_recall(candidates, gold)
        row = {
            'tier': name,
            'n_candidates': len(candidates),
            'precision_vs_P': round(precision, 4),
            'recall_vs_P': round(recall, 4),
            'jaccard_vs_P': round(jaccard(candidates, gold), 4),
        }
        if baseline is not None and name != args.baseline_tier:
            row.update({f'{k}_vs_{args.baseline_tier}': v
                       for k, v in count_inflation(baseline, candidates).items()})
        rows.append(row)

    fieldnames = sorted({k for row in rows for k in row})
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```
pixi run pytest tests/test_genome_wide_concordance.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add bin/genome_wide_concordance.py tests/test_genome_wide_concordance.py
git commit -m "$(cat <<'EOF'
Add genome_wide_concordance.py for Analysis 2

Precision/recall/Jaccard of each tier's novelty (or loss) candidate list
against Tier P's, plus Tier R's candidate-count-inflation check -- the
over-splitting signature the 16/2-control sets in Analysis 1 cannot detect.
EOF
)"
```

- [ ] **Step 6: Run for real (gains, `pezizo_set1`)**

Tier R's own candidate list first needs building from its refined families —
reuse `profile_to_matrix.py`'s keep-rule directly via a tiny inline filter
(Tier R's "candidates" are: family present in `>= ingroup_min_frac` of
ingroup AND `<= other_max_frac` of outgroup, computed from the refined
cluster membership the same way `build_cluster_membership_presence` does).
Since Task 5 already produced `pezizo_set1.R.controls_scored.tsv` (per
-control, not genome-wide), generate the genome-wide version now:

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
python3 - <<'EOF'
import sys
sys.path.insert(0, '/bigdata/stajichlab/jstajich/projects/NovInvenio/bin')
sys.path.insert(0, '/bigdata/stajichlab/jstajich/projects/NovInvenio/lib')
import pandas as pd
from config_parser import INGROUP_ROLES, OUTGROUP_ROLES, parse_config
import score_controls as sc

matrix = pd.read_csv('results/pezizo_set1_cluster/presence_matrix.tsv', sep='\t')
protein_to_proteome = dict(zip(matrix['protein_id'], matrix['source_proteome']))
samples = parse_config('studies/fungi/pezizo_set1/config.csv')
ingroup = {s.short for s in samples if s.group in INGROUP_ROLES}
outgroup = {s.short for s in samples if s.group in OUTGROUP_ROLES}
proteome_cols = [c for c in matrix.columns if c not in ('protein_id', 'source_proteome')]

rep_to_members = {}
with open('studies/fungi/pezizo_set1_cluster/tier_comparison/gain.refined_cluster.tsv') as fh:
    for line in fh:
        rep, member = line.rstrip('\n').split('\t')[:2]
        rep_to_members.setdefault(rep, []).append(member)

presence = sc.build_cluster_membership_presence(rep_to_members, protein_to_proteome, proteome_cols)
with open('studies/fungi/pezizo_set1_cluster/tier_comparison/gain.R.candidates.txt', 'w') as out:
    for rep, members in rep_to_members.items():
        pres = presence[rep]
        ing_frac = sum(pres.get(s, 0) for s in ingroup) / len(ingroup)
        out_frac = sum(pres.get(s, 0) for s in outgroup) / len(outgroup) if outgroup else 0.0
        if ing_frac >= 0.75 and out_frac <= 0.0:
            for m in members:
                out.write(f'{protein_to_proteome[m]}::{m}\n')
EOF

python3 bin/genome_wide_concordance.py \
  --gold results/pezizo_set1/candidates.txt \
  --tier "C+H" results/pezizo_set1_cluster/candidates.txt \
  --tier C studies/fungi/pezizo_set1_cluster/tier_comparison/gain.C.candidates.txt \
  --tier R studies/fungi/pezizo_set1_cluster/tier_comparison/gain.R.candidates.txt \
  --output studies/fungi/pezizo_set1_cluster/tier_comparison/pezizo_set1.gain.concordance.tsv
```

`--tier C`'s candidates file doesn't exist yet either — build it the same way
as the inline script above, but pointing `rep_to_members` at the
**unrefined** `results/pezizo_set1_cluster/families/families_cluster.tsv`
instead of the refined one (same loop body, different input file — write it
as a second small inline block or a `--cluster-tsv` CLI arg on the same
script; either is fine, but do not skip generating it — Task 7's
`--baseline-tier C` default needs it to compute Tier R's inflation numbers).

Repeat the same two commands for losses (`--gold
results/pezizo_set1/loss_candidates.txt`, Tier C/R candidates built from
`loss_families/families_cluster.tsv` + `loss_presence_matrix.tsv` with
`query_group=OUT` — ingroup/outgroup roles swapped relative to the gain
direction, per the spec's Analysis 2 loss-side construction) into
`pezizo_set1.loss.concordance.tsv`.

- [ ] **Step 7: Commit the real run's output**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add studies/fungi/pezizo_set1_cluster/tier_comparison/
git commit -m "$(cat <<'EOF'
Run Analysis 2 (genome-wide concordance, gains and losses) on pezizo_set1

Precision/recall/Jaccard of Tiers C+H/C/R against Tier P's full candidate
lists, plus Tier R's candidate-count-inflation check against Tier C, for
both directions.
EOF
)"
```

---

### Task 8: `trace_cost_report.py` — Analysis 3 (real compute cost)

**Files:**
- Create: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/bin/trace_cost_report.py`
- Test: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/tests/test_trace_cost_report.py`

**Interfaces:**
- Produces: `parse_duration(s: str) -> float` (nextflow's `1h 2m 3s`-style
  duration strings → seconds; also handles bare `7m 34s`, `346ms`, `0`);
  `load_trace(path: Path) -> list[dict]` (nextflow trace TSV → list of `{name,
  status, realtime_seconds}`); `cpu_hours_by_process_group(rows: list[dict], groups: dict[str, list[str]]) -> dict[str, float]`
  (`group_name -> total CPU-hours`, summing `realtime_seconds` across every
  row whose `name` contains any of that group's substrings, `COMPLETED`/
  `CACHED` only — `FAILED` retries don't count as real cost since they don't
  contribute to the final result).

- [ ] **Step 1: Write the failing tests**

```python
# NovInvenio_Investigations/tests/test_trace_cost_report.py
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
import trace_cost_report as tcr  # noqa: E402


def test_parse_duration_hours_minutes_seconds():
    assert tcr.parse_duration('2h 8m 50s') == 2 * 3600 + 8 * 60 + 50


def test_parse_duration_minutes_seconds():
    assert tcr.parse_duration('7m 34s') == 7 * 60 + 34


def test_parse_duration_milliseconds():
    assert tcr.parse_duration('346ms') == 0.346


def test_parse_duration_zero():
    assert tcr.parse_duration('0') == 0.0


def test_load_trace(tmp_path):
    p = tmp_path / 'trace.txt'
    p.write_text(
        'task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n'
        '1\ta\tb\tDIAMOND_SEARCH (x)\tCOMPLETED\t0\t2026-01-01\t1h\t45m\t100%\t1\t1\t1\t1\n'
        '2\ta\tb\tHMMSEARCH (y)\tFAILED\t1\t2026-01-01\t1h\t10m\t100%\t1\t1\t1\t1\n'
    )
    rows = tcr.load_trace(p)
    assert rows == [
        {'name': 'DIAMOND_SEARCH (x)', 'status': 'COMPLETED', 'realtime_seconds': 45 * 60.0},
        {'name': 'HMMSEARCH (y)', 'status': 'FAILED', 'realtime_seconds': 10 * 60.0},
    ]


def test_cpu_hours_by_process_group_excludes_failed():
    rows = [
        {'name': 'DIAMOND_SEARCH (x)', 'status': 'COMPLETED', 'realtime_seconds': 3600.0},
        {'name': 'DIAMOND_SEARCH (y)', 'status': 'FAILED', 'realtime_seconds': 3600.0},
        {'name': 'HMMSEARCH (z)', 'status': 'CACHED', 'realtime_seconds': 1800.0},
    ]
    groups = {'diamond': ['DIAMOND_SEARCH'], 'hmm': ['HMMSEARCH']}
    result = tcr.cpu_hours_by_process_group(rows, groups)
    assert result == {'diamond': 1.0, 'hmm': 0.5}
```

- [ ] **Step 2: Run tests to verify they fail**

```
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
pixi run pytest tests/test_trace_cost_report.py -v
```
Expected: FAIL — module does not exist.

- [ ] **Step 3: Implement**

```python
#!/usr/bin/env python3
"""Analysis 3: real per-tier compute cost from nextflow trace files already published
in each run's nextflow_log/ -- see notes/superpowers/specs/2026-09-13-cluster-vs-
pairwise-sensitivity-design.md.
"""
import argparse
import csv
import glob
import re


DURATION_RE = re.compile(r'(?:(\d+)h)?\s*(?:(\d+)m(?!s))?\s*(?:(\d+(?:\.\d+)?)s)?\s*(?:(\d+)ms)?')


def parse_duration(s):
    s = s.strip()
    if s in ('', '-', '0'):
        return 0.0
    h, m, sec, ms = DURATION_RE.match(s).groups()
    total = 0.0
    if h:
        total += int(h) * 3600
    if m:
        total += int(m) * 60
    if sec:
        total += float(sec)
    if ms:
        total += int(ms) / 1000
    return total


def load_trace(path):
    rows = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            rows.append({
                'name': row['name'],
                'status': row['status'],
                'realtime_seconds': parse_duration(row['realtime']),
            })
    return rows


def cpu_hours_by_process_group(rows, groups):
    totals = {name: 0.0 for name in groups}
    for row in rows:
        if row['status'] not in ('COMPLETED', 'CACHED'):
            continue
        for group_name, substrings in groups.items():
            if any(sub in row['name'] for sub in substrings):
                totals[group_name] += row['realtime_seconds'] / 3600
                break
    return {k: round(v, 3) for k, v in totals.items()}


def latest_trace(nextflow_log_dir):
    traces = sorted(glob.glob(f'{nextflow_log_dir}/*-trace.txt'))
    if not traces:
        raise FileNotFoundError(f'no *-trace.txt files under {nextflow_log_dir}')
    return traces[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pairwise-trace-dir', required=True, dest='pairwise_trace_dir',
                    help="pairwise run's nextflow_log/ directory")
    ap.add_argument('--cluster-trace-dir', required=True, dest='cluster_trace_dir',
                    help="cluster+HMM run's nextflow_log/ directory (Tier R's own new "
                         "within-family diamond step is timed separately, not read "
                         "from this trace -- see --tier-r-diamond-seconds)")
    ap.add_argument('--tier-r-diamond-seconds', type=float, default=None,
                    dest='tier_r_diamond_seconds',
                    help='Wall-clock of the refine_ambiguous_families.py diamond step '
                         '(time it separately with `time` when running Task 5/6 Step 6 '
                         '-- there is no trace file for a plain subprocess call)')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    pairwise_rows = load_trace(latest_trace(args.pairwise_trace_dir))
    cluster_rows = load_trace(latest_trace(args.cluster_trace_dir))

    pairwise_groups = {'search': ['DIAMOND_SEARCH', 'DIAMOND_SELF', 'TBLASTN']}
    cluster_groups = {'cluster_and_hmm': ['MMSEQS', 'BUILD_FAMILY_PROFILES', 'HMMSEARCH',
                                          'BUILD_CHUNK']}

    p_cost = cpu_hours_by_process_group(pairwise_rows, pairwise_groups)
    ch_cost = cpu_hours_by_process_group(cluster_rows, cluster_groups)

    rows = [
        {'tier': 'P', 'cpu_hours': p_cost['search']},
        {'tier': 'C+H', 'cpu_hours': ch_cost['cluster_and_hmm']},
        {'tier': 'C', 'cpu_hours': 0.0},  # free byproduct of Tier C+H's own clustering
        {'tier': 'R', 'cpu_hours': round((args.tier_r_diamond_seconds or 0.0) / 3600, 3)},
    ]
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['tier', 'cpu_hours'], delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```
pixi run pytest tests/test_trace_cost_report.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add bin/trace_cost_report.py tests/test_trace_cost_report.py
git commit -m "$(cat <<'EOF'
Add trace_cost_report.py for Analysis 3

Real per-tier CPU-hours from each run's own nextflow trace files (Tier P's
diamond/tblastn search cost vs Tier C+H's clustering+HMM cost; Tier C is
free, Tier R's own small diamond step is timed separately and passed in).
EOF
)"
```

- [ ] **Step 6: Run for real**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
# Re-time Tier R's diamond step from Task 5 Step 6b so Analysis 3 has a real number:
time python3 /bigdata/stajichlab/jstajich/projects/NovInvenio/bin/refine_ambiguous_families.py \
  --cluster-tsv results/pezizo_set1_cluster/families/families_cluster.tsv \
  --families results/pezizo_set1_cluster/families/families.tsv \
  --matrix results/pezizo_set1_cluster/presence_matrix.tsv \
  --oversized-families results/pezizo_set1_cluster/families/oversized_families.tsv \
  --config studies/fungi/pezizo_set1/config.csv \
  --pep Afum=studies/fungi/pezizo_set1/data_dir/pep/Afum.pep.fa \
        Amega=studies/fungi/pezizo_set1/data_dir/pep/Amega.pep.fa \
        Ztri=studies/fungi/pezizo_set1/data_dir/pep/Ztri.pep.fa \
        Ncra=studies/fungi/pezizo_set1/data_dir/pep/Ncra.pep.fa \
        Cimm=studies/fungi/pezizo_set1/data_dir/pep/Cimm.pep.fa \
  --self-hits results/pezizo_set1/self_hits/*.paralog_cutoffs.tsv \
  --output-cluster-tsv /tmp/gain.refined_cluster.tsv \
  --output-families /tmp/gain.refined_families.tsv
# Note the "real" wall-clock seconds from `time`'s output, then:

python3 bin/trace_cost_report.py \
  --pairwise-trace-dir results/pezizo_set1/nextflow_log \
  --cluster-trace-dir results/pezizo_set1_cluster/nextflow_log \
  --tier-r-diamond-seconds <seconds from `time`, e.g. 12.4> \
  --output studies/fungi/pezizo_set1_cluster/tier_comparison/pezizo_set1.cost.tsv
```

- [ ] **Step 7: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add studies/fungi/pezizo_set1_cluster/tier_comparison/pezizo_set1.cost.tsv
git commit -m "$(cat <<'EOF'
Run Analysis 3 (real compute-cost accounting) on pezizo_set1

CPU-hours per tier: Tier P's diamond+tblastn search, Tier C+H's clustering+
HMM, Tier C (free), Tier R's own small within-family diamond step (measured
directly, not read from a trace file).
EOF
)"
```

---

### Task 9: Assemble the final report (Analysis 4 + everything above)

**Files:**
- Create: `/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/bin/render_tier_comparison_report.py`
- Output: `studies/fungi/pezizo_set1_cluster/tier_comparison/report.md`

**Interfaces:**
- Consumes: every TSV produced by Tasks 5-8 (`*.tier_summary.tsv`, `*.per_control.tsv`,
  `*.gain.concordance.tsv`, `*.loss.concordance.tsv`, `*.cost.tsv`).
- Produces: one markdown file, no new pure functions to unit-test (this is a
  straight read-TSVs-render-markdown script — the individual pieces it reads
  were already tested and run for real in Tasks 5-8).

- [ ] **Step 1: Implement**

```python
#!/usr/bin/env python3
"""Assemble Tasks 5-8's TSV outputs into one markdown report with the Analysis 4
decision-framework recommendation -- see notes/superpowers/specs/2026-09-13-
cluster-vs-pairwise-sensitivity-design.md.
"""
import argparse
import csv
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def render_table(rows):
    if not rows:
        return '_(no rows)_\n'
    cols = list(rows[0].keys())
    lines = ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join(['---'] * len(cols)) + '|']
    for row in rows:
        lines.append('| ' + ' | '.join(str(row.get(c, '')) for c in cols) + ' |')
    return '\n'.join(lines) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tier-comparison-dir', required=True, dest='dir')
    ap.add_argument('--label', default='pezizo_set1')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    d = Path(args.dir)
    sections = []
    sections.append(f'# Cluster-vs-pairwise sensitivity: {args.label}\n')

    sections.append('## Analysis 1: controls recall/FP by tier\n')
    sections.append(render_table(read_tsv(d / f'{args.label}.tier_summary.tsv')))
    sections.append('\n### Per-control outcomes\n')
    sections.append(render_table(read_tsv(d / f'{args.label}.per_control.tsv')))

    gain_conc = d / f'{args.label}.gain.concordance.tsv'
    if gain_conc.exists():
        sections.append('\n## Analysis 2: genome-wide concordance (gains)\n')
        sections.append(render_table(read_tsv(gain_conc)))
    loss_conc = d / f'{args.label}.loss.concordance.tsv'
    if loss_conc.exists():
        sections.append('\n## Analysis 2: genome-wide concordance (losses)\n')
        sections.append(render_table(read_tsv(loss_conc)))

    cost = d / f'{args.label}.cost.tsv'
    if cost.exists():
        sections.append('\n## Analysis 3: real compute cost (CPU-hours)\n')
        sections.append(render_table(read_tsv(cost)))

    sections.append('\n## Analysis 4: decision framework\n')
    sections.append(
        '_Fill in by hand after reading the tables above — see the spec\'s Analysis 4 '
        'for the exact framework (escalate Tier R near-misses to Tier P; treat '
        'candidate-count inflation and HEX1-style clustering-inherent misses as the '
        'two conditions that would invalidate the framework)._\n'
    )

    Path(args.output).write_text('\n'.join(sections))


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Run for real and read the result**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
python3 bin/render_tier_comparison_report.py \
  --tier-comparison-dir studies/fungi/pezizo_set1_cluster/tier_comparison \
  --label pezizo_set1 \
  --output studies/fungi/pezizo_set1_cluster/tier_comparison/report.md
```

Read `report.md`. Manually write the Analysis 4 section's actual
recommendation (replacing the placeholder paragraph the script writes) based
on what the real Tasks 5-8 numbers actually show — this is a judgment call
the spec's Analysis 4 describes but cannot be templated, and it is the
investigation's actual deliverable. Do not leave the placeholder paragraph in
the committed file.

- [ ] **Step 3: Commit**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
git add bin/render_tier_comparison_report.py studies/fungi/pezizo_set1_cluster/tier_comparison/report.md
git commit -m "$(cat <<'EOF'
Assemble the cluster-vs-pairwise sensitivity investigation's final report

Renders Tasks 5-8's tier-summary, per-control, concordance, and cost tables
into one markdown report, with a hand-written Analysis 4 recommendation
based on the real numbers.
EOF
)"
```

---

## Self-Review

**Spec coverage:**
- Ground-truth verification (already done pre-plan, no task needed).
- Tier P/C+H/C/R definitions → Task 1 (`--presence-mode`), Task 2-3
  (`refine_ambiguous_families.py`).
- Analysis 1 (16-control + Agaricales) → Tasks 4-6.
- Analysis 2 (genome-wide concordance, gains + losses, inflation check) → Task 7.
- Analysis 3 (real compute cost) → Task 8.
- Analysis 4 (decision framework) → Task 9.
- "No `nf_NovInvenio` Nextflow changes" constraint → honored throughout (every
  new/modified file is a `bin/` script; Global Constraints states it explicitly).
- Promotion path (port Tier R into a DSL2 process) → explicitly out of scope
  per the spec, not a task in this plan.

**Placeholder scan:** the only non-code placeholder is Task 9's
`report.md`'s Analysis 4 paragraph, which Step 2 explicitly requires
replacing with a real judgment call before commit — not left as a TBD.

**Type/interface consistency:** `build_cluster_membership_presence` (Task 1)
is called identically by Tier C and Tier R (Task 5's `compare_cluster_tiers.py`
only varies which `--cluster-tsv`/`--families` it points at); `load_candidates`'s
id format (`source_proteome::protein_id`) matches the format
`build_presence_matrix.py`/`profile_to_matrix.py` already emit, verified
against real `candidates.txt` content, not assumed.
