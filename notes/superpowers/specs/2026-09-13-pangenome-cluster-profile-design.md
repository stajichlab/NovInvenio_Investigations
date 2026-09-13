# Pangenome gene-family profiling: core/shell/cloud binning, co-loss/co-gain, and Starship-driven physical clustering

Date: 2026-09-13

## Problem

`fungi/Afumigatus_pangenome` (293 ingroup *A. fumigatus* strains, 2 outgroup) is the
first NII study built for an intraspecific pangenome question rather than the
ingroup-vs-outgroup novelty question NovInvenio's `--cluster_tool mmseqs`/`pairwise`/
`novelty_discovery` pathways answer. NovInvenio's existing family-clustering machinery
(`lib/clusters.py`'s `FamilyIndex`, `modules/mmseqs_cluster.nf`) and its `core.html`
report (`--core_min_frac`) both stop at "is this candidate/protein present broadly
enough" — neither bins accessory gene families by frequency, tests for co-occurring
gain/loss between families, nor asks whether co-varying families are physically
clustered in the genome (the Starship/giant-transposon mobility signature this
species is known for). This design adds that layer as a standalone analysis (not a
new NovInvenio Nextflow workflow — see "Why standalone" below), specific to this
study but written so the method generalizes to a future pangenome study.

Motivating biological question: given a pangenome's presence/absence matrix, (1)
lump genes found in essentially every strain as core and set them aside, (2) bin the
rest by how many strains carry them, excluding true singletons, (3) find pairs/groups
of accessory families whose gain/loss is correlated across strains, and (4) check
whether correlated families are also physically adjacent in the genome — the
signature of a mobile element (Starship) carrying a cargo cluster in or out of a
strain as a unit, as already characterized for several *A. fumigatus* secondary
metabolite clusters. A first concrete target: screen all 293 strains for the HAC
(hrmA-Associated Cluster) gene family and the independent `hacA` locus.

## Why standalone (not a new NovInvenio `--cluster_tool`)

Considered and rejected for now: a new `--cluster_tool pangenome_profile` Nextflow
workflow in NovInvenio itself, with its own modules/tests/report page. Rejected
because the method (binning cutoffs, co-occurrence statistic, synteny window size)
is still being worked out against real data — committing to Nextflow process
boundaries and a containerized tool matrix before the analysis itself is validated
would slow iteration for no present benefit. Once the method is settled and there's
a second pangenome study that would reuse it verbatim, revisit promoting the
validated pieces into NovInvenio proper (mirroring how `bin/build_study_config.py`
itself was promoted out of study-specific one-offs — see
`2026-09-11-study-onboarding-design.md`).

Consequence: this lives as scripts under
`studies/fungi/Afumigatus_pangenome/bin/` (study-specific, per NII's own
"study-specific vs. shared scripts" convention), reusing NovInvenio's `lib/`
where the contract already fits (`lib/clusters.py`'s cluster-TSV parsing) rather
than re-implementing it.

## Non-goals

- **Not** attempting to resolve expanded/contracted multigene family copy number
  from Illumina-based assemblies — acknowledged as unreliable at this read length/
  assembly quality, and out of scope even where PacBio/Nanopore assemblies exist in
  this set (kept consistent across the whole strain panel rather than special-cased
  per assembly technology).
- **Not** running OrthoFinder. mmseqs2/diamond's fast, near-identical-sequence
  clustering is a deliberate choice for a same-species pangenome (see "Clustering
  backend" below) — the multigene-family-merging behavior OrthoFinder is good at is
  not the goal here.
- **Not** a new HTML report in the `novelties.html`/`core.html` family. Outputs are
  TSVs/plots for exploratory analysis; a shareable report is a possible future step
  once the method stabilizes.

## Architecture

```
config.csv (293 IN / 2 OUT, existing)
        |
        v
 [1] cluster backend (mmseqs2 easy-cluster | diamond cluster/deepclust)
        |  -- ALL ingroup protein sequences, not just novelty candidates
        v
 rep/member TSV  (same contract as lib/clusters.py's read_cluster_tsv)
        |
        v
 [2] family x strain presence matrix  (+ optional per-strain copy count)
        |
        v
 [3] frequency binning: core (>=~95%, exact cutoff TBD from real histogram)
     | soft-core | shell | cloud | singleton (tracked, excluded from 4/5)
        |
        v
 [4] co-occurrence (Jaccard / phi) between shell+cloud families
     -> annotate each correlated pair/group with DAPC clade (TaxonGroup) composition
        |
        v
 [5] per-strain GFF3 gene-order synteny window
     -> flag N-gene windows whose family presence/absence correlates across strains
        |
        v
 [6] validation / benchmark suite against mbio.01092-25 supplement
     -> scorecard: does clustering/binning/synteny recover known Starships?
        |
        v
 [7] targeted screens: HAC (hrmA-associated family) + hacA (independent locus)
```

## Components

### 1. Clustering backend (swappable, mmseqs2 default)

One function per backend, both reduced to the same `(rep_id, member_id)` TSV shape
`lib/clusters.py::read_cluster_tsv()` already parses, so `build_families()` /
`FamilyIndex` work unchanged regardless of backend:

- `mmseqs2 easy-cluster` (default) — same tool/defaults NovInvenio's
  `modules/mmseqs_cluster.nf` uses (`--min-seq-id 0.3 -c 0.8 --cov-mode 0`), tuned
  for near-identical-sequence grouping rather than distant homology.
- `diamond cluster` / `diamond deepclust` — alternative backend, same output
  contract. Comparative performance is exactly what the benchmark suite (component
  6) scores, rather than assumed up front.

Input: all ingroup proteomes' sequences (293 strains), not the novelty-candidate
subset `workflows/cluster.nf`'s `CLUSTER` operates on — this is full pangenome
clustering.

**Open tuning question (bioinformatics review target):** mmseqs2/diamond identity
and coverage thresholds were chosen for NovInvenio's novelty-candidate use case;
whether `--min-seq-id 0.3 -c 0.8` (or diamond's equivalent) is the right operating
point for *cleanly separating* closely related but functionally distinct paralogs
within one species — as opposed to correctly merging alleles/near-identical
orthologs across strains of the same species — is unverified. This is precisely
what the benchmark suite in component 6 is for: known cargo-gene sets (Table
S12/S13) give a ground truth for "should these particular IDs cluster together,"
which a bare identity/coverage default choice does not guarantee.

### 2. Presence/frequency matrix

Family x strain matrix (0/1; optional per-strain member count for copy number,
non-authoritative per the Non-goals section) built from the cluster TSV plus
`config.csv`'s Short-to-proteome mapping — the same join `lib/config_parser.py`'s
`Sample` records already support.

### 3. Frequency binning

- **Core**: present in >= ~95% of strains (Roary/Panaroo-style convention) — lumped
  aside per the original ask; exact cutoff to be set from the real frequency
  histogram of this dataset rather than committed here, since 293 strains at
  varying assembly/annotation completeness may show a different natural break than
  a round number.
- **Soft-core / shell / cloud**: named bands below the core cutoff, boundaries
  likewise set from the observed histogram (candidate reference points raised
  during brainstorming: ~90-95% soft-core, ~15-90% shell, >1 strain but <~15% cloud)
  rather than fixed sight-unseen.
- **Singleton** (exactly 1 strain): tracked in the matrix, excluded from components
  4 and 5 (no co-occurrence or synteny signal is possible from n=1).

### 4. Co-occurrence (co-loss/co-gain)

Pairwise Jaccard (or phi coefficient) across shell+cloud family presence vectors.
Every correlated pair/group is annotated with its DAPC clade breakdown (`TaxonGroup`
column, already in `config.csv`) — **not** a formal phylogenetic correction, but
enough visibility to distinguish "these two families are gained/lost together
across independent clades" from "these two families both happen to be markers of
one clonal clade," per the explicit ask to flag rather than silently ignore that
confound.

### 5. Physical clustering (synteny)

Per-strain GFF3 gene order (`data_dir/gff3/`, already resolved per-strain by the
existing config) rather than a shared reference coordinate system — robust to
indels/rearrangement between strains, and does not require every strain to share
one strain's coordinate frame. Slide a window of N adjacent genes (start with a
couple of window sizes, e.g. 3 and 10, rather than committing to one — sensitivity
noted as an open control) and flag windows whose family presence/absence pattern is
correlated across strains: the physical-linkage signature of a block moving
together (candidate Starship cargo).

### 6. Validation / benchmark suite

Uses the newly ingested `mbio.01092-25-s0002.xlsx` (Gluck-Thaler, Forsythe, Puerner,
Gutierrez-Perez, Stajich, Croll, mBio 2025, doi:10.1128/mbio.01092-25 — provenance
recorded in this study's `DATA_MANIFEST.yaml`) as a battery of **known positive
controls**, not a single spot check:

| Table | Content | Used for |
|---|---|---|
| S6 | 20 high-confidence Starships + population frequency | Pick a frequency spread (common/rare/intermediate) so the benchmark covers easy and hard cases, not just the most obvious one |
| S21 | Per-strain presence/absence genotyping for those 20 | Ground-truth presence vector to compare component 3's matrix-derived calls against, for whichever of the 293 strains overlap the paper's population |
| S12 / S13 | BLAST-recovered cargo gene lists per Starship | Ground-truth "these gene IDs belong together" sets — tests whether component 1's clustering actually groups them into one family rather than splitting or over-merging |
| S7 | Manual annotation in AF293/A1163/CEA10 | Small, very clean 3-strain cross-check if those three (or their equivalent assemblies) are in the 293 |
| S14-S16 | Segregating insertion regions with genotyping | Boundary ground truth for component 5's synteny-window recovery |

Scorecard, per clustering backend (mmseqs vs. diamond): for each control Starship,
does clustering recover the correct cargo-gene grouping, does the frequency/
presence call match the paper's genotyped presence/absence, and does the synteny
window catch the correct region boundaries. This is the empirical answer to the
mmseqs-vs-diamond question, and validates components 1-5 as a whole before either
backend is trusted on unpublished/novel candidate clusters.

**Required control, not yet resolved:** the paper's gene/strain identifiers
(`AFUB_*`, `Afu*g*`, strain-specific `g#` locus tags, `XP_*` accessions; isolate/
genome-code strain names) do not string-match this pangenome's own annotation IDs
or `config.csv` `Short` values. Recovery must go through a sequence-based crosswalk
(diamond/mmseqs against the paper's cited protein accessions), not literal ID
matching, and needs an explicit strain-overlap check between the 293-strain set
here and the paper's ~519 (population) / 13 (reference-quality) / 3 (AF293/A1163/
CEA10) strain sets before assuming any given control strain is actually present in
this study.

### 7. Targeted screens: HAC + hacA

Per the resolved definition: `hacA` (Afu3g04070, chromosome 3) and `hrmA`
(Afu5g14900, chromosome 5) are **not** one contiguous physical cluster — confirmed
from Table S19 (both listed as independently characterized virulence loci, no
shared cluster tag) and basic locus-tag chromosome placement. "HAC" (hrmA-
Associated Cluster) targets the paralogous gene family Table S5 annotates as
*"Subtelomeric hrmA-associated cluster protein AFUB_079030/YDR124W-like"* (Pfam
PF11001 / IPR047092), documented there riding inside several different named
Starships across strains (`Nebuchadnezzar-h1`, `Osiris-h3`, `Logos-h1`, `Gnosis-h2`,
`Logos-h2`, `navis10-var35`) — a real, paper-documented case of the same gene
family moving between distinct mobile elements. Screen: for all 293 strains,
which HAC-family paralog(s) (if any) are present, and — where the strain overlaps
the paper's annotated set — which Starship (if any) it currently rides in. `hacA`
is screened separately as its own independent presence/absence locus.

## Open controls (recorded, not resolved by this design)

- **ID crosswalk** between this pangenome's own protein/gene IDs and the paper's
  mixed ID schemes (component 6 and 7 both depend on this).
- **Strain overlap** between the 293-strain set and the paper's population(s).
- **Assembly/annotation completeness per strain** — a "missing" shell/cloud gene
  could be a genuine loss or an assembly gap/annotation miss. Recommend a
  completeness column (e.g. BUSCO) per strain before over-interpreting rare/cloud
  calls, particularly for strains not among the reference-quality set.
- **Synteny window-size sensitivity** — no single N is committed; try at least two
  sizes and compare.
- **Core/shell/cloud band cutoffs** — set from the real frequency histogram once
  computed, not fixed in this design.

## Testing / validation plan

No NovInvenio-style pytest suite (this is exploratory, standalone scripts per "Why
standalone" above). The benchmark suite in component 6 **is** the validation: it
must run and produce a scorecard before component 5's synteny method or either
clustering backend is applied to any novel (non-paper-validated) candidate region,
per the explicit ask that this test how well mmseqs/diamond and the synteny method
actually perform before trusting them on new findings.

## Study-specific notes

Strain-level open questions, the HAC/hacA gene ID list, and controls specific to
running this against `Afumigatus_pangenome`'s actual 293 strains are tracked in
`studies/fungi/Afumigatus_pangenome/PANGENOME_CLUSTER_PROFILE_NOTES.md`, kept
separate from this general design so the method described here stays reusable by
a future pangenome study.
