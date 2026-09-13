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
 [1] cluster backend, two-tier (mmseqs2 easy-cluster | diamond cluster)
        |  -- ALL ingroup protein sequences, not just novelty candidates
        |  -- tier 1: allele/ortholog groups (~90-95% id) = the unit for 2-5
        |  -- tier 2: re-cluster tier-1 reps loosely (~30-50% id), superfamily label only
        v
 rep/member TSV  (same contract as lib/clusters.py's read_cluster_tsv)
        |
        v
 [1b] genome-level rescue pass (tblastn/miniprot) for coverage-failed calls
      + isoform collapse + strain dereplication (Mash/ANI)
        |
        v
 [2] family x strain presence matrix  (0/1/genome-only three-state, + copy count)
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

### 1. Clustering backend (swappable, mmseqs2 default) — two-tier, per Fable review finding 1

**Revised after bioinformatics review (2026-09-13, see "Review disposition" below).**
NovInvenio's existing default (`--min-seq-id 0.3 -c 0.8 --cov-mode 0`) was tuned for
cross-species novelty-candidate detection and is far too loose here: at 30% identity
with greedy/transitive clustering, every PF11001 paralog in a strain (the HAC-family
notes list seven distinct loci in CEA10 alone) would collapse into one family,
making component 7's "which paralog, riding which Starship" question unanswerable,
and deflating shell/cloud counts by saturating paralog-family presence. *A. fumigatus*
orthologs across strains are typically >98% identical (the same reasoning bacterial
pangenome tools like Roary/Panaroo build on), so intraspecific family resolution
needs a much tighter operating point than cross-species novelty detection does.

Two-tier scheme:
- **Tier 1** (the unit for components 2-5 — allele/ortholog groups): mmseqs2
  `--min-seq-id 0.9 -c 0.8 --cov-mode 0 --cluster-reassign`; diamond equivalent
  `diamond cluster --approx-id 90 --member-cover 80`. `deepclust` targets remote
  homology and is not appropriate for this step.
- **Tier 2** (superfamily label only, not used for presence/frequency): re-cluster
  tier-1 representatives at a looser 30-50% identity, purely for annotation (e.g.
  grouping the PF11001 paralogs under one HAC-family label without merging their
  per-strain presence calls).

Both tiers reduced to the same `(rep_id, member_id)` TSV shape
`lib/clusters.py::read_cluster_tsv()` already parses, so `build_families()` /
`FamilyIndex` work unchanged. The benchmark suite (component 6) tunes the tier-1
identity empirically (try 0.85/0.9/0.95 against Table S12/S13's known cargo-gene
groupings) rather than treating 0.9 as final.

Input: all ingroup proteomes' sequences (293 strains), not the novelty-candidate
subset `workflows/cluster.nf`'s `CLUSTER` operates on — this is full pangenome
clustering.

### 1b. Genome-level rescue + strain dereplication — new, per Fable review findings 1-2

**Added after bioinformatics review.** Two failure modes the original design didn't
address:

- **Coverage-based false absence.** `--cov-mode 0` requiring 80% coverage of both
  sequences means a truncated model at a contig edge, a split gene, or a
  fusion drops out of its family — becoming a spurious singleton or a false
  "loss" that is actually an annotation artifact, not biology. Since strains here
  are annotated via more than one upstream pipeline/source (NCBI- and
  UniProt-derived), annotation-source batch effects are a real risk of
  contaminating co-occurrence patterns (component 4) with pipeline signal, not
  strain biology. Rescue: for every family, search its representative against
  each strain's own genome sequence (tblastn — already in NovInvenio's toolchain
  — or miniprot) and record a three-state presence call: protein-model-present /
  genome-only-hit / absent. A completeness column (e.g. BUSCO) does not substitute
  for this — it flags gross incompleteness, not this specific failure mode.
- **Isoform/transcript inflation and strain duplication.** Collapse isoforms
  before clustering (`bin/collapse_isoforms.py` already exists in NovInvenio) or
  per-strain "copy number" will include alternative transcripts, not paralogs.
  Separately, dereplicate strains by Mash/ANI before any frequency count — public
  strain sets can contain the same isolate under two different names/accessions,
  which would inflate every downstream presence-fraction call.

### 2. Presence/frequency matrix

Family (tier-1) x strain matrix — three-state (present / genome-only / absent) per
component 1b, plus per-strain copy count (isoform-collapsed) — built from the
cluster TSV plus `config.csv`'s Short-to-proteome mapping — the same join
`lib/config_parser.py`'s `Sample` records already support. Frequency counts use
dereplicated strains (component 1b), not the raw 293.

### 3. Frequency binning

- **Core**: present in >= ~95% of strains (Roary/Panaroo-style convention) — lumped
  aside per the original ask; exact cutoff to be set from the real frequency
  histogram of this dataset rather than committed here, since 293 strains at
  varying assembly/annotation completeness may show a different natural break than
  a round number. **Per Fable review finding 7:** compute the histogram (and the
  final cutoff) after excluding low-completeness strains, or report the cutoff both
  ways — an unfiltered histogram will be biased downward by incomplete assemblies
  showing spurious absences.
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

**Revised after bioinformatics review (findings 3-4):**
- **Statistical rigor, not just a raw Jaccard threshold.** With ~10^4 accessory
  families there are ~10^7-10^8 pairs, and two cloud families each present in the
  same handful of strains score Jaccard 1.0 trivially. Require a minimum frequency
  floor (e.g. present in >=5 dereplicated strains) before testing a pair at all,
  use Fisher's exact / hypergeometric p-values with Benjamini-Hochberg FDR rather
  than a bare similarity cutoff, and report effect size alongside significance.
  Jaccard is also symmetric — it cannot itself distinguish co-gain from co-loss;
  polarize direction using the two outgroups (`Aslen_ref`, `Neofi_ref`): present in
  the outgroup(s) and most ingroup clades implies loss in the strains lacking it,
  absent from the outgroup(s) implies gain. Don't label a pair "co-loss" without
  this polarization step.
- **A cheap stratified null, not full phylogenetic comparative methods (yet).**
  Clade annotation shows a confound without measuring it. Two low-cost additions:
  (a) collapse near-clonal strains (already required for dereplication in 1b)
  before computing statistics; (b) a permutation null that shuffles each family's
  presence vector *within* `TaxonGroup` (not globally), so a pair significant
  under the stratified null co-varies beyond clade membership. Requiring
  co-variation in >=2 independent clades is a crude but useful independent-
  contrasts proxy. Coinfinder (Whelan et al. 2020) implements lineage-aware
  co-occurrence testing and is worth running as an external cross-check rather
  than reimplementing from scratch. Full phylogenetic independent contrasts are
  not warranted at this exploratory stage, but every reported pair should be
  labeled a *hypothesis*, not a confirmed association, until it clears the
  stratified null.

### 5. Physical clustering (synteny)

Per-strain GFF3 gene order (`data_dir/gff3/`, already resolved per-strain by the
existing config) rather than a shared reference coordinate system — robust to
indels/rearrangement between strains, and does not require every strain to share
one strain's coordinate frame. Slide a window of N adjacent genes (start with a
couple of window sizes, e.g. 3 and 10, rather than committing to one — sensitivity
noted as an open control) and flag windows whose family presence/absence pattern is
correlated across strains: the physical-linkage signature of a block moving
together (candidate Starship cargo).

**Refined after bioinformatics review (finding 5) — the frame was right, the
procedure was underspecified:**
- Windows must be enumerated **per strain** and reduced to family-tuples, not
  anchored to one reference strain's gene order — otherwise the test only sees
  windows that happen to exist in whichever strain was picked as reference.
- Exclude windows that cross a contig boundary — a draft assembly will otherwise
  produce truncated/spurious "windows" at contig edges.
- A fixed N=3/N=10 window only tests local adjacency; Starships carry tens to
  hundreds of genes. Add a complementary **accessory-island** test: define maximal
  runs of consecutive non-core genes per strain, then test island-level presence/
  co-occurrence, rather than relying solely on fixed-size windows.
- For each correlated family pair/group, report the fraction of co-carrying
  strains where both members actually lie within k genes of each other on the same
  contig — a direct linkage statistic, not just a correlation coefficient.
- Use the Starship captain gene (DUF3435 tyrosine recombinase) family as a
  positive marker where present: a correlated block adjacent to a captain is much
  stronger mobility evidence than adjacency alone.

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

**Added after bioinformatics review (finding 6): negative controls, not positive
only.** The battery above is entirely positive controls (known Starships that
should be recovered). Add: (a) a handful of conserved, non-mobile secondary
metabolite gene clusters (present in nearly all strains) to estimate the synteny
step's false-positive rate on clusters that should *not* look variably mobile, and
(b) a set of random accessory-family pairs (matched for frequency) to estimate
component 4's co-occurrence false-positive rate under the actual data. Also note:
Table S12/S13's cargo-gene lists are themselves BLAST-derived in the source paper,
i.e. a method-dependent reference rather than absolute ground truth — treat
disagreement with them as "investigate," not automatically "our method is wrong."

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
- **Fraction of draft (Illumina/short-read) vs. long-read (Nanopore/PacBio)
  assemblies among the 293 strains** (Fable review finding 8) — not yet
  determined; the severity of the coverage-based-false-absence and contig-edge
  synteny issues (component 1b, component 5) scales directly with this fraction,
  so it should be computed and recorded before running the full analysis.
- **Strain dereplication (Mash/ANI)** — required before any frequency count
  (component 1b); not yet run.

## Review disposition (Fable bioinformatics review, 2026-09-13)

An independent bioinformatics review (Claude Fable model, full report retained in
session history) returned three must-fix findings and three should-consider
findings; all are folded into the components above rather than listed separately,
so implementation follows this document directly:

- **Must-fix, addressed above:** (1) clustering identity/coverage regime — see
  component 1's two-tier revision; (2) fragmented/heterogeneous gene models
  causing spurious absences — see component 1b's rescue pass; (3) co-occurrence
  needing a frequency floor + FDR correction rather than a raw Jaccard threshold —
  see component 4's revision.
- **Should-consider, addressed above:** stratified permutation null for
  co-occurrence (component 4), underspecified synteny procedure (component 5),
  positive-only benchmark suite (component 6).
- **Minor, addressed above:** core cutoff sensitivity to strain completeness
  (component 3); draft-vs-long-read assembly fraction (Open controls).

Net effect: component 1's clustering identity moved from mmseqs2's inherited
30%/80%-coverage novelty-detection default to a two-tier ~90%-identity
allele/ortholog scheme with a separate loose superfamily tier, and component 2
gained a genome-level rescue pass — both directly load-bearing for whether the
HAC-family screen (component 7) can distinguish individual paralogs at all.

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
