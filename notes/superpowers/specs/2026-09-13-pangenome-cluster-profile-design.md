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
- ~~**Not** a new HTML report in the `novelties.html`/`core.html` family. Outputs
  are TSVs/plots for exploratory analysis; a shareable report is a possible future
  step once the method stabilizes.~~ **Superseded 2026-09-13** (Opus review M9):
  component 9 below now specs exactly this report. Kept struck through rather than
  deleted so the reversal is visible, not silently contradicted.

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

### 8. Pair classification: physical (Starship-explained or not) vs. trans — addendum, 2026-09-13

Motivated by real per-strain results (HAC screen, hacA/hrmA reference screen) that
sharpened a distinction the original design left implicit. Components 4 and 5 were
built deliberately separate — co-occurrence (4) tests **every** shell/cloud family
pair regardless of genomic position; synteny (5) tests **only** physically adjacent
families — but nothing joins their outputs, so a pair significant in component 4
currently has no answer to "is this because they're physically linked, or because
they covary despite being unlinked?" That join is the missing piece, and it is the
direct answer to two distinct real questions this pathway needs to serve:

1. **Physical gene clusters that gain/lost together** — most are expected to be
   Starship-explainable (this population's dominant known mobility mechanism), but
   the design must actively look for the exceptions: a physically-clustered
   co-gain/co-loss block with **no** Starship/captain-gene evidence nearby is
   itself a finding (a novel mobile element, an unannotated Starship variant, or a
   non-Starship mechanism entirely — recombination hotspot, segmental duplication/
   loss, etc.) and should be surfaced, not discarded as "expected."
2. **Single, physically unlinked genes that covary** — the classic case is a
   toxin-antitoxin-style functional pair, but the textbook toxin-antitoxin
   architecture is usually *itself* clustered (often inside a secondary-metabolite
   cluster), which is exactly why it's the wrong template for the genuinely novel
   class this component targets: two genes that co-gain/co-lose across the
   population **without** being genomic neighbors at all. That pattern is only
   explicable by some other correlated selective/functional pressure (a real
   candidate interaction/co-evolution signal, or a shared regulatory/ecological
   driver), which makes it the more interesting — and more error-prone — output of
   the two, so it needs the stricter labeling below.

**Revised after independent review (Claude Opus model, 2026-09-13) — the frame was
right, the classifier as first drafted would have degenerated to "everything is
unexplained_physical" and let clonal/completeness confounds straight into `trans`.
Must-fix findings and how each is now addressed:**

- **A real evidence source for "Starship-explained" (must-fix M1).** As first
  drafted, nothing in this pipeline actually detects a Starship or a captain gene —
  component 5 only *recommends* using one as a marker, and the paper's own
  Table S5/S21 coordinates only cover 254/293 strains by name-match, don't
  cross-walk to this study's own assembly versions (see `Asfu_A1163` in the study
  notes), and aren't a general per-strain Starship-finder. Without a real source,
  every physical pair falls to `unexplained_physical` by default, inverting its
  intended meaning ("we looked and found nothing" vs. "we never looked"). **Fix:**
  add an explicit prerequisite step — `hmmsearch` a DUF3435 (tyrosine recombinase,
  the Starship captain gene) profile against all 293 proteomes, producing
  `captain_loci.tsv` (strain, protein_id, contig, position) — as a hard dependency
  of this component, run alongside component 1b rather than assumed. A pair's
  physical block counts as `starship_explained` only if a captain hit falls within
  the same accessory island/window on the same contig in at least one carrying
  strain, or (for the 254 name-matched strains only) a paper Starship coordinate
  overlaps. **Add a fourth label, `physical_unclassified`**: syntenic, no captain
  hit, AND no paper coordinate available for any carrying strain — distinct from
  `unexplained_physical` (syntenic, captain search *was* run, genuinely came back
  empty). Absence of evidence must stay distinguishable from evidence of absence.
- **Non-overlapping, non-exhaustive bands, fixed (must-fix M2).** `trans` was
  defined as `linkage_fraction == 0` and `unexplained_physical`/`starship_explained`
  as "syntenic," an unquantified threshold — pairs in between (e.g. linked in 2 of
  40 co-carrying strains, exactly what a partial/relic insertion looks like) fell
  in neither bucket. **Fix:** band `linkage_fraction` with named cut points —
  `>= 0.5` physical (eligible for `starship_explained`/`unexplained_physical`/
  `physical_unclassified`), `<= 0.05` `trans`-eligible, in between
  `ambiguous_linkage` (reported, not discarded) — and require a minimum of 5
  co-carrying strains before any classification is attempted (fewer than that,
  label `insufficient_data`). `k` (the within-k-genes rank distance
  `linkage_fraction` uses) needs the same sensitivity-analysis treatment component
  5 already gives window size `N` — they are different parameters and were
  previously conflated.
- **`trans`'s exclusion criterion was smaller than the mechanism it excludes
  (must-fix M3) — the single most damaging fix.** Component 5 itself states
  "Starships carry tens to hundreds of genes," which is why it added the
  accessory-island test; the first draft then excluded physical linkage for
  `trans` using the *same* within-k-genes (`k` on the order of 10) rank-distance
  statistic. Two genes 40 genes apart inside one Starship, or in the same
  subtelomeric block, would score `linkage_fraction == 0` and get classified
  `trans` — the worst possible false positive for the bucket the user cares about
  most. **Fix:** `trans` requires ALL of: (a) not in the same accessory island in
  any carrying strain, (b) not within a physical bp window sized to a real Starship
  footprint (not a gene-rank count) on the same contig in any carrying strain, and
  (c) on different contigs, OR separated by more than that bp window, in every
  strain where both are on an assembled (not fragmented) contig.
- **Multi-copy families broke the linkage statistic (must-fix M4).** This study's
  own headline result — PF11001 present at 3-9 copies in every one of 293 strains,
  Table S5 documenting 6 distinct loci on 6 different Starships in CEA10 alone —
  means `linkage_fraction`'s single-position-per-family-per-strain assumption is
  wrong for exactly the families most likely to be mobile. **Fix (code change to
  `synteny_windows.py`'s `linkage_fraction`, not just this component):** for a
  multi-copy family, take the MINIMUM rank/bp distance over all copy-pairs (i.e.
  "is *any* copy of family A near *any* copy of family B"), not whichever single
  position the code happens to have stored.
- **`trans` must actually be gated by the permutation null, not just FDR (must-fix
  M5).** `find_cooccurring_pairs` only filters on `fdr_q`; `permutation_p` was
  computed but never used to gate anything, even though this component's own text
  admitted `trans` is the bucket with no physical fallback if the null misses a
  confound. **Fix:** a pair may only be labeled `trans` (not merely reported) if it
  ALSO clears `permutation_p < 0.05` AND co-occurs across >= 2 distinct clades
  under the SAME clade-labeling scheme (operationalizing component 4's existing
  "independent-contrasts proxy" prose, which had no code or component-8 reference
  before this fix) — see the stratification-scheme fix immediately below for why
  "same scheme" matters.
- **The stratification column is a mixture of three incompatible clustering
  schemes (must-fix M5, continued).** `TaxonGroup` mixes `Barber_clusterN` (252
  strains), `DAPC_Clade_N` (10), and `Mash_clade_N` (31) — three different
  clustering methods in one column, 14 strata total. Shuffling within a
  `Mash_clade_N` stratum that actually spans several real Barber clusters under-
  controls; treating a `Mash_clade_0` strain and a real-same-clade
  `Barber_cluster6` strain as different strata over-controls. Both directions of
  error, simultaneously. **Fix:** add a `taxon_scheme` column (`barber` / `dapc` /
  `mash`) alongside `TaxonGroup`; restrict the clade-stratified permutation null
  (and the `>= 2 clades` gate above) to the 252 single-scheme (`barber`) strains
  unless/until the schemes are reconciled, and state this restriction on every
  `trans` row rather than silently mixing schemes.
- **Strain quality/annotation-source confound must be checked, not just cited
  (must-fix M6).** Component 1b already flags that strains come from more than one
  upstream annotation source and calls this "a real risk of contaminating
  co-occurrence patterns with pipeline signal" — component 8 didn't carry that
  forward. A low-completeness or differently-annotated assembly misses many genes
  at once, so everything it misses "co-occurs" with everything else it misses;
  that signal is orthogonal to clade (the permutation null won't catch it) and has
  no physical linkage (it lands straight in `trans`). The 17-strain smoke test
  already shows real exposure here: 9,043 of 20,257 families (45%) are singletons.
  **Fix:** before reporting any `trans` pair, check it is not enriched among
  low-completeness (BUSCO, once computed per the Open controls) or single-source
  strains; report a completeness/source composition column alongside
  `clade_composition`, analogous to it.
- **Outgroup-based polarization breaks under 90%-identity tier-1 clustering
  (must-fix M7).** Component 4 polarizes gain/loss using `Aslen_ref`/`Neofi_ref`
  presence read from the same tier-1 (90% identity) matrix (per the study notes'
  Pipeline conventions item 3) — but *A. lentulus*/*A. fischeri* orthologs of
  genuinely conserved genes will often fall below a 90%-identity, 80%-coverage
  threshold against *A. fumigatus*, scoring falsely absent from both outgroups and
  defaulting every such gene to `"gain"`. With only 2 outgroups and a hard cutoff,
  this is not a rare edge case. **Fix:** call outgroup presence via a more
  permissive route before polarizing — tier-2 superfamily membership (already
  computed, component 1) or the genome-level tblastn rescue pass (component 1b) —
  never tier-1 membership alone. This is a direct, previously-unnoticed
  consequence of the original Fable review's must-fix #1 (the two-tier identity
  scheme) interacting with component 4, not a new problem introduced here.
- **Output table, corrected column set (was missing `direction_b`/`jaccard`,
  couldn't represent a multi-Starship family, had no block grouping — must-fix
  m1-m3 in the review):** `pair_classification.tsv` columns: `family_a, family_b,
  classification` (`starship_explained` / `unexplained_physical` /
  `physical_unclassified` / `ambiguous_linkage` / `trans` / `insufficient_data`),
  `linkage_fraction, block_id` (connected components over the syntenic pair graph
  -- the physical unit is a block, not a pair), `jaccard, fisher_p, fdr_q,
  permutation_p, direction_a, direction_b` (a pair where A gains and B loses is not
  a co-gain/co-loss pair — both directions must be reported, not just one),
  `taxon_scheme, clade_composition, completeness_composition,
  tier2_sibling_flag` (true if `family_a`/`family_b` are tier-2 superfamily
  siblings — a possible clustering-artifact flag, component 1's tier-2 output was
  otherwise unused here). Starship identity moves to its own long-format table,
  `pair_starship_evidence.tsv` (family, strain, starship_id, captain_protein_id),
  since a single family can ride several different named Starships across strains
  (a scalar `starship_id_or_null` column cannot represent that, per the study's own
  HAC finding).
- **`genome_only` calls have no gene-order position (must-fix m4 in the review).**
  A tblastn-rescued (`GENOME_ONLY`) call has no GFF3 gene model and therefore no
  position for `linkage_fraction`/accessory-island membership to use. Rescued
  calls must be excluded from the synteny side of this join entirely (not silently
  treated as "present but position-less," which would systematically depress
  `linkage_fraction` and bias pairs toward `trans`) — state this exclusion
  explicitly wherever `genome_only` calls are folded into "present" for other
  purposes (e.g. `PresenceMatrix.is_present`).
- **Do not finalize classification thresholds against the smoke test.** The
  17-strain run's own `cooccurring_pairs.tsv` is header-only (zero pairs at
  `min_strain_count=5` out of 17 strains) — component 8 has never seen a single
  real input row. Component 9's "don't build the report against synthetic/smoke
  data" rule (below) applies equally here: tune bands/thresholds only against the
  real 293-strain run.

### 9. Report design — addendum, 2026-09-13

Reuse NovInvenio's existing report architecture (`lib/report_common.py`'s shared
CSS/JS chrome, `lib/skins.py`'s theming, the self-contained-no-network-fetch
constraint, the canvas-heatmap-plus-accessible-table-tab pattern) rather than
building a new page design from scratch — this pathway was always meant to
graduate into NovInvenio (see the "Why standalone" section), and building its
report as a structural twin of `novelties.html`/`core.html`/`losses.html` from the
start is what makes that graduation cheap later instead of a rewrite.

**Revised after independent review (Opus, 2026-09-13) — three concrete collisions
with this codebase's own documented report conventions, plus an underspecified
network panel:**

- **A new `pangenome.html`**: the family×strain presence/copy-number heatmap
  (canvas, virtualized rows), with:
  - **Payload encoding specified now, not deferred (should-consider S1).** 10-20k
    families × 293-700 strains is 6-14M cells — a 30-70x jump over the existing
    novelty report's ~20k rows x a handful of proteome columns, not the "modest
    increase" first claimed. A naive per-cell JSON encoding is tens of MB against
    this codebase's single-self-contained-file constraint. Encode the three-state
    matrix as a packed bitset (2 bits/state) per family, base64 in the JSON
    payload, not one JSON value per cell; the copy-number sidecar stays a sparse
    `{(family, strain): n}` map (most cells are 0/1, not >1) rather than a dense
    array.
  - **Column virtualization, not just row virtualization (should-consider S2).**
    The existing pages' `drawGrid()` only virtualizes rows because they've never
    had more than a handful of proteome columns; 293-700 strain columns needs the
    same treatment on the column axis, or a per-clade-collapsed default view with
    drill-down to per-strain columns on demand — pick one explicitly rather than
    assuming the existing row-only virtualization already covers this.
  - **Strain grouping by clade, WITHOUT new hue tokens (should-consider S3).**
    `lib/skins.py`'s own contract: hue carries evidence *type*, not group —
    ingroup/outgroup is deliberately carried by column position, not color, and
    there is no categorical N-color palette machinery to reuse (only 2 data-mark
    tokens exist, `--series-1`/`--series-2`, each WCAG-checked across every shipped
    skin). Coloring 14 `TaxonGroup` strata would mean adding 14 new tokens that
    `tests/test_skins.py` enforces contrast/CVD-separation on for every skin — real
    new work, not reuse, and a second data dimension competing for hue with the
    existing evidence-type encoding. **Do this by column position/ordering plus a
    labeled band strip instead** (grouped columns, a thin text-labeled row above
    the heatmap marking clade boundaries) — consistent with the existing
    position-not-hue convention rather than fighting it.
  - Row ordering/blocking by frequency bin (component 3) first, then by
    `pair_classification`'s `block_id` (not raw pair membership — the physical
    unit is a block, per component 8's fix) within shell/cloud.
  - A per-family detail panel surfacing: frequency bin, copy-number distribution,
    any `pair_classification` rows it participates in (with `classification` shown,
    not just "is in a pair"), and per-strain Starship identity from
    `pair_starship_evidence.tsv` where applicable (a family can ride several named
    Starships across strains, not one).
- **The `trans`-pair network panel, right-sized (should-consider S5 — was
  underspecified in the three ways that decide whether it works):** a fixed
  (non-simulated) deterministic layout, not an in-page force-directed sim with an
  unbounded perf profile; a hard cap on displayed edges (top-N by effect size, N
  fixed, e.g. 100); and a **mandatory edge-list table alongside the graph** — per
  this codebase's own rule that "the table tab is the accessibility twin of the
  heatmap and must keep every value reachable without hovering," a graph with no
  tabular equivalent breaks that rule. If the real `trans` set turns out large
  (component 8's confound fixes may shrink it substantially, or may not — this was
  an unbounded assumption, not a measured one), fall back to table-only above a
  fixed node-count threshold rather than rendering an unreadable hairball.
- **Test enforcement, not just conventions (should-consider S4).** The existing
  three report pages are covered by `tests/test_report_templates.py` (raw-hex
  detection, `node --check` on script bodies) and `tests/test_report_js_behaviour.py`
  (jsdom-driven DOM behavior); `tests/test_skins.py` enforces the WCAG floors on
  every skin. `pangenome.html` must be added to those same test globs from the
  start, not left to inherit the constraints without the guardrails that enforce
  them — a report "living outside the suite" is exactly how the hardcoded-color
  test's own rationale gets silently violated.
- **Do not build this against synthetic or 17-strain smoke-test data.** Per the
  2026-09-13 conversation that raised this addendum: the report and any
  "discoveries" it surfaces only mean something once a real presence matrix exists
  (293-strain run at minimum) — sequence the real run first, the report second.

### 10. Next steps and scaling to ~500-700 genomes — addendum, 2026-09-13

**Immediate next steps on the current 293-strain set** (in dependency order):

**Revised after independent review (Opus, 2026-09-13) — must-fix M8 caught that
the original list asserted "confirmed fast at 293-strain scale" for steps that had
NOT actually been run at that scale yet (only Mash+PCoA and a 2-query/1-profile HAC
screen had); must-fix M10 caught that the list dropped three gates this document
elsewhere calls mandatory. Both fixed below. One step (2) is now backed by a real
measurement taken after the review, which is included rather than discarded:**

1. Isoform-collapse all 293 proteomes (blocked on a per-strain gene<->protein
   mapping — `NovInvenio/bin/collapse_isoforms.py` needs an NCBI feature_table.txt
   per strain; these custom assemblies may not have one, so this may need a GFF3-
   based collapse instead — an open question, not yet resolved). **Interim
   decision, recorded not silently assumed:** a same-day spot check of 3 strains'
   GFF3s found gene:mRNA ratios of 10188:10503, 10109:9929, and 8825:8825 — a
   real, low (~3% at most) isoform rate consistent with fungal biology, not the
   near-1:1-isoform-per-gene pattern that would make skipping this step
   dangerous. The full 293/295-strain run below proceeded WITHOUT isoform
   collapse on this basis — a checked, not a blind, simplification — but this
   remains an open item to resolve properly (a GFF3-based collapse) rather than a
   closed decision.
2. Build the Short-prefixed, tier-1-clustered, full presence matrix
   (`cluster_backend.py mmseqs-tier1` -> `build_presence_matrix.py`). **Now backed
   by a real measurement, not an extrapolation:** run against all 295 strains
   (293 IN + 2 OUT, 2,788,402 proteins) on 2026-09-13, tier-1 clustering finished
   in ~14 minutes wall-clock, producing 47,983 families. This directly answers
   must-fix M8's challenge for THIS specific step — the claim was previously
   unsupported, it no longer is, for clustering specifically. (Presence-matrix
   building, frequency binning and co-occurrence at this scale remain unmeasured
   until steps 4-5 below actually run — do not extend this one real data point to
   the rest of the pipeline, which is exactly the error M8 caught.)
3. Genome-level tblastn rescue pass (component 1b) — now known to be a *real*,
   not theoretical, need: the HAC screen caught one fragmented gene model (
   `Asfu_E165L4`'s hrmA split across two protein records) by accident on a single
   locus; a genome-wide rescue pass will find more.
4. **Gate (must-fix M10, previously missing): strain dereplication (Mash/ANI,
   component 1b) and assembly-quality determination (BUSCO completeness,
   draft-vs-long-read fraction — Open controls) BEFORE any frequency count.**
   This document states dereplication is required "before any frequency count"
   elsewhere; the step list didn't actually place it before step 5 below. Making
   it explicit here, not just in the notes' unchecked-item list.
5. Frequency binning (component 3) on the real, dereplicated data — this is when
   the core/soft-core/shell/cloud cutoffs stop being provisional.
6. **Gate (must-fix M10, previously missing): the component 6 benchmark suite
   MUST run and produce a scorecard before this step.** The "Testing / validation
   plan" section already says synteny/either clustering backend must not be
   applied to any novel candidate region before this; the step list previously
   ran straight from frequency binning into co-occurrence/synteny with no
   benchmark step anywhere.
7. Co-occurrence (component 4) — **with the O(n^2) pair-enumeration cliff (must-
   fix M8) mitigated before running at this scale, not after seeing it hang.**
   `find_cooccurring_pairs` currently enumerates every shell/cloud family pair via
   `itertools.combinations` before BH correction. At 17 strains that's ~3,485
   shell+cloud families (a few million pairs); at 293 strains, if the smoke
   test's 45% singleton rate is any guide, a conservative 20-40k eligible
   families gives 2x10^8-8x10^8 pairs — each a materialized `(vec_a, vec_b)` pair
   of 293-element lists, hundreds of GB resident, plus ~10^8 `scipy.fisher_exact`
   calls, plus 1,000 permutations per FDR-survivor. This is genuinely quadratic in
   family count and family count itself grows with strain count in an open
   pangenome — **the real scaling cliff in this whole pathway**, previously
   unflagged. Required before running: streaming two-pass BH (don't retain all
   vector pairs simultaneously), bitset/popcount-encoded presence vectors instead
   of Python bool lists, an analytic prefilter (e.g. a fast chi-square screen)
   with exact Fisher only on survivors, and permutation testing only on
   FDR-survivors with early stopping.
8. Synteny (component 5) + the new pair classification (component 8) — in that
   order and only after step 7's mitigations are in place, since 8 depends on
   both 5 and 7's (properly-run) output.
9. The `pangenome.html` report (component 9) once step 8 produces real,
   non-synthetic output.
10. Separately, and not blocking 1-9: the full 293-strain ParSNP run
    (`run_parsnp_ingroup293.sh`), deferred per the 2026-09-13 conversation decision
    (Mash+PCoA's real concordance against 262 published labels — 5/7 Barber
    clusters 99-100% pure — was judged sufficient for now). If revisited: see the
    corrected ParSNP scaling discussion below before assuming the 17-strain
    prototype's regime still applies.

**Scaling to ~500-700 genomes later** — what changes and what doesn't (corrected
per Opus review must-fix M8; the previous version's blanket "no scaling cliff"
claim was wrong for the one component that actually has one):

- **Mash + PCoA + k-means (component "TaxonGroup fill")**: scales trivially: Mash
  sketching/distance is near-linear and the 293-strain run already finished in
  well under a minute. No changes needed. (This is the one component whose real
  293-strain timing was already measured before this addendum was written.)
- **mmseqs2/diamond tier-1 clustering**: now backed by a real 295-strain
  measurement (step 2 above, ~14 min) — expect roughly the same order of runtime
  at 500-700 genomes, not a step change, though this is one data point, not a
  scaling curve; re-measure rather than assume linearity.
- **Presence-matrix building and the id_crosswalk/HAC screen's diamond+hmmsearch
  approach**: the HAC screen's own diamond/hmmsearch runs at 293-strain scale
  (minutes) are real evidence, but that was a 2-query/1-profile search — a
  fundamentally smaller workload than an all-vs-all family search. Don't
  extrapolate from it to the full pipeline's heavier steps.
- **Co-occurrence (component 4) is NOT confirmed fast and has a real, quadratic
  scaling cliff (must-fix M8) — this is the correction to the previous version's
  central claim.** `find_cooccurring_pairs`'s all-pairs enumeration over
  shell+cloud families is quadratic in family count, and family count itself
  grows with strain count in an open pangenome (more strains -> more accessory
  families discovered, not a fixed set) — so 293->700 is a plausibly superlinear
  step, not "the same order of runtime." This must be measured at the 293-strain
  scale (step 7 above, with its mitigations) before any claim is made about
  500-700; do not repeat the unmeasured-optimism error here.
- **ParSNP will NOT scale the same way as the 17-strain prototype implies, and
  the 293-strain regime is ALREADY different (should-consider S6) — not only a
  500-700 concern.** The 17-strain prototype's core alignment ran strictly
  single-threaded (parsnp only parallelizes across partitions, which don't form
  below ~50 genomes) and took 1.46h wall-clock (matching CPU time for a
  single-threaded step; the study notes' earlier "45+ CPU-minutes" was a
  still-running estimate taken before the job finished and has been corrected to
  this same 1.46h figure — should-consider S7's cross-document inconsistency is
  resolved). **293 genomes already exceeds the ~50-genome partition threshold**,
  so a 293-strain ParSNP run activates a qualitatively different
  (partition-parallel) code path than the prototype measured — the 48h/8c/32gb
  ask in `run_parsnp_ingroup293.sh` is therefore a reasoned guess extrapolated
  across a regime change, not a timing "informed by" a directly-comparable
  measurement, and should be labeled as such. 500-700 genomes stays on the same
  partitioned code path as 293, so the 293-strain run (if/when performed) is the
  right benchmark for 500-700 — not the 17-strain prototype.
- **Synteny/accessory-island detection (component 5)**: the per-strain,
  contig-boundary-excluded design was chosen specifically to be assembly-count-
  agnostic (no shared reference coordinate frame), so it should scale linearly in
  the number of strains with no algorithmic change needed. Component 8's new
  captain-gene hmmsearch step (must-fix M1) is one profile against N proteomes —
  same near-linear shape as the HAC screen's PF11001 search.
- **The report (component 9)**: no longer assumed a "modest increase" — see
  component 9's corrected payload-encoding discussion (bitset packing specified,
  not deferred) for why 293->700 columns needed a real answer, not an assumption.
- **More strains is NOT automatically more power (should-consider S10, corrects
  the previous version's one-sided framing).** More strains buys more independent
  clades to stratify the permutation null against and a better-resolved frequency
  histogram — real benefits. But it also grows the accessory family count, and
  therefore the pair count component 4/8 must test, roughly quadratically (see
  the corrected co-occurrence scaling discussion above) — so the multiple-testing
  burden can grow faster than the per-pair power gain. Whether net discoverability
  of real `trans` pairs actually improves at 500-700 vs. 293 is an open empirical
  question this design should state as such, not assume as a payoff.

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

## Review disposition (Opus review, 2026-09-13)

An independent bioinformatics review (Claude Opus model, full report retained in
session history) covered components 8-10 (the same-day addenda above) and returned
ten must-fix findings, ten should-consider findings, and seven minor findings — all
folded into components 8/9/10 above (each fix is inline where the finding applies,
tagged with its finding ID), not listed separately, so implementation follows this
document directly, same convention as the Fable disposition above.

- **Must-fix, addressed above:** (M1) no real evidence source for
  `starship_explained`, degenerating the classifier — component 8 gained a captain-
  gene hmmsearch prerequisite and a fourth `physical_unclassified` label; (M2)
  non-exhaustive/overlapping bands — explicit `linkage_fraction` cut points and an
  `ambiguous_linkage`/`insufficient_data` label added; (M3) `trans`'s exclusion
  criterion (`k`-genes) smaller than the mechanism it excludes (Starships span
  tens-hundreds of genes) — `trans` now requires accessory-island and bp-distance
  exclusion, not gene-rank distance alone; (M4) `linkage_fraction` assumed one
  copy per family per strain, broken by this study's own 3-9-copies-per-strain HAC
  finding — fixed to take the minimum over all copy-pairs; (M5) `trans` was gated
  only by FDR, not the permutation null the component's own text said it needed,
  and `TaxonGroup` mixes three incompatible clustering schemes — both fixed
  (permutation + >=2-clade gate, `taxon_scheme` column, restrict to the 252
  single-scheme strains); (M6) strain completeness/annotation-source confound
  (cited in component 1b, not carried into component 8) — completeness/source
  composition check added; (M7) outgroup-based gain/loss polarization breaks under
  90%-identity tier-1 clustering (a previously-unnoticed interaction with the
  Fable review's own must-fix #1) — outgroup presence must now use tier-2/rescue,
  not tier-1 alone; (M8) the scaling claim in component 10 was unsupported for
  most of the pipeline and missed the one real quadratic cliff
  (`find_cooccurring_pairs`'s all-pairs enumeration) — claim corrected, mitigations
  specified, and partially re-supported by a real 295-strain tier-1 clustering
  measurement (~14 min) taken during this same revision; (M9) component 9
  contradicted the Non-goals section without amending it — Non-goals struck
  through with a pointer, not left silently contradictory; (M10) the next-steps
  ordering dropped three gates this document calls mandatory elsewhere (benchmark
  suite before synteny/novel-region application, dereplication before frequency
  counting, assembly-quality checks before the full run) — all three added back as
  explicit numbered steps.
- **Should-consider, addressed above:** report payload/cell-count understated
  (S1, now a specified bitset encoding); no column-axis virtualization plan (S2);
  clade coloring collides with the skins system's hue-carries-evidence-type
  contract (S3, now column position/ordering + a labeled band strip instead);
  new report page would inherit no test enforcement (S4, now added to the same
  test globs as the existing three pages); `trans` network panel had no layout/
  cap/accessibility-twin plan (S5, all three specified); ParSNP's 293-strain
  regime already differs from the 17-strain prototype's (S6, and cross-document
  timing inconsistency S7 resolved to the single measured 1.46h figure — see the
  study notes); BH-FDR's independence assumption is questionable under clonal
  population structure (S9, noted as a real limitation — not fixed by a specific
  code change, flagged for the real run to reconsider, e.g. BY correction or a
  permutation-derived threshold); more-strains-is-more-power was one-sided (S10,
  corrected to note the countervailing pair-count growth).
- **Should-consider, deliberately not changed (recorded, not silently dropped):**
  mutually-exclusive/anti-associated pairs are out of scope by this component's
  one-sided Fisher test (S8) — noted here as a known, accepted scope limit rather
  than fixed, since adding a two-sided pass is a real scope expansion best decided
  against real data, not preemptively; component 8 has never seen a real input row
  (S11, the 17-strain smoke test's `cooccurring_pairs.tsv` is header-only) — not a
  fix, a standing reminder folded into component 8's own "do not finalize
  thresholds against the smoke test" line.
- **Minor, addressed above:** `pair_classification.tsv`'s column set was missing
  `direction_b`/`jaccard` and couldn't represent a multi-Starship family or a
  physical block (m1-m3, all fixed in the corrected column list, plus a new
  `pair_starship_evidence.tsv` and `block_id`); `genome_only` calls have no gene-
  order position and would otherwise bias pairs toward `trans` (m4, now an
  explicit exclusion); tier-2 superfamily siblings were never flagged as a
  possible clustering-artifact signal (m5, added as `tier2_sibling_flag`).
- **Minor, noted not separately fixed:** an uncited Mash timing claim (m6, already
  correctable to a real number — see component 10's step 2 rewrite, which cites
  the real measurement taken during this revision) and the accessory-island
  >10%-missing warning not being an explicit checkpoint in the step list (m7,
  covered by the same "assembly-quality determination" gate added for M10).

Net effect: component 8 changed from a two-input join with an unfilled evidence
gap into a procedure with an explicit new upstream data source (captain-gene
search) and machine-checkable gates before a pair is ever labeled `trans` — the
bucket the user most wants to trust. Component 10's central claim reversed from
"nothing here scales badly" to "one specific, previously-unflagged step
(co-occurrence's pair enumeration) is the real bottleneck, and it must be
mitigated before the pathway is run at 293 scale, let alone 500-700."

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
