# `bin/ni discover`: bulk species-population enumeration for pangenome studies

Date: 2026-09-11

## Problem

`bin/ni resolve` (shipped earlier today) answers "what's the accession for
*this* species/strain" — one row in, one accession out. It has no way to
answer "give me every annotated genome NCBI has for *Fusarium oxysporum*, and
let me pick which become ingroup vs. outgroup" — the actual first step of
building a pangenome study. Doing that by hand today means manually paging
through NCBI's website for a species that can have 70+ genome assemblies.

## Scope

This is a new, separate capability from `resolve` — it doesn't resolve a named
species/strain, it enumerates an entire population and helps partition it.
Two related capabilities from the same conversation are explicitly **not**
part of this spec:

- **Cross-species outgroup** — already possible today: run `bin/ni resolve`
  a second time for whatever separate outgroup species you want. No new code.
- **Population-structure-based partitioning (PCoA/DAPC/phylogeny)** —
  deferred to its own future spec. It requires actual genomic distance data
  (ANI, core-gene SNPs, or a phylogeny) computed from the genomes first, then
  a clustering step on top — a real analysis pipeline stage, not an
  accession-lookup tool. This repo's `nf_phyling` skill (BUSCO-marker
  phylogenomics) is the natural building block for the "phylogeny" option
  when that's designed; `bin/ni discover` never touches sequence data, only
  NCBI's own metadata.

## What NCBI's metadata gives you — and its real limits

**2026-09-11 revision — this section was rewritten after an independent
bioinformatics review (Fable model, with its own live NCBI queries) found the
original version's central claim was factually wrong for its own worked
example. Every point below reflects that review's live findings, not the
original draft.**

Live-verified (`datasets summary genome taxon "Fusarium oxysporum"
--annotated --as-json-lines`): 70 annotated assemblies. NCBI's own
`organism.organism_name` field frequently encodes host specialization ("forma
specialis") directly, e.g. `"Fusarium oxysporum f. sp. lycopersici 4287"` vs.
plain `"Fusarium oxysporum Fo47"`. Grouping by the `f. sp. <word>` token splits
the 70 into 30 with no forma specialis in the name, plus 16 named groups
(6 vasinfectum, 5 cubense, 4 conglutinans, 4 mori, 3 each of albedinis/cepae/
lycopersici, 2 each of rapae/raphani, 9 singletons).

**This is a real, usable grouping signal — but four things about the raw data
make it unsafe to use naively, all confirmed by tracing individual records,
not just counting them:**

1. **The "no forma specialis" bucket is not a coherent reference population.**
   It contains the *same strains* that also appear under named formae
   speciales, just registered without the pathotype in their name string:
   `Fo4287` (the reference *lycopersici* strain, also present as
   `f. sp. lycopersici 4287`), `Fo5176` (the reference *conglutinans* strain,
   also present as `f. sp. conglutinans Fo5176`), `II5` (the reference
   *cubense* TR4 strain), and a clinical human isolate (`NRRL 32931`).
   Whether a submitter included `f. sp. X` in the registered name is a
   registration choice, not a biological fact about the strain — the original
   draft of this spec called this bucket "typically generic/reference-type
   strains," which is wrong and has been removed.
2. **The same strain can appear as multiple separate assembly records.**
   `Fo5176` appears 3 times, `Fo47` 3 times (2 after GCA/GCF pair-collapse),
   `4287`/`Fo4287` under two different spellings. Treating each as a distinct
   pangenome member (which a naive per-record row + Short-suffix scheme
   would) inflates apparent core-genome size and can suppress real novelty
   calls.
3. **NCBI Taxonomy has already split parts of this species complex out under
   different species names.** `Fusarium oxysporum` (taxid 5507) sits under a
   parent `Fusarium oxysporum species complex` (taxid 171631, rank
   SPECIES_GROUP, 40 children). Querying the complex node directly returns 4
   more annotated genomes a plain `"Fusarium oxysporum"` name query never
   sees — including `Fusarium odoratissimum` `NRRL 54006` (= `II5`), the
   actual reference strain for cubense race TR4. This will recur for any
   species complex (Colletotrichum, Aspergillus section Nigri, the *F.
   solani*/*F. fujikuroi* complexes), not just this one.
4. **Annotation quality/method varies far more than forma specialis does, and
   nothing in the raw listing surfaces it.** The 70 records come from 26
   different annotation providers; annotated protein-coding gene counts range
   from 14,256 to 31,923 across genomes of the same nominal species — a
   spread large enough to be the dominant confounder in a novelty/loss
   pipeline, dwarfing whatever the forma-specialis grouping contributes.

**Separately, forma specialis itself is a host-range pathotype label, not a
clade** — several formae speciales are documented as polyphyletic in the
mycology literature (O'Donnell et al. 1998; multiple later studies
distinguishing cubense race 1 from TR4 as separate lineages), and horizontal
transfer of pathogenicity chromosomes between lineages is a documented
phenomenon specific to *F. oxysporum* (Ma et al. 2010, PMID 20237561).
`TaxonGroup` carrying the forma-specialis label is fine as metadata; treating
it as a phylogenetically coherent `IN`/`OUT` split is not something this tool
should do automatically (see `--auto`'s redefinition below).

## Command

```
bin/ni discover --species "<name>" --study-dir studies/<domain>/<set>
    [--ingroup-groups <comma-separated group labels>]
    [--outgroup-groups <comma-separated group labels>]
    [--auto]
    [--include-species-complex]
```

- Queries `datasets summary genome taxon "<name>" --annotated --as-json-lines`
  (restricting to annotated is not a flag — it's always on, since an
  unannotated assembly can't produce a protein FASTA, same rule `resolve`
  already enforces for its own NCBI protein path).
- **Groups by NCBI Taxonomy's own `FORMA_SPECIALIS` rank, not string-parsing
  `organism_name`.** Live-verified mechanism (re-review, all 74 records):
  each genome record's `organism.tax_id` (usually strain-rank, e.g. `1229664`
  "cubense race 1") needs a **separate, batched** `datasets summary taxonomy
  taxon <id>` lookup — the genome-summary record itself carries no lineage,
  and the taxonomy report's own `parents` field is a bare list of IDs with no
  rank annotations, so a **second batched taxonomy call on those parent IDs**
  is needed to learn which one is rank `FORMA_SPECIALIS` (confirmed live: 2
  batched calls total resolve all 74 records' groups). The FORMA_SPECIALIS
  node's own name is the full binomial-prefixed string (`"Fusarium oxysporum
  f. sp. cubense"`) — strip the species prefix to get the bare label
  (`cubense`) for `TaxonGroup`. This two-call mechanism is structurally
  correct and immune to tokenization edge cases (multi-word forma names,
  race/pathotype suffixes) the original string-parse draft was only
  accidentally safe against — confirmed to agree with the regex on 74/74
  live records. Fall back to the `f. sp. <token>` regex on `organism_name`
  only for the (rare) case a taxon has no `FORMA_SPECIALIS` ancestor at all,
  clearly logging that the fallback path was used for that row.
- **Checks whether the resolved taxon has a parent species-complex
  (`SPECIES_GROUP`-rank) node, and reports it regardless of flags.** If one
  exists, print how many additional annotated genomes live under the
  complex but not under the literal species name queried (live example: 4
  more for *F. oxysporum*, including the actual TR4 reference strain
  registered as *F. odoratissimum*). Without `--include-species-complex`,
  `discover` still only enumerates the literal species — but the report
  makes the gap visible rather than silent. With the flag, it queries the
  complex node instead and includes every child species' genomes, each row's
  `Species` reflecting whatever species name that record actually carries
  (do not force everything under the queried name).
- **Detects same-strain duplicate registrations** before assigning Shorts.
  **Correction from re-review (live-verified against the real dataset):**
  the criteria must be **BioSample match OR normalized-strain match**,
  applied as independent, unioned rules — not "BioSample first, strain only
  as a fallback when BioSample is absent." Every record in the live
  *F. oxysporum* + *F. odoratissimum* set (74/74) carries a BioSample, so a
  fallback-only rule never fires; the real duplicates this feature exists to
  catch (`Fo5176` re-registered 3 times, `Fo47` re-registered independently
  of its own GCA/GCF pair, `4287`/`Fo4287` under two spellings) each have
  **different** BioSample accessions — they're independent resubmissions of
  the same physical strain, not one BioSample with multiple assembly
  versions. So: group records whose `assembly_info.biosample.accession`
  matches, **union** with records whose normalized strain/isolate string
  matches case-insensitively after stripping a leading
  `Fo`/species-abbreviation prefix, treat each resulting group as **one**
  candidate row. When a duplicate group's members disagree on annotation
  pipeline/assembly level, prefer the one selected by the GCA/GCF
  annotation-preference rule below; report which records were merged and
  which was kept, never merge silently without a report line. **Known
  residual gap, report-only:** a duplicate that's both a different strain
  *string* and registered under a *different species name* in the
  species-complex case (live example: `II5` under *F. oxysporum* vs. `NRRL
  54006` under *F. odoratissimum` — the same physical strain) cannot be
  caught by either rule; `--auto`'s report should say so explicitly rather
  than implying duplicate-detection is exhaustive.
- **Always prints the group/count table** to stdout, regardless of which
  flags are given, now also carrying **provider, protein-coding gene count,
  assembly level, and contig N50** per group (min/max or a flagged outlier
  note) so annotation-quality variation — the dominant confounder — is
  visible before you assign groups, not just the forma-specialis counts.
- Refuses to run if `studies/<domain>/<set>/species.csv` already exists (this
  command seeds a *new* study; use `resolve`/hand-editing for an existing
  one) — `sys.exit` with a clear message, matching this repo's existing
  don't-silently-clobber conventions.

### Selection modes

1. **Neither `--ingroup-groups`/`--outgroup-groups` nor `--auto` given
   (default):** every candidate genome (after duplicate-collapsing) is
   written to `species.csv` with `Group` left **blank**. You assign `IN`/`OUT`
   (or delete rows you don't want) by hand afterward, informed by the
   printed table. This is the safe default — matches `resolve`'s own
   "never guess" posture, just applied to grouping instead of accession-
   picking.
2. **`--ingroup-groups`/`--outgroup-groups` given explicitly:** only genomes
   in a named group get written, with `Group` set accordingly. A group label
   named in either flag that doesn't exist in the actual data is an error
   (`sys.exit`, lists the real group labels found) — never silently ignored.
3. **`--auto` given — redefined after review, no longer writes `Group` at
   all.** The original design had `--auto` pick "largest named group = OUT,
   ungrouped bucket = IN" and write that directly into `species.csv`. Review
   found this produces a scientifically incoherent split on its own worked
   example (the "ungrouped" bucket is not a comparable baseline — see above)
   — an automatic `IN`/`OUT` inference from forma-specialis name presence is
   not a safe default under any framing. `--auto` now **prints an extended
   proposal section** to the group/count table (which group looks largest/
   best-sampled, which strains are duplicated across groups, which sibling
   species exist in the complex) but writes every row with `Group` blank,
   identically to mode 1 — `--auto` differs from the default only in how
   much analysis goes into the printed report, never in what gets written.

`--ingroup-groups`/`--outgroup-groups` and `--auto` are still mutually
exclusive (`sys.exit` if both given) since combining an explicit assignment
with a proposal report doesn't mean anything.

## Output rows

Since the `datasets` query already returns everything needed, `discover`
writes fully-resolved rows directly — it does not call `resolve_row`:

```
Short          — sanitized strain/isolate identifier, species-prefixed (e.g.
                 "Foxy_TR4", not a bare "TR4" or a bare numeric strain like
                 "9") so Shorts are never ambiguous or non-identifier-safe;
                 when strain looks like a race/pathotype token (e.g. "TR4",
                 "race 1") and a more specific isolate value exists, prefer
                 the isolate for the Short instead of the strain
Species        — the species-level name from NCBI Taxonomy's own
                 classification (walk the lineage to the SPECIES-rank node),
                 not string-surgery on organism_name -- organism_name varies
                 in shape even for records with no "f. sp." suffix at all
                 (e.g. "Fusarium oxysporum Fo47" vs "Fusarium oxysporum
                 NRRL 32931" both need the same Species value, which naive
                 suffix-stripping cannot guarantee)
Strain         — organism.infraspecific_names.strain (or .isolate if strain
                 is absent)
Group          — IN / OUT / blank, per the selection mode above
TaxonGroup     — the forma-specialis label from the FORMA_SPECIALIS rank
                 lookup above, or "no-fsp-in-name" when none exists (not
                 "(ungrouped)" -- that phrasing reads as "no information,"
                 when what it actually means is "this submitter didn't
                 register a pathotype name," which per the findings above is
                 not the same as "this is a baseline/reference strain")
Protein_Source — ncbi
Protein_Accession — the assembly accession (same GCA the genome uses --
                 this record's own annotation_info confirms a protein set
                 exists in this same NCBI Datasets package)
Taxon_ID       — organism.tax_id
Genome_Source  — ncbi
Genome_Accession — the assembly accession
GFF3_Source/GFF3_Accession — blank (auto -- Genome_Source=ncbi already
                 pulls whatever GFF3 that package includes, same as
                 resolve's existing behavior)
```

**GCA/GCF pair-collapse** (reusing `_collapse_paired` from `resolve`'s own
code) applies here too — a species-level listing can include both members of
a pair, and they must not become two rows for one genome.

**Correction from re-review (live-verified against all 4 real pairs in the
*F. oxysporum* dataset):** the first draft of this refinement described
`GCF_013085055.1` (Fo47's RefSeq member) as "NCBI's own independent RefSeq
annotation" versus the GCA member's "submitter's own annotation," and
proposed keying preference on `annotation_info`'s **provider** field. Neither
holds up: all 4 real pairs (1 Xi'an Jiaotong, 3 Broad Institute) have the
**same `provider` on both members** — RefSeq is *propagating* the
submitter's own annotation, via
`annotation_info.pipeline: "NCBI Eukaryotic Annotation Propagation
Pipeline"`, not producing an independent one. The actual protein-count
difference between pair members in all 4 real cases is tiny (0-5 genes,
e.g. Fo47's GCA=16,202 vs. GCF=16,197). Given that, **the pair-preference
rule stays the same as `resolve`'s existing GCA-preference rule — no new
rule needed.** What *is* worth capturing: the GCF/RefSeq-propagated member
uniquely carries a `busco` completeness block in these cases (Fo47's GCF
member: 97.8%) that the GCA member doesn't — record that BUSCO score in the
printed table (feeding the same per-group quality columns as the annotation-
heterogeneity reporting above) even though the GCA member is still the one
kept. Report which member of each pair was kept and its BUSCO score (if the
other member has one) in the printed table.

## Testing

Same pattern as `resolve`'s test suite: `query_ncbi_assemblies`-style parsing
is already tested; `discover`'s new pieces (taxonomy-rank grouping, the
group-table formatting with its new provider/gene-count/N50 columns, the
selection modes, `--auto`'s report-only behavior, duplicate-strain
collapsing, species-complex detection/inclusion, the annotation-aware
GCA/GCF preference rule, the existing-`species.csv` refusal) get unit tests
against canned `datasets`-shaped fixtures, no live calls in the test suite —
mirroring `tests/test_ni_resolve.py`'s existing `monkeypatch`-the-subprocess
pattern.

Given this session's own experience twice now (`resolve`'s 36 "passing"
tests that mocked past two real defects; this spec's own first draft making
a false biological claim its own worked example contradicts), fixtures must
be drawn from actual live captures, not hand-written shapes:

- The full 70-record `Fusarium oxysporum --annotated` capture, frozen as a
  JSON file under `tests/fixtures/` — covers the messy `organism_name`
  strings, the duplicate-strain cases (`Fo5176`, `Fo47`, `4287`/`Fo4287`),
  and the annotation-provider/gene-count spread.
- The 4 additional `Fusarium odoratissimum` records from the species-complex
  node — required so the species-complex-detection test is checked against
  the real record that motivated this feature (the actual TR4 reference
  strain, `II5`/`NRRL 54006`). This set deliberately includes `race 4`
  (`GCF_000350365.1`), which has no `isolate` field and no `f. sp.` in its
  name at all — a real edge case for both the Short-derivation rule (I4) and
  the FORMA_SPECIALIS-fallback path, not just a duplicate-detection case.
- At least one live-captured GCA/GCF pair (the `Fo47` pair) to confirm the
  pair-collapse keeps the GCA member as `resolve`'s existing rule already
  does — this spec's re-review found no separate annotation-preference rule
  is actually needed (see the Output Rows correction above), so this test
  just needs to confirm the existing behavior still applies at
  species-listing scale, not exercise new logic.

## Out of scope

- Population-structure-based partitioning (PCoA/DAPC/phylogeny) — future spec,
  a real analysis pipeline stage. `--auto`'s report may name candidate groups,
  but never infers or writes an `IN`/`OUT` split from anything short of the
  user's own explicit `--ingroup-groups`/`--outgroup-groups` choice.
  Given forma specialis's documented polyphyly, a real fix here is exactly
  what a future phylogeny/ANI-based spec should provide.
- Cross-species outgroup automation — already covered by re-running `resolve`.
- Automatically resolving BioSample `host`/isolation-source metadata into a
  grouping signal for species that don't encode forma specialis in their
  taxonomic name at all — a real future enhancement (the review noted `host`
  is populated for 40/70 *F. oxysporum* records and is often the better
  signal), but this spec's grouping stays scoped to the FORMA_SPECIALIS
  taxonomy rank (with the `organism_name` regex fallback) for now.
- Full species-complex generality: `--include-species-complex` is
  implemented specifically because it was needed for *F. oxysporum*'s own
  worked example, and the mechanism (query the `SPECIES_GROUP`-rank parent
  instead of the literal species) should generalize, but this spec doesn't
  attempt to verify it against every genus this tool might ever be pointed
  at — implementation should live-verify the mechanism against at least one
  other species complex before considering it broadly trustworthy.
