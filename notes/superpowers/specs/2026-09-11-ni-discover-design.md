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

## What NCBI's metadata already gives you for free

Live-verified 2026-09-11 (`datasets summary genome taxon "Fusarium
oxysporum" --annotated --as-json-lines`): 70 annotated assemblies. NCBI's own
`organism.organism_name` field frequently already encodes host specialization
("forma specialis") directly, e.g. `"Fusarium oxysporum f. sp. lycopersici
4287"` vs. plain `"Fusarium oxysporum Fo47"`. Grouping by the `f. sp. <word>`
token (when present) splits the 70 into:

```
 30  (no forma specialis in organism_name)
  6  f. sp. vasinfectum
  5  f. sp. cubense
  4  f. sp. conglutinans
  4  f. sp. mori
  3  f. sp. albedinis / f. sp. cepae / f. sp. lycopersici
  2  f. sp. rapae / f. sp. raphani
  1  each of 9 more named formae speciales
```

This is real, structured, already-curated grouping information — no sequence
analysis needed to get a usable first partition for many plant-pathogenic
species. Genomes with no forma specialis in their name (the 30-strong bucket
above) are typically generic/reference-type strains (e.g. `Fo47`).

## Command

```
bin/ni discover --species "<name>" --study-dir studies/<domain>/<set>
    [--ingroup-groups <comma-separated group labels>]
    [--outgroup-groups <comma-separated group labels>]
    [--auto]
```

- Queries `datasets summary genome taxon "<name>" --annotated --as-json-lines`
  (restricting to annotated is not a flag — it's always on, since an
  unannotated assembly can't produce a protein FASTA, same rule `resolve`
  already enforces for its own NCBI protein path).
- Parses each record's `organism.organism_name` for an `f. sp. <token>`
  substring; records without one go in a `(ungrouped)` bucket.
- **Always prints the group/count table** (like the one above) to stdout,
  regardless of which flags are given — you should see the population shape
  before committing to a partition.
- Refuses to run if `studies/<domain>/<set>/species.csv` already exists (this
  command seeds a *new* study; use `resolve`/hand-editing for an existing
  one) — `sys.exit` with a clear message, matching this repo's existing
  don't-silently-clobber conventions.

### Selection modes (all three from the original conversation, as one flag's behavior)

1. **Neither `--ingroup-groups`/`--outgroup-groups` nor `--auto` given
   (default):** every annotated genome is written to `species.csv` with
   `Group` left **blank**. You assign `IN`/`OUT` (or delete rows you don't
   want) by hand afterward, informed by the printed table. This is the safe
   default — matches `resolve`'s own "never guess" posture, just applied to
   grouping instead of accession-picking.
2. **`--ingroup-groups`/`--outgroup-groups` given explicitly:** only genomes
   in a named group get written, with `Group` set accordingly. A group label
   named in either flag that doesn't exist in the actual data is an error
   (`sys.exit`, lists the real group labels found) — never silently ignored.
3. **`--auto` given (opt-in heuristic, not the default):** proposes a split —
   the single largest named forma-specialis group becomes `OUT`, the
   `(ungrouped)` bucket becomes `IN` (reference-type strains vs. the
   best-sampled specialized lineage — a simple, defensible starting point,
   not a scientific claim). If there's a tie for the largest named group,
   `--auto` refuses (`sys.exit`, lists the tied groups) rather than picking
   arbitrarily. **The written `species.csv` is a proposal to review, not a
   trusted result** — say so in the command's own stderr output, since
   nothing downstream re-checks this choice.

`--ingroup-groups`/`--outgroup-groups` and `--auto` are mutually exclusive
(`sys.exit` if both given).

## Output rows

Since the `datasets` query already returns everything needed, `discover`
writes fully-resolved rows directly — it does not call `resolve_row`:

```
Short          — derived from strain/isolate name (sanitized to
                 alphanumeric+underscore, truncated), de-duplicated with a
                 numeric suffix on collision
Species        — organism.organism_name, with any "f. sp. ..." suffix
                 stripped back to the plain binomial (the forma specialis
                 becomes TaxonGroup, not part of Species)
Strain         — organism.infraspecific_names.strain (or .isolate if strain
                 is absent)
Group          — IN / OUT / blank, per the selection mode above
TaxonGroup     — the forma-specialis label ("f. sp. lycopersici") or
                 "(ungrouped)"
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

GCA/GCF pair-collapse (reusing `_collapse_paired` from `resolve`'s own code)
applies here too — a species-level listing can include both members of a
pair, and they must not become two rows for one genome.

## Testing

Same pattern as `resolve`'s test suite: `query_ncbi_assemblies`-style parsing
is already tested; `discover`'s new pieces (forma-specialis grouping, the
group-table formatting, the three selection modes, the `--auto` tie-refusal,
the existing-`species.csv` refusal) get unit tests against canned
`datasets`-shaped fixtures, no live calls in the test suite — mirroring
`tests/test_ni_resolve.py`'s existing `monkeypatch`-the-subprocess pattern.
Given this session's own experience (36 "passing" tests that turned out to
mock past two real defects), at least one fixture must be drawn from an
actual live `datasets summary genome taxon "Fusarium oxysporum" --annotated`
capture (frozen as a JSON file under `tests/fixtures/`, not hand-written),
so the grouping-parser test is checked against real organism_name strings
(including the messy ones: `"...f. sp. conglutinans race 2 54008"`,
`"...f. sp. cubense race 1"`) rather than a clean invented shape.

## Out of scope

- Population-structure-based partitioning (PCoA/DAPC/phylogeny) — future spec,
  a real analysis pipeline stage.
- Cross-species outgroup automation — already covered by re-running `resolve`.
- Any grouping signal beyond `organism_name`'s `f. sp.` token (e.g. parsing
  BioSample host/isolation-source metadata for species that don't encode
  forma specialis in their taxonomic name at all) — a real future
  enhancement, but `organism_name` already covers `F. oxysporum`'s well-known
  case and this spec doesn't try to solve the general problem.
