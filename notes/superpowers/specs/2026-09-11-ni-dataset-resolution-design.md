# `bin/ni`: dataset resolution and config assembly

Date: 2026-09-11

## Problem

Building a new study today means the user (or an AI session standing in for them)
already knows, for every species, whether it has a UniProt reference proteome and
what its NCBI assembly accession is — or has to find out by hand, one species at a
time, via UniProt's website and NCBI Datasets. This session hit exactly that
friction twice: `Scom` (*Schizophyllum commune*) turned out to have no UniProt
reference proteome at all (confirmed by checking `rest.uniprot.org/uniparc` live),
and `UHM_lachnoNovelclade`'s four named ingroup species needed their NCBI
accessions supplied by hand because nothing in this repo could look them up.

`bin/ni` automates that lookup: given a `species.csv` with `Short`/`Species`/
`Strain`/`Group`/`TaxonGroup` already filled in and the `Protein_Source`/
`Genome_Source`/`*_Accession` columns blank, it resolves each row to a real,
justified accession — or leaves it blank and says exactly why, never guessing.

## Scope: Phase 1 only

This spec covers **accession resolution for already-chosen species** — the user
supplies species names, `bin/ni resolve` finds their accessions. Open-ended
discovery ("find me good Ascomycota outgroups") is a deliberately separate,
later phase — nothing here blocks building it, but it needs its own design
(candidate ranking against phylogenetic/taxonomic criteria is a different problem
than "resolve this exact species," with different failure modes and no
comparably crisp non-ambiguous stopping condition). Noted in `DESIGN.md`'s
"Explicitly deferred" section once this ships.

## Command surface

```
bin/ni resolve --study-dir studies/<domain>/<set>
bin/ni fetch    --study-dir studies/<domain>/<set>   # thin pass-through
```

`resolve` is the new capability this spec describes. `fetch` does no new work —
it execs `bin/build_study_config.py --study-dir <same dir>` unchanged, so the
already-tested download/materialize logic (Task 1 of the 2026-09-11 onboarding-
dispatch plan) isn't duplicated. `bin/ni` itself is a thin argparse dispatcher;
all resolution logic lives in `lib/ni_resolve.py`, unit-testable against mocked
API responses the same way `tests/test_build_study_config.py` mocks `run()`.

## Resolution algorithm, per row

Only rows with **both** `Protein_Source` and `Genome_Source` blank are touched —
a row with either already filled (e.g. a `local_faa` row a human already set) is
left alone, so `resolve` is safe to re-run after partial manual edits.

**2026-09-11 revision — this section was corrected after an independent
bioinformatics review (Fable model) verified the original algorithm against
live UniProt REST and `datasets` queries and found it would silently blank or
mis-resolve real species, including this project's own flagship ingroup
species. See the four "Critical" fixes below; each is load-bearing, not
stylistic.**

### Step 0 — resolve `Species` to an NCBI taxon ID first

`GET https://rest.uniprot.org/taxonomy/search?query=<Species>&fields=taxon_id` (or
equivalent) resolves the species name to a taxon ID before either Step 1 or
Step 2 runs. This is required, not optional — see Step 1's first fix below.
Multiple taxon-ID matches (e.g. a name that's ambiguous at the species level)
→ leave the row blank, report the candidates, same "never guess" posture as
every other ambiguity in this design. A name that resolves to zero taxa →
leave blank, report "not a resolvable taxon name" (this covers placeholder
names like `Xylaria sp.` that don't correspond to a single NCBI taxon).

### Step 1 — UniProt reference-proteome search

Query `GET https://rest.uniprot.org/proteomes/search?query=taxonomy_id:<taxon
ID from Step 0> AND proteome_type:REFERENCE` — **not** a free-text
`organism_name:"<Species>"` query. Live testing during review found
`organism_name` matches mycoviruses named after their host and marked
`Reference proteome` in their own right (querying "Neurospora crassa" returns
both the real fungal proteome UP000001805 *and* "Neurospora crassa fusarivirus
1", UP000831571; *A. fumigatus* returns four viral hits alongside the fungal
one; *F. oxysporum* returns six). Scoping to the resolved taxon ID's subtree
eliminates this — confirmed live, `taxonomy_id:5334 AND
proteome_type:REFERENCE` (*N. crassa*'s taxon ID) returns exactly one record.
`proteomeType` response strings are `Reference proteome` / `Non Reference
proteome` / `Excluded` (the search enum value is `REFERENCE`) — confirmed live;
earlier drafts of this spec guessed a legacy string ("Reference and
representative proteome") that is not what the API actually returns.

- **Zero candidates** → go to Step 2 (NCBI fallback).
- **Exactly one candidate** → accept it. A `Strain` mismatch or absent strain
  metadata does NOT block acceptance — most reference proteomes are
  one-per-species, and requiring exact strain-string agreement would
  manufacture false negatives for the common case where UniProt's own strain
  string is just formatted differently from what a user typed. **Caveat added
  after review:** this is a real risk, not just a formatting nuisance, for
  species with documented accessory-genome/strain-specific biology —
  *F. oxysporum* (formae speciales differ by hundreds to thousands of
  lineage-specific-chromosome genes), *Z. tritici* (dispensable chromosomes;
  `IPO323`, already this repo's `pezizo_set1` strain, is one of the isolates
  that carries them), *A. fumigatus* (strain-specific gene islands between
  Af293 and A1163). Accepting a single UniProt candidate regardless of strain
  stays the default (a per-row strain check would need its own
  species-complex-aware logic this spec doesn't attempt), but `resolve` MUST
  emit a `Strain mismatch` warning in its report (see Output) whenever the
  accepted record's own strain string doesn't contain `species.csv`'s
  `Strain` value, so a human reviewing the summary sees it rather than it
  passing silently. This is a report-only change — it does not block
  acceptance or add a new blank-row case.
- **Multiple candidates** → try the `Strain` tiebreak: if exactly one
  candidate's UniProt strain string contains/matches `species.csv`'s `Strain`
  value (case-insensitive substring), accept that one. Otherwise — no match
  disambiguates, or more than one still does — leave the row blank and report
  all candidates (proteome ID, strain string, BUSCO completeness score if
  present, `genomeAssembly.assemblyId` if present) for a human decision. Never
  pick "the first one" or "the most recently modified one" as a silent
  tiebreaker; those aren't biologically meaningful signals. With Step 0's
  taxon-ID scoping in place, true multi-candidate cases should now mean
  "multiple strain-level reference proteomes genuinely exist for this species"
  (e.g. *C. neoformans* vs. *C. deneoformans* if a search isn't scoped tightly
  enough), not viral contamination.

Accepted UniProt match → `Protein_Source=uniprot`, `Protein_Accession=<proteome
ID>`, `Taxon_ID=<taxon ID from the record>`. Record the response's
`proteomeCompletenessReport.buscoReport.score` in the report output (not used
for any automatic decision — see Out of scope — but free to collect and
answers `DESIGN.md` Sec 9's open BUSCO-path question). If the record's own
`genomeAssembly.assemblyId` field is present (the same field `DESIGN.md` Sec 5
already established as the version-matched genome/protein pairing mechanism),
also set `Genome_Source=ncbi`, `Genome_Accession=<that assemblyId>` — one query
resolves both columns, and the genome stays guaranteed version-matched to the
protein set. **Fix from review:** this field is reliably a GCA accession
(confirmed live across every fungal record checked), sometimes a *superseded*
GCA version — e.g. UP000007431 (*S. commune* H4-8) points at GCA_000143185.1
while NCBI's current assembly is .2. `resolve` must call `datasets` on the
recorded accession to fetch its `assembly_info.assembly_status` and include a
`(superseded — current is <accession>)` note in the report when it isn't
`current`; it still writes the assemblyId UniProt actually references
(matching the protein set is the point), it just doesn't do so silently. If
the UniProt record has no linked assembly, `Genome_Source` stays blank and
Step 2's NCBI search runs for the genome only (still recording
`Protein_Source=uniprot` from this step).

**Fix from review — FTP availability check before committing to
`Protein_Source=uniprot`.** `resolve` reads the live REST API; the actual
download (`fetch_uniprot_proteome.py`) pulls from a per-release FTP snapshot
(`ftp.uniprot.org/.../reference_proteomes/`). These can disagree — confirmed
live during review: `Scom` (*S. commune*), the exact species that motivated
this whole design because it had no UniProt reference proteome as of
2026-09-10, now has one in UniProt's current release (added between
2026-09-10 and 2026-09-11). Before writing `Protein_Source=uniprot`, `resolve`
issues a `HEAD` request against that proteome's FTP directory and only accepts
the match if it resolves; otherwise it falls through to Step 2 and reports "
`REFERENCE` proteome found in REST but not yet in the FTP release — treating as
absent" so the row doesn't resolve to something the fetch step can't
download. Record the live release string this check ran against (not a
hardcoded constant — `DESIGN.md` Sec 5 already lists 2026_02; the live FTP
README read 2026_03 during review, one release later).

### Step 2 — NCBI fallback (no UniProt reference proteome, or no linked assembly)

Query NCBI Datasets (`datasets summary genome taxon <taxon ID from Step 0>
--as-json-lines`, the same CLI `bin/fetch_genome_assembly.py` already depends
on) for assemblies of this organism. **Fix from review:** query by taxon ID,
not organism-name string, for the same precision reason as Step 0/1; and use
`--as-json-lines` with server-side filters (`--reference`, `--annotated`)
rather than pulling every assembly and filtering client-side — a
species-complex query can return hundreds (*F. oxysporum*: 852 assemblies,
*A. fumigatus*: 421, confirmed live), and the report's "list all tied
candidates" output (see Output) must stay readable, not enumerate hundreds of
rows.

**Fix from review — collapse GCA/GCF pairs before ranking.** NCBI's RefSeq
assemblies are typically submitted alongside a paired GenBank assembly of the
identical sequence, and `datasets` marks **both** members of the pair
`refseq_category: reference genome` (confirmed live: *N. crassa*'s
GCA_000182925.2/GCF_000182925.2 pair, and the same pattern for *S. commune*,
*A. fumigatus*, *C. neoformans*). Naively ranking would see two "reference
genome" hits and misreport every such species as an unresolvable tie. `resolve`
groups candidates by `assembly_info.paired_assembly` before ranking (one group
per unique underlying assembly, however many accessions of it exist), then
ranks *groups*, and within a chosen group records the GCA (GenBank/INSDC)
accession — matching UniProt's own model (Step 1 always yields a GCA) and this
repo's existing `species.csv` convention, so a study's `Genome_Accession`
column doesn't mix GCA and GCF arbitrarily depending on which path resolved
it.

Rank surviving assembly groups by NCBI's own metadata:

1. RefSeq `reference genome` (NCBI's own top designation — one per species,
   when it exists).
2. RefSeq `representative genome`. **Caveat from review:** not observed as a
   populated value on any of the ~10 fungal species spot-checked live during
   review (each had either `reference genome` or no `refseq_category` at
   all) — kept in the ranking as documented NCBI Datasets vocabulary, but
   expect rank 3 to be the common path for a taxon Step 1 didn't already
   resolve.
3. Best available complete/chromosome-level assembly (RefSeq or GenBank),
   preferred over scaffold/contig level — standard assembly-QC hygiene: a
   fragmented, low-N50 assembly manufactures spurious "gene absent" calls from
   an incomplete rather than truly missing region, which is exactly the
   failure mode a novelty/loss pipeline is least equipped to tell apart from a
   real biological absence.

**Fix from review — require annotation when this rank feeds the protein
path.** If Step 1 found no UniProt proteome at all (so this genome's own
NCBI-packaged protein set becomes `Protein_Source=ncbi`), the ranking above
must additionally require `annotation_info` to be present on the candidate —
confirmed live that most assemblies for a given species are unannotated (54
total for *S. commune*, only 9 carrying `annotation_info`), and
`fetch_genome_assembly.py --include-protein` only WARNS (doesn't error) on a
missing protein file, so an unfiltered rank-3 pick can silently produce a row
with no usable protein data. When Step 1 already resolved the protein (from
UniProt) and Step 2 is only filling in the genome, this annotation
requirement does not apply.

- **Exactly one candidate group survives ranking** → accept it:
  `Genome_Source=ncbi`, `Genome_Accession=<GCA accession>`. If Step 1 found no
  UniProt proteome, also set `Protein_Source=ncbi`,
  `Protein_Accession=<same accession>`, `Taxon_ID=<taxon ID from `datasets`'
  own `organism.tax_id`>` (recorded for the `ncbi` path too — Step 1 already
  records it for `uniprot`; both paths should leave `Taxon_ID` populated for
  downstream `TaxonGroup`/lineage consistency checks).
- **Multiple candidate groups tie within the top surviving rank** (e.g. two
  distinct chromosome-level assemblies from different submitters, neither
  marked reference/representative) → leave blank, report all tied candidates
  (accession, assembly level, submitter, submission date) for a human
  decision, capped at a readable count with a note when more were found.
  Assembly choice at that point is a real judgment call (recency vs. an
  existing analysis' precedent vs. submitter reputation) this tool shouldn't
  make silently.
- **Zero assemblies found anywhere** → leave the row fully blank, report "no
  UniProt reference proteome and no NCBI assembly found for `<Species>`" —
  this is the correct terminal state for a species that turns out not to have
  usable public genomic data yet, not a bug to work around.

## Strain-targeted resolution (pangenome / multi-strain studies)

**Added 2026-09-11, same day as the spec, in response to a real use case:**
resolving several *named strains of the same species* against each other (a
within-species pangenome novelty/loss screen — which genes are
strain-specific vs. core to the species), not just one reference-quality
proteome per species. This directly leans into the exact biology Step 1's
strain-mismatch caveat above warns about (*F. oxysporum* formae speciales,
*Z. tritici* dispensable chromosomes, *A. fumigatus* strain-specific islands)
— for a pangenome screen those *are* the signal being looked for, not noise
to average away.

This is not a new command — it's a widening of Steps 1 and 2 that only
activates when the default (reference-only) search comes up empty **and**
the row's `Strain` field is non-blank (a blank `Strain` keeps today's
reference-only behavior; there's no way to auto-pick "some other strain"
without the user naming which one).

- **Step 1b (UniProt, widened):** if Step 1's `proteome_type:REFERENCE`
  search found no match (zero candidates, or a multi-candidate tie the
  `Strain` field didn't settle) and `Strain` is set, repeat the query without
  the `proteome_type` filter — `taxonomy_id:<id>` alone, so `Non Reference`
  proteomes are included — and match candidates' strain strings against
  `Strain` the same way as Step 1's tiebreak (case-insensitive substring).
  Exactly one match → accept it (`Protein_Source=uniprot`, same as Step 1,
  just not proteome-type-restricted). Zero or multiple matches → fall through
  to Step 2b.
- **Step 2b (NCBI, widened):** if Step 2's reference/representative/best-
  complete ranking found no match for this species and `Strain` is set,
  repeat the `datasets summary genome taxon <id>` query without requiring
  `--reference`, collapse GCA/GCF pairs as before, and match each surviving
  group's own strain metadata (`organism.infraspecific_names.strain` /
  the assembly's BioSample strain attribute) against `Strain`. Exactly one
  match → accept it (`Genome_Source=ncbi`, and `Protein_Source=ncbi` too if
  Step 1b also found nothing — same annotation-required rule as Step 2 when
  feeding the protein path). Zero or multiple matches → blank and report,
  same as every other unresolved case.

**Practical effect:** a `species.csv` with five rows all `Species=Fusarium
oxysporum` but five different `Strain` values (e.g. one f. sp. per row) now
resolves to five distinct genome/protein pairs instead of all five silently
collapsing onto whichever single row UniProt's reference proteome happens to
match — which is exactly the failure mode a naive implementation of Step
1/Step 2 alone would have produced for this use case. `bin/build_study_config.py`'s
existing duplicate-stem check (Task 1 of the onboarding-dispatch plan) is the
safety net if two rows ever *do* resolve to the identical accession by
mistake (e.g. the same strain listed twice) — it already hard-errors on that,
unchanged by this addition.

## Output

`species.csv` is rewritten in place (every already-filled row byte-for-byte
unchanged — `resolve` only ever adds to blank cells, consistent with being safe
to re-run). Stdout ends with a two-part summary:

```
Resolved 7 of 9 species (UniProt release 2026_03):
  Lsph  Lacrimispora sphenoides       -> uniprot UP0000XXXXX (BUSCO 98.2) + ncbi GCA_... (current)
  Easp  Enterocloster asparagiformis  -> ncbi GCF_...  (no UniProt reference proteome)
  Scom  Schizophyllum commune         -> uniprot UP001497681 (BUSCO 97.5) + ncbi GCA_023508785.1
    WARNING strain mismatch: accepted record's strain is "Tattone D", species.csv says "H4-8"
  Foxy1 Fusarium oxysporum            -> uniprot UP000xxxxxx (strain-targeted match: "Fo47")
  ...

2 species need a decision -- species.csv left blank for these rows:
  Xsp   Xylaria sp.
    Not a resolvable NCBI taxon name.
  Ysp   Yarrowia sp.
    No UniProt reference proteome.
    3 NCBI assembly groups tied at chromosome level, none marked reference/representative:
      GCA_111111111.1  (submitted 2024-03, DOE JGI)
      GCA_222222222.1  (submitted 2025-01, university consortium)
      GCA_333333333.1  (submitted 2025-06, university consortium)
```

Every accepted row's report line names which fix path produced it (default
reference match, strain-targeted widened match) and flags a strain-string
mismatch when one exists, so a human skimming the summary sees exactly what
was assumed versus what the row asked for — nothing is silently "close
enough."

## Testing

`lib/ni_resolve.py`'s UniProt/NCBI query functions are thin, mockable wrappers
(same shape as `fetch_uniprot_proteome.py`'s `fetch_json()`); the resolution
*logic* (single-match acceptance, strain tiebreak, strain-targeted widening,
GCA/GCF pair collapsing, NCBI rank ordering, the FTP-availability check, the
blank-and-report path) is unit-tested against canned API-response fixtures, no
live network calls in the test suite — mirroring `tests/test_build_study_config.py`'s
`monkeypatch`-the-network-call pattern. Fixtures must include at least one
mycovirus-contaminated `organism_name`-style response (to prove the taxon-ID
scoping actually excludes it) and one paired GCA/GCF response (to prove pair
collapsing works), since both were the mechanism of a review-caught defect,
not just a hypothetical.

**Pagination:** UniProt's search endpoint caps at 500 results per page with
cursor-based continuation; `lib/ni_resolve.py`'s UniProt query helper must
follow the cursor rather than assume a single page contains every candidate —
a species-complex query could plausibly exceed 500 non-reference proteomes,
and an ambiguity decision based on a truncated first page would be wrong in a
way this design otherwise works hard to avoid.

## Out of scope (this spec)

- Phase 2 discovery/browsing (open-ended clade queries) — future spec.
- BUSCO/completeness scoring as an automatic *decision* input (rank-breaking or
  filtering) — collected and reported (see Output) for both the default and
  strain-targeted paths, but never used to silently pick between two
  candidates; NCBI's own reference/representative designation and the
  strain-name match already provide the decision signal, and the ambiguous-tie
  path stops for a human rather than needing this tool to adjudicate assembly
  quality itself.
- Ranking or scoring *which* strains make a good pangenome comparison set
  (e.g. "pick 10 maximally-diverse *F. oxysporum* strains") — strain-targeted
  resolution above requires the user to already have the specific strain names
  in mind, matching Phase 1's overall "resolve already-chosen species" scope;
  choosing *which* strains to compare is closer to the deferred Phase 2
  discovery problem than to resolution.
- Any UI/interactive prompt loop — `resolve` is a single non-interactive pass;
  the human (or an AI session) edits `species.csv` for whatever it left blank
  and reruns.
