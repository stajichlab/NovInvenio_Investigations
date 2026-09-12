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

### Step 1 — UniProt reference-proteome search

Query `GET https://rest.uniprot.org/proteomes/search?query=organism_name:"<Species>"`
(the same `UNIPROT_REST` host `bin/fetch_uniprot_proteome.py` already uses),
filtered to `proteomeType` values UniProt marks as reference-quality ("Reference
and representative proteome" / "Reference proteome" — the exact field/enum values
get pinned during implementation against a live query, not guessed here).

- **Zero candidates** → go to Step 2 (NCBI fallback).
- **Exactly one candidate** → accept it. Per this session's confirmed decision
  (Q2), a `Strain` mismatch or absent strain metadata does NOT block acceptance —
  most reference proteomes are one-per-species, and requiring exact strain-string
  agreement would manufacture false negatives for the common case where UniProt's
  own strain string is just formatted differently (e.g. `"OR74A"` vs
  `"OR74A / FGSC 9010"`). Bioinformatically, a single UniProt-designated reference
  proteome for a species is authoritative regardless of exactly which deposited
  strain record it's attached to, for the purposes of ingroup/outgroup presence
  calls (strain-level polymorphism essentially never changes a gene-family
  presence/absence call at the phylogenetic distances these comparisons operate
  at).
- **Multiple candidates** → try the `Strain`-tiebreak (Q3): if exactly one
  candidate's UniProt strain string contains/matches `species.csv`'s `Strain`
  value (case-insensitive substring, the same tolerance as the single-candidate
  case), accept that one. Otherwise — no match disambiguates, or more than one
  still does — leave the row blank and report all candidates (proteome ID, strain
  string, `genomeAssembly.assemblyId` if present) for a human decision. Never
  pick "the first one" or "the most recently modified one" as a silent
  tiebreaker; those aren't biologically meaningful signals.

Accepted UniProt match → `Protein_Source=uniprot`, `Protein_Accession=<proteome
ID>`, `Taxon_ID=<taxon ID from the record>`. If the record's own
`genomeAssembly.assemblyId` field is present (the same field `DESIGN.md` Sec 5
already established as the version-matched genome/protein pairing mechanism),
also set `Genome_Source=ncbi`, `Genome_Accession=<that assemblyId>` — one query
resolves both columns, and the genome stays guaranteed version-matched to the
protein set (the whole reason Sec 5 uses this field rather than an independent
NCBI search). If the UniProt record has no linked assembly, `Genome_Source`
stays blank and Step 2's NCBI search runs for the genome only (still recording
`Protein_Source=uniprot` from this step).

### Step 2 — NCBI fallback (no UniProt reference proteome, or no linked assembly)

Query NCBI Datasets (`datasets summary genome taxon "<Species>"`, the same CLI
`bin/fetch_genome_assembly.py` already depends on) for assemblies of this
organism. Rank by NCBI's own `assembly_level`/`refseq_category` metadata (Q4):

1. RefSeq `reference genome` (NCBI's own top designation — one per species, when
   it exists)
2. RefSeq `representative genome`
3. Best available complete/chromosome-level assembly (RefSeq or GenBank),
   preferred over scaffold/contig level — standard assembly-QC hygiene:  a
   fragmented, low-N50 assembly manufactures spurious "gene absent" calls from
   an incomplete rather than truly missing region, which is exactly the failure
   mode a novelty/loss pipeline is least equipped to tell apart from a real
   biological absence.

- **Exactly one candidate survives this ranking** (i.e., ranking produces a
  unique top choice, not a tie within the same rank) → accept it:
  `Genome_Source=ncbi`, `Genome_Accession=<that accession>`. If Step 1 already
  found no UniProt proteome at all, also set `Protein_Source=ncbi`,
  `Protein_Accession=<same accession>` (this session's `Protein_Source=ncbi`
  path — the protein comes from the same NCBI Datasets genome package,
  `--include-protein`).
- **Multiple assemblies tie within the top surviving rank** (e.g. two
  chromosome-level GenBank assemblies, neither marked reference/representative)
  → leave blank, report all tied candidates (accession, assembly level,
  submitter, submission date) for a human decision. Assembly choice at that
  point is a real judgment call (recency vs. an existing analysis' precedent vs.
  submitter reputation) this tool shouldn't make silently.
- **Zero assemblies found anywhere** → leave the row fully blank, report "no
  UniProt reference proteome and no NCBI assembly found for `<Species>`" — this
  is the correct terminal state for a species that turns out not to have usable
  public genomic data yet, not a bug to work around.

## Output

`species.csv` is rewritten in place (every already-filled row byte-for-byte
unchanged — `resolve` only ever adds to blank cells, consistent with being safe
to re-run). Stdout ends with a two-part summary:

```
Resolved 7 of 9 species:
  Lsph  Lacrimispora sphenoides       -> uniprot UP0000XXXXX (+ ncbi GCA_...)
  Easp  Enterocloster asparagiformis  -> ncbi GCF_...  (no UniProt reference proteome)
  ...

2 species need a decision -- species.csv left blank for these rows:
  Xsp   Xylaria sp.
    No UniProt reference proteome.
    3 NCBI assemblies tied at chromosome level, none marked reference/representative:
      GCA_111111111.1  (submitted 2024-03, DOE JGI)
      GCA_222222222.1  (submitted 2025-01, university consortium)
      GCA_333333333.1  (submitted 2025-06, university consortium)
```

## Testing

`lib/ni_resolve.py`'s UniProt/NCBI query functions are thin, mockable wrappers
(same shape as `fetch_uniprot_proteome.py`'s `fetch_json()`); the resolution
*logic* (single-match acceptance, strain tiebreak, NCBI rank ordering, the
blank-and-report path) is unit-tested against canned API-response fixtures, no
live network calls in the test suite — mirroring `tests/test_build_study_config.py`'s
`monkeypatch`-the-network-call pattern.

## Out of scope (this spec)

- Phase 2 discovery/browsing (open-ended clade queries) — future spec.
- BUSCO/completeness scoring of candidate proteomes or assemblies beyond NCBI's
  own `assembly_level`/`refseq_category` fields — a real quality signal, but a
  second axis of judgment this spec doesn't fold in; NCBI's own reference/
  representative designation already encodes a curated quality judgment for the
  common case, and the ambiguous-tie path already stops for a human rather than
  needing this tool to adjudicate assembly quality itself.
- Any UI/interactive prompt loop — `resolve` is a single non-interactive pass;
  the human (or an AI session) edits `species.csv` for whatever it left blank
  and reruns.
