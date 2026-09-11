# Study onboarding: unified source dispatch in `build_study_config.py`

Date: 2026-09-11

## Problem

Bringing a new study's genome/proteome data into `studies/<domain>/<set>/` currently
has two disconnected paths:

1. **UniProt+NCBI studies** (pezizo_set1, mushrooms_tremella, pezizo_set1_cluster,
   yeast_filamentous, zoosporic_dikarya, agaricomycetes_mmseqs,
   agaricomycetes_novelty_discovery, agaricomycetes_pairwise) use the shared
   `bin/build_study_config.py`, keyed off a UniProt proteome ID + NCBI GCA accession
   per species.csv row. This script has no branch for data that already exists as a
   local file.
2. Any study needing local, already-downloaded FASTA (a collaborator's MAGs, a
   dereplicated genome set, etc.) has had to write its own one-off
   `studies/<domain>/<set>/bin/build_<study>_config.py` — `UHM_Koxytoca` and
   `UHM_Akkermansia` each invented their own `species.csv` `Source`/`Accession`
   vocabulary, and neither is reusable by a third study.

`UHM_lachnoNovelclade` (this design's motivating case) needs a *third* variant: some
species' genomes are local files, others must be fetched from NCBI by GCF accession,
and this can differ row-by-row independent of where the protein FASTA comes from.
Three studies independently needing "local file instead of a fetch" is the signal
`CLAUDE.md`'s own "Study-specific vs. shared scripts" rule uses to promote a pattern
into shared `bin/` — this design does that promotion for the *reusable* parts, while
explicitly declining to force-fit the one case (`UHM_Akkermansia`) whose complexity
goes beyond "copy a file from a path."

## Unified `species.csv` schema

```
Short,Species,Strain,Group,TaxonGroup,
Protein_Source,Protein_Accession,
Genome_Source,Genome_Accession,
GFF3_Source,GFF3_Accession        (optional triplet -- may be blank/omitted)
```

`Protein_Source` and `Genome_Source` are independent per row — a species' protein and
genome can come from unrelated places (this is exactly `UHM_lachnoNovelclade`'s GCF
case: local protein FASTA already pulled, but the genome itself still needs an NCBI
fetch).

| Source value | Accession means | Handler |
|---|---|---|
| `uniprot` | UniProt proteome ID (Protein only; needs a companion `Taxon_ID` column, kept from today's schema) | `fetch_uniprot_proteome.py` |
| `ncbi` | GCA/GCF assembly accession (Genome/GFF3 only — NCBI Datasets returns both together) | `fetch_genome_assembly.py` |
| `local_faa` | Path (or glob) to an existing protein FASTA | plain copy + provenance record |
| `local_genome` | Path (or glob) to an existing genome FASTA | plain copy + provenance record |
| `local_gff3` | Path (or glob) to an existing GFF3 | plain copy + provenance record |
| `none` / blank | — | leaves that `config.csv` cell empty, same as today's already-optional GFF3 |

## `build_study_config.py` dispatch

For each `species.csv` row, independently for the Protein/Genome/GFF3 triplet:

```
if Source == "uniprot":      run fetch_uniprot_proteome.py --proteome-id <Accession> --taxid <Taxon_ID> ...
elif Source == "ncbi":       run fetch_genome_assembly.py --accession <Accession> ...
elif Source.startswith("local_"): shutil.copyfile(Accession, data_dir/<kind>/<stem><ext>); build_record(...)
elif Source in ("none", ""): leave config.csv cell empty
else: sys.exit(f"[{short}] ERROR: unknown {kind}_Source {Source!r}")
```

- **Stem naming**: rows with `Protein_Source=uniprot` keep today's
  `{Proteome_ID}_{Taxon_ID}` stem (self-documenting, collision-proof across strains).
  Rows with any other `Protein_Source` use `Short` as the stem (mirrors
  `UHM_Koxytoca`'s existing convention) — `Short` is already required to be unique per
  study, so no new uniqueness invariant is introduced.
- **Duplicate detection**: unchanged — still a `seen_stems` dict, now keyed by
  whatever stem that row resolved to.
- **Annotation extraction** (`extract_dat_annotations.py`, GO/Pfam/InterPro/gene
  name/description): only runs when `Protein_Source == "uniprot"` (unchanged from
  today — no `.dat` file exists for local/NCBI-only proteins). Studies with local
  proteins simply get no GO/Pfam enrichment input for those rows, same as today's
  `UHM_Koxytoca`/`UHM_Akkermansia`.
- **`--skip-fetch`**: unchanged semantics — skips the `uniprot`/`ncbi` fetch
  subprocess calls (re-copies from cache), always re-copies `local_*` files (they're
  not "fetched" so there's nothing to skip) and always re-runs annotation extraction.
- **Provenance**: `local_*` copies get a `build_record(...)` call analogous to
  `UHM_Koxytoca`'s script today (`source_url="(internal -- UCR HPCC, <dir>)"`,
  `license` supplied via a new `--local-license` flag defaulting to
  `"Internal / unpublished"`, `derived_by` naming the source path and this script).

## Migration scope

| Study | Action |
|---|---|
| pezizo_set1, mushrooms_tremella, pezizo_set1_cluster, yeast_filamentous, zoosporic_dikarya, agaricomycetes_mmseqs, agaricomycetes_novelty_discovery, agaricomycetes_pairwise | Mechanical `species.csv` column rename: `UniProt_Proteome_ID`→`Protein_Source=uniprot`+`Protein_Accession`, `GCA_Accession`→`Genome_Source=ncbi`+`Genome_Accession`. `Taxon_ID` column kept as-is. No `config.csv`/`data_dir` rebuild needed (already correct); a `--skip-fetch` dry run against the existing cache validates the migrated `species.csv` still resolves to the same files. |
| `UHM_Koxytoca` | Migrate onto the unified script. Both groups get `Protein_Source=local_faa`; `Genome_Source=ncbi` for OUT (replaces the separate `fetch_koxytoca_outgroup_dna.sh` shim — `build_study_config.py` now calls `fetch_genome_assembly.py` directly instead of assuming a pre-populated cache), `Genome_Source=local_genome` for IN. Delete `studies/bacteria/UHM_Koxytoca/bin/build_koxytoca_config.py` and `fetch_koxytoca_outgroup_dna.sh`. |
| `UHM_Akkermansia` | **Left alone.** Its `genome_map.tsv` indirection (proteome→source_genome lookup with NCBI-accession overrides), GFF3 matched by a *different* id than the proteome accession, and a byte-level prodigal header-artifact rewrite are real per-file transforms, not "copy a file from a resolved path" — forcing them through a generic dispatcher would either lose functionality or require the dispatcher to grow study-specific hooks, defeating the point of unifying. `build_akkermansia_config.py` stays as the study-specific script `CLAUDE.md`'s own rule anticipates; it already emits the same `config.csv`/`DATA_MANIFEST.yaml` shape the unified script produces, so nothing downstream needs to know the difference. |

## Onboarding skill

A project skill (`.claude/skills/new-study/`) wrapping the *judgment* step the
scripts can't do themselves — deciding, per candidate species, which `Source` value
applies:

1. Pick `<domain>/<set_name>`, create `studies/<domain>/<set_name>/`.
2. Given one or more candidate data directories (protein FASTA dirs, genome FASTA
   dirs, or bare accession lists), inventory filenames and classify each species:
   - A protein/genome file already present under a given local directory →
     `local_faa`/`local_genome` with `Accession` = that file's path.
   - A bare `GCF_*`/`GCA_*` id with no matching local file → `Genome_Source=ncbi`.
   - A bare UniProt proteome ID → `Protein_Source=uniprot` (+ `Taxon_ID` lookup).
3. Write `species.csv`.
4. Run `bin/build_study_config.py --study-dir studies/<domain>/<set_name>`.
5. Hand off to `bin/run_study.sh <domain>/<set_name> [pipeline args]`.

This is deliberately the part that stays AI-assisted rather than fully automated:
classifying "does this species already have a local file, or does it need a fetch"
requires directory inventory and matching heuristics (as done manually for
`UHM_lachnoNovelclade` below) that aren't worth hard-coding into `bin/` today. As
the classification heuristics stabilize (e.g. "always check the genome dir for an
exact-stem match before falling back to `ncbi`"), they can move from the skill's
prompt into `build_study_config.py --study-dir ... --infer-sources <dir>` itself —
out of scope for this design.

## Worked example: `UHM_lachnoNovelclade`

Inventory of the three source directories (done manually for this design; the
onboarding skill above automates this step for future studies):

- `/bigdata/stajichlab/jpere468/unknown_tree/clade1_faa/` — 4 outgroup protein FASTA,
  all `UHM*.bin.*` MAG stems.
- `/bigdata/stajichlab/jpere468/unknown_tree/lachno_ingroup_faa/` — 5 ingroup protein
  FASTA: 4 named by GCF accession (`GCF_000687555.1.faa` etc.), 1 `UHM*.bin.*` MAG
  stem.
- `/bigdata/stajichlab/jpere468/drep_herptile_95/high_quality_genomes/` — 1658 genome
  FASTA (dereplicated, lab-wide). Exact-stem match found for all 5 `UHM*.bin.*`
  species (outgroup and ingroup); no match for any `GCF_*` stem (expected — those are
  reference genomes, not local MAGs).

Resulting `species.csv` rows:

| Species stem | Group | Protein_Source / Accession | Genome_Source / Accession |
|---|---|---|---|
| `UHM1088.41098__UHM1088.41098_R.bin.56` | OUT | local_faa / clade1_faa/\<stem\>.faa | local_genome / high_quality_genomes/\<stem\>.fa |
| `UHM1210.23070__UHM1210.23070_R.bin.85` | OUT | local_faa | local_genome |
| `UHM896.23052__UHM896.23052_R.bin.22` | OUT | local_faa | local_genome |
| `UHM904.23055__UHM904.23055_R.bin.119` | OUT | local_faa | local_genome |
| `UHM207.23041__UHM207.23041_R.bin.71` | IN | local_faa | local_genome |
| `GCF_000687555.1` | IN | local_faa (already-pulled protein) | ncbi / GCF_000687555.1 |
| `GCF_002797975.1` | IN | local_faa | ncbi / GCF_002797975.1 |
| `GCF_025149125.1` | IN | local_faa | ncbi / GCF_025149125.1 |
| `GCF_900112885.1` | IN | local_faa | ncbi / GCF_900112885.1 |

9 species total (4 outgroup, 5 ingroup). `TaxonGroup`/`Species`/`Strain` display
values still need to be supplied by the user (not derivable from filenames alone) —
left as a TODO in the implementation plan, not this design.

## Out of scope

- Automating the `Source` classification heuristic itself into `bin/` (noted above).
- Any `bin/ni` CLI entrypoint — no such thing exists in this repo yet; this design
  only touches the existing `bin/build_study_config.py`/`species.csv` contract that
  such a future CLI would itself call.
- `UHM_Akkermansia` migration (explicitly declined above).
