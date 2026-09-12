---
name: new-study
description: Use when onboarding a new study into studies/<domain>/<set_name>/ -- inventories candidate protein/genome directories, classifies each species' Protein_Source/Genome_Source/GFF3_Source, writes species.csv, and hands off to bin/build_study_config.py + bin/run_study.sh. Use when the user wants to start a new NII study, add a new species set, or bring in a new dataset (local FASTA directories, NCBI accessions, or UniProt proteome IDs).
---

# Onboarding a new study

Reference: `notes/superpowers/specs/2026-09-11-study-onboarding-design.md` for the
full schema rationale.

**Alternative: pangenome discovery.** If you are starting a pangenome/multi-strain
study and don't already know which strains or forma-specialis groups to include, use
`bin/ni discover` instead. It queries NCBI to enumerate all annotated genomes of a
species, groups them by pathotype, and writes `species.csv` ready for manual Group
assignment or explicit selection via flags. Pass `--include-species-complex` when the
species has sibling populations registered under a different species name within the
same NCBI species-group/complex (e.g. *Fusarium odoratissimum*/TR4 genomes that a
plain *Fusarium oxysporum* query misses) -- see `bin/ni discover --help` and the
`README.md` section on `bin/ni discover` for details; the sections below assume you
already have a species list ready to resolve.

## Steps

1. **Pick `<domain>/<set_name>`.** Domain is one of `conf/domains.yaml`'s slugs
   (`fungi`, `animal`, `plant`, `bacteria`, `other`). `set_name` is a short,
   descriptive slug (matches existing studies like `pezizo_set1`,
   `UHM_lachnoNovelclade`). Create `studies/<domain>/<set_name>/`.

2. **Inventory every candidate data source the user gives you.** For each
   directory or accession list:
   - List files (`ls <dir>`). A directory of protein FASTA is a candidate
     `Protein_Source=local_faa` source; a directory of genome FASTA is a
     candidate `Genome_Source=local_genome` source.
   - For every species stem you find in a protein directory, check whether the
     *same stem* exists in a genome directory (exact filename match, not
     substring — a bin numbered `22` is not a substring match for `229`). If it
     does, that species is `local_genome`. If it doesn't, and the stem looks
     like an NCBI accession (`GCF_*`/`GCA_*`), that species is `Genome_Source=ncbi`
     with `Genome_Accession` = that accession. If neither, stop and ask the user
     where that species' genome comes from — do not guess.
   - A bare UniProt proteome ID (no local file at all) is `Protein_Source=uniprot`
     (needs a `Taxon_ID` too — look it up via UniProt if the user hasn't given
     one).
   - A species with no UniProt reference proteome, whose protein comes from the
     same NCBI Datasets genome package as its genome, is `Protein_Source=ncbi`
     with `Protein_Accession` = the same `GCF_*`/`GCA_*` accession as
     `Genome_Accession`.

3. **Get real `Species`/`Strain`/`TaxonGroup` values from the user.** Filenames
   (MAG bin IDs, accessions) are not species names. Do not fabricate a
   Latin binomial or taxon group — ask, or point to whatever taxonomy/metadata
   file the user's source project already has (e.g. a GTDB-tk summary, a
   `groups.tsv`) and confirm your reading of it with the user before writing
   `species.csv`.

4. **Write `species.csv`** with the header:
   `Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID,Genome_Source,Genome_Accession,GFF3_Source,GFF3_Accession`
   `Group` is `IN` or `OUT`. Leave `GFF3_Source`/`GFF3_Accession` blank unless the
   user has an explicit local GFF3 to attach. Setting `GFF3_Source=none` forces
   the GFF3 cell to stay empty even when `Genome_Source=ncbi` would otherwise
   supply a GFF3 from that same NCBI Datasets package — use it when a study
   wants to explicitly decline an NCBI-provided GFF3.

5. **Resolve sources for species needing UniProt/NCBI lookups.** Once `species.csv`
   has `Short`/`Species`/`Strain`/`Group`/`TaxonGroup` filled in, run:
   ```
   pixi run python bin/ni resolve --study-dir studies/<domain>/<set_name>
   ```
   This auto-fills the source columns (`Protein_Source`, `Protein_Accession`,
   `Taxon_ID`, `Genome_Source`, `Genome_Accession`) for any species classified
   as needing a UniProt or NCBI lookup (as opposed to `local_faa` or `local_genome`,
   which you classified by hand in step 2). `bin/ni resolve` only touches blank
   rows — species with `Protein_Source=local_faa` or `Genome_Source=local_genome`
   are left alone. Any rows `bin/ni resolve` reports as needing a decision (e.g.,
   multiple competing UniProt proteome entries for the same species, or no NCBI
   genome found) **must be resolved by hand** — check the reported candidates,
   pick one, or determine the species genuinely has no usable public data yet —
   before proceeding to the next step.

6. **Run** `pixi run python bin/ni fetch --study-dir studies/<domain>/<set_name>`
   (preferred -- the same `bin/ni` entrypoint used in step 5; it is a thin
   pass-through to `bin/build_study_config.py --study-dir ...`, which also
   still works directly). Fix any `ERROR:` it reports (missing file, unknown
   Source value) before moving on.

7. **Launch**: `pixi run python bin/ni run --study-dir studies/<domain>/<set_name> [-- extra nextflow args]`
   (preferred -- same `bin/ni` entrypoint as steps 5-6; a thin pass-through to
   `bin/run_study.sh <domain>/<set_name> [nextflow args]`, which also still works
   directly). See that script's own header comment for `NII_PIPELINE`/
   `NOVINVENIO_ROOT` local-checkout requirements.

## What this skill does NOT automate

Classifying `Source` values is judgment, not a fixed algorithm — directory
naming conventions vary per source project. This skill's job is to do that
classification carefully and show its reasoning, not to guess silently.
