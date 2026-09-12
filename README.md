# NovInvenio_Investigations (NII)

The science half of the NovInvenio split: trait data, curated species/config setup,
per-study results, and analysis extensions built on top of the
[`nf_NovInvenio`](https://github.com/stajichlab/nf_NovInvenio) pipeline. NII holds
data and configs; the pipeline's code lives in `nf_NovInvenio`, invoked by reference.

See `DESIGN.md` for the full design record and `CLAUDE.md` for the operational
ruleset (data provenance rules, where new code goes). This file is a practical
quick start for onboarding a new study.

## Quick start

Every study goes through the same four stages: **generate** a `species.csv`
(the table the pipeline is built from) → **download** the data it points to →
**manually tweak** anything the automatic step couldn't resolve → **launch**
the pipeline. `bin/ni` is the one entrypoint for all four stages except the
manual tweak, which is just editing a CSV in your editor.

```
generate table        download data        manual tweaks       launch
------------------    -----------------    ------------------   ----------------
ni resolve   (1a)      ni fetch             edit species.csv     ni run
  or                     |                  by hand, then          |
ni discover  (1b) -----> |                  re-run resolve/    -> (build_study_config.py
  or                     |                  fetch if you            runs automatically
hand-write species.csv   |                  changed anything)       if not already built)
                         v                                          |
                    config.csv + data_dir/                          v
                                                                nextflow run ...
```

1. Create `studies/<domain>/<set_name>/` (`<domain>` is one of `conf/domains.yaml`'s
   slugs: `fungi`, `animal`, `plant`, `bacteria`, `other`).
2. Generate `species.csv` one of three ways:
   - **You already know your species** (a standard ingroup/outgroup comparison,
     or a pangenome set where you already know every strain): hand-write it with
     identity columns filled in and source columns blank, then run
     `bin/ni resolve` to fill in the accessions (Use case 1/2 below).
   - **You want every annotated genome of one species** (a pangenome study where
     you don't know the strain population upfront): run `bin/ni discover`
     instead — it queries NCBI directly and writes `species.csv` for you
     (Use case 3 below).
3. Download the data: `pixi run python bin/ni fetch --study-dir studies/<domain>/<set_name>`.
   This builds `config.csv` + `data_dir/` from `species.csv`.
4. **Manual tweaks, if needed** — a row `resolve`/`discover` couldn't resolve
   unambiguously is left blank with the reason printed; fix it by hand (a local
   FASTA path, a corrected strain name, a manually chosen accession), or assign
   `Group` (IN/OUT) on `discover`'s output. Re-run `ni fetch` after any edit —
   it's safe to re-run, and only touches what changed.
5. Launch the pipeline: `pixi run python bin/ni run --study-dir studies/<domain>/<set_name> [-- extra nextflow args]`.
   `ni run` is a thin pass-through to `bin/run_study.sh` — same effect, same
   script, just callable from the same `ni` entrypoint as the earlier steps.
   It builds `config.csv`/`data_dir` first if step 3 wasn't run yet, so steps
   3-5 can also be collapsed into this one command for a study with no rows
   needing manual attention.

Concrete examples of steps 2-5 for each use case are in "Use case 1/2/3" below.
The `.claude/skills/new-study/SKILL.md` skill walks through the hand-written
`species.csv` path (step 2's first bullet) with more detail on classifying
`Protein_Source`/`Genome_Source` when a species' data is already a local file
rather than something to look up.

## `bin/ni`: automatic dataset resolution

`bin/ni` fills in `species.csv`'s source columns automatically (querying UniProt and
NCBI) instead of you looking up each accession by hand, then hands off to the
existing fetch/build step. Full design: `notes/superpowers/specs/2026-09-11-ni-dataset-resolution-design.md`.

```bash
pixi run python bin/ni resolve --study-dir studies/<domain>/<set_name>
pixi run python bin/ni fetch    --study-dir studies/<domain>/<set_name>
pixi run python bin/ni run      --study-dir studies/<domain>/<set_name>
```

- `resolve` only ever touches rows where **both** `Protein_Source` and
  `Genome_Source` are blank — a row you've already filled in by hand (e.g. a
  `local_faa` species) is left alone, so it's safe to re-run after partial manual
  edits.
- It never guesses: a species it can't resolve unambiguously (multiple candidate
  proteomes/assemblies, no data found, an unresolvable name) is left blank, with the
  reason printed in the summary. Fix those rows by hand, then re-run `resolve`.
- `fetch` is a thin pass-through to `bin/build_study_config.py` — identical behavior,
  same script, just accessible from the same `ni` entrypoint.
- `run` is a thin pass-through to `bin/run_study.sh` — it builds `config.csv`/`data_dir`
  first if `fetch` wasn't run yet, then launches `nextflow run stajichlab/nf_NovInvenio`
  against the study. Any extra arguments after `--` go straight to `nextflow` (e.g.
  `-- -resume` or `-- --pfam_hmm /path/to/Pfam-A.hmm`); a study's own
  `run_params.txt`, if present, supplies its committed params first.

**First real use of this tool should be supervised** — check the resolved
accessions against what you expect before trusting a run unattended. (An earlier
review of this tool, before it shipped, caught it silently resolving nothing on
real species names and writing a `Taxon_ID` that would 404 on download — both
fixed and re-verified live, but it's a reminder this is new.)

### Use case 1: species comparison (standard ingroup/outgroup novelty/loss)

You already know the species you want; `bin/ni resolve` just looks up each one's
accessions. Leave `Strain` blank — that gets you each species' single UniProt/NCBI
reference:

```csv
Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID,Genome_Source,Genome_Accession,GFF3_Source,GFF3_Accession
Ccin,Coprinopsis cinerea,,IN,Agaricales,,,,,,,
Scom,Schizophyllum commune,,IN,Agaricales,,,,,,,
Cneo,Cryptococcus neoformans,,OUT,Tremellomycetes,,,,,,,
```

```bash
pixi run python bin/ni resolve --study-dir studies/fungi/my_new_comparison
pixi run python bin/ni fetch    --study-dir studies/fungi/my_new_comparison
pixi run python bin/ni run      --study-dir studies/fungi/my_new_comparison
```

`resolve` prefers UniProt's reference proteome per species, falling back to NCBI
when a species has no UniProt reference (this is what `Scom`, *Schizophyllum
commune*, needed).

### Use case 2: pangenome / within-species strain comparison

Same `species.csv` shape, but every row shares one `Species` and you **do** give a
real `Strain` per row. A non-blank `Strain` triggers strain-targeted widening:
`resolve` searches non-reference UniProt proteomes / non-reference NCBI assemblies
by strain name, instead of collapsing every row onto the one official reference:

```csv
Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID,Genome_Source,Genome_Accession,GFF3_Source,GFF3_Accession
Fo47,Fusarium oxysporum,Fo47,IN,Fusarium,,,,,,,
FoR4,Fusarium oxysporum,race 4,IN,Fusarium,,,,,,,
FoRad,Fusarium oxysporum,f. sp. raphani 54005,OUT,Fusarium,,,,,,,
```

```bash
pixi run python bin/ni resolve --study-dir studies/fungi/my_pangenome_study
pixi run python bin/ni fetch    --study-dir studies/fungi/my_pangenome_study
pixi run python bin/ni run      --study-dir studies/fungi/my_pangenome_study
```

Notes specific to this mode:
- `Strain` needs to be an exact-ish match (case-insensitive substring) to what
  UniProt/NCBI actually call that isolate — a vague value like `"strain 1"` won't
  resolve.
- Expect more rows to come back blank/ambiguous than in a normal cross-species run,
  since there's no "reference" designation to fall back on for a non-default strain
  — that's the tool refusing to guess which strain you meant, not a bug.
- Whether a row is `IN` or `OUT` is entirely up to how you're framing the pangenome
  question (e.g. one clinical lineage vs. everything else); `bin/ni` only resolves
  accessions, it doesn't know or care about your grouping.

### Use case 3: pangenome from population discovery

Start a pangenome study when you do not know the species' annotated-strain population
upfront. `bin/ni discover` queries NCBI for all annotated genomes of a species, groups
them by forma specialis (a taxonomic pathotype label), reports population statistics,
and writes `species.csv` ready for the pipeline:

```bash
pixi run python bin/ni discover --species "Fusarium oxysporum" --study-dir studies/fungi/my_pangenome_study
```

The tool prints a table of forma-specialis groups (e.g. `lycopersici`, `cubense`,
`no-fsp-in-name` for genomes registered without a pathotype label) with genome counts
and per-group quality metadata: annotation provider, protein-coding gene-count range,
assembly level, contig N50, and (for a collapsed GCA/GCF pair) the BUSCO score of the
discarded higher-quality member — annotation-quality variation across providers is
often a bigger confounder than the forma-specialis grouping itself. Review the table,
then assign `Group` (IN/OUT) to each strain by hand in `species.csv`, or use
`--ingroup-groups` and `--outgroup-groups` if you already know which groups to
separate (label values are whatever the printed table shows for this species — a
genome with no forma-specialis label in its name always groups under the literal
label `no-fsp-in-name`, not `others`):

```bash
pixi run python bin/ni discover --species "Fusarium oxysporum" --study-dir studies/fungi/my_pangenome_study \
  --ingroup-groups "lycopersici,cucurbitacearum" --outgroup-groups "no-fsp-in-name"
```

Pass `--auto` to request an extended report with a largest-named-group suggestion
(the `no-fsp-in-name` fallback bucket is deliberately excluded from that comparison —
it is not a coherent population, see the design spec) and cross-group duplicate
warnings — this proposal is informational only and does NOT set Group in the written
file. Nota bene: forma specialis is a host-specificity label used for plant pathogens,
not a phylogenetic split; see `notes/superpowers/specs/2026-09-11-ni-discover-design.md`
for grouping rationale and caveats.

Pass `--include-species-complex` when NCBI Taxonomy has split part of the species'
population out under a different species name within the same species-group/complex
node — the real motivating case is *Fusarium oxysporum*: querying the literal species
name misses `Fusarium odoratissimum` genomes (including the reference TR4 strain),
which are only visible by querying the parent `Fusarium oxysporum species complex`
node. Without the flag, `discover` still reports how many additional genomes exist in
the complex; with it, every child species' genomes are enumerated and each row's
`Species` column reflects whatever species name that record actually carries (not
forced under the name you queried).

Once `species.csv`'s `Group` column is filled in (by hand, or already set via
`--ingroup-groups`/`--outgroup-groups` above), continue exactly like the other
two use cases:

```bash
pixi run python bin/ni fetch --study-dir studies/fungi/my_pangenome_study
pixi run python bin/ni run   --study-dir studies/fungi/my_pangenome_study
```

## Data provenance

Every tracked data file needs a provenance record (source URL, release/version,
access date, checksum, license) — see `CLAUDE.md`'s "Data provenance & tracking"
section before committing anything new under `studies/` or `config_support/`.
