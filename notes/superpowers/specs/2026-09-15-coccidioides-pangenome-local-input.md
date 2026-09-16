# Coccidioides pangenome study: local-input onboarding + `--pipeline pangenome` fix

Date: 2026-09-15 (v2 — revised after Fable second-opinion review)
Study: `studies/fungi/coccidioides_pangenome/` (NII repo)
Pipeline: `nf_NovInvenio` (local checkout: `/bigdata/stajichlab/jstajich/projects/NovInvenio`, entry `pangenome.nf`)

## Revision note (v2)

This spec was reviewed by a second model (Fable) before being turned into an
execution plan. Verdict: **approve with changes**. All five of its findings are
folded in below (parser edge cases, threshold denominator, script-placement
rationale, test-small pass criterion, and two new blocking open questions —
`asm_stats.tsv` join key, copy-vs-symlink data volume). Also independently
verified during revision: `species.csv`'s real header (per the in-repo
`new-study` skill) includes `Taxon_ID`, which the v1 draft omitted; and the
`annotation_freeze` `pep/DNA/GFF` directories are themselves symlinks into
`unscaffolded_annotation/`, not real files — relevant to the copy-vs-symlink
question. Taxon IDs for both species (5501 / 199306) were supplied by the user
and cross-checked with `taxonkit lineage` — open question 8 is now resolved.

## Problem

`studies/fungi/coccidioides_pangenome/README.md` asks to run the recently-added,
**untested** `--pipeline pangenome` mode of `nf_NovInvenio` on 535 Coccidioides
strains (171 *C. immitis* / 364 *C. posadasii*), whole-set and per-species-split, using
data that already exists on disk as a **funannotate**-style annotation freeze — not a
UniProt or NCBI Datasets fetch, which is the only path `bin/build_study_config.py`'s
existing dispatch and every prior NII study have exercised.

Two independent problems block a straight run:

1. **Input dispatch**: no `species.csv` exists yet for this study, and the data isn't
   shaped like any prior study's fetch-driven input.
2. **Pipeline bug**: the pangenome workflow's gene-position step assumes NCBI-style
   GFF3 (`protein_id=` on CDS lines). Coccidioides' funannotate GFF3 has no
   `protein_id=` attribute at all — the pipeline will not error, it will silently
   produce empty position tables for every strain.

This spec covers both, plus the quality-filtering and species-group logistics the
README asks for. It intentionally stops short of being an execution plan — per the
user's ask, this spec goes to a second-opinion review (Fable) before being turned into
a `writing-plans`-style task list and before any pipeline code or full run.

## Fact base (verified 2026-09-15, not assumed)

### Local data shape

- `/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/annotation_freeze/20260112/{pep,DNA,GFF,CDS}/`
  — 535 files per subdir, uniform naming
  `Coccidioides_{immitis|posadasii}_{Strain}.{proteins.fa,scaffolds.fa,gff3,cds-transcripts.fa}`.
  Species and strain are both parseable from the filename alone; each strain's files
  are internally self-consistent (no cross-strain ID collisions to worry about).
  171 immitis / 364 posadasii by filename count.
- `Assembly/samples.csv`: 559 rows (one per sequencing run, not per strain), columns
  `RunAcc,Strain,BioSample,Center,Experiment,Project,Organism,FileBase,Notes,LocusTag`.
  `LocusTag` (e.g. `C1196DB`) is the internal gene-ID prefix baked into that strain's
  own GFF3/protein IDs — useful for QC cross-referencing, not needed for file parsing.
  **559 samples vs. 535 annotated strains: 24 samples have no matching
  annotation_freeze output.** Not yet root-caused (dropped at assembly or annotation
  stage — unclear which).
- Quality data already exists, no fresh BUSCO run needed:
  `Assembly/BUSCO/<strain>.AAFTF/short_summary.specific.ascomycota_odb10.<strain>.AAFTF.{txt,json}`
  (505 strains have a BUSCO dir — another shortfall vs. 535, not yet reconciled) and
  `Assembly/asm_stats.tsv` (506 rows) which already has `BUSCO_Complete`,
  `BUSCO_Single`, `BUSCO_Duplicate`, `BUSCO_Fragmented`, `BUSCO_Missing`,
  `BUSCO_NumGenes` columns alongside assembly contiguity stats (`N50`, `L50`,
  `T2T_SCAFFOLDS`, etc.) keyed by `SampleID`.
- Prior OrthoFinder pangenome comparison: `.../Pangenome/{OrthoFinder_diamond,results,
  Functional_ann,...}` — exists for later benchmarking, not inspected in depth here.
- Trees for later figure work: ASTRAL consensus
  (`Phylogeny/results/tree_cds_ascomycota/astral/consensus_aster2.nw`) and an ML tree
  (`Phylogeny/results/msa_filter_cds_ascomycota-buildtree/Cocci_cds.488taxa_ascomycota.fa.raxml.rba.raxml.bestTree`,
  488 taxa — not all 535 strains present).

### funannotate GFF3 vs. the pipeline's assumption

Confirmed by direct grep against a real file
(`annotation_freeze/20260112/GFF/Coccidioides_immitis_1M0.gff3`):

```
CDS  893 1495 . + 0  ID=C1196DB_000001-T1.cds;Parent=C1196DB_000001-T1;
```

No `protein_id=` attribute anywhere on CDS lines. The matching protein FASTA header is
`>C1196DB_000001-T1 C1196DB_000001` — i.e. the ID that actually matches the protein
FASTA is the CDS's `Parent=` (the mRNA/transcript ID), not any `protein_id`.

**Verified at full scale during Fable review**: 0 of 535 GFF3s carry `protein_id=`;
all carry `Parent=`. For strain 1M0, the CDS `Parent=` ID set equals the protein
FASTA header ID set exactly (8655/8655, zero mismatches). 40/40 strains spot-checked
have mRNA count == protein count, and no strain has a `-T2`-style multi-isoform
transcript — so isoform collapsing is not a live concern for *this* dataset, but the
fix must not assume that in general.

`nf_NovInvenio`'s `bin/pangenome_build_gene_positions.py`:
- Line 30: `_PROTEIN_ID_RE = re.compile(r"(?:^|;)protein_id=([^;\n]+)")`
- Lines 34–53: builds `{protein_id: (contig, start, end)}` purely from that regex, per
  strain, from CDS rows.
- No existence/fraction check — a GFF3 with zero regex matches produces an empty
  positions dict and the function returns normally (`gff3_path.exists()` is true, so no
  error path fires).

This is exactly the gap already flagged, but not fixed, in
`studies/fungi/Afumigatus_pangenome/NEXTFLOW_MIGRATION_NOTES.md` (§ around line 161):
> "The GFF3 `protein_id=` CDS attribute is assumed (NCBI style) by the gene-positions
> step — should add a hard error/warning when more than some threshold fraction of
> proteins lack resolvable positions, rather than letting those calls silently
> collapse to `insufficient_data` downstream."

Downstream blast radius if unfixed: `FAMILY_POSITIONS` (built from `GENE_POSITIONS`)
gets empty inputs for every Coccidioides strain, and `PAIR_CLASSIFICATION` collapses
every candidate co-gain/co-loss pair to `insufficient_data` — silently producing a
pipeline run that "succeeds" but delivers none of README goals 3–5 (co-gain/co-loss
family detection, cluster-of-genes patterns, trans-acting correlated gain/loss). This
is the single highest-risk item in this plan: it fails silently, not loudly, on a
535-strain / possibly-hours-long run if not caught first in a small test.

### `--pipeline pangenome` cluster backend

`pangenome.nf` (lines ~98–109): `--pangenome_cluster_backend` accepts `mmseqs`
(default, validated) or `diamond`, but the `diamond` branch is a **hard `error()`** —
implemented but explicitly disabled pending validation against real data (diamond
cluster-ID restoration assumption is unverified). This matches the README's own ask
("we are using mmseqs as primary... name folders so diamond could be swapped in
later") — diamond is a documented future option, not something to attempt now.

### Study onboarding: local-file dispatch already exists, nothing new needed there

`bin/build_study_config.py` (per
`notes/superpowers/specs/2026-09-11-study-onboarding-design.md`) already supports
`Protein_Source=local_faa` / `Genome_Source=local_genome` / `GFF3_Source=local_gff3`
— "path (or glob) to an existing file → plain copy + provenance record." This is
exactly the Coccidioides case for every one of the 535 strains (no fetch of any kind).
The design doc (line ~68) also confirms `Short` is required unique per row and used as
the file-stem for local-source rows — strain names like `485B-1_L_OLD_CPA0023`
already appear verbatim in filenames, so no character-set conflict exists.

**What's missing is only the generator**: nothing in `bin/` or elsewhere walks
`annotation_freeze/20260112/{pep,DNA,GFF}` and emits the `species.csv` rows. This is
correctly a *study-specific* script per this repo's CLAUDE.md ("Study-specific vs.
shared scripts" — it hardcodes the `Coccidioides_{species}_{Strain}.*` naming
convention and a specific external directory), so it belongs in
`studies/fungi/coccidioides_pangenome/bin/`, not shared `bin/`.

### No reusable precedent elsewhere in-repo

`studies/bacteria/cyanobacteria/` and `studies/fungi/fusarium_FOXY/` (both currently
untracked/in-progress) don't demonstrate a local-file `species.csv` generator pattern
beyond what `build_study_config.py`'s existing `local_*` source types already cover —
nothing to reuse from there.

## Proposed changes

### 1. Pipeline fix (in `nf_NovInvenio`, not NII — this is pipeline code)

`bin/pangenome_build_gene_positions.py`:
- Prefer `protein_id=` when present; else fall back to CDS `Parent=`. Do not mix
  both keys within one strain (pick one dialect per strain, based on which attribute
  is actually present on that strain's CDS rows).
- **Split `Parent=` on comma** before use — GFF3 allows `Parent=A,B` for CDS rows
  shared by multiple transcripts (alternative splicing at the exon level); the
  current regex would otherwise emit the literal string `"A,B"` as one bogus ID.
  Not observed in the Coccidioides data (0/40 strains checked have this), but it's
  valid GFF3 and cheap to handle correctly now.
- **Cross-check resolved IDs against the protein FASTA**, not just presence of an
  attribute. A `Parent=` (or `protein_id=`) match proves the GFF3 syntax parses; it
  does not prove the resolved ID is the one the protein FASTA actually uses (JGI/
  Ensembl-style dialects can carry a GFF3 transcript ID that differs from the
  FASTA's protein ID). The parser needs the strain's protein FASTA as an input
  (it currently only reads the GFF3) so it can compute the match fraction below
  against a real denominator, not GFF3-internal bookkeeping alone.
- Add the hard-error/warning threshold check the migration notes already call for,
  using that FASTA-based match fraction as the denominator:
  **hard error if a strain's resolved-and-FASTA-matching fraction is below 50%**
  (this is the actual "wrong dialect entirely" signal — real Coccidioides data
  resolves at ~100%, so this only fires on genuine structural mismatch), **warning
  above a 2% unresolved rate** (catches a handful of stray annotation-tool quirks
  without hard-failing a mostly-fine run). These replace the v1 draft's unverified
  "10%" guess — no such threshold convention exists elsewhere in this pipeline
  (checked via grep during Fable review), so there was nothing to match against;
  50%/2% are chosen directly from what actually distinguishes "wrong dialect" from
  "a few odd genes," not copied from a precedent that doesn't exist.
- Add/extend `nf_NovInvenio`'s existing pytest suite
  (`tests/test_pangenome_*.py` pattern already established) with **two new
  fixtures**: (a) a funannotate-style GFF3 (`Parent=` only, no `protein_id=`,
  including one `Parent=A,B` row) plus its matching protein FASTA, and (b) a
  mismatched-dialect fixture (GFF3 IDs that don't appear in the paired FASTA) to
  exercise the new hard-error path. This test file doesn't exist yet for this
  module — this is new test coverage, not an extension of something already there.

### 2. New study-specific script

`studies/fungi/coccidioides_pangenome/bin/build_coccidioides_species_csv.py`:
- Walks `annotation_freeze/20260112/{pep,DNA,GFF}`, parses
  `Coccidioides_{species}_{strain}.*` from filenames.
- Emits `species.csv` using the repo's actual header (per the in-repo `new-study`
  skill, not the abbreviated one in the v1 draft of this spec):
  `Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID,Genome_Source,Genome_Accession,GFF3_Source,GFF3_Accession`
  — `Short=<strain>`, `Species=Coccidioides {species}`, `Strain=<strain>`,
  `Group=<immitis|posadasii>` (used as the `IN`/`OUT` split point between the two
  per-species reruns, not a literal `IN`/`OUT` value — see open question 4),
  `TaxonGroup` filled from the species name parsed from the filename, `Taxon_ID`
  filled from a fixed lookup table of two entries (species-level, not per-strain):
  *Coccidioides immitis* = **5501**, *Coccidioides posadasii* = **199306** (genus
  *Coccidioides* = 5500, for reference/future use if a genus-level row is ever
  needed). **Source**: provided directly by the user (2026-09-15) and independently
  verified in this session via `module load taxonkit; taxonkit lineage` against
  the installed NCBI taxonomy dump on HPCC — `taxonkit lineage` confirms all three
  IDs resolve to the expected `...Onygenaceae;Coccidioides;Coccidioides immitis`/
  `...posadasii` lineages. This replaces the v1/v2-draft placeholder that flagged
  these as unverified. `Protein_Source=local_faa` / `Genome_Source=local_genome` /
  `GFF3_Source=local_gff3` with the corresponding freeze paths as
  `*_Accession`. Feeds straight into the existing `build_study_config.py` dispatch,
  no changes needed there.
- **Placement rationale, sharpened after Fable review**: the in-repo `new-study`
  skill is explicit that species.csv classification is *judgment*, not a fixed
  algorithm, because source-project directory naming conventions vary and are
  often ambiguous (e.g. MAG bin IDs). Coccidioides is the opposite case — one
  fully deterministic naming convention (`Coccidioides_{species}_{strain}.ext`),
  verified consistent across all 535 files with a scriptable regex, no ambiguity
  to adjudicate. That's what justifies writing a script here at all rather than
  doing this by hand per the skill's normal per-species judgment flow. It still
  goes in the study's own `bin/`, not shared `bin/`, per CLAUDE.md's "study-specific
  vs. shared" rule (this is the first study needing this exact directory shape —
  promoting a pattern to shared `bin/` is this repo's own stated response to a
  *third* study needing the same thing, not the first). To make future reuse cheap
  without pretending this is already a generalized tool: the root path and the
  filename regex are pulled to two named constants at the top of the file with a
  docstring pointing at this spec, so a future study with the same freeze-directory
  shape is a copy-and-edit-two-constants job, not a rewrite — genuine
  parameterization into a shared CLI tool is deferred until a second study actually
  needs it.
- **Correction (caught while drafting the execution plan, not by Fable):** the v2
  draft above said the three runs (whole-set/immitis/posadasii) would be handled by
  "one shared `species.csv` with a `--group` filter at run time." That's wrong —
  `pangenome.nf`'s `GROUP` column (lines ~132-140, 162-163) is compared against
  `--pangenome_ingroup_label`/`--pangenome_outgroup_label` (default `IN`/`OUT`); it's
  the pipeline's **ingroup/outgroup** split for its own internal analysis, not a
  species-subset selector, and every row with any other `GROUP` value is silently
  ignored by `pangenome.nf`'s samplesheet filter. This study has no ingroup/outgroup
  design (README doesn't call for one) — every row gets `Group=IN`, and `--group`
  is not a real flag anywhere in this pipeline. Corrected design: generate **one**
  `species.csv`/`config.csv`/`data_dir` for all 535 strains (single ~20.7 GB copy,
  `Group=IN` throughout), then produce the immitis-only and posadasii-only samplesheets
  as **filtered copies of `config.csv`** (by the `TaxonGroup` column, which does carry
  species identity) that point at the *same* `data_dir` — `pangenome.nf` resolves
  FASTA by basename against `--pangenome_data_dir`, so a smaller samplesheet
  referencing the same data_dir needs no extra copying. This needs a second small
  script: `studies/fungi/coccidioides_pangenome/bin/filter_config_by_taxon.py
  --config config.csv --taxon-group "Coccidioides immitis" --out
  config_immitis.csv` (and same for posadasii) — trivial (pandas/csv filter), but
  worth naming explicitly here since it's the actual mechanism, not a `--group` flag.
  `Short` is confirmed unique across all 535 strains, so no collision risk from
  reusing one `data_dir` across all three samplesheets. **Not yet verified**: whether
  `pangenome_profile.nf` tolerates an empty outgroup channel (no `OUT` rows at all) —
  no hard `error()`/emptiness check was found on a quick read, but this should be
  confirmed by the small-subset test run (§4), not assumed.
- Optional `--min-busco-complete <pct>` / `--min-n50 <bp>` filter, joining against
  `Assembly/asm_stats.tsv`. **Join-key fix (Fable review caught this): `asm_stats.tsv`'s
  `SampleID` column is `<strain>.AAFTF` (e.g. `M192.AAFTF`), not the bare strain
  name** — the v1 draft's proposed `SampleID==strain` join would silently match
  zero rows. The generator must strip/append `.AAFTF` (or match on the strain
  substring) to actually join. Exact quality thresholds still need a decision —
  not assumed here.
- Documented in a header docstring + this spec (see placement rationale above) —
  satisfies the "document scripts used to enable reuse" requirement from the task,
  scoped honestly to "copy-and-adapt," not "drop-in generic tool."
- Per CLAUDE.md's study-specific-script convention: hardcodes
  `NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")`
  rather than deriving paths from `__file__`.
- Writes a provenance note (source dir + freeze date `20260112`) into the study's
  `DATA_MANIFEST.yaml`/provenance sidecar per this repo's hard provenance rule — no
  raw data is fetched (all `local_*` sources), so this is a lightweight pointer
  record, not a full checksum/version record for a downloaded file.
- **Copy vs. symlink — resolved: copy** (user decision, 2026-09-15).
  `build_study_config.py` materializes `local_*` sources with
  `shutil.copyfile`/`copyfileobj`, which dereferences symlinks and writes real
  bytes; the freeze directory's `pep/DNA/GFF` entries are themselves symlinks into
  `unscaffolded_annotation/`. **Real measured total** (summed directly across all
  535 files per directory via `stat -L`, not extrapolated from one strain):
  `pep/` 2.27 GB + `DNA/` 14.95 GB + `GFF/` 3.50 GB = **~20.7 GB** copied into the
  study's `data_dir/`. `data_dir/` is gitignored (ephemeral, class 1 per this
  repo's provenance rules), so this never hits git — it's local `/bigdata` disk
  use only, and the user confirmed that's acceptable at this size. No symlink-mode
  addition to `build_study_config.py` needed for this study.

### 3. Data-quality audit (before any run, small or full)

- Reconcile the 559-samples vs. 535-annotated-strains gap: which 24 `samples.csv`
  strains have no `annotation_freeze` output, and why (dropped at assembly QC? still
  in progress? a naming mismatch the generator script would otherwise silently skip?).
  This should be a `species.csv`-generation-time report (strains found vs. strains in
  samples.csv but missing), not a separate manual step.
- Reconcile the 535-vs-505/506 BUSCO/asm_stats gap the same way (report, don't
  silently exclude). Quality bar itself is now settled: **exclude
  `BUSCO_Complete` < 90%** (see Open Questions item 2 for the full distribution
  analysis behind this number).

### 4. Test-small-first path

**Pass criterion corrected after Fable review.** The v1 draft implied a small
Nextflow subset run's `PAIR_CLASSIFICATION` output was the signal to watch — it
isn't: with only 5-10 strains, `PAIR_CLASSIFICATION`'s own statistical minimums
(`min_co_carrying`, `min_clades`) will legitimately return `insufficient_data`
regardless of whether the GFF3 fix works, so it's not a valid pass/fail check for
this bug. **The correct pass criterion is `gene_positions.tsv` row count per strain
== that strain's protein count** (i.e. the parser resolved a position for
(near-)every protein) — checked directly, not inferred from a downstream stage
several steps removed from the fix.

**Cheaper and more complete step, added per Fable review, done first**: run the
fixed `pangenome_build_gene_positions.py` directly (no Nextflow, no SLURM) against
all 535 real GFF3+FASTA pairs — this takes seconds, and directly validates the fix
against every strain in the actual dataset before any pipeline run at all. The
5-10-strain Nextflow subset run still follows, but as a pipeline-wiring/staging
smoke test, not the primary correctness check for the parser fix.

The subset run also does **not** exercise mmseqs memory/time behavior at the real
scale (~4.6M total proteins across 535 strains) or real SLURM job-count/queue
behavior — that needs a separate resource-sizing step before the full run, sized
per the `nextflow-hpcc` skill's scatter-gather guidance, run under
`NII_PIPELINE`/`stajichlab` account only once the small-scale checks above pass.

### 5. Result-folder naming

Per the README's own ask, name pangenome result directories to make a future
`diamond`-backend rerun distinguishable (e.g. `results/mmseqs_<runtag>/` rather than
bare `results/<runtag>/`) — cosmetic, no pipeline change needed since diamond is
already hard-disabled.

## Open questions (need a decision before this becomes an execution plan)

1. ~~Threshold for the GFF3 gene-position hard-error/warning~~ — **resolved**: 50%
   hard-error / 2% warning, based on the FASTA-cross-check denominator (see §1).
2. ~~BUSCO/assembly-quality cutoff for excluding strains~~ — **resolved**:
   **hard-exclude `BUSCO_Complete` < 90%.** Computed directly from
   `asm_stats.tsv` (505 of 535 strains have a row; 497 of those have a
   non-blank `BUSCO_Complete` value): median 97.2%, p10 95.5%, min 0.2%. There is
   a clean, well-separated gap: exactly 8 strains fall below 90% (range
   0.2%–83.8%; e.g. `NM_9861.AAFTF` at 0.2% complete, 707 contigs, N50=6168 —
   an evidently failed/fragmented assembly), and every strain from 90.3% up
   through 98.8% forms one smooth, continuous distribution with no comparable
   gap. No separate N50 hard cutoff: N50 correlates strongly with the same
   failure mode already caught by the BUSCO cutoff, and several strains sit at
   92–94% BUSCO completeness with excellent N50 (200–330 kb, e.g. `B3224`:
   92.6%/206 kb) — a standalone N50 filter would wrongly exclude these
   legitimately-lower-completeness-but-well-assembled genomes. **The 30 strains
   with no `asm_stats.tsv` row at all (535 annotated − 505 in the QC table) are
   included by default with a flag/warning in the generation-time report**, not
   silently excluded — absence of a QC record is not evidence of bad quality,
   and conflating the two would hide a distinct data-completeness issue (see
   item 7 below) inside a quality-filtering decision.
3. ~~Whether to investigate the 24 missing / ~30 BUSCO-gap strains before or as
   part of generating `species.csv`~~ — **resolved**: as part of, via a
   generation-time report, so it's never silently dropped.
4. ~~Whether the three runs should be three separate `species.csv` trees or one
   shared file with per-run filtering~~ — **resolved**: one shared `species.csv`,
   `--group` filter at run time (`Short` verified unique across all 535 strains).
5. Scope of pipeline-side pytest fixture work — **resolved**: two new fixtures
   (funannotate-dialect happy path incl. a `Parent=A,B` row, and a
   mismatched-dialect failure case), since this test file doesn't exist yet for
   this module.
6. ~~Copy vs. symlink for `local_*` sources~~ — **resolved**: copy, ~20.7 GB total
   (measured, see §2). User confirmed acceptable.
7. **New — `asm_stats.tsv` join key**: must join on `SampleID` == `<strain>.AAFTF`,
   not bare strain name (verified during Fable review; v1 draft had this wrong).
   The user has since added a README entry directly in the `BUSCO/` results
   directory explaining the `<strainID>.AAFTF` naming (I could not locate/read
   this file myself when checking — path or permissions may need confirming — but
   the naming pattern it documents is consistent with what's already used here).
8. ~~NCBI Taxon_ID for the two species~~ — **resolved**: *C. immitis* = 5501,
   *C. posadasii* = 199306 (genus *Coccidioides* = 5500), user-supplied and
   verified via `taxonkit lineage` on 2026-09-15 (see §2).

## Explicitly not in scope for this spec

- Figures/reports (open/closed pangenome plots, co-occurrence stats) — README goals
  2-5, downstream of getting a correct run.
- Starship/starfish cluster analysis (README goal 4) — pipeline doesn't support this
  yet per the README itself ("though no starship/starfish analysis has been run yet").
- Any diamond-backend validation work.
