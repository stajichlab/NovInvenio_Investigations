# Pangenome cluster-profile analysis — study-specific notes

General method/design: see NII's
`notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md`. This file
holds the parts specific to running that method against this study's actual data:
open questions, source data for the targeted screens, and controls to check before
trusting a result.

## Dataset

- `config.csv`: 293 ingroup *A. fumigatus* strains, 2 outgroup (*A. lentulus*
  `Aslen_ref`, *A. fischeri* `Neofi_ref`). Fixed 2026-09-13 (GROUP column had
  `INGROUP` instead of `IN`, plus 10 rows with no GROUP at all — see NovInvenio
  session history; backup at `config.csv.bak`).
- Reference/validation data: `mbio.01092-25-s0002.xlsx` (Gluck-Thaler et al. 2025,
  mBio, doi:10.1128/mbio.01092-25 — provenance recorded in `DATA_MANIFEST.yaml`).
  25 supplementary tables covering Starship (giant transposon) annotation across a
  519-strain *A. fumigatus* population plus 13 reference-quality and 3 fully
  manually-annotated (AF293/A1163/CEA10) strains.

## HAC / hacA targeted screen — source data

From `mbio.01092-25-s0002.xlsx`:

- **Table S19** (virulence/stress-resistance gene catalog): `HacA` / `Afu3g04070`
  (XP_748727.1) and `hrmA` / `Afu5g14900` — both listed as independently
  characterized loci, no shared `cluster` tag between them. Confirms `hacA` and
  `hrmA` are **not** one physical cluster (different chromosomes: Afu3g.. vs.
  Afu5g..).
- **Table S5** (Starship-predicted feature sequences, 519-strain population): the
  actual "HAC" target — genes annotated *"Subtelomeric hrmA-associated cluster
  protein AFUB_079030/YDR124W-like"* (EMAP=ENOG503P78B, IPR=IPR047092/IPR021264,
  PFAM=PF11001). Documented instances riding in distinct named Starships:
  `Nebuchadnezzar-h1` (AF293, contains `hrmA` itself — AF293_XP-753167.1,
  contig `AF293_NC-007198.1:3851274-3852152`), `Osiris-h3` (CEA10_g3677.t1),
  `Logos-h1` (CEA10_g7203.t1), `Gnosis-h2` (CEA10_g7270.t1), `Logos-h2`
  (CEA10_g5575.t1), `Nebuchadnezzar-h1` again in CEA10 (CEA10_g5946.t1), and one
  unnamed-navis instance `navis10-var35` (CEA10_g4899.t1).
- Screen target: presence/copy-number of the PF11001/IPR047092 family across all
  293 strains, plus (where an overlapping strain has Starship annotation in the
  paper) which named Starship it currently rides in. `hacA` (Afu3g04070) tracked
  as an independent single-locus presence/absence check, not assumed linked.

## Benchmark-suite controls specific to this run

- **ID crosswalk (required, not yet done).** This study's own protein/gene IDs
  (from the NCBI/UniProt-sourced annotations in `data_dir/pep`, `annotations/`)
  need a sequence-based crosswalk (diamond/mmseqs, not string matching) against
  the paper's mixed ID schemes (`AFUB_*`, `Afu*g*`, per-strain `g#` locus tags,
  `XP_*` accessions) before any Table S12/S13/S19 gene list can be checked against
  this study's own clustering output.
- **Strain overlap (required, not yet done).** Cross-reference `config.csv`'s 293
  `Short`/`Strain` values against the paper's isolate/genome-code strain names
  (Table S2 "Metadata for publicly available *A. fumigatus* strains", Table S1 for
  the 11 newly Nanopore-sequenced strains) to determine which of the 293 strains
  here actually appear in the paper's population — the benchmark suite (Tables
  S6/S21/S7) can only be scored against strains present in both.
- **Assembly/annotation completeness.** No completeness metric (e.g. BUSCO) is
  currently attached per strain in `config.csv`/`DATA_MANIFEST.yaml`. Before
  interpreting a "missing" shell/cloud family call (including a missing HAC-family
  member) as a real loss, check whether that strain's assembly/annotation is
  markedly less complete than the panel median.
- **Draft vs. long-read assembly mix.** Not yet determined how many of the 293 are
  Illumina/short-read drafts vs. Nanopore/PacBio long-read assemblies (this
  study's own 11 newly Nanopore-sequenced strains from Table S1 are one known
  long-read subset). This directly affects how much the genome-level rescue pass
  (fragmented/split gene models) and synteny contig-edge exclusion matter here.
- **Strain dereplication.** Not yet run — check for the same isolate appearing
  under two names/accessions (Mash/ANI) before computing any frequency.

*(The three items above were raised by an independent bioinformatics review of the
general design on 2026-09-13; see that design doc's "Review disposition" section
for the full must-fix/should-consider list — most of it changes the method itself,
not just this study's data, so it's recorded there rather than duplicated here.)*

## Pipeline conventions (established by `bin/build_presence_matrix.py`)

The scripts here form a chain; these are the contracts that hold it together.

1. **Short-prefixed protein IDs.** The clustering input FASTA concatenates every
   strain's isoform-collapsed proteome with headers rewritten to
   `><Short>|<original_protein_id>`. The cluster TSV carries only sequence IDs,
   so this prefix is the only way a family member can be traced back to a
   strain. `bin/build_presence_matrix.py`'s docstring holds the prefixing
   one-liner; `--id_sep` overrides the `|`.
2. **Family ID = tier-1 cluster representative ID, verbatim.** Everything
   downstream (`rescue_pass.py --matrix`, `cooccurrence.py`,
   `hac_screen.py --hac_family_id/--haca_family_id`) keys on it.
3. **Matrix columns.** `build_presence_matrix.py --groups` selects them and now
   **defaults to `IN,OUT`**, since `cooccurrence.py`'s outgroup gain/loss
   polarization reads the outgroup columns out of this same matrix. If a
   matrix is ever built ingroup-only (`--groups IN`), `cooccurrence.py`
   detects the missing outgroup columns, warns loudly on stderr, and reports
   every family's direction as `ambiguous` rather than silently mislabeling
   everything `gain` (a real bug caught by the branch's final review and
   fixed — see `notes/superpowers/plans/2026-09-13-pangenome-cluster-profile-plan.md`'s
   SDD ledger). Frequency and co-occurrence statistics themselves are always
   computed over the ingroup only.
4. **Copy numbers live in a sidecar**, `<matrix>.copy_number.tsv`, written and
   read automatically by `PresenceMatrix.to_tsv`/`from_tsv`. The matrix file
   itself stays a pure three-state table.
5. **Dereplication is opt-in per run.** `frequency_bins.py --inventory` and
   `cooccurrence.py --inventory` take `dereplicate_strains.py`'s
   `strain_inventory.tsv` and count only `is_representative == 1` strains —
   the spec's "frequency counts use dereplicated strains, not the raw 293".
   Without `--inventory` both fall back to all ingroup strains.
6. **`TaxonGroup` is currently empty for 285 of the 295 config rows**, so
   `cooccurrence.py`'s clade-stratified permutation null degenerates to an
   unstratified shuffle. The script now warns loudly on stderr when that is the
   case; `permutation_p` is not a phylogenetic control until DAPC clade labels
   are filled in.

Run order:

```
collapse_isoforms (NovInvenio) -> Short-prefixed all_ingroup.fa
  -> cluster_backend.py mmseqs-tier1|diamond-tier1
  -> build_presence_matrix.py        (presence_matrix.tsv [+ .copy_number.tsv])
  -> rescue_pass.py                  (presence_matrix.rescued.tsv)
  -> frequency_bins.py               (frequency_table.tsv)
  -> cooccurrence.py / synteny_windows.py / hac_screen.py
```

Scripts that read NovInvenio's `config_parser` locate the sibling checkout via
`lib/novinvenio_path.py`; set `NOVINVENIO_ROOT` if it is not a sibling of this
repo's parent directory.

## Open items

- [ ] Fill in `config.csv`'s `TaxonGroup` (DAPC clades) so the co-occurrence
      permutation null is genuinely clade-stratified.
- [ ] Determine draft-vs-long-read assembly mix across the 293 strains.
- [ ] Dereplicate strains (Mash/ANI) before any frequency count.
- [ ] Run tier-1 clustering (~90% identity, per the revised design) plus the
      genome-level tblastn/miniprot rescue pass; isoform-collapse first.
- [ ] Compute the real family-frequency histogram (component 1-2 of the general
      design) before fixing core/soft-core/shell/cloud cutoffs — after excluding
      low-completeness strains.
- [ ] Build the ID crosswalk (paper IDs <-> this study's IDs).
- [ ] Build the strain-overlap table (this study's 293 vs. the paper's populations).
- [ ] Run the benchmark suite (Tables S6/S21/S12/S13/S7/S14-16 as positive
      controls, plus conserved non-mobile SM clusters and random family pairs as
      negative controls) against whichever strains overlap, scoring mmseqs vs.
      diamond before trusting either on novel candidate clusters.
- [ ] Run the HAC/hacA targeted screen once the crosswalk exists.
