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

## Open items

- [ ] Compute the real family-frequency histogram (component 1-2 of the general
      design) before fixing core/soft-core/shell/cloud cutoffs.
- [ ] Build the ID crosswalk (paper IDs <-> this study's IDs).
- [ ] Build the strain-overlap table (this study's 293 vs. the paper's populations).
- [ ] Run the benchmark suite (Tables S6/S21/S12/S13/S7/S14-16) against whichever
      strains overlap, scoring mmseqs vs. diamond before trusting either on novel
      candidate clusters.
- [ ] Run the HAC/hacA targeted screen once the crosswalk exists.
