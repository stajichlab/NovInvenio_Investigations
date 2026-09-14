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
  session history; `config.csv` is git-tracked, so prior states are recoverable
  via `git log`/`git show` rather than a standalone `.bak` file). `TaxonGroup`
  filled 2026-09-13 -- see "Pipeline conventions" item 6 below.
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

## HAC/hacA reference-sequence screen -- results (2026-09-13)

Ran directly against all 293 strains' own proteomes (no need to wait for the
full tier-1 clustering pipeline): `diamond blastp` of the AF293 reference
hacA/hrmA proteins (pulled straight from the local
`NovInvenio/db/modelorgs/FungiDB-68_AfumigatusAf293_AnnotatedProteins.fasta`,
no external fetch needed) against a Short-prefixed concatenation of all 293
proteomes, plus `hmmsearch` of the PF11001 profile (extracted from the local
`NovInvenio/db/pfam/Pfam-A.hmm` via `hmmfetch AFUB_07903_YDR124W_hel`,
its actual NAME field -- `PF11001` is only the accession) against the same
set. New code: `bin/id_crosswalk.py`'s `parse_diamond_blastp_hits_per_strain`/
`parse_hmmsearch_tblout_hits`, and `bin/hac_reference_screen.py` (tested,
`results/hac_crosswalk/`).

**Real, striking finding -- PF11001 is NOT a usable presence/absence marker
for "the HAC locus."** hmmsearch found 3-9 strongly-scoring (>170 bitscore,
tight distribution) PF11001 hits in every single one of the 293 strains: it's
a broader, ancient, non-mobile paralogous domain family ("YDR124W-like,
helical bundle domain"), not something unique to the Starship-mobile
subtelomeric copy the paper calls "HAC." Screening on raw PF11001
presence would (wrongly) call every strain positive.

**The actual signal is the specific hrmA ortholog**, found by diamond blastp
identity: a clean bimodal split at ~29-39% identity (61 strains -- hits to
OTHER PF11001-family paralogs, not the true ortholog) vs. 98.6-100% identity
(9 hits, 8 strains after dedup -- see below). At a `>=50% identity, >=50%
query coverage` threshold (`hac_reference_screen.py` defaults), **hacA is
present in all 293 strains** (mostly ~100% identity -- a core, essential,
non-mobile gene, as expected for a UPR transcription factor) while **the true
hrmA ortholog is present in only 8 of 293 strains**: `Asfu_W72310, Asfu_UD1,
Asfu_C180s1, Asfu_ATCC_46645, Asfu_A1163, Asfu_Af293, Asfu_niveus,
Asfu_E165L4`. `Asfu_A1163` and `Asfu_Af293` both carrying it matches the
paper's own Table S5 (both listed with named Starship instances) -- a real,
independent validation of the method.

**Bonus catch: a real fragmented gene model**, exactly the failure mode
`rescue_pass.py` was designed for. `Asfu_E165L4` hit hrmA with TWO separate
100%-identity protein records covering complementary, non-overlapping query
ranges (residues 93-353 and 1-52) -- a gene split across an assembly/
annotation boundary, correctly still called present because the
higher-scoring fragment alone clears the coverage floor.

**Starship cross-reference (2026-09-13, same session).** Of the 8 hrmA-
positive strains, 4 matched the paper's own population by name in Table S21
(`Asfu_C180s1`, `Asfu_A1163`, `Asfu_Af293`, `Asfu_E165L4`) -- the other 4
(`W72310`, `UD1`, `ATCC_46645`, `niveus`) aren't in the paper's dataset at
all, so no Starship data exists for them. Results in
`results/hac_crosswalk/hac_starship_crosswalk.tsv` and
`table_s5_hac_rows.tsv`:

- **`Asfu_Af293` -- resolved exactly.** Table S5 tags `AF293_XP-753167.1`
  (contig `AF293_NC-007198.1:3851274-3852152`) as literally the `hrmA` gene,
  riding in `Nebuchadnezzar-h1`. Our own AF293 GFF3 coordinate for hrmA
  (`Chr5:3851274-3852415`, from the local FungiDB AF293 proteome) starts at
  the exact same base -- direct, independent confirmation, no ambiguity.
- **`Asfu_A1163` -- narrowed to a candidate set, not resolved to one.**
  Table S5 lists 6 distinct PF11001-family gene loci scattered across CEA10
  (5 named + 1 unnamed navis), each on its own Starship:
  `Osiris-h3`/`Logos-h1`/`Gnosis-h2`/`Logos-h2`/`Nebuchadnezzar-h1`
  (+ `navis10-var35`, unnamed). Table S21 shows A1163 genotyped POSITIVE for
  exactly those same 5 named Starships -- an exact match confirming A1163
  really does carry this whole multi-Starship paralog set, as the paper
  describes. But our own A1163 annotation uses the OLD Broad Institute
  assembly (`DS499599.1` contigs, `AFUB_*` locus tags), while Table S5's
  CEA10 coordinates are from a NEWER complete-genome assembly (`CP097565-568`
  contigs) -- different assembly versions of the same strain, so the exact
  genomic coordinates don't cross-walk directly (would need a liftover/
  realignment between assemblies, not attempted). What we DO have: our
  positive hit (98.6% identity, `EDP49868.1`/`AFUB_079010`) sits ~6kb from
  `AFUB_079030` -- literally the gene Pfam's own PF11001 family description
  names as the reference (`"...AFUB_079030/YDR124W-like..."`) -- confirming
  our hit is a real member of the same physical subtelomeric paralog block,
  even without pinning down which of the 5 specific Starships carries it.
- **`Asfu_E165L4` -- a genuine, unresolved discrepancy, not a bug.**
  Sequence-level screen found a true (100% identity) hrmA ortholog present
  (split across two fragments -- the same fragmented-gene-model case noted
  above). But Table S21 genotypes E165L4 NEGATIVE for all 5 known
  HAC-bearing Starships (positive only for `Osiris-h4`, not one of the
  HAC-bearing ones). Two real possibilities, not distinguished by data on
  hand: (1) a relic/solo hrmA copy left behind after its carrying Starship
  excised, or (2) an intact but different/unlisted element (Table S21 only
  genotypes 20 "high-confidence" Starships) -- either way, this shows gene-
  level presence and Starship-genotype presence can genuinely disagree, a
  real nuance for interpreting any future HAC screen result.
- **`Asfu_C180s1`**: Table S21 shows it positive for 1 of the 5 known
  HAC-bearing Starships (`Nebuchadnezzar-h1`; its other positive,
  `Osiris-h4`, is not one of the HAC-bearing ones) -- consistent with, but
  not proof of, that being the source of its hrmA copy.

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
6. **`TaxonGroup` is now filled for all 293 IN strains** (2026-09-13; see
   `bin/assign_clades.py`). 10 rows kept the paper's own hand-entered
   `DAPC_Clade_N` labels untouched; the other 283 were filled with
   `Mash_clade_N` -- Mash whole-genome distance -> classical PCoA -> k-means
   (k chosen by silhouette sweep). **Run ingroup-only** (`--groups IN`):
   running it across IN+OUT first found "k=2" as the silhouette-best split,
   but that was 294-vs-1 -- the single *A. lentulus* outgroup strain's
   Mash distance to every *A. fumigatus* strain dwarfs any within-species
   distance and swamps the clustering. Ingroup-only gave a much more
   plausible k=4 (silhouette 0.543; sizes 129/19/62/83).
   **Concordance check** against the 10 known `DAPC_Clade_N` labels: 7/10
   land with their true-label peers (`DAPC_1` -> mostly `Mash_clade_3`,
   `DAPC_3` -> mostly `Mash_clade_2`); 2 outliers (`Asfu_AF1009B`,
   `Asfu_H237`) plus `DAPC_2`'s single example landing in the same cluster as
   most `DAPC_3` strains (can't tell from n=1 whether that's a real
   Mash-indistinguishable DAPC_2/DAPC_3 boundary or just one ambiguous
   strain). Good enough to use as a working clade-stratification variable now
   -- `cooccurrence.py`'s permutation null is genuinely stratified as of this
   fill -- but Mash is a k-mer sketch distance, coarser than a real
   SNP/allele-based DAPC; treat `Mash_clade_N` as a working label, not a
   validated population-structure result, until the ParSNP path below
   confirms or revises it.

   **Better validation set found and applied (2026-09-13, same session).**
   `mbio.01092-25-s0002.xlsx` Table S21 ("Presence/absence genotyping...")
   carries a `groupID` column with real published population-cluster labels
   for 479 isolates, drawn from TWO separate prior studies bundled into that
   one table (see its `dataset` column): `clade1/2/3` from **Lofgren** et al.
   (the same DAPC scheme as the 10 `DAPC_Clade_N` rows -- consistent
   nomenclature) and `cluster1-7` from **Barber** et al. (a distinct
   population-genetic clustering, different method/study entirely). Matching
   `config.csv`'s `Strain`/`Short` against Table S21's `isolateID`/
   `originalID` (normalized alnum-only comparison) matched **254 of 293 IN
   strains** -- far more than the 10 already hand-entered. `TaxonGroup` was
   updated with priority DAPC_Clade_N (kept) > Table S21 match
   (`Barber_clusterN`/`Lofgren_cladeN`) > `Mash_clade_N` fallback: final
   counts are 10 DAPC + 252 Table-S21 + 31 Mash-only (2 OUT rows stay blank).
   Concordance against the *k=4* Mash clustering (the one actually written to
   `TaxonGroup`) is strong: 5 of 7 Barber clusters are 99-100% pure against a
   single Mash clade; only `cluster2` (82%) and `cluster4` (73%) split across
   two.

   **k=7 experiment (did not improve things).** Tried forcing
   `assign_clades.py --k_fixed 7` to match Barber's 7-cluster count directly,
   hoping to resolve the cluster2/cluster4 split. Result was NOT a clean
   improvement -- k-means at a different k is an unconstrained
   re-partition, not a refinement of the k=4 solution, so some clusters got
   purer (`cluster3`: from merged-into-clade_0 to 94% pure) while others got
   markedly worse (`cluster5`: 100% -> 71%; `cluster2`: 82% -> 46%). Net: k=4
   (the silhouette-selected value) remains the better overall match to the
   real population structure, so it's what stays in `config.csv`. Lesson for
   later: don't assume "more k = finer, strictly-better resolution" for
   Mash+k-means without checking against ground truth each time.

   **ParSNP prototype** (17-strain subset, `results/test_smoke_17strain/
   parsnp_run/`): confirmed parsnp's `-p THREADS` only parallelizes across
   partitions (`--max-concurrent-partitions`), which don't form below ~50
   genomes ("Too few genomes to run partitions of size >50. Running all
   genomes at once.") -- the actual MUM-search/LCB-alignment core ran
   single-threaded regardless of `-p`. **Final measured timing (job completed
   2026-09-13 15:42, corrected from an earlier still-running 45+ CPU-minute
   estimate):** 1.46h wall-clock (14:14:50 -> 15:42:26) for the 18-genome
   alignment (17 strains + reference), ~1M core-genome SNPs called
   (`parsnp.vcf`). `run_parsnp_ingroup293.sh` (48h/8c/32gb, `highclock` queue,
   `$SCRATCH`-staged) is written and ready for the full 293-strain run
   (outgroups excluded, same rationale as the Mash run above; reference =
   `Asfu_Af293`) but **deliberately not submitted** -- 2026-09-13 decision:
   Mash+PCoA's real concordance against 262 published labels (5/7 Barber
   clusters 99-100% pure, see below) was judged sufficient for now, and 293
   genomes already crosses the ~50-genome partition threshold the 17-strain
   prototype never reached, so that prototype's timing does not directly
   predict the 293-strain run's cost (see the design spec's Opus-review-revised
   component 10 for the full scaling discussion). Revisit only if a specific
   downstream result needs finer resolution than Mash provides.

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

- [x] Fill in `config.csv`'s `TaxonGroup` (DAPC clades) so the co-occurrence
      permutation null is genuinely clade-stratified. Done 2026-09-13 via
      Mash+PCoA (`bin/assign_clades.py`); see the "Pipeline conventions"
      section item 6 above for the method and concordance check. Revisit if/
      when the full 293-strain ParSNP run (`run_parsnp_ingroup293.sh`)
      produces a real SNP-based clade assignment to compare against.
- [ ] Determine draft-vs-long-read assembly mix across the 293 strains.
- [ ] Dereplicate strains (Mash/ANI) before any frequency count.
- [x] Run tier-1 clustering (~90% identity, per the revised design). Done
      2026-09-13: all 295 strains (293 IN + 2 OUT), no isoform collapse (checked
      assumption -- spot-checked gene:mRNA ratios across 3 strains showed ~3% or
      less excess, consistent with fungal biology's low alt-splicing rate; see
      design spec's Opus-reviewed component 10 step 1), 2,788,402 proteins ->
      47,983 tier-1 families in ~14 min wall-clock
      (`results/full_293run/tier1_cluster.tsv`).
- [ ] Genome-level tblastn/miniprot rescue pass -- still pending. Now known to be
      a real, not theoretical, need: the HAC reference screen caught one
      fragmented gene model (`Asfu_E165L4`'s hrmA split across two protein
      records) by accident on a single locus.
- [ ] Compute the real family-frequency histogram (component 1-2 of the general
      design) before fixing core/soft-core/shell/cloud cutoffs — after excluding
      low-completeness strains.
- [ ] Build the ID crosswalk (paper IDs <-> this study's IDs).
- [x] Build the strain-overlap table (this study's 293 vs. the paper's
      populations). Done 2026-09-13 as a side effect of the TaxonGroup fill:
      254/293 IN strains matched Table S21's isolateID/originalID (name
      normalization, no sequence-based check) -- see item 6 above and
      `results/clade_assignment/s21_matches.tsv`. This is a NAME-based
      overlap only, not the sequence-based ID crosswalk the next item still
      needs for actual gene-list lookups (Tables S12/S13/S19).
- [ ] Run the benchmark suite (Tables S6/S21/S12/S13/S7/S14-16 as positive
      controls, plus conserved non-mobile SM clusters and random family pairs as
      negative controls) against whichever strains overlap, scoring mmseqs vs.
      diamond before trusting either on novel candidate clusters.
- [x] Run the HAC/hacA targeted screen. Done 2026-09-13 via direct
      reference-sequence diamond/hmmsearch against all 293 proteomes
      (`bin/hac_reference_screen.py`), not the family-ID/matrix route
      `hac_screen.py` assumes -- that still needs the full tier-1
      clustering run (a separate open item above). See "HAC/hacA
      reference-sequence screen -- results" section above for the actual
      finding (hacA universal, true hrmA ortholog in only 8/293 strains).

## Component 8 (pair classification) -- real results (2026-09-14)

Full pipeline run against real 295-strain data, in dependency order:

1. Tier-1 clustering (2026-09-13): 2,788,402 proteins -> 47,983 families
   (~14 min).
2. Presence matrix, strain dereplication (123/295 reps), frequency binning:
   6,600 core / 582 soft-core / 3,943 shell / 9,238 cloud / 27,620 singleton.
3. Co-occurrence (`bin/cooccurrence.py`, rewritten for scale -- bitset
   presence vectors, an analytic prefilter, sparse BH-FDR; see the design
   spec's component 10 step 7): 25,765,431 candidate pairs screened ->
   1,259,432 clear FDR (0.05). First attempt was OOM-killed by the
   interactive session's 8GB job memory cgroup after ~10h (confirmed via
   `dmesg`, not a host-wide OOM or an algorithm bug); fixed (smaller
   survivor tuples, permutation-phase progress logging) and resubmitted as
   a real SLURM batch job (`run_cooccurrence_full293.sh`, `-p stajichlab`,
   32gb) -- finished in well under an hour once given real memory.
   **Sanity check on the 1.26M "significant" pairs**: only 547,553 (43.5%)
   also clear the clade-stratified permutation null -- the majority get
   correctly knocked back down, exactly the behavior the design's
   permutation-null safeguard was built for, not a sign of trouble.
4. `bin/build_gene_positions.py` (new): per-strain protein_id -> genomic
   position from each strain's own GFF3 CDS `protein_id=` attribute (not
   gene ID -- family membership is already keyed on protein ID, so this
   sidesteps the gene-ID-vs-protein-ID crosswalk problem entirely).
   2,788,389 positions across all 295 strains.
5. `bin/build_family_positions.py` (new): joins gene_positions.tsv with
   `tier1_cluster.tsv` into per-strain, per-family RANK positions (multi-
   copy families get a list of positions, not one -- fixes
   `synteny_windows.py`'s `linkage_fraction`, which previously assumed a
   single position per family per strain and would have silently
   mislabeled physical linkage for exactly the multi-copy families this
   study's own HAC finding showed matter most). 2,758,017 (strain, family,
   copy) positions.
6. DUF3435 (Starship captain gene, the real evidence source component 8's
   M1 fix required) `hmmsearch` against all 295 strains: 2,248 hits, every
   one of the 295 strains has at least one (median ~6-8 per strain) --
   biologically plausible for this Starship-rich species.
7. `bin/pair_classification.py` (new): classified all 1,259,432 FDR-
   significant pairs (~4.5 min real runtime, no SLURM job needed).

**Final classification** (all 1,259,432 pairs):

| Label | Count | |
|---|---|---|
| `trans_unconfirmed` | 649,996 | non-physical, fails the permutation/clade gate |
| `trans` | 522,933 | non-physical, statistically robust -- interaction/co-evolution candidates |
| `insufficient_data` | 75,816 | <5 co-carrying strains with a resolvable position |
| `unexplained_physical` | **5,242** | physically clustered, no Starship captain nearby -- the "other patterns" this component exists to surface |
| `ambiguous_linkage` | 3,029 | partial physical linkage |
| `starship_explained` | 2,416 | physically clustered AND near a captain gene -- the expected mechanism |

Physical linkage (`starship_explained` + `unexplained_physical` +
`ambiguous_linkage` = 10,687, 0.85%) is genuinely rare among FDR-significant
pairs, as expected -- most co-occurrence signal is not physical, which is
exactly the distinction component 8 was built to make visible rather than
leave conflated.

**Known partial-implementation caveats** (documented in
`pair_classification.py`'s own module docstring, not hidden): rank-window
linkage only (component 5's original "within k genes" statistic; the fuller
must-fix M3 ask -- an additional same-accessory-island check and a real
bp-distance window sized to an actual Starship footprint -- is not yet
implemented); no `taxon_scheme` restriction (must-fix M5's fix for
`TaxonGroup` mixing three incompatible clustering schemes -- `trans` calls
should be treated as provisional until that's added); no Table S5/S21
paper-coordinate cross-reference (DUF3435-hit-adjacency only).

**Genome-level tblastn rescue pass** (component 1b): submitted as a
separate SLURM batch job (`run_rescue_pass_tblastn.sh`, job 28371098,
`-p stajichlab`, 32c/64gb/4 days) -- one indexed tblastn run (47,983 tier-1
family reps vs. an 8GB combined genome database covering all 295 strains,
260,375 contigs) rather than 295 separate per-strain searches. Still
running as of this writing; not yet folded back into the presence matrix
(so nothing above yet reflects the rescue pass's fixes for fragmented gene
models like the `Asfu_E165L4` hrmA case).
