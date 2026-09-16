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
- **Assembly/annotation completeness -- done 2026-09-15, see results section below.**
  Genome-level BUSCO (fungi_odb12) pulled from the lab's precomputed BFD archive
  for all 295 strains, not re-run. Completeness is uniformly high (96.7-99.4%
  complete, median 98.8%) -- the panel is NOT a completeness confound; a "missing"
  shell/cloud family call can't be waved away as a low-quality-assembly artifact
  for any of these 295 strains.
- **Draft vs. long-read assembly mix -- partial answer, byproduct of the BUSCO pull.**
  The same BFD BUSCO summaries carry scaffold/contig counts: 279/295 strains have
  >100 contigs (consistent with short-read/draft assemblies), only 13 have <=20
  (consistent with long-read/near-chromosome-level). This is a contig-count
  heuristic, not a confirmed sequencing-technology determination (no read-type
  metadata cross-referenced) -- good enough to flag which strains warrant more
  caution for the rescue pass's fragmented-gene-model logic and synteny's
  contig-edge exclusion, not a final answer to the open item.
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
- [x] Assembly/annotation completeness (BUSCO). Done 2026-09-15: matched all 295
      strains' `Genome_Accession` (2 needed a manual GCA->GCF crosswalk -- see
      results section below) against the lab's precomputed BFD archive
      (`/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs`, genome-level
      `BUSCO_genome`, fungi_odb12 lineage) instead of re-running BUSCO. See
      "Assembly/annotation completeness (BUSCO) -- results" section below.
- [~] Determine draft-vs-long-read assembly mix across the 293 strains. Partial:
      a contig-count heuristic came free from the BUSCO pull above (279/295
      >100 contigs, 13 <=20) -- not a confirmed sequencing-technology call, real
      read-type metadata still needed to close this out.
- [ ] Dereplicate strains (Mash/ANI) before any frequency count. Note: component
      8's real run (below) already used a 123/295-representative dereplication
      via `dereplicate_strains.py` -- this checkbox describes the general open
      item, not evidence that step never happened; reconcile wording next time
      this file is touched.
- [x] Run tier-1 clustering (~90% identity, per the revised design). Done
      2026-09-13: all 295 strains (293 IN + 2 OUT), no isoform collapse (checked
      assumption -- spot-checked gene:mRNA ratios across 3 strains showed ~3% or
      less excess, consistent with fungal biology's low alt-splicing rate; see
      design spec's Opus-reviewed component 10 step 1), 2,788,402 proteins ->
      47,983 tier-1 families in ~14 min wall-clock
      (`results/full_293run/tier1_cluster.tsv`).
- [x] Genome-level tblastn/miniprot rescue pass -- a real, not theoretical,
      need (the HAC reference screen caught one fragmented gene model,
      `Asfu_E165L4`'s hrmA split across two protein records, by accident on
      a single locus). Done 2026-09-15, after two real bugs found and fixed
      the same day: (1) `-max_target_seqs 5` truncated across strains, not
      per strain (fixed with a per-strain-database redesign, job 28399806);
      (2) rescue-pass presence calls had no resolvable genomic position for
      `pair_classification.py` (fixed with `extract_rescue_positions.py`).
      See "Rescue pass redo COMPLETE" below for the full corrected numbers.
- [x] Compute the real family-frequency histogram (component 1-2 of the general
      design) before fixing core/soft-core/shell/cloud cutoffs — after excluding
      low-completeness strains.
- [x] Build the ID crosswalk (paper IDs <-> this study's IDs). Done
      2026-09-15 -- see "ID crosswalk + benchmark scorecard -- real results"
      below. 614/616 requested AF293 RefSeq accessions (Tables S19/S5/S13)
      resolved via NCBI efetch + diamond blastp against this study's own
      2.79M-protein proteome.
- [x] Build the strain-overlap table (this study's 293 vs. the paper's
      populations). Done 2026-09-13 as a side effect of the TaxonGroup fill:
      254/293 IN strains matched Table S21's isolateID/originalID (name
      normalization, no sequence-based check) -- see item 6 above and
      `results/clade_assignment/s21_matches.tsv`. This is a NAME-based
      overlap only, not the sequence-based ID crosswalk the next item still
      needs for actual gene-list lookups (Tables S12/S13/S19).
- [~] Run the benchmark suite (Tables S6/S21/S12/S13/S7/S14-16 as positive
      controls, plus conserved non-mobile SM clusters and random family pairs as
      negative controls) against whichever strains overlap, scoring mmseqs vs.
      diamond before trusting either on novel candidate clusters. Partial:
      done 2026-09-15 for positive controls (S6/S21/S13) and one negative
      control (S19 conserved SM/virulence genes), scored against BOTH a real
      mmseqs run (rescued matrix) and a real diamond-tier1 run (raw,
      unrescued matrix) -- see "ID crosswalk + benchmark scorecard -- real
      results" below. NOT done: S12/S7/S14-16 (13-reference-strain and
      manual-annotation tables -- most of those strains' gene IDs are the
      paper's own internal locus numbering with no public accession, so they
      are not sequence-crosswalkable without the paper's own genome
      annotations, which are not in this supplement) and the random-family-
      pair negative control (deferred -- `cooccurring_pairs.rescued.tsv` is
      1.3GB, sampling it for a fair matched-frequency negative-control draw
      needs more than a quick grep and was not attempted this session).
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

*(Note, 2026-09-15: job 28371098 was subsequently cancelled and replaced by a
chunked/compressed rescue-pass implementation run in a separate session
-- see `git log` for `Afumigatus_pangenome: chunked/compressed tblastn rescue
pass` and `Afumigatus_pangenome: rescue pass folded in, real result +
re-running downstream`. This paragraph is left as-written for its own
history; the current rescue-pass method/result is whatever that later
session documented, not this one.)*

**Rescue pass CORRECTNESS BUG found and fixed, full redo in progress
(2026-09-15, later same day).** The "100,952 ABSENT -> GENOME_ONLY,
singleton fraction 57.6%->30.7%" result above (the `bf15230` commit) is
**known wrong (an undercount) and superseded** -- do not cite those numbers.

- **Root cause**: `run_rescue_pass_tblastn.sh` / `_chunked.sh` used
  `-max_target_seqs 5` against ONE combined genome database covering all
  295 strains. That flag caps distinct subject sequences **across the
  whole combined database, per query** -- not 5 per strain. Verified
  directly against the completed run's own output: of 23,747 queries with
  any hit, 21,223 (89.4%) hit **exactly** 5 distinct strains, the
  signature of systematic truncation. Any family genuinely present in more
  than a handful of strains (i.e. most shell/core families -- the common
  case, not an edge case) had its true rescue candidates silently dropped
  once the first 5 strains it happened to encounter filled the budget.
  The comment justifying the flag ("capping hits per query keeps output
  size sane without changing the presence/absence call itself") was
  incorrect.
- **Fix, not a patch**: redesigned from one combined-genome-DB search to a
  per-strain search (`bin/extract_absent_family_queries.py` +
  `run_rescue_pass_per_strain.sh`, replacing
  `run_rescue_pass_tblastn.sh`/`_chunked.sh`) -- each strain gets its own
  tiny genome database and is searched only against the families currently
  `ABSENT` in that strain (11,372,509 (family, strain) entries across 295
  strains, ~80% of all cells, matching the earlier-measured ABSENT
  fraction). This removes the cross-strain competition for hits entirely
  (no more shared budget to exhaust), and `-max_target_seqs 200` per
  strain is now cheap since each database is ~1/295th the size -- 200 is
  generous headroom over any observed intra-genome copy number (e.g. the
  DUF3435 captain gene's median 6-8 copies/strain). 24-way SLURM array job
  (strains dealt round-robin by descending query count for load balance),
  each task's database build + tblastn run staged on `$SCRATCH`
  (node-local) with output zstd-compressed as it's produced, only the
  final `.tsv.zst` copied back to shared storage. Submitted as job
  28399806. The prior (wrong) `results/rescue_pass/chunks/*.tsv.zst`
  outputs were moved to `results/rescue_pass/chunks_pre_fix_backup/`, not
  deleted.
- **Downstream consequence**: job 28392739 (the co-occurrence re-run
  against the wrong rescued matrix) was killed before completion once this
  was found. `presence_matrix.rescued.tsv`, `frequency_table.rescued.tsv`,
  and everything computed from them (the 30.7% singleton figure) must be
  regenerated once job 28399806's per-strain tblastn results are folded
  back in via `bin/rescue_pass.py` (unchanged -- it already accepts
  repeated `--tblastn_tsv` regardless of how the search was partitioned).
- Found by an independent Fable-model pipeline review requested in
  preparation for the eventual Nextflow migration (see the migration
  section below), then verified directly against the run's own tblastn
  output before acting on it.

**Rescue pass redo COMPLETE and folded through the full chain (2026-09-15,
later same day).** Job 28399806 (per-strain tblastn) finished cleanly (295/295
strains, all array tasks COMPLETED). `run_post_rescue_pipeline.sh` (job
28409589) then folded it back through frequency binning, co-occurrence, pair
classification, and figures -- **47m54s total wall-clock** for the whole
chain (vs. the ~30-90h the old Monte Carlo co-occurrence step alone would
have needed).

**Real, now-trustworthy composition result:**

| Bin | Pre-rescue (wrong, bf15230) | Post-rescue (corrected) |
|---|---:|---:|
| core | 13.8% (6,600) | **53.0% (25,436)** |
| soft_core | 1.2% (582) | 1.3% (643) |
| shell | 8.2% (3,943) | 14.6% (6,998) |
| cloud | 19.3% (9,238) | 8.4% (4,054) |
| singleton | 57.6% (27,620) | **22.6% (10,852)** |

Singleton fraction more than halved, core more than tripled -- a far
stronger correction than the earlier (also-buggy) rescue attempt's
57.6%->30.7% claim suggested, confirming the annotation/coverage-artifact
concern (Opus review must-fix M6) was real and substantial, not marginal.

**Second bug found and fixed the same day: `pair_classification.py` couldn't
resolve most of the newly-significant pairs.** `insufficient_data` jumped
from 6.0% (pre-rescue) to 96.2% of FDR-significant pairs immediately after
the corrected co-occurrence rerun -- traced to a real cause: `classify_pair()`
only counts a strain as "co-carrying with a resolvable position" via
`family_positions.tsv`, which was built purely from GFF3-annotated
`protein_id=` positions. A rescue-pass `GENOME_ONLY` call is a raw tblastn
genomic hit with no annotated protein, so it had no position at all --
correctly detected as present by the rescue fix, but invisible to the
physical-linkage classifier. Fixed with a new script,
`bin/extract_rescue_positions.py` (re-parses the same per-strain tblastn
output already produced, keeping the single highest-bitscore qualifying hit's
own genomic coordinates -- contig + start/end -- for every (family, strain)
pair the final matrix calls GENOME_ONLY; 5,871,769 positions recovered this
way), plus a `build_family_positions.py` change (`--rescue_positions`,
merges these directly into the SAME per-strain rank ordering by genomic
coordinate, interleaved with real annotated genes, not appended after them).
Rerunning `pair_classification.py` against the corrected
`family_positions.rescued.tsv` (8,629,786 total positions, up from
2,758,017) dropped `insufficient_data` to **0.3% (13,643 pairs)** and grew
every physically-resolvable category in absolute terms too:

| Classification | Pre-rescue | Post-rescue, positions bug | Post-rescue, FULLY corrected |
|---|---:|---:|---:|
| trans | 522,933 | 81,500 | **2,613,303** |
| trans_unconfirmed | 649,996 | 71,950 | 1,533,331 |
| unexplained_physical | 5,242 | 2,221 | 44,576 |
| ambiguous_linkage | 3,029 | 4,208 | 29,949 |
| starship_explained | 2,416 | 3,327 | 8,890 |
| insufficient_data | 75,816 | 4,080,486 | 13,643 |
| **Total FDR-significant pairs** | 1,259,432 | 4,243,692 | 4,243,692 |

Canonical, current files (in `results/full_293run/`):
`presence_matrix.rescued.tsv`, `frequency_table.rescued.tsv`,
`cooccurring_pairs.rescued.tsv`, `family_positions.rescued.tsv`,
`pair_classification.rescued.tsv` (this is the position-complete v2 --
the position-incomplete intermediate was deleted after its numbers were
recorded in the table above, not kept as a standing file; `rescue_positions.tsv`
is likewise deleted once merged into `family_positions.rescued.tsv`, since
both are fully regenerable from the tblastn chunks in
`results/rescue_pass/per_strain_chunks/`). Full aggregated report:
`results/SUMMARY.md` (`bin/build_summary_report.py`). Figures: `results/figures/`.

*(Both fixes found by checking real numbers against a specific, verifiable
hypothesis before trusting them -- not by assuming a plausible-looking
result was correct. See `run_rescue_pass_per_strain.sh` and
`extract_rescue_positions.py`'s own module docstrings for the full technical
detail of each.)*

**Third finding, same day: `extract_rescue_positions.py`'s 68-minute runtime
was an algorithmic bug, not a workload that genuinely needed
multiprocessing.** Initially parallelized the per-file tblastn parsing
(`--processes`, multiprocessing.Pool) expecting that to be the win. Real
profiling showed parsing all 295 files takes only ~17s of CPU time even
sequentially -- the actual 68 minutes came from `family not in
matrix.families`, an O(47,983) linear scan on a plain Python list, executed
once per entry in the hit-position map (hundreds of thousands to millions of
times at real scale). Converting `matrix.families`/`matrix.strains` to sets
once before that loop gave a measured **~7.7x speedup on a 20-file subset
(290s -> 38s), identical output**. The `--processes` multiprocessing was
kept on top of that fix (a genuine, smaller further win -- ~1.2x on the same
subset -- once the real bottleneck was gone), but the lesson for next time:
profile before parallelizing, since the obvious "big loop over many files"
shape can hide an unrelated O(N*F) bug that parallelism only ever masks a
fraction of, never fixes.

**Separately, `cooccurrence.py`'s within-clade permutation null was
replaced with an exact closed-form test (2026-09-15).** The Monte Carlo
implementation (`permutation_null_pvalue`, 200 permutations/pair) was
projected to take ~30-90h single-core at the real rescued-matrix scale
(3.3M FDR-significant pairs). An external review found that the
within-clade shuffle keeps every margin fixed, so the null distribution of
the overlap count is exactly a sum of independent per-clade
Hypergeometric draws -- no sampling needed. `exact_stratified_pvalue`
(same file) computes this via cached hypergeometric PMFs + `np.convolve`;
benchmarked ~460x faster (0.19h vs ~90h projected for 3.3M pairs) and
verified against brute-force enumeration and independent scipy references
by a separate Opus-model review (six findings, all addressed: bounded the
PMF cache size, clamped the return value to <=1.0, and added edge-case
tests for saturated/degenerate clades and >2 unequal clades). The
`permutation_p` output column name is unchanged -- same statistical
target (a within-clade-stratified association test), computed exactly
instead of by sampling. `--n_perms` is still accepted on the CLI for
backward compatibility but no longer consumed. All 150 study tests pass.
This change is independent of the rescue-pass bug above and stays correct
regardless of which rescued matrix eventually feeds it.

## Assembly/annotation completeness (BUSCO) -- results (2026-09-15)

Ran against the design doc's open "Assembly/annotation completeness" control
(component 10 gate, "assembly-quality checks before the full run"). Per-strain
BUSCO was **not re-run** -- pulled from the lab's existing precomputed archive
instead, since it already covers this exact strain set.

**Source and method:**
- `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/results/genome_stats/BUSCO_genome/`
  (genome-mode BUSCO, `fungi_odb12` lineage, 1,122 BUSCOs, BUSCO v6.0.0,
  `miniprot` gene predictor) -- sharded by hash prefix, files named by the
  full NCBI assembly accession (`<GCA|GCF>_<version>_<assembly_name>`).
- Matched by `species.csv`'s `Genome_Accession` column against the BFD run's
  own `samples.csv` `ASMID` column (accession prefix match, ignoring the
  assembly-name suffix): **293/295 matched directly on the `GCA_*` accession
  this study uses.**
- **2 strains needed a manual crosswalk**, both reference genomes where BFD's
  archive holds the RefSeq (`GCF_*`) record instead of the GenBank (`GCA_*`)
  one this study's `species.csv` points to, and (for `Asfu_Af293`) the version
  suffix also differs:
  - `Neofi_ref` (*A. fischeri* outgroup): study uses `GCA_000149645.4`,
    matched to BFD's `GCF_000149645.3_ASM14964v4` -- same underlying assembly
    (`ASM14964v4`), GenBank/RefSeq pair, no real version conflict.
  - `Asfu_Af293`: study uses `GCA_000002655.1`, matched to BFD's
    `GCF_000002655.1_ASM265v1` -- same assembly, GenBank/RefSeq pair.
  Both are simple GCA<->GCF crosswalks of the identical assembly (confirmed
  by matching assembly-name suffix, e.g. `ASM265v1`), not a different genome
  version -- not a red flag, just BFD's own accession-source convention for
  these two well-known reference strains.
- **Result: 295/295 strains matched, zero left needing a fresh BUSCO run.**
  Raw summary files copied verbatim into
  `results/busco_genome/<Short>.BUSCO_summary.fungi_odb12.txt` (gitignored,
  per `results/` convention -- regenerable from BFD + this method) and parsed
  into a consolidated table (Complete/Single/Duplicated/Fragmented/Missing %,
  n_BUSCOs, scaffold/contig counts, total length, percent gaps -- one row per
  strain), also written to `results/busco_genome/busco_completeness_summary.tsv`.
  **The consolidated table is additionally committed at the study root**,
  `busco_completeness_summary.tsv` (same tier as `config.csv`/
  `DATA_MANIFEST.yaml` -- small, doesn't scale with candidate/sequence count,
  so it's kept tracked outside the `results/` ignore rule). Provenance: derived
  2026-09-15 by matching `species.csv`'s `Genome_Accession` against
  `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs`'s `samples.csv`
  `ASMID` + `results/genome_stats/BUSCO_genome/` (fungi_odb12, BUSCO v6.0.0,
  miniprot) exactly as described above -- no checked-in script yet, this
  section is the derivation record until one exists.

**Finding: completeness is uniformly high, not a confound for this panel.**
Complete BUSCO% across all 295 strains: min 96.7%, median 98.8%, max 99.4% --
a narrow, high range. **This directly answers the open control's question**:
a "missing" shell/cloud family call (including a missing HAC-family member)
for any of these 295 strains is very unlikely to be an assembly-completeness
artifact -- the panel doesn't contain any markedly-incomplete outlier
assemblies. Lowest-completeness strains, for reference: `Asfu_E175s2`
(96.7%), `Asfu_C172L1` (97.3%), `Asfu_CNMCM8714`/`Asfu_CNMCM8812` (97.4-97.7%)
-- still high in absolute terms, not a quality gate concern.

**Bonus, not the primary ask: a contig-count heuristic for draft-vs-long-read.**
The same summaries carry scaffold/contig counts, which partially answers the
still-open "draft vs. long-read assembly mix" item: 279/295 strains have >100
contigs (consistent with short-read/draft assemblies, median 534 contigs
across the panel), only 13 have <=20 contigs (consistent with long-read/
near-chromosome-level assembly), 3 in between. This is a contig-count proxy,
not a confirmed sequencing-technology determination (no read-type metadata
cross-referenced) -- real signal, but the open item isn't fully closed by
this alone.

## Leiden module rerun on the corrected trans network (2026-09-15)

Rerun of `bin/detect_trans_modules.py` against `pair_classification.rescued.tsv`
(the fully-corrected file, both rescue bugs fixed) -- the previous run
(commit `d021c58`) used the pre-rescue 522,933-edge network and is now stale;
its output was moved to `results/full_293run/modules_preview_prerescue_stale/`
for the record, not for use.

**New network is ~5x larger**: 2,613,303 `trans` edges (up from 522,933),
8,902 nodes (families with at least one `trans`-classified edge, up from
7,064). Same resolution sweep (0.5/1.0/2.0/5.0/10.0), same seed=0, all
finishing in 20-27s each (igraph+leidenalg scales fine at this size,
comfortably interactive -- no SLURM job needed):

| Resolution | Modules | Largest module | Singleton modules | Singleton % |
|---|---:|---:|---:|---:|
| 0.5 | 6 | 3,288 | 0 | 0% |
| 1.0 | 10 | 2,105 | 0 | 0% |
| 2.0 | 27 | 1,595 | 5 | 18.5% |
| 5.0 | 724 | 1,405 | 569 | 78.6% |
| 10.0 | 1,537 | 1,140 | 1,359 | 88.4% |

Same qualitative pattern as the pre-rescue run (module count and singleton
fraction both grow steeply with resolution; the pre-rescue sweep found
5->806 modules and 0%->74% singletons) -- now confirmed on a network with
5x the edges and real, corrected classifications, not just the earlier
smaller/wrong one. Strengthens rather than changes the earlier
interpretation: this looks like a genuinely persistent, densely-
interconnected core trans-co-occurrence structure across resolutions,
not an artifact of the smaller stale network. The same open question
remains unresolved: whether that persistent core reflects real biology
(e.g. a shared regulatory/selective-pressure network) or a residual
population-structure confound not fully caught by the clade-stratified
permutation null -- still no obviously "right" resolution value to pick,
and no published pangenome-specific precedent to set one from. Not
resolved this session; flagged for whoever writes the final interpretation.

## Presence/absence matrix outlier: Asfu_H1106 (2026-09-15)

The regenerated `presence_absence_matrix.png` shows one visible vertical
streak of extra presence calls, distinct from the general core/shell/cloud/
singleton banding. Investigated before deciding whether to exclude the
strain from the profile.

**Identity, for the record**: `Asfu_H1106` (this study's Short name) =
strain `H-1-10-6` (species.csv `Strain` column), assembly accession
`GCA_020501995.1` (both `Genome_Accession` and `Protein_Accession` --
NCBI-sourced for both), NCBI Taxon ID `746128`, `TaxonGroup` =
`Barber_cluster5`.

**Finding**: `Asfu_H1106` has 6,197 shell+cloud+singleton family presences
vs. a panel median of 3,591 (stdev 856 -- roughly 3 standard deviations
above median), the single most extreme strain by this measure. Breaking
that down by call type: 2,195 of those are protein-level `PRESENT` calls
(vs. a panel median of 640 -- ~3.4x median) while only 4,002 are
rescue-pass `GENOME_ONLY` calls (vs. a panel median of 2,926 -- ~1.4x
median, much less extreme). **The elevation is concentrated in the
original protein-level annotation, not the rescue pass** -- ruling out the
tblastn rescue redesign as the cause.

**Checked and ruled out as explanations:**
- BUSCO completeness: 98.7% (Complete) -- normal, not a low-quality
  assembly.
- Contig count: 678 -- above the panel median (534) but not in the
  extreme/fragmented tail (279/295 strains already have >100 contigs).
- Total annotated protein count: 10,961 -- within the normal range for
  this dataset (compared strains range ~8,800-11,100), not an
  over-permissive/over-fragmented gene-calling artifact by this measure.
- Annotation source: `ncbi` (species.csv `Protein_Source`), the same
  source used by 292/295 strains (only 3 use `uniprot`) -- not an
  annotation-pipeline-difference explanation either.

**Decision: kept in the panel, not excluded.** No positive evidence of a
technical artifact was found despite checking the obvious candidates
(assembly completeness, fragmentation, gene-calling volume, annotation
source) -- excluding a real strain on an unexplained-but-not-clearly-wrong
pattern would be an unjustified data-driven exclusion. Flagged as an open
QC item for further investigation if it recurs or matters to a specific
downstream claim (e.g. check allele-level divergence/clade placement for
this strain specifically, or whether it carries an unusually high real
count of divergent gene copies that split into shell/cloud families at the
tier-1 clustering identity threshold rather than staying in their true
core family).

## Dereplication threshold sensitivity check (2026-09-15)

Ran `bin/dereplicate_strains.py` at 4 Mash thresholds (0.0005, 0.001
[default, already used for the real 123/295-representative run], 0.002,
0.005) to check how sensitive the representative-strain count is to this
choice -- previously only the default had been tried.

| Mash threshold | Representative strains (of 295) |
|---|---:|
| 0.0005 | 207 |
| 0.001 (default, used in the real run) | 123 |
| 0.002 | 14 |
| 0.005 | 3 |

**Real finding: extreme sensitivity, and the default sits in the steep
part of the curve, not a stable plateau.** Halving the threshold
(0.001->0.0005) nearly doubles the representative count (123->207);
doubling it (0.001->0.002) collapses it by an order of magnitude
(123->14). This is not "the default is wrong" -- 123/295 (41.7%) is a
plausible representative fraction for a species with this much clonal/
outbreak-cluster structure -- but it does mean the dereplicated-strain
count (and therefore `frequency_bins.py --inventory`'s dereplicated
frequency counts) should be reported as threshold-dependent, not as a
fixed, robust number, and any downstream claim that specifically leans on
"123 representative strains" should note this sensitivity rather than
treat 123 as self-evidently the right cutoff. Not resolved further this
session (would need an independent criterion -- e.g. a known outbreak/
clonal-cluster ground truth from the reference paper's own population
structure -- to pick a threshold on grounds other than "the value used so
far").

## ID crosswalk + benchmark scorecard -- real results (2026-09-15)

Closes out the two remaining component-6 open items above, with real data
end to end -- no synthetic/placeholder inputs. New scripts:
`bin/fetch_paper_reference_proteins.py`, `bin/build_ground_truth_tables.py`;
`bin/id_crosswalk.py` and `bin/benchmark_scorecard.py` were already written
and unit-tested (the latter's `main()` was a deliberate stub -- now wired
up for real, see below). All outputs in `results/id_crosswalk/`
(gitignored, regenerable).

### What was actually crosswalkable, and what was not

Inspected the real columns of Tables S6/S7/S12/S13/S19/S21 before assuming
the general design spec's description matched (it mostly did, with one
correction): S12/S13's `geneID` column is **not** a public accession for
most of the paper's 13 "reference-quality" strains -- it's the paper's own
internal per-strain locus numbering (e.g. `47-10_000766`, `F7763_000019`),
tied to genome annotations that exist only in the paper's own (unpublished
alongside this supplement) analysis, not NCBI. Of the 552 Table S13
cargo-gene rows (17 named Starships), only the 40 rows for strain `AF293`
use a real NCBI accession (spelled `AF293_XP-<digits>.<version>` in the
table, one dash-for-underscore edit away from the real `XP_<digits>.
<version>` RefSeq accession). Table S19 (the virulence/SM-cluster gene
catalog) is much more usable: 647 of 802 rows carry a real `XP_*`
accession (798 loci are `Afu*g*` AF293 gene IDs; the other 4 are `AFUB_*`
A1163 ones). Table S5's `featureID` column has some `AF293_XP-*`-style
rows too. **616 unique real accessions total** were pulled from these
three tables' AF293-attributable rows; 2 of them turned out to be literal
placeholder text (`"AspGD only"` / `"Only available on AspGD"`) in Table
S19's `accession` column, not real accessions -- a data artifact in the
paper's own table, not a fetch failure.

**Consequence for benchmark coverage, stated plainly rather than papered
over**: the crosswalk (and everything built on it) is effectively an
**AF293-anchored crosswalk**. Of the 20 high-confidence Starships (Table
S6/S21), only 2 (`Nebuchadnezzar-h1`, `Nebuchadnezzar-h2`, both carrying
the previously-characterized hrmA/HAC paralog) have any crosswalked cargo
gene at all, because those are the only ones with an AF293 instance in
Table S13. The other 15 named Starships' cargo lists are 100% non-AF293
paper-internal locus tags with no path to a real sequence from this
supplement alone -- genuinely not resolvable without either (a) the
paper's own raw genome annotations (not published in this supplement) or
(b) running `starfish` (see below) directly against this study's own
assemblies for those strains, which sidesteps the paper's IDs entirely by
detecting Starships independently rather than crosswalking to the paper's
naming.

### starfish (considered, not run)

The user flagged mid-task that `starfish` (Gluck-Thaler's own Starship
caller, github.com/egluckthaler/starfish) is not installed here and could
independently call Starships on this study's own genomes rather than
relying on the paper's supplement tables. Not used this session: the
AF293-anchored crosswalk above already gave a real, resolvable,
sequence-verified positive control (Nebuchadnezzar-h1/h2, scored below),
and installing + validating a new tool + running it on even a handful of
genomes is a materially bigger addition than what was needed to produce a
real scorecard result. It would be the right next step specifically to
extend cargo-grouping/presence-recovery coverage to the other 15 named
Starships (whose paper-internal gene IDs are otherwise a dead end) -- left
as an explicit, real, not-yet-attempted option, not silently dropped.

### Fetch, crosswalk, and clustering — what ran

- **NCBI efetch**: this environment DOES have outbound internet access
  (verified with a live test fetch before committing to this approach) --
  614/616 accessions fetched cleanly in 4 batches
  (`results/id_crosswalk/paper_reference_proteins.fa`), provenance recorded
  in `DATA_MANIFEST.yaml`.
- **Diamond blastp crosswalk**: built a diamond db from the same
  Short-prefixed `results/full_293run/all_strains.fa` (already existed,
  2,788,402 proteins, 295 strains -- no rebuild needed) and blasted the 614
  paper reference sequences against it (`--more-sensitive
  --max-target-seqs 25`, ~9 min wall-clock on 4 threads). **614/614 queries
  got at least one hit** (`results/id_crosswalk/paper_vs_study.tsv`, 14,905
  hit lines); `id_crosswalk.py`'s existing (already-tested)
  `parse_diamond_blastp_besthits` reduced this to one best-hit crosswalk
  row per paper accession (`results/id_crosswalk/crosswalk.tsv`). No bug
  found in `id_crosswalk.py` -- used as-is, per the task's own instruction
  not to add code there without a real reason. Sanity check: median hit
  identity across all 14,905 lines is 100% (mostly near-identical
  cross-strain orthologs, as expected for this species), and 199/614 best
  hits land on this study's own `Asfu_Af293` (the same strain the query
  sequences came from) -- the rest landing on other strains' equally-close
  paralogs/orthologs is plausible at >98%-identical intraspecific distance
  and not a red flag.
- **Diamond-backend tier-1 clustering, attempted and completed for real**
  (the task's "possibly multi-hour" open question): submitted as SLURM job
  28425696 (`run_diamond_tier1_cluster.sh`, `-p stajichlab`, 16c/64gb) against
  the same `all_strains.fa` used for the real mmseqs-tier1 run. **Finished
  in 3 minutes wall-clock** -- far faster than mmseqs's own ~14 minutes on
  the identical input, not slower as the "could be slower" framing in the
  task worried. Result: **66,108 tier-1 families** (vs. mmseqs's 47,983 at
  the same nominal 90%-identity/80%-coverage operating point) --
  `results/full_293run/tier1_diamond_cluster.tsv`. Diamond's clustering
  splits families ~38% more finely than mmseqs at the same nominal
  threshold; not investigated further here (a real, reportable backend
  difference, but explaining *why* -- different seed/extension heuristics,
  different effective coverage handling -- is future work, not needed for
  the scorecard itself). One minor code change made in support of this run:
  `bin/cluster_backend.py`'s `run_diamond_cluster`/CLI gained an optional
  `--threads` passthrough (previously silently used diamond's own thread
  auto-detection, which does not respect a SLURM cgroup's CPU allocation)
  -- covered by the existing `test_cluster_backend.py` (still passes
  unchanged; that test only covers `two_tier_families`, not this function).
- **Diamond presence matrix**: built from `tier1_diamond_cluster.tsv` via
  the existing `build_presence_matrix.py` (unchanged, ~1 min) ->
  `results/full_293run/presence_matrix.diamond.tsv`. **This is the RAW,
  UNRESCUED matrix** -- no genome-level tblastn rescue pass was run for the
  diamond backend (that is a genuinely separate, multi-hour, 295-strain
  SLURM array job per strain, not a rerunnable-in-minutes step like
  clustering turned out to be) -- stated explicitly here and in every
  comparison below so the mmseqs-vs-diamond scorecard numbers are not
  mistaken for an apples-to-apples clustering-only comparison; the mmseqs
  side benefits from the rescue pass's ~2x core-fraction correction and the
  diamond side does not.

### `benchmark_scorecard.py` -- wired up for real, no longer a stub

Rewrote `main()` (the stub deliberately exited non-zero) with real I/O:
`load_crosswalk`/`load_cargo_by_name`/`load_starships_by_short`/
`build_protein_to_family` (thin wrappers, mostly reusing
`lib/pangenome_matrix.py`'s existing `read_cluster_tsv`), and
`score_backend`/`score_negative_control`, which call the pre-existing,
already-unit-tested `score_presence_recovery`/`score_cargo_grouping`
functions unchanged. New CLI args beyond the original stub's placeholder
shape: `--negative_control_genes`, `--tier1_cluster_mmseqs`/
`--tier1_cluster_diamond` (needed to map a crosswalked protein ID to its
tier-1 family -- the matrix alone doesn't carry per-protein family
membership). `--tier1_cluster_diamond`/`--matrix_diamond` are optional: if
omitted, the scorecard runs mmseqs-only rather than refusing to run, per
the task's explicit "score mmseqs alone rather than block" guidance. 8 new
unit tests added (`tests/test_benchmark_scorecard.py`,
`tests/test_build_ground_truth_tables.py`); all 179 study tests pass
(`pixi run pytest studies/fungi/Afumigatus_pangenome/tests/ -q`), no
regressions.

Also found and fixed one real ID-format mismatch while wiring this up (not
a pre-existing bug in tested code, a bug in the same-session ground-truth
builder): Table S13's AF293 geneIDs are spelled `AF293_XP-<digits>.<v>`,
but the crosswalk is keyed on the plain `XP_<digits>.<v>` accession
actually fetched/blasted -- without normalizing, EVERY cargo-grouping
control silently scored 0/0 (looked like "no crosswalk," not a format
bug). Added `normalize_geneid_for_crosswalk` to
`build_ground_truth_tables.py` (converts only the AF293-prefixed form;
other strains' internal locus tags are left untouched since there is
nothing to normalize them to) plus 2 unit tests, then reran.

### Real scorecard numbers

**Positive control -- cargo grouping** (Table S13, nameID-grouped): of the
20 named Starships, only `Nebuchadnezzar-h1`/`Nebuchadnezzar-h2` have any
crosswalked cargo gene (20/20 genes crosswalked for each -- see coverage
caveat above). Identical for both backends: **purity 1.0, completeness
0.05**. This is NOT a clustering failure -- purity 1.0 means every member
of the predicted family that captured the best overlap is a true cargo
gene (no contamination), and completeness 0.05 (1 of 20) means the 20
PF11001-family paralogs correctly split across ~20 different tier-1
families rather than collapsing into one. That is exactly the intended
behavior of the two-tier 90%-identity clustering scheme (component 1's
whole reason for existing was to keep this study's own previously-found
PF11001 paralogs from over-merging) -- a real, if unintuitive-sounding,
validation rather than a negative result. The other 18 Starships:
`n_crosswalked_genes = 0`, correctly reported as "no crosswalk hit," not
silently scored as a failure.

**Positive control -- presence recovery** (Table S6/S21, re-keyed to this
study's own `Short` via a 261/295-strain normalized-name match against
Table S21's `isolateID`/`originalID` -- slightly more than the
2026-09-13 TaxonGroup fill's 254/293 because this includes the 2 OUT
strains and uses `config.csv`'s current, corrected content): only
`Nebuchadnezzar-h1`/`Nebuchadnezzar-h2` have a resolvable diagnostic
family (same coverage limit as above). Real numbers, both backends nearly
identical despite picking different diagnostic strains as their top
crosswalk hit (`Asfu_niveus`/`Asfu_W72310` for mmseqs/diamond on h1,
`Asfu_Af293`/`Asfu_ATCC_46645` on h2 -- expected, since many strains carry
near-identical copies and the "best hit" is a coin flip among near-ties):

| Starship | Backend | n strains scored | TP | FP | TN | FN | Accuracy | Jaccard |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Nebuchadnezzar-h1 | mmseqs | 261 | 3 | 1 | 256 | 1 | 0.992 | 0.600 |
| Nebuchadnezzar-h1 | diamond | 261 | 3 | 1 | 256 | 1 | 0.992 | 0.600 |
| Nebuchadnezzar-h2 | mmseqs | 261 | 1 | 0 | 260 | 0 | 1.000 | 1.000 |
| Nebuchadnezzar-h2 | diamond | 261 | 1 | 0 | 260 | 0 | 1.000 | 1.000 |

Small numerator (this Starship is genuinely rare in the overlapping
population, matching the earlier hrmA reference screen's own finding of
only 8/293 strains hrmA-positive) but a real, correctly-recovered signal
either way -- both backends agree exactly on this control, giving no
basis (from this one control) to prefer one over the other.

**Negative control -- conserved SM/virulence genes** (Table S19, 645
genes with a real accession, expected to be broadly core/soft-core rather
than showing Starship-style presence/absence variability): **100%
crosswalk coverage** (645/645, since these accessions came from the same
fetch/blast pool).

| Backend | Matrix | core | soft_core | shell | cloud | singleton | mean freq |
|---|---|---:|---:|---:|---:|---:|---:|
| mmseqs | rescued | 554 (85.9%) | 9 (1.4%) | 54 (8.4%) | 25 (3.9%) | 3 (0.5%) | 0.916 |
| diamond | raw/unrescued | 513 (79.5%) | 14 (2.2%) | 57 (8.8%) | 48 (7.4%) | 13 (2.0%) | 0.867 |

Both backends correctly call the large majority of conserved SM/virulence
genes core/soft-core (as expected -- this is the sanity check the negative
control exists for), and both do markedly better on this gene set than the
whole-genome baseline (mmseqs's own genome-wide rescued composition is
53.0% core -- these 645 conserved genes land core at 85.9%, a real,
expected enrichment). mmseqs's higher core fraction and lower cloud/
singleton tail here is very plausibly the tblastn rescue pass at work
(exactly what it was built to fix -- coverage-based false absences turning
a truly-core gene into an apparent cloud/singleton), not necessarily a
mmseqs-vs-diamond clustering-quality difference -- **this comparison is
confounded by the rescue-pass gap, not a clean backend comparison**, as
flagged above.

### What's still open after this session

- **S12/S7/S14-16 as additional positive controls**: not attempted --
  S12 (13-reference-strain BLAST recovery) has the same paper-internal-ID
  problem as S13 for 12 of its 13 strains; S7 (manual 3-strain annotation)
  and S14-16 (segregating insertion regions) were not inspected this
  session at all. Real, available next step if broader benchmark coverage
  is wanted.
- **Random-family-pair negative control**: not attempted (see Open items
  checklist above) -- `cooccurring_pairs.rescued.tsv` is 1.3GB and a fair
  matched-frequency sampling needs more care than a quick pass.
- **mmseqs-vs-diamond comparison is not yet apples-to-apples**: would need
  a full tblastn rescue pass run against the diamond clustering (a
  genuinely separate multi-hour SLURM array job, not attempted this
  session) before the negative-control frequency-bin numbers above are a
  clean backend comparison rather than a backend-plus-rescue-status
  comparison.
- **starfish**: not installed/run (see above) -- would be the real path to
  extending positive-control coverage to the 15 non-AF293 named
  Starships.

## Accessory islands: real physical clustering beyond Starships (2026-09-15)

Answers a direct question: are the 44,576 `unexplained_physical` +
8,890 `starship_explained` + 29,949 `ambiguous_linkage` co-occurring pairs
(83,415 total, all statistically robust FDR<0.05 physically-linked calls)
part of REAL, larger multi-gene genomic islands, or just isolated
two-gene adjacencies? Ran `bin/synteny_windows.py`'s `accessory_islands()`
(built 2026-09-13, tested, never previously run against real data) against
all 295 strains, using `family_positions.rescued.tsv`'s already-computed
per-strain gene order directly (no GFF3 re-parsing needed) via a new
script, `bin/find_accessory_islands.py`.

**Real performance finding, fixed before trusting any output**: the first
implementation looped over all 83,415 significant pairs for every island
found (`accessory_islands()` returns thousands of non-core runs across
295 strains) -- an O(islands x pairs) blowup, killed after 12+ minutes
with zero output. Fixed with a `family -> [(pair, classification), ...]`
index (`build_pair_index`) so each island only checks pairs touching its
own member families. Fable-model review confirmed the fix is correct
(verified no misses/double-counts against the naive semantics) and
empirically safe at this dataset's actual scale (max hub-family degree
103, median 23, out of 7,645 families touching any significant pair --
no pathological hub case here). Real run after the fix: **2 minutes**
(down from 12+ minutes and not finished).

**Real results**: 70,665 significant islands found across all 295
strains, deduplicated by (member-family-set, classification-set) to
**12,861 distinct islands**. Island size ranges 2-837 genes (median 7,
p90=56, p99=221) -- the largest ones are very likely extended
subtelomeric/repeat-rich accessory regions, not compact biosynthetic gene
clusters, and should not be over-interpreted as single functional units
without further checking; realistic candidate secondary-metabolite gene
clusters are more plausibly found restricting to the 3-30-gene range
(typical fungal BGC size).

Cross-referenced each island's member families against the DUF3435
(Starship captain) hmmsearch results already on hand, plus a NEW
hmmsearch of two secondary-metabolite backbone Pfam domains (PKS
ketosynthase, `PF00109`/`ketoacyl-synt`; NRPS condensation domain,
`PF00668`/`Condensation`) run against the same all-strains proteome
(`results/sm_backbone/SM_backbone_vs_study.tblout`, 11,342 raw hits,
92s runtime):

| has_captain_gene | has_sm_backbone_gene | Distinct islands |
|---|---|---:|
| N | N | 11,741 (91.3%) |
| Y | N | 904 (7.0%) |
| N | Y | 151 (1.2%) |
| Y | Y | 65 (0.5%) |

**The large majority (91.3%) of significant physical-linkage islands have
NEITHER a Starship captain gene NOR a PKS/NRPS backbone gene anywhere in
the cohort** -- real, adjacency-confirmed multi-gene co-loss/co-gain
blocks with no mechanism identified by either check this study has run.
This is a genuinely open target for follow-up (candidate: a different
mobile-element family, or a distinct secondary-metabolite backbone class
not covered by the two Pfam domains checked here).

**A concrete, checkable example** (restricting to the 3-30-gene range):
a 30-family island in strain `Asfu_G2141`, supported by 246 significant
pairs (`ambiguous_linkage`+`unexplained_physical`), contains the reviewed
Swiss-Prot entry `Asfu_Af293|sp|Q4WKX2|FGND_ASPFU`, has a PKS/NRPS
backbone hit (`has_sm_backbone_gene=Y`) but no captain gene -- a real
candidate secondary-metabolite-associated accessory island independent
of Starship mobilization, worth targeted follow-up (full gene list:
`results/accessory_islands/significant_islands.tsv`).

**Caveats (both from the Fable review, real not hypothetical)**:
1. `has_captain_gene`/`has_sm_backbone_gene` are cohort-level flags (any
   copy of that family, in any of the 295 strains, anywhere) -- a "Y"
   does NOT mean the specific copy inside that specific island in that
   specific strain carries the domain. Treat as "this family is known to
   carry this domain somewhere in the population," not a per-instance
   confirmation.
2. The dedup-by-member-set step can report a nested/partial island (one
   strain's larger island containing another strain's smaller one as a
   proper subset) as two separate "distinct islands" -- read the
   12,861 count as "distinct exact member sets," not a claim that there
   are 12,861 non-overlapping biological loci.

## Functional (Pfam domain) enrichment of significant islands (2026-09-15)

Direct follow-up: what functions are enriched among the significant-island
member families, vs. the correct background (all 11,052 families actually
eligible for co-occurrence testing -- shell+cloud bins, `cooccurrence.py`'s
own selection -- never the whole genome, since core genes were never
eligible for testing and would spuriously inflate "accessory-typical"
domain enrichment regardless of which island is tested)?

**New script, real run**: `bin/summarize_island_functions.py`. Ran a real
Pfam-A hmmscan (30,134 profiles, pressed database) against the 9,224
unique island-member family representative sequences
(`run_island_pfam_scan.sh`, job 28438315, 1h44m/16c) plus the 2,620
additional eligible-background families not already covered
(`run_background_pfam_scan.sh`, job 28438374, 25min/8c) -- together, full
Pfam coverage of all 11,052 eligible families, not just the 9,224 subset
in islands. 5,944 of 11,052 families (53.8%) have >=1 Pfam domain hit
(E<=1e-3 at both the sequence and domain level).

**Real bug caught by a test before running on real data**: island-member
families include singleton-bin families (`accessory_islands()` merges any
non-core run: shell, cloud, OR singleton), but the eligible background is
shell+cloud only -- a singleton family was never eligible for
co-occurrence testing and can't sensibly count as "in" or "out" of a
population it was never part of. A test constructing exactly this case
(`test_domain_enrichment_warns_when_island_members_not_subset_of_
background`) caught scipy crashing on a resulting negative 2x2-table cell
before the real run. Fixed by intersecting island-member families with the
eligible background before testing (792 excluded on the real data, 8,432
remain) -- documented as an expected, reported occurrence, not silently
dropped.

**Real result: 901 domains tested, 71 significant at FDR<0.05** (one-sided
Fisher's exact, BH-corrected). Top hits, all biologically plausible, not
noise:

| Domain | In islands / background | Notes |
|---|---:|---|
| DUF3435 | 199/205 | The Starship captain-gene domain itself -- expected, an internal-consistency check (this domain directly defines `starship_explained`) |
| DUF3723 | 109/109 | Increasingly documented in the literature as part of the broader Starship "backbone" beyond the captain gene alone -- suggests some of the "unexplained" islands below may still be Starship-associated via a signal this study's narrower DUF3435-only screen didn't check |
| Ank / Ank_2-5 (ankyrin repeat) | up to 215/231 | Documented Starship-cargo-associated domain family in other fungi |
| DDE_1 (DDE transposase) | 68/68 | **Independent evidence of a different mobile-element family** driving some physically-linked co-occurrence -- direct support for the "different TE, not Starships" hypothesis raised when the 91.3%-unexplained finding was first reported |
| NACHT | 70/73 | Fungal heterokaryon-incompatibility / innate-immune-like gene family -- classic accessory-genome functional category |
| Glyco_hydro_71, Patatin | ~55 each | Secreted, cell-wall-modifying/lipase enzymes -- another classic horizontally-variable fungal accessory-genome category |

Full table: `results/accessory_islands/domain_enrichment.tsv`. Per-island
domain annotations (for browsing, sorted by island size descending):
`results/accessory_islands/significant_islands.with_domains.tsv`.

**A real methodological nuance found while spot-checking, not a bug**: an
island's own annotated Pfam domains and the specific reason one of its
member pairs was classified `starship_explained` by `pair_classification
.py` don't always overlap -- e.g. a real 15-family island
(`Asfu_C169L1`, domains `Beta-prop_ATRN-LZTR1,Beta-prop_FBX42,DDE_1,
Glyco_transf_90,Kelch_HCF,Kelch_KLHDC2_KLHL20_DRC7`) contains a
`starship_explained` pair despite none of its own annotated members being
DUF3435. This is because `pair_classification.py`'s captain-gene check
uses a fixed rank-window AROUND THE PAIR (k=10 genes), which is not the
same region as `accessory_islands()`'s merged non-core run -- the captain
gene can sit just outside the island boundary while still being "nearby"
by the pair-level window definition. Not a contradiction, but a reminder
that "island" and "pair-adjacency window" are two different, only
partially overlapping notions of physical proximity in this analysis.
