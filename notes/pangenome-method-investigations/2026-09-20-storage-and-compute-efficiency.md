# NovInvenio pangenome path: storage and compute efficiency review

Date: 2026-09-20. Reviewer: read-only pass over
`/bigdata/stajichlab/jstajich/projects/NovInvenio` (HEAD `8344725`) and the real
`coccidioides_pangenome/results/mmseqs_genus_vs_ureesii` outputs (529 ingroup + 1
outgroup strain, 54,421 families). Every number below was measured on that run or on
this login node unless it says "estimate". Scratch files are in
`/rhome/jstajich/.claude/jobs/3354f7d2/tmp/`.

## 0. The biggest surprise first: the rescue pass did nothing, and nobody noticed

`presence_matrix.tsv` and `presence_matrix.rescued.tsv` are byte-identical
(`cmp` clean, 207,934,530 B each). The matrix holds 4,509,470 `present` cells,
24,333,660 `absent`, and **zero `genome_only`**. `rescue_positions.tsv` is a header
line (26 B).

Cause, from the task log
(`.nf_launch/genus_vs_ureesii/work/16/145a917f.../.command.log`): 530 lines of
`Warning : <Short>.tblastn.tsv.zst is a symbolic link, ignoring`, then
`Rescue: 0 ABSENT->GENOME_ONLY, 0 skipped`. Nextflow stages inputs as symlinks;
`zstd -dc` without `-f` refuses a symlink and exits 1 (verified here with zstd 1.5.6).
The `compressed_io.py` of 2026-09-17 ignored that exit code. `RESCUE_PASS` only errors
when `hits and skipped == len(hits)`; an empty `hits` set passes silently.
`EXTRACT_RESCUE_POSITIONS` hit the same wall (its log ends "0 GENOME_ONLY positions
written") after peaking at 21.6 GB RSS.

Status at HEAD: fixed. `55069ab` (2026-09-20) makes the reader raise on a non-zero zstd
exit, and `e086c78` resolves the symlink first; `tests/test_pangenome_compression.py`
covers both. But the flagship outputs on disk predate the fix.

What the rescue would have done. I decompressed all 530 tblastn files
(131,145,913 hit rows) and counted distinct `(family, strain)` pairs passing the
pipeline's own thresholds (`pident >= 90`, `qcovs >= 80`):
**12,877,646 pairs qualify = 53% of the 24.3 M `absent` cells.** Per strain the count
ranges 2,909 (Uree, the outgroup) to 25,678 (UTAH_9443); median 24,507. With rescue
applied the matrix density goes from 15.6% present to about 60%. Every downstream
table (frequency bins, 14,896 eligible shell/cloud families, 13.07 M co-occurring
pairs, islands, Leiden modules, `report.md`) was computed on the unrescued matrix.

Cost of the dead computation: `TBLASTN_PER_STRAIN` = 530 tasks, 18,097 task-minutes
(301.6 task-hours, `med_cpu` label) plus 7.5 GB retained under `rescue/`
(4.8 GB BLAST DBs in `storeDir`, 2.7 GB tblastn `.zst`).

Two things follow, neither of them optional:

1. `bin/pangenome_rescue_pass.py` needs a guard: if at least one `--tblastn_tsv` input
   has non-zero size and zero hits parsed, exit non-zero. Same for
   `pangenome_extract_rescue_positions.py`. The HEAD reader fix catches this specific
   failure, but a parse-level regression (column shift, `qcovs` missing) would again be
   silent.
2. The 53% figure is a scientific question for the author, not for me. 24,080 of the
   54,421 families (44%) are single-strain; a tblastn hit at >= 90% identity and
   >= 80% query coverage from a 90%-identity cluster representative will fold most of
   them into near-core. Whether that is the intended behaviour of "rescue" or a sign
   that the thresholds need to be at or above the clustering identity is a design
   decision I do not have data to make. The pipeline should be rerun with rescue working
   before the current numbers are used further.

## 1. On-disk encoding

### Measured on `presence_matrix.rescued.tsv` (54,421 x 530, 15.6% present)

| encoding | bytes | zstd -3 | zstd -19 | gzip -6 |
|---|---|---|---|---|
| current TSV, literal words | 207,934,530 | 2,788,882 (0.15 s) | 1,461,721 (12.8 s, 548 MB RSS) | 1,968,767 |
| TSV, one `0`/`1` char per cell | 30,420,701 | 1,309,629 | 912,523 | 1,224,719 |
| packed bits, 67 B per row, no names | 3,646,207 | 652,100 | 538,395 | 648,851 |
| information content (1 bit/cell) | 3,605,391 | | | |

Decompress of the zstd -3 TSV: 0.36 s, 6 MB RSS.

Reading: the 58x-vs-information figure is real but it is a *plain-text* problem, not an
*encoding* problem. `zstd -3` on the file as it is today gets 74.6x. The best packed
format I can build gets a further 4.3x, which is **2.2 MB per run**. The run directory
is 13 GB. There is no storage argument for any format beyond `zstd` on the existing TSV.

So: **"zstd the existing TSV" is the right answer for storage, and the rest is
gold-plating.** PR #113's "smallest share of the win" judgement was correct in bytes
saved per run; it was wrong in bytes saved per second of work (0.15 s to compress, one
rename in two modules). Revisit it.

What breaks when the matrix becomes `presence_matrix.tsv.zst`:

- `modules/pangenome/presence_matrix.nf`, `modules/pangenome/rescue.nf`
  (`RESCUE_PASS` output name), and every process that names the file in its `input:`.
- `bin/pangenome_report_tables.py::per_strain_summary` (line 187) opens the matrix with
  plain `open()`, not `open_maybe_compressed`. It is the only consumer that would break;
  the other seven go through `PresenceMatrix.from_tsv`, which already handles `.zst`.
- `PresenceMatrix.to_tsv` already writes `.zst` when the path says so
  (`_ZstdWriter`), and `copy_number_path()` already keeps the suffix last. No lib change.

Level: use `-3` (default). `-19` costs 12.8 s and 548 MB RSS for 1.3 MB. The
`_ZstdWriter` is single-threaded `zstd -q -f`; at 208 MB input that is 0.15 s. Fine.

### Copy-number sidecar

`presence_matrix.copy_number.tsv` is written only when a cell has copies (the builder
always passes `copies`, so it is written every run; it is not in the published listing
here because `PRESENCE_MATRIX` marks it `optional: true` and the run predates it, or it
was not published). It is a sparse `(family, strain, copies)` list, one row per present
cell: 4.5 M rows here. It is only consumed by `read_copy_number_sidecar`, which no
`bin/` script calls except through `from_tsv(read_copy_number=True)`, and no consumer
then reads `.copy_number` (grep: only `tests/test_pangenome_matrix.py:94-96`). Two
options: (a) leave it, compressed alongside the matrix; (b) drop it from the default
`from_tsv` load (`read_copy_number=False` default) so nobody pays 4.5 M dict entries for
a field nothing reads. (b) is a one-line change; do it with the numpy refactor below.

### Columnar / bitset / Parquet / HDF5 / npz

Not here. Reasons, concretely:

- The consumers want three access patterns: whole-matrix statistics (`FREQUENCY_BINS`,
  `COOCCURRENCE`, `REPORT_RENDER`), per-strain column sums (`per_strain_summary`), and a
  small family subset (`ISLAND_SYNTENY`, 298 of 54,421). All three are served by one
  sequential scan of a 2.8 MB file in about 1 s (measured below). A random-access
  container solves a problem that does not exist at 1 s.
- Parquet/Arrow adds `pyarrow` (about 100 MB in the pixi env, not currently a dep);
  HDF5 adds `h5py`. Neither can be inspected with `zstd -dc | cut -f1,5`, which is how
  this group debugs a matrix. `npz` is free (numpy is already there via pandas/scipy)
  but is opaque to non-Python tools and freezes the row/column order into two side
  arrays. A `.npz` of the bool array would be 3.6 MB uncompressed / 0.65 MB zstd, the
  same as the packed row above, for a format change touching all eight consumers.
- The sparse alternative already exists. `family_positions.tsv` (4,512,101 rows) is one
  row per present `(strain, family, copy)`, so it *is* the COO form of the matrix plus a
  rank. Do not create a second sparse form. Also, once rescue works, density goes to
  about 60% and a sparse list loses to a dense bitset anyway.

If a binary form is ever wanted, the packed-bits row (67 B per family for 530 strains,
`int(bits, 2).to_bytes()`) plus a family-name column is 3.6 MB, needs no dependency,
and is exactly what `pangenome_cooccurrence.py::presence_bitset` already builds in
memory. That is the only "bitset format" I would consider, and only as an in-memory
representation, not a file.

### Other large intermediates, measured zstd ratios (200 MB head samples)

| file | plain | zstd -3 ratio | zstd -19 ratio | note |
|---|---|---|---|---|
| `pair_classification.tsv` | 2.13 GB | 11.5x (about 185 MB) | 18.6x | #113 compresses it |
| `cooccurring_pairs.tsv` | 1.82 GB | similar | | **redundant, see 4.3** |
| `family_positions.tsv` | 260 MB | 6.8x | 11.6x | #113 |
| `gene_positions.tsv` | 232 MB | 4.3x | 7.6x | #113 |
| `cluster/tier1_cluster.tsv` | 265 MB | 10.1x | 16.5x | not in #113 |
| `rescue/db/*` (530 BLAST DBs) | 4.8 GB | n/a | | `storeDir`, see 4.5 |
| `rescue/per_strain_chunks/*.zst` | 2.7 GB | already zstd | | see 4.6 |

Run directory today: 13 GB. After #113 plus the four items in section 4 (drop DB
`storeDir`, drop the duplicate pairs file, compress matrix and cluster TSV): about
3.2 GB, of which 2.7 GB is raw tblastn output.

## 2. Read-side access

Measured loaders on the real 208 MB matrix, Python 3.12 from the pixi env:

| loader | wall | peak RSS |
|---|---|---|
| `PresenceMatrix.from_tsv` (dict keyed by `(family, strain)`, 28.8 M entries) | 57.7 s | 4.52 GB |
| same file into a `(54421, 530)` numpy bool array, line by line | 11.2 s | 0.13 GB |
| filtered stream, keep 298 named families' lines | 0.92 s | 13.5 MB |
| `zstd -dc` of the -3 file to `/dev/null` | 0.36 s | 6 MB |

The trace confirms the dict cost at every consumer: `PRESENCE_MATRIX` 2.2 GB,
`FREQUENCY_BINS` 4.3 GB, `RESCUE_PASS` 4.3 GB, `EXTRACT_ABSENT_QUERIES` 4.6 GB,
`EXTRACT_RESCUE_POSITIONS` 21.6 GB, `COOCCURRENCE` 16.2 GB (about 4.3 GB of it the
matrix). The island-synteny script's 5.54 GB -> 2.32 GB after filtering is the same
dict, and the remaining 2.3 GB is `load_positions` doing the same thing to 4.5 M rows.

Answer to the interface question:

- **Filtered streaming is correct and sufficient.** A sequential scan is 0.9 s
  uncompressed and would be about 1.3 s from `.zst`. An offset index would save under a
  second and add a file to keep in sync. Do not build one.
- **The thing to fix is the in-memory representation, not the file.** Keep the
  `PresenceMatrix` API (`is_present`, `call`, `presence_vector`, `strain_count`,
  `frequency`, `to_tsv`, `from_tsv`, `set_call`) and back it with a `uint8` array of
  shape `(F, S)` with values 0/1/2 for absent/present/genome_only, plus
  `{family: row}` and `{strain: col}` index dicts. All 12 external call sites use the
  methods; nothing outside `tests/test_pangenome_matrix.py:94-96` touches `.calls` or
  `.copy_number`. Add two things while there: `from_tsv(path, families=<set>)` that
  skips non-matching lines before splitting them (the 0.92 s path), and a
  `presence_array(strains)` accessor returning the bool sub-array so
  `FREQUENCY_BINS`, `REPORT_RENDER.accumulation_curve` and `COOCCURRENCE` stop calling
  `is_present` 28.8 M times in Python (about 30 s per pass; `REPORT_RENDER` does it
  twice).
- Expected effect: every matrix consumer drops from 2-5 GB to under 0.2 GB and from
  about 60 s to about 11 s load. `COOCCURRENCE` loses 4.3 GB of its 16.2 GB.

Same treatment for `load_positions` in `bin/pangenome_island_synteny.py` and
`load_family_positions` in `bin/pangenome_pair_classification.py`: filter by the needed
`(strain, family)` set while streaming (I did this to build the payload below: 500
islands' positions load in a few seconds, tens of MB).

## 3. The embedded web payload (`island_synteny.html`)

I rebuilt the real payload with `lib/island_synteny.build_payload` against the actual
islands table and matrix (`report_tables/islands_with_domains.tsv`: 27,666 located
islands, 10,602 with >= 2 strains).

| top_n | payload | `haplotypes` block | of which `pattern` strings | of which `strains` name lists | families+domains | haplotypes total | max / median per island |
|---|---|---|---|---|---|---|---|
| 50 | 0.55 MB | 0.48 MB | 0.072 MB | 0.29 MB | 0.04 MB | 3,520 | 225 / 53 |
| 500 | 4.49 MB | 4.06 MB | 0.350 MB | 2.88 MB | 0.25 MB | 25,020 | 225 / 44 |

Verdict on the `"0110"` strings: **a non-problem.** They are 13% of the payload at
top-50 and 8% at top-500. Bitpacking or base64 cannot save more than about 60 KB at
top-50 and would add a decoder to hand-written JS. The spec's bitpacking item should be
closed as not needed, with these numbers as the reason.

The real cost is the per-haplotype `strains` name arrays: 53% of the payload at
top-50, 64% at top-500. Every strain name is repeated once per island (530 names x
about 11.5 B x 50 islands). Two lossless replacements, both measured:

| encoding | top-50 | top-500 |
|---|---|---|
| current: `haplotypes[].strains` as name lists | 0.55 MB | 4.49 MB |
| `haplotypes[].strains` as integer indices into `payload.strains` | 0.36 MB | 2.64 MB |
| per island one `strain_hap` array (length = n_strains, value = haplotype index) and no per-haplotype strain list | **0.28 MB** | **2.03 MB** |

`strain_hap` is the right one: `payload.strains` already exists in the payload and is
already sorted, so `strain_hap[i]` maps `payload.strains[i]` to its row. The JS needs
one loop to rebuild `strains` per haplotype for the hover label at line 438 of
`lib/island_synteny_template.py`. `collapse_haplotypes` stays as is (tests keep passing);
`build_payload` adds the array and drops the list. This halves the file at any scale.

Scaling model (estimates; the haplotype count is bounded by the strain count and by
`2^size`, and in this data the median island has about 10% as many haplotypes as
strains):

- 2,000 strains, top-50, current encoding: names about 1.2 MB, patterns about
  0.2 MB, total about 1.5 MB. With `strain_hap`: about 0.6 MB.
- 2,000 strains, top-500, current encoding: about 14 MB. With `strain_hap`: about 5 MB.

Threshold: browsers parse a 5-10 MB inline JSON block in well under a second, and the
page is opened from `file://` so there is no transfer cost. The existing
`novelties.html` precedent is about 5 MB. The current encoding crosses 10 MB only at
top-500 with about 1,500 strains; with `strain_hap` it does not cross 10 MB below about
4,000 strains at top-500. Rendering is canvas (`fillRect` per cell; max 225 x 56 per
island), so there is no DOM cliff at top-500. The 500-button sidebar is the first thing
that would feel slow, not the data.

Encodings to avoid here: run-length or sparse index lists per pattern (the patterns are
17-56 chars; the overhead of the list syntax exceeds the saving), base64 bitpacking (a
decoder for 60 KB), and gzip-in-JS (the page has no dependencies and `file://` has no
`DecompressionStream` guarantee across browsers).

## 4. Computation

Trace for the run (`results/mmseqs_genus_vs_ureesii/trace.txt`, 2,137 tasks). Every
process not listed finished in under 5 min.

| process | wall | peak RSS | what it did |
|---|---|---|---|
| `TBLASTN_PER_STRAIN` x 530 | 301.6 task-hours | 0.44-0.58 GB each | produced 131 M hit rows; 0 applied (section 0) |
| `PAIR_CLASSIFICATION` | 1 h 33 m | 2.3 GB | 13,073,012 pairs |
| `COOCCURRENCE` | 1 h 23 m | 16.2 GB | 110.9 M pairs screened, 28.3 M exact Fisher, 13.07 M stratified tests |
| `CLUSTER_TIER1` | 3 m 41 s | 3.9 GB | mmseqs; fine |
| `GENE_POSITIONS` | 4 m 55 s | 44 MB | fine |
| `EXTRACT_RESCUE_POSITIONS` | 44 s | **21.6 GB** | wrote 26 bytes |

Per-call costs measured in the pixi Python 3.12 (n = 411 ingroup representatives):

| operation | per call |
|---|---|
| `scipy.stats.fisher_exact(alternative="greater")` | 271 us |
| `exact_stratified_pvalue` (4 clades, `np.convolve`) | 197 us |
| `quick_screen_pvalue` | 1.4 us |
| `fisher_counts_from_bits` | 1.25 us |
| `linkage_fraction`, 530 strains, both families in 30 | 174 us |
| `linkage_fraction`, 530 strains, both families in all | 1,070 us |
| `ast.literal_eval` of a `clade_composition` repr | 29 us |

### 4.1 `COOCCURRENCE`: 28.3 M `fisher_exact` calls that are one vectorised `hypergeom.sf`

Log: 14,896 eligible families, 110,937,960 candidate pairs, 28,321,259 clear the
analytic prefilter, 13,073,012 clear FDR. 28.3 M x 271 us = 2.1 CPU-hours in
`fisher_exact`, plus 13.07 M x 197 us = 43 min in the stratified test.

The one-sided Fisher p-value for a 2x2 table is exactly
`hypergeom.sf(both - 1, N, row1, col1)`. Verified against `fisher_exact` on 3,000
random tables: max relative difference 0.0. `hypergeom.sf` takes arrays. On this
loaded login node it ran at 0.2 M pairs/s single-threaded, so all 110.9 M pairs exact
in about 10 min worst case, with no prefilter and no `benjamini_hochberg_sparse`
padding trick needed (110.9 M float64 = 887 MB, or a two-pass streaming BH). The
`both` counts for all pairs come from one `A @ A.T` on a `float32` `(14896, 411)`
matrix (use float32, not int32: numpy integer matmul is not BLAS-backed; the int32
version took 161 s here, a float32 one is BLAS). Expected: 83 min -> about 10-15 min,
and the 16.2 GB drops by the 4.3 GB dict matrix plus the 28.3 M-tuple `survivors` list
(about 2.5 GB) plus the 13 M-dict `results` list (about 5 GB) if rows are written as
they are produced instead of collected. Compatibility: output file identical; `fdr_q`
values change only because the prefilter's "p = 1.0 for screened-out pairs" assumption
is replaced by the true p-values, which can only make BH slightly *more* conservative
where the screen was loose. Say so in the PR.

The stratified test: 13.07 M calls at 197 us is 43 min. Many pairs share the same
`(both_observed, (n_c, k_c, m_c) per clade)` key; the PMF is already cached, the
convolution is not. Cache the final p on the full key and measure the hit rate on this
run before deciding whether to vectorise. Deferring the test to after linkage
classification does not help: 12.7 M of the 13.07 M pairs are `trans*` and need it.

Unverified observation: the trace records `%cpu = 992` for `COOCCURRENCE`, whose loop is
single-threaded Python. I did not determine the cause. BLAS worker threads spinning
during the many tiny `np.convolve`/`hypergeom.pmf` calls is a common cause; a rerun
with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` in `beforeScript` would confirm or
rule it out, and costs nothing if wrong.

### 4.2 `PAIR_CLASSIFICATION`: O(pairs x strains) where O(pairs x co-carriers) is available

For each of 13.07 M pairs, `classify_pair` scans all 530 strains for `n_co_carrying`
and `linkage_fraction` scans them again (174-1,070 us), then `ast.literal_eval` parses a
Python repr (29 us; 6.3 min total). About 400 us x 13 M = 1.4 h, which matches the 1 h
33 m.

Fix: build `carriers[family] = set(strains with a position)` once from
`family_positions` (54,421 sets, cheap), then `co = carriers[a] & carriers[b]` and pass
only `co` to `linkage_fraction`. Intersection cost is proportional to the smaller set
(median family is in 1-5 strains here). Expected 10x. Write `clade_composition` as JSON
(or emit `n_clades` as its own integer column) so the consumer does not `literal_eval`
13 M strings. Output unchanged.

### 4.3 `cooccurring_pairs.tsv` is a strict subset of `pair_classification.tsv`

Same 13,073,012 rows; `pair_classification.tsv` carries every column of
`cooccurring_pairs.tsv` plus `classification` and `linkage_fraction`. Both are
published: 1.82 GB + 2.13 GB plain, about 350 MB after #113. Stop publishing
`cooccurring_pairs` (keep the channel; it is `PAIR_CLASSIFICATION`'s input). Saves
1.8 GB plain / about 160 MB zst per run for a one-line `publishDir` change.

### 4.4 92% of `pair_classification.tsv` is read three times and used by nobody

Classification counts: `trans_unconfirmed` 11,773,437 (90.1%), `trans` 966,865,
`insufficient_data` 281,947, `unexplained_physical` 34,124, `ambiguous_linkage` 16,639.
`BUILD_ISLANDS` keeps the physical classes (50,763 rows), `LEIDEN_MODULES` keeps
`trans` (966,865), `REPORT_TABLES` only counts. Each decompresses and `csv.DictReader`s
the whole 2.1 GB. Emit a second, filtered file (`pair_classification.significant.tsv`,
about 1 M rows, about 15 MB zst) from the same process and point the three consumers at
it; keep the full table as the archival record. This also makes the 12 M
`trans_unconfirmed` rows a candidate for not being written at all if the author decides
they are not worth 160 MB per run.

### 4.5 `MAKE_STRAIN_GENOME_DB` in `storeDir`: 4.8 GB kept forever for 1.5-2.2 s of work

530 `makeblastdb` tasks, each 1.5-2.2 s realtime, stored permanently under
`rescue/db/` (4.8 GB, 37% of the run directory) and submitted as 530 separate SLURM
jobs. Build the DB inside `TBLASTN_PER_STRAIN`'s task directory (or `$SCRATCH`) and do
not publish it. Saves 4.8 GB per study and 530 jobs; loses 18 min of `low_cpu` reuse on
a `-resume`. This also matches the author's own 1-1.5 h job-sizing rule: 530 jobs of 2 s
is the anti-pattern.

### 4.6 Rescue outputs are parsed twice and mostly discarded

`RESCUE_PASS` and `EXTRACT_RESCUE_POSITIONS` each read all 530 tblastn files (2.7 GB
zst, 131 M rows) and each keep only qualifying `(family, strain)` pairs; the second one
also loads the matrix again and holds 530 per-file best-hit dicts before merging
(21.6 GB peak). Options, in order of value/effort: (a) in
`EXTRACT_RESCUE_POSITIONS`, return early when the matrix has no `GENOME_ONLY` cell, and
merge per file instead of collecting 530 dicts; (b) fold position extraction into
`RESCUE_PASS` so the tblastn set is parsed once and both outputs are written from one
pass; (c) have `TBLASTN_PER_STRAIN` also emit a per-strain best-qualifying-hit table
(`awk` on the tblastn stream: family, contig, sstart, send, pident, qcovs, bitscore;
about 12.9 M rows total, tens of MB zst) so the two consumers never touch the raw
131 M rows. (c) changes what is published; keep the raw `.zst` as optional output for
re-thresholding.

### 4.7 The two search directions

Pairwise mode (`SEARCH` / `LOSS_SEARCH`): no duplication. `SEARCH` queries ingroup
proteomes against all others; `LOSS_SEARCH` queries outgroup proteomes against all
others. Reverse pairs are distinct searches with non-symmetric e-values, and `storeDir`
already caches every one. Do not try to derive one direction from the other.

Profile mode (`PROFILE_SEARCH` / `PROFILE_LOSS_SEARCH`): real duplication.
Each alias clusters its own seed group, runs famsa + hmmbuild per family, and
hmmsearches every proteome. Every family conserved across both groups gets two HMMs and
two full scans. A joint clustering of ingroup + outgroup, one HMM set, one scan, and
both directions derived from one matrix would roughly halve profile-mode compute on a
balanced study. Cost: it changes what a "family" is (seeded from both groups, not one),
which is an ADR-0002 semantic decision, not a refactor. Ticket, with the author.

### 4.8 Things I checked that are not waste

- `accumulation_curve`: 20 permutations x 530 strains of numpy `|=`/`&=` over 54 k
  families is under a second. The 30 s around it is the `is_present` loop (2.4), not the
  rarefaction.
- `CLUSTER_TIER1` 3 m 41 s, `GENE_POSITIONS` 4 m 55 s, Leiden on 966,865 edges: fine.
- `build_pair_index` in `BUILD_ISLANDS` already avoids the O(islands x pairs) scan.

## 5. Ranked lists

### (a) Do now

1. Rerun the flagship study with the HEAD `compressed_io.py`, after adding the
   zero-hits guard to `pangenome_rescue_pass.py` and `pangenome_extract_rescue_positions.py`
   (section 0). Everything downstream of the matrix is currently computed on
   silently-unrescued data.
2. `presence_matrix*.tsv` -> `.tsv.zst` (two module output names, one `open()` in
   `pangenome_report_tables.py:187`). 416 MB -> 5.6 MB per run, 0.15 s.
3. Drop `MAKE_STRAIN_GENOME_DB`'s `storeDir`; build the DB inside the tblastn task.
   4.8 GB and 530 jobs per study.
4. Stop publishing `cooccurring_pairs.tsv`. 1.8 GB plain per run, one line.
5. Island payload: `strain_hap` index array instead of per-haplotype name lists. Halves
   the page (0.55 -> 0.28 MB at top-50, 4.49 -> 2.03 MB at top-500). Close the
   bitpacking item as not needed with the section 3 numbers.

### (b) Ticket

6. `PresenceMatrix` on a `uint8` array behind the same API, plus a filtered
   `from_tsv(families=...)`. 4.5 GB / 58 s -> 0.13 GB / 11 s at every consumer;
   `EXTRACT_RESCUE_POSITIONS` 21.6 GB and `COOCCURRENCE` 16.2 GB both shrink.
7. `COOCCURRENCE`: vectorised `hypergeom.sf` for all pairs, stream results. 83 min ->
   about 10-15 min; memory under 3 GB. Then measure the stratified-test cache hit rate.
8. `PAIR_CLASSIFICATION`: carrier-set intersection; JSON `clade_composition`. About 10x.
9. Filtered `pair_classification.significant.tsv` for the three consumers (4.4).
10. Rescue parse-once and per-strain best-hit tables (4.6).
11. Profile-mode joint clustering across the two directions (4.7) - needs an ADR
    decision first.
12. `tier1_cluster.tsv` (265 MB, 10x compressible) into the #113 set.
13. Try `OPENBLAS_NUM_THREADS=1` on the Python stats steps and read the `%cpu` column.

### (c) Do not do

- Parquet, Arrow, HDF5, or npz for the presence matrix. Saves at most 2.2 MB per run
  over zstd-TSV, adds a dependency or an opaque file, touches eight consumers.
- Any bitset, RLE, base64 or sparse-index encoding of the island page's pattern strings.
  They are 8-13% of the payload; the names are 53-64%.
- A row-offset index or random-access container for the matrix. A filtered sequential
  scan is 0.9 s.
- Changing the three-state cell words to digits in the TSV. 2.8 MB -> 1.3 MB zst is not
  worth breaking the documented cell contract.
- A second sparse (family, strain) file. `family_positions.tsv` already is one, and a
  working rescue pushes density to about 60% where sparse loses.
- Deriving `LOSS_SEARCH` from `SEARCH`'s reciprocal hits in pairwise mode. E-values are
  not symmetric and `storeDir` already caches.
- `zstd -19` anywhere in the pipeline. 12.8 s and 548 MB RSS for 1.3 MB on the matrix;
  `-3` is the right default and `-T0` only matters for the multi-GB pair tables.
