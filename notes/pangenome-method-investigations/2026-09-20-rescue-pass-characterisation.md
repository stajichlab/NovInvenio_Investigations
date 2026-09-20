# What the RESCUE_PASS would actually rescue — Coccidioides 529-strain pangenome

Date: 2026-09-20. Read-only investigation; no repo files changed.
Study: `NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/results/mmseqs_genus_vs_ureesii/output/pangenome`

## Method

Re-implemented the shipped rescue predicate (`bin/pangenome_rescue_pass.py`:
`pident >= 90` on any HSP **and** `qcovs >= 80` for that query-subject pair) over
**all 530 per-strain TBLASTN files, all 130M HSP lines** (not a sample). For every
qualifying (family, strain) cell the best-bitscore subject's HSP span was taken and
intersected with that strain's predicted genes (`gene_positions.tsv`), each gene
labelled with its tier-1 family (`cluster/tier1_cluster.tsv`).
Scripts and raw per-cell output: `/rhome/jstajich/.claude/jobs/3354f7d2/tmp/rescue/`.

Sanity check of the headline figure: **12,877,646** cells pass — exact match to the
reported number. 24,333,660 absent cells → **52.92%**. Density 15.63% → 60.28%.
`genome_only` in the shipped matrix = 0, confirming the zero-cells bug.

## 1. Three-way overlap breakdown (all 12,877,646 cells, exhaustive)

| Class | Cells | % |
|---|---|---|
| Overlaps a predicted gene of a **different** family | 9,988,267 | **77.56%** |
| **Intergenic** (no predicted gene overlap) | 2,889,379 | **22.44%** |
| Overlaps a predicted gene of the **same** family | 0 | 0.00% |

Zero same-family hits — the presence matrix is internally consistent; no cell was
already-present-but-called-absent.

**77.6% of the rescue adds no new gene locus to the strain.** The DNA it finds is
already annotated, already in the matrix, under a different family row.

## 2. What the other_family hits really are (paralog vs clustering artifact)

Co-occurrence test over the 451,912 distinct (query family, overlapped family) pairs,
using the presence matrix: overlap coefficient = |both present| / min(|A|,|B|).

- **94.1%** of other_family cells come from family pairs that **essentially never
  co-occur** (overlap coeff < 0.01). A paralog pair co-occurs; a split of one locus
  into two family rows does not.
- Only **4.3% of other_family cells (429,579; 3.3% of all rescued)** have both families
  present in ≥100 strains — the genuine paralog/multi-copy candidates.
- The top recurring pairs are unmistakable: query family present in **1** strain,
  overlapped family present in **529**, co-presence **0**. e.g.
  `CA15|CEB9A78_001363-T1` (1 strain) → `4545-MICE|CB4077B_006100-T1` (529 strains), 529 cells.

Mechanism, measured: **query-rep length / overlapped-family-rep length is < 0.8 in
71.9%** of other_family cells, with a sharp cliff exactly at 0.8 (12.4% in 0.7–0.8,
4.8% in 0.8–0.9) — **mmseqs' own `-c 0.8` coverage cutoff**. Spot checks: 276 vs 387 aa,
263 vs 338 aa, 411 vs 618 aa. These are **truncated/fragmented gene models that tier-1
clustering refused to merge with their full-length orthologue**, becoming singleton
families that then "rescue" into all 529 genomes.

Query-family prevalence of rescued cells: **44.7% come from families present in exactly
one strain**; 2,925,354 cells are a singleton family hitting a **core** (500–530 strain)
family's gene.

## 3. Copy number (Q2)

Copy number is not the story. Mean copies/strain is 1.00 in every stratum
(heavily rescued ≥265 strains: 1.001; lightly rescued: 1.001; never rescued: 1.000).
Families multi-copy in ≥1 strain: heavily rescued 0.7%, never rescued 3.5% — heavily
rescued families are *less* multi-copy than background. Not a paralog-expansion signal.

## 4. Pfam (Q3)

Mild promiscuous-domain enrichment in the ≥80%-other_family families (n=23,274; 36.5%
Pfam-annotated) vs never-rescued (9.5%): Pkinase 408, APH 404, MFS_1 328,
PK_Tyr_Ser-Thr 267, DUF3435 155 (a *Starship* transposon captain domain), AAA_22/AAA_16.
Real but second-order — these counts are ~1-2% of the 23k families. Intergenic-dominated
families are barely annotated (18.6%), top terms in single digits, consistent with short
dubious ORFs rather than a TE class.

## 5. Alignment statistics and where the thresholds sit (Q4)

Pass-set distributions: pident 40% at exactly 100, median 99; qcov 60% at exactly 100,
median 100; 99.6% of cells have a single qualifying contig. This is a genus of ~99%-identical
genomes — identity carries almost no discriminating information.

**Tightening pident/qcov does not change what is rescued:**

| threshold (pident, qcov) | cells kept | % intergenic | % other_family |
|---|---|---|---|
| 90, 80 (shipped) | 12,877,646 (100%) | 22.4 | 77.6 |
| 95, 95 | 8,965,152 (69.6%) | 19.7 | 80.3 |
| 98, 95 | 7,159,089 (55.6%) | 18.8 | 81.2 |
| 99, 99 | 4,967,281 (38.6%) | 20.4 | 79.6 |
| 99.9, 99 | 3,652,243 (28.4%) | 22.4 | 77.6 |

The class mix is **flat (19–23% intergenic) across the entire grid**. There is no natural
cut point in identity or coverage — because the artifact hits are 100%-identical real DNA.

## 6. Are the intergenic hits genuine missed genes?

15 random strains re-scanned with coordinates (84,280 intergenic cells):

- **83.2% come from family reps < 150 aa** (60% in the 100–149 aa bin, 23% < 100 aa);
  median rep length across all families is 201 aa. The intergenic bucket is overwhelmingly
  short proteins.
- **47.5%** fall in a 2 kb window hit by **≥3 different families** in the same strain —
  repeat/low-complexity hotspots.
- 89.9% lie within 2 kb of a predicted gene; median 548 bp.
- Contig-edge proximity *does* track assembly quality (4.1% of hits <2 kb from a contig end
  in the best-assembled strain vs 33.5% in the worst) — a real but small fragmentation signal.

Conservative genuine-dropout filter (intergenic AND rep ≥150 aa AND not in a repeat hotspot):
**5.6% of intergenic = ~163,000 cells = 1.27% of the rescuable set = 0.67% of absent cells.**
A looser filter (rep ≥150 aa only) gives 484,900 cells = 3.8% / 2.0%.

## 7. Assembly quality (Q5) — does not support world (a)

Across all 530 strains: Spearman(n_contigs, rescue rate) = **−0.122** (p=0.005);
Spearman(n50, rescue rate) = **+0.297** (p=2.9e-12); Spearman(n_contigs, intergenic
fraction) = −0.030 (p=0.5). Best-assembled quartile rescue rate 0.527 vs worst 0.523.
Better-assembled genomes get rescued slightly *more*. The prior near-clone finding is not
contradicted (Spearman(n_contigs, n_absent) = −0.43), but it does not drive the rescue.

## 8. Concentration (Q6) — diffuse, not a block

41,060 of 54,421 families (75%) would gain ≥1 rescued cell; 30,439 families would be
rescued in ≥100 strains; top 1% of families account for only 1.7% of rescued cells,
top 10% for 16.8%. This is a long, pervasive tail — it would move essentially the whole
accessory genome into the core, not one suspicious block.

## Verdict

| Bucket | Cells | % of 12.88M |
|---|---|---|
| Redundant locus — already-annotated gene, different family row (clustering/gene-model fragmentation) | 9,558,688 | **74.2%** |
| True paralog / multi-copy hit | 429,579 | **3.3%** |
| Intergenic, short (<150 aa) or repeat-hotspot — unresolved, mostly dubious ORFs | 2,726,261 | **21.2%** |
| Intergenic, ≥150 aa, not a repeat hotspot — **plausible genuine annotation dropout** | ~163,000 | **1.3%** |

**Mostly spurious — but not in the way the question framed it.** It is not paralogy (3.3%).
It is redundancy: three quarters of the rescue re-asserts genes the matrix already contains
under another family, because tier-1 mmseqs clustering at `-c 0.8` split truncated gene
models off from their full-length orthologues. Applying the rescue as shipped would take
density 15.6% → 60.3% while adding **at most ~2% real new gene content**, and would make
24,080 singleton families look near-core.

## Recommendation

1. **Do not apply the rescue as configured.** It would corrupt the core/accessory split.
2. **Do not bother tightening `pangenome_rescue_min_pident` / `min_qcov`.** Measured: the
   artifact fraction is invariant from (90,80) to (99.9,99). Raising them only discards
   cells uniformly.
3. **Add the decisive structural criterion**: reject a rescue hit whose genomic span
   overlaps a predicted gene already assigned to a *different* family in that strain.
   Removes 77.6% of cells, is cheap (`gene_positions.tsv` is already produced), and is
   exactly the "is this a new locus?" question the rescue is supposed to answer.
4. **Add two cheap guards on what remains**: minimum query-rep length (≥150 aa; removes 83%
   of intergenic noise) and rejection of hits in 2 kb windows hit by ≥3 distinct families
   (repeat hotspots). Result: ~163k rescued cells, density 15.63% → 16.2%.
5. **Fix the real problem upstream.** 44.7% of rescuable cells come from singleton families
   whose rep is a truncated copy of a core gene. The rescue pass is a symptom detector for a
   clustering-parameter problem. Consider `--cov-mode 1` (target coverage) or a post-cluster
   merge of families whose reps align at ≥90% identity over the shorter sequence, and
   re-measure absence afterwards — the absent-cell count itself will drop before any rescue.
