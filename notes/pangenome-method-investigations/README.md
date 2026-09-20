# Pangenome method investigations — 2026-09-20

Four investigations into the pangenome subworkflow's assumptions, run against the
529-strain Coccidioides `genus_vs_ureesii` study and (for contrast) the 293-strain
A. fumigatus study. Prose and measurements only — no candidate or sequence data, so
these are class-2 and committed normally per `CLAUDE.md`.

Each fed one or more NovInvenio issues; the issue numbers are the actionable record,
these notes are the evidence behind them.

| note | what it answers |
|---|---|
| `2026-09-20-storage-and-compute-efficiency.md` | Where storage and compute are wasted, and what is gold-plating not worth doing |
| `2026-09-20-phylogrouping-and-relatedness.md` | Whether sub-species structure exists; what relatedness metric to use and what it should feed |
| `2026-09-20-rescue-pass-characterisation.md` | What the TBLASTN rescue pass is actually rescuing |
| `2026-09-20-design-decisions.md` | The seven design decisions agreed with the PI, with the evidence for each |
| `2026-09-20-fragmentation-vs-content.png` | Assembly quality vs pangenome content, 529 strains |

## The finding that ties them together

Three symptoms, one root cause: **`pangenome_tier1_cov = 0.8` fragments truncated gene
models into singleton families.**

- **71.9%** of cells the rescue pass would have flipped come from query proteins shorter
  than 80% of the gene they land on — with the frequency dropping off a cliff at
  **exactly 0.8**, mmseqs' own coverage threshold. The rescue was reading back a
  clustering parameter, not biology.
- Fragmented assemblies carry **2.7x more strain-private families** (median 54 vs 20) and
  **+11% accessory content** on slightly *less* DNA.
- Coccidioides' gain:loss came out **143:1 with 0% ambiguous**, against A. fumigatus'
  **7.1:1 with 11.6%** — where the rescue pass actually worked (5,871,769 cells applied
  vs 0).

That one unvalidated parameter propagates into the accessory genome, the islands, the
Leiden modules, and the gain/loss direction calls.

## Other headline numbers

- **Sub-species structure does not exist** in this dataset. Split-half ARI (split each
  genome's contigs odd/even, sketch and cluster each half independently): species k=2
  gives **1.000**; every within-species partition at k=2..10 gives **0.02-0.26**. Two
  real groups plus a continuum. DAPC would not help.
- **All co-occurrence confounding lives on the species axis.** Measured FPR at
  alpha=0.05: unstratified Fisher **0.327**, the existing k=2 stratified test **0.050**,
  k=10/20/40/80 all 0.033-0.040. So a phylogenetically-aware or mixed-model null is
  explicitly *not* worth building.
- **Two strains are mislabelled**, not hybrids: `485B-1_L_OLD_CPA0023` (labelled immitis,
  0.000072 from a posadasii-labelled sibling in the same isolate series) and `B3476`
  (labelled posadasii, sitting centrally in the immitis cloud). Across all 529 strains,
  **zero** fall in the 0.9-1.5 hybrid window; only these two, at 0.16 and 0.22. The error
  is upstream in the annotation-freeze filenames.
- **The presence matrix is 58x larger than its information content** — it stores the
  words "present"/"absent" for what is one bit; 198 MB carrying 3.4 MB, compressing to
  1.4 MB under `zstd -19`.
- **RESCUE_PASS silently applied zero cells** on the flagship run, wasting 301.6 task-hours
  of TBLASTN, because `zstd -dc` refuses Nextflow's staged symlinks and the exit code was
  ignored. Fixed; guard added.

## Issues these produced

#130 assembly-quality QC · #131 assumptions register · #132 cutoff sensitivity analysis ·
#133 structural rescue criterion · #134 diagnostics surfacing with `--pangenome_strict`

## Caveat

Every number here was measured on an **unrescued** matrix (density 15.63%), because the
rescue pass had silently failed. Findings that depend on presence/absence density — the
gain/loss skew especially — should be re-measured after the clustering question (#132)
is settled and rescue is re-run. The DNA-only findings (relatedness, clade structure,
mislabelled strains) are unaffected.
