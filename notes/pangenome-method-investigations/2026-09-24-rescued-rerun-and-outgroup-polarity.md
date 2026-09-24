# Rescued re-run results and the outgroup-polarity problem — 2026-09-24

Follows `2026-09-23-session-handoff.md`. The #133 validation re-run (SLURM job
29026388, pipeline pinned at `3277bdb`) completed after 11 h 36 min. This note
records what it shows, and the review of `2026-09-21-state-of-the-analysis.md`
section 5 that the handoff asked for if the gain:loss ratio did not move.

Run: `studies/fungi/coccidioides_pangenome/results/rescue_structural_genus_vs_ureesii/`
(529 ingroup strains + 1 outgroup, *U. reesii* UAMH1704, Short `Uree`).

## 1. Gain:loss did not move

Counted from `cooccurring_pairs` column `direction_a`. The "rows" method (every
significant pair, labelled by family_a's direction) reproduces the handoff's
99.3% / 0.7% on the old run, so it is the matched comparison.

| run | method | gain | loss | gain:loss | n |
|---|---|---|---|---|---|
| old, unrescued (`mmseqs_genus_vs_ureesii`) | rows | 99.31% | 0.69% | 144:1 | 13,073,012 pairs |
| new, rescued (`rescue_structural_genus_vs_ureesii`) | rows | 99.28% | 0.72% | 138:1 | 18,412,084 pairs |
| old | unique family_a | 99.19% | 0.81% | 122:1 | 14,820 |
| new | unique family_a | 98.83% | 1.17% | 85:1 | 16,969 |

Section 5's success test was "gain:loss moves from 143:1 toward 7:1 with a
non-zero ambiguous fraction". It failed on the first half. The second half
cannot pass for this study at all: see section 3.

## 2. #167 diagnostics on this run (hand run)

NovInvenio `origin/main` `63c81f2` (includes merged PR #167), run by hand because
the pinned run predates `ASSEMBLY_QUALITY_QC`/`DIAGNOSTICS`. Outputs and the job
script: `results/rescue_structural_genus_vs_ureesii/manual_qc_pr167/` (gitignored).
The rescue funnel file was rebuilt from `RESCUE_PASS`'s `.command.err` counts
(the pinned run did not write `--funnel_tsv`).

| diagnostic | status | value |
|---|---|---|
| `assembly_quality_confound` | triggered | max \|rho\| 0.42 (accessory_present vs n_contigs, partial on total_length), n = 529 |
| `rescue_redundancy` | triggered | 15,604,028 / 20,005,269 (78.0%) rescuable hits overlap a different family's gene |

Spearman rho, n = 529:

| comparison | raw | partial (total_length) |
|---|---|---|
| accessory_present vs n_contigs | 0.362 | 0.424 |
| accessory_present vs N50 | -0.328 | -0.383 |
| core_missing vs n_contigs | 0.034 | 0.015 |
| core_missing vs N50 | -0.198 | -0.192 |

Proteins in strain-private families lie within 1 kb of a contig end 39.01% of
the time (n = 14,871), against 5.91% of all 4,504,703 proteins (6.6x).
The assembly-fragmentation confound persists after rescue. These |rho| values
are lower than the ±0.53 in section 2 of the state-of-the-analysis note; that
figure's source matrix is not recorded, so the change is not attributed to the
rescue fix here.

## 3. Why "0% ambiguous" means nothing with one outgroup

`polarize_direction` (NovInvenio `bin/pangenome_cooccurrence.py`) returns
`ambiguous` only when `outgroup_total == 0` or the outgroup strains disagree
(`0 < present < total`). With exactly one outgroup strain, disagreement is
impossible: every family is `gain` or `loss`. `config_genus_vs_ureesii.csv` has
529 IN and 1 OUT. A. fumigatus (`config.csv`) has 293 IN and 2 OUT
(`Aslen_ref`, `Neofi_ref`), which is why it can have an ambiguous fraction.

## 4. Section-5 review: what else is wrong

### The outgroup is too divergent for the pipeline's identity cutoffs

Tier-1 clustering needs `pangenome_tier1_min_id = 0.9`, `pangenome_tier1_cov = 0.8`;
the rescue pass needs `pangenome_rescue_min_pident = 90.0`,
`pangenome_rescue_min_qcov = 80.0` (values at `3277bdb`). A *U. reesii* ortholog
of a Coccidioides gene mostly fails both, so it never joins the Coccidioides
family, and the family reads as "absent in the outgroup", i.e. "gain".

Direct test, DIAMOND blastp best hit (`-k 1`) against *C. immitis* RS:

| query | proteins | with a hit | median identity | pass ≥90% id and ≥80% qcov |
|---|---|---|---|---|
| *U. reesii* UAMH1704 | 7,397 | 6,828 (92.3%) | 77.4% | 823 (11.1%) |
| *C. posadasii* C735-TKO (control) | 8,489 | 7,763 (91.4%) | 98.8% | 7,070 (83.3%) |

The same effect in the pangenome itself (new run, `per_strain_summary.tsv` and
the rescued matrix):

- `Uree` has 7,470 families, of which **6,518 (87.3%) are singletons**. The
  529 ingroup strains average 9,593 families and 40.4 singletons.
- `Uree` is present in only **11.7% of core families** (5,759 core families).
  This matches the 11.1% of its proteins that pass the cutoffs.
- The rescue pass added only **74** `genome_only` cells to `Uree` (same 90/80
  cutoffs, TBLASTN).

Per frequency bin (new run; "present" = present or genome_only, as in
`PresenceMatrix.is_present`):

| bin | families | outgroup present | polarised (family_a) | gain | loss |
|---|---|---|---|---|---|
| core | 5,759 | 11.7% | 0 | – | – |
| soft_core | 507 | 7.7% | 0 | – | – |
| shell | 6,445 | 2.0% | 6,437 | 6,306 | 131 |
| cloud | 17,764 | 0.6% | 10,532 | 10,465 | 67 |
| singleton | 23,843 | 27.3% | 0 | – | – |

(Only shell and cloud families enter the co-occurrence screen, so only they are
polarised. The 27.3% of singletons "present" in the outgroup are `Uree`'s own
private families.)

Conclusion: **the 138:1 skew is set by outgroup divergence against a 90%
identity cutoff, not by absence calling in the ingroup.** A shell or cloud
family is called "loss" only when the outgroup carries it, and the outgroup
cannot be seen to carry it at 90% identity. The rescue fix could not change
this, because the rescue uses the same 90% cutoff.

### Comparison with A. fumigatus

Same script, `Afumigatus_pangenome/results/full_293run/*.rescued.tsv`:

| bin | families | outgroup: all present | any present | any via genome_only |
|---|---|---|---|---|
| core | 25,436 | 53.7% | 72.4% | 45.2% |
| shell | 6,998 | 10.0% | 25.7% | 20.4% |
| cloud | 4,054 | 7.3% | 21.0% | 15.0% |

Its outgroups are seen in far more families, and much of that comes from the
rescue pass (9,616 and 12,446 `genome_only` cells in `Aslen_ref` / `Neofi_ref`).

**Caveat:** this method does not reproduce the A. fumigatus reference figures
(77.5% / 10.9% / 11.6%). It gives 83.92% / 3.17% / 12.91% by pair rows, and
75.4% / 9.3% / 15.3% by unique family_a (6,690 / 827 / 1,358). The source of the
reference figures is not recorded in these notes. So the A. fumigatus "7.1:1"
is not a method-matched benchmark for Coccidioides; only the old-vs-new
Coccidioides comparison above is.

## 5. What this changes

- Gain/loss direction for `genus_vs_ureesii` is not interpretable as it is
  computed now. The `gain_loss_skew` diagnostic (#134, not yet wired) would be
  measuring outgroup divergence here.
- Islands, Leiden modules and co-occurrence p-values do not use polarity; this
  finding does not by itself invalidate them.
- The two mislabelled strains, the relatedness results and the assembly-quality
  confound are unaffected.

Options (not decided; each needs its own test):

1. **Polarise the outgroup at a lower identity.** Keep ingroup clustering at
   90%, but call outgroup presence per family from a separate, lower-identity
   search (e.g. the family representative against the outgroup proteome/genome
   at the ~77% level measured above). Cheapest change. Needs a chosen cutoff and
   a check against known orthologs.
2. **A closer or second outgroup.** A second outgroup would also make
   `ambiguous` possible. Depends on which genomes exist; not checked.
3. **Tree-based polarity.** `pangenome_cooccurrence.py` already has
   `direction_a_tree` (Dollo, `lib/ancestral_states.py`, issue #140). Needs a
   tree; the phylogrouping note says a tree must come from markers/SNPs, not mash.
4. **Report polarity as not computed for this study** until one of the above is
   done.

## Reproduce

- Gain:loss counts: `awk` over `cooccurring_pairs.tsv[.zst]`, column
  `direction_a`, counting rows and first occurrence per `family_a`.
- Per-bin table: `bin/pangenome_outgroup_polarity.py --matrix <run>/presence_matrix.rescued.tsv
  --freq <run>/frequency_table.tsv --pairs <run>/cooccurring_pairs.tsv.zst --outgroup Uree`
  (A. fumigatus: `--outgroup Aslen_ref Neofi_ref` on `full_293run/*.rescued.tsv`).
- Identity test: `diamond makedb --in RS.pep.fa`; `diamond blastp -k 1
  --outfmt 6 qseqid sseqid pident qcovhsp scovhsp` for each query; first hit per
  query; pass = pident ≥ 90 and qcovhsp ≥ 80.
