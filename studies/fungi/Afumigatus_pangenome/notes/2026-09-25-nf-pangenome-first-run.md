# A. fumigatus: first nf_NovInvenio pangenome.nf run (2026-09-25)

## Run

- Pipeline: NovInvenio `af6fd68b69311f75417e06a962c906cfe736af41` (main).
  The run used a `git archive` snapshot at
  `.nf_launch/nf_latest/pipeline_af6fd68/`, so edits to the live checkout
  could not change it mid-run.
- Launcher: `run_pangenome.sh` (head job via sbatch). Study config
  overrides: `stajichlab_queue.config`.
- Input: `config.csv`, 295 rows (293 IN, 2 OUT). Backend mmseqs, tier-1
  0.9 identity / 0.8 coverage (pipeline defaults). Rescue on, with the
  structural filter (default). Island Pfam HMM: `NovInvenio/db/pfam/Pfam-A.hmm`
  (Pfam 38.2).
- Output: `results/nf_latest/output/pangenome/`. Trace:
  `results/nf_latest/trace.29078035.txt`.
- Head job 29078035 (exfab). One attempt. No resume was needed.
- Result: 1,220 tasks completed, 0 failed, 0 cached.

## GFF3 workaround

3 of 295 proteomes are UniProt FASTAs: Aslen_ref, Neofi_ref (the
outgroups) and Asfu_Af293 (IN). Their headers are `tr|ACC|ENTRY ... GN=<locus_tag>`.
Their GFF3s are NCBI/DDBJ annotation keyed by `protein_id=` GenBank accessions.
`pangenome_build_gene_positions.py` matched 0 IDs for these strains and
stopped with a hard error. The other 292 strains (NCBI FASTAs) matched 100%.

`run_pangenome.sh` passes a GFF3 override directory (`--pangenome_gff3_dir`).
The override has symlinks for 292 strains. For the 3 UniProt strains it has a
copy of the GFF3. In that copy, `protein_id=` is set to the FASTA ID whose
`GN=` equals the CDS `locus_tag` (1:1 matches only). Resolved positions:

| Strain | Resolved / FASTA proteins |
|---|---|
| Aslen_ref | 10,248 / 10,351 (99.0%) |
| Neofi_ref | 10,010 / 10,387 (96.4%) |
| Asfu_Af293 | 8,686 / 9,647 (90.0%) |

The old `full_293run/gene_positions.tsv` also has rows for these 3
strains. Those rows use GenBank IDs, which do not match the FASTA IDs. So
the old run most likely had no usable positions for these 3 strains. I did
not check this against its downstream tables.

## Runtime (from the trace)

- Wall clock: 07:59:33 to 12:45:21, which is 4 h 46 min.
- CPU time (realtime x %cpu, summed): 709 cpu-h. Summed task realtime: 185 h.
- TBLASTN_PER_STRAIN (rescue): 295 tasks, 700 cpu-h. The longest task
  was 0.73 h. Before the run I estimated about 650 cpu-h. That estimate
  scaled the 529-strain Coccidioides run (1,182 cpu-h) by strain count.
- All other processes together: about 9 cpu-h. FAMILY_PFAM_SCAN used
  6.0 cpu-h and CLUSTER_TIER1 used 1.9 cpu-h.

## Results

Tier-1 families: 48,001. Frequency bins (rescued matrix) are computed over
121 dedup-representative ingroup strains (293 strains to 123 representatives
at mash distance 0.001; 2 of the 123 are the outgroups):

| Bin | Families |
|---|---:|
| core | 7,166 (14.9%) |
| soft_core | 660 (1.4%) |
| shell | 4,623 (9.6%) |
| cloud | 10,569 (22.0%) |
| singleton | 24,983 (52.0%) |

Rescue funnel: 10,036,064 cells passed the identity/coverage threshold.
8,816,848 were rejected because they overlap a gene in another family.
575,015 were rejected because the representative was short. 277,879 were
rejected as hotspot hits. 288,340 cells were applied. The `rescue_redundancy`
diagnostic triggered (87.9% overlap). The filter is what removes these
cells.

Co-occurring pairs: 1,356,621. Pair classes:

| Class | Pairs |
|---|---:|
| trans_unconfirmed | 707,252 |
| trans | 561,492 |
| insufficient_data | 72,813 |
| unexplained_physical | 10,519 |
| ambiguous_linkage | 4,545 |

- Significant islands: 11,754 (size 2-159 families, median 5). 6,806 of
  them are single-strain islands.
- Leiden modules (resolution 1.0): 11 modules, 8,293 families. The module
  sizes are 2,937, 2,064, 1,486, 885, 575, 178, 152, 8, 4, 2 and 2.
- Island Pfam enrichment: top domains by FDR are APH (q=1.8e-11), DUF3435
  (Starship captain, q=1.9e-9), GloB_C, DDE_1 and DUF3723.
- Assembly-quality confound: max |rho| = 0.08 (threshold 0.3). Status OK.
- Heaps' law gamma = 0.317 (open pangenome).

## Comparison with results/full_293run (bespoke scripts)

Numbers from the old run are taken from its files and from `results/REPORT.md`.

| Metric | full_293run | nf_latest |
|---|---:|---:|
| Tier-1 families | 47,983 | 48,001 |
| core (pre-rescue / rescued) | 6,600 / 25,436 | 7,166 (rescued) |
| soft_core | 582 / 643 | 660 |
| shell | 3,943 / 6,998 | 4,623 |
| cloud | 9,238 / 4,054 | 10,569 |
| singleton | 27,620 / 10,852 | 24,983 |
| trans pairs (rescued) | 2,613,303 | 561,492 |
| co-occurring pairs (rescued) | 4,243,692 | 1,356,621 |
| distinct significant islands | 12,861 (2-837, median 7) | 11,754 (2-159, median 5) |
| Leiden modules at r=1.0 | 10 (largest 2,105) | 11 (largest 2,937) |

- Family counts are close (difference: 18 families). The frequency
  denominator is 121 strains in both runs. I did not check that the old
  clustering parameters were identical.
- The rescued numbers are not directly comparable. The old rescue had no
  structural filter, so it moved about 18,800 families into core. The new
  run's rescued bins are close to the old pre-rescue bins. The old
  pre-rescue pair set has 522,933 trans pairs. That is the nearer
  comparison to the new 561,492.
- Island counts are not directly comparable. The old islands used the old
  rescued pair set, a different island builder and deduplication by exact
  member set. They also had captain/SM-backbone flags. The new run did not
  set a captain HMM (`EMPTY_CAPTAIN_STUB` ran).
- The old module run used the 2.6 M rescued trans pairs. The new run uses
  561 k trans pairs. So the module counts are not directly comparable.
- The old run had a `starship_explained` class (2,416 pre-rescue pairs).
  The new run has no such class because no captain HMM was given.

## Open items

- Pipeline fix needed (not made here): `pangenome_build_gene_positions.py`
  needs a crosswalk for UniProt FASTAs paired with NCBI GFF3s. It should
  map FASTA `GN=` to GFF3 CDS `locus_tag`.
- No captain-gene HMM was set. A rerun with `--pangenome_captain_hmm`
  (DUF3435) would add Starship-aware pair classes.
