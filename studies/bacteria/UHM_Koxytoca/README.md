# UHM_Koxytoca

*Klebsiella oxytoca* novelty/loss study: 18 UHM metagenome-assembled ingroup
genomes vs. 16 NCBI RefSeq outgroup references, with a *K. pneumoniae* MGH 78578
model-organism gene-name lookup for the ingroup.

## Datasets

- **Ingroup (`Group=IN`, `TaxonGroup=UHM_koxytoca`)** — 18 *K. oxytoca*
  metagenome-assembled genomes from the UHM combined MAG set (binned with
  metashot, genes called with prodigal). Protein
  (`data_dir/pep/<Short>.pep.fa`, `Protein_Source=local_faa`) and genome DNA
  (`data_dir/dna/<Short>.dna.fa`, `Genome_Source=local_genome`) are both
  pulled — no GFF3 for these (no NCBI annotation package for local MAGs); see
  "Why no GFF3 for the ingroup" below.
- **Outgroup (`Group=OUT`, `TaxonGroup=Mammal_koxytoca`)** — 16 NCBI RefSeq
  *K. oxytoca* assemblies (`GCF_*`). Protein (`data_dir/pep/<Short>.pep.fa`,
  copied from the source `koxytoca_outgroup_faa/` directory via
  `Protein_Source=local_faa` in `species.csv`) and genome DNA
  (`data_dir/dna/<Short>.dna.fa`, fetched live from NCBI via
  `Genome_Source=ncbi`) are both pulled by the unified
  `bin/build_study_config.py --study-dir studies/bacteria/UHM_Koxytoca` — DNA is
  needed here because `cluster_tool=pairwise`'s VALIDATE workflow only ever runs
  TBLASTN against the outgroup.
- **Model organism** — *K. pneumoniae* MGH 78578 (`Kpn78578`), fetched from both
  NCBI (`GCF_000016305.1`) and UniProt (`UP000000265`). Full provenance and
  regeneration recipe: `../../../config_support/MODELORG_KPN78578_PROVENANCE.md`.

Full per-file provenance (source URL, release, checksum, license): `DATA_MANIFEST.yaml`.

## Why no GFF3 for the ingroup

Since this study's 2026-09-11 migration onto the unified `bin/build_study_config.py`
dispatcher, the 16 outgroup rows DO get a real GFF3 — `Genome_Source=ncbi` pulls
the NCBI Datasets package's `genomic.gff` alongside the `.fna`, and
`build_study_config.py` copies it to `data_dir/gff3/<Short>.gff3` and sets
`config.csv`'s `GFF3` column automatically. Only the 18 ingroup rows
(`Genome_Source=local_genome`, prodigal-called MAGs with no NCBI annotation
package) have no GFF3, and that has no consequence beyond the report tables:

nf_NovInvenio's GFF3 column is **report-only**: it feeds `lib/gff3_genes.py`'s
gene/mRNA position lookup, which adds a `chrom`/`start` column to the interactive
novelties/core/losses HTML reports (`bin/make_report.py` and friends). It plays no
role in clustering, TBLASTN validation, or Pfam/SwissProt annotation, and it's
optional per row — a missing GFF3 just means blank `chrom`/`start` in the report
tables, never an error (`lib/gff3_genes.py:158-159`).

Since that's purely cosmetic, the ingroup's missing GFF3 was left as-is. If
genomic coordinates for the ingroup MAGs ever become useful for this study, a
prodigal/prokka annotation pass producing a `.gff` per MAG, added as a
`local_gff3` row in `species.csv`, is what `build_study_config.py` expects.

## Running

```bash
sbatch studies/bacteria/UHM_Koxytoca/run_Koxytoca.sh
```

Since this study's 2026-09-11 migration onto the unified dispatcher,
`bin/run_study.sh`'s generic "config.csv/data_dir missing -> rebuild" fallback
(`bin/build_study_config.py --study-dir studies/bacteria/UHM_Koxytoca`) is the
correct rebuild path for this study too -- no special-case caveat needed anymore.

## Gotcha: editing modelorgs.yaml does NOT bust `-resume`'s cache

`nf_NovInvenio/workflows/annotate.nf` declares `morgs_config` as a `val` (plain
path string), not a `path` input — deliberate for the (huge) `--pfam_hmm`/
`--swissprot_dmnd` files, an unwanted side effect for this small YAML: Nextflow's
resume cache hashes the *string value* of a `val` input, not the file it points
to, so editing `modelorgs.yaml`'s content while keeping the same path leaves a
`-resume` run fully cached with the OLD annotation, silently (no error, no
warning — `nextflow log <run> -f process,hash,workdir` still shows
`ANNOTATE_MATRIX`/`LOSS_ANNOTATE:ANNOTATE_MATRIX` as `CACHED`).

This actually happened once (species.csv's original `modelorgs.yaml` only had a
`short: Kpn78578` entry, which never matched any real Short code -- see
`modelorgs.yaml`'s own header comment -- fixed to 18 per-ingroup-Short entries,
but a `-resume` run after that fix showed `completed=0 cached=1209`, i.e. it
didn't even try).

**After any modelorgs.yaml content change**, before your next `-resume` run:
```bash
source /etc/profile.d/modules.sh && module load nextflow
cd .nf_launch/bacteria/UHM_Koxytoca
nextflow log <last-run-name> -f process,hash,workdir | grep ANNOTATE_MATRIX
rm -rf <those two workdir paths>   # ANNOTATE:ANNOTATE_MATRIX + LOSS_ANNOTATE:ANNOTATE_MATRIX
```
Everything downstream (MAKE_NOVELTIES/MAKE_REPORT/MAKE_CORE_REPORT/
MAKE_LOSSES_REPORT/MAKE_PDF_REPORT/COLLATE_REPORTS) will correctly cascade-
invalidate on its own once those two tasks produce different output.
