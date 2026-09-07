# UHM_Koxytoca

*Klebsiella oxytoca* novelty/loss study: 18 UHM metagenome-assembled ingroup
genomes vs. 16 NCBI RefSeq outgroup references, with a *K. pneumoniae* MGH 78578
model-organism gene-name lookup for the ingroup.

## Datasets

- **Ingroup (`Group=IN`, `TaxonGroup=UHM_koxytoca`)** — 18 *K. oxytoca*
  metagenome-assembled genomes from the UHM combined MAG set (binned with
  metashot, genes called with prodigal). Proteins only
  (`data_dir/pep/<Short>.pep.fa`) — no genome DNA or GFF3 is pulled for these;
  see "Why no GFF3" below.
- **Outgroup (`Group=OUT`, `TaxonGroup=Mammal_koxytoca`)** — 16 NCBI RefSeq
  *K. oxytoca* assemblies (`GCF_*`). Protein (`data_dir/pep/<GCF_accession>.faa`,
  same filename as the source `koxytoca_outgroup_faa/` directory) and genome DNA
  (`data_dir/dna/<GCF_accession>.fna`) are both pulled, via
  `studies/bacteria/UHM_Koxytoca/bin/fetch_koxytoca_outgroup_dna.sh` + `studies/bacteria/UHM_Koxytoca/bin/build_koxytoca_config.py` — DNA is
  needed here because `cluster_tool=pairwise`'s VALIDATE workflow only ever runs
  TBLASTN against the outgroup.
- **Model organism** — *K. pneumoniae* MGH 78578 (`Kpn78578`), fetched from both
  NCBI (`GCF_000016305.1`) and UniProt (`UP000000265`). Full provenance and
  regeneration recipe: `../../../config_support/MODELORG_KPN78578_PROVENANCE.md`.

Full per-file provenance (source URL, release, checksum, license): `DATA_MANIFEST.yaml`.

## Why no GFF3

The NCBI Datasets package each outgroup genome comes from (via
`studies/bacteria/UHM_Koxytoca/bin/fetch_koxytoca_outgroup_dna.sh` → `bin/fetch_genome_assembly.py`) includes a
`genomic.gff` alongside the `.fna` we do use — it just isn't pulled into
`config.csv`'s `GFF3` column for this study.

nf_NovInvenio's GFF3 column is **report-only**: it feeds `lib/gff3_genes.py`'s
gene/mRNA position lookup, which adds a `chrom`/`start` column to the interactive
novelties/core/losses HTML reports (`bin/make_report.py` and friends). It plays no
role in clustering, TBLASTN validation, or Pfam/SwissProt annotation, and it's
optional per row — a missing GFF3 just means blank `chrom`/`start` in the report
tables, never an error (`lib/gff3_genes.py:158-159`).

Since that's purely cosmetic here, it was left out. If genomic coordinates in the
report ever become useful for this study, the `.gff` files are already sitting
next to each `.fna` under `data/ncbi/<accession>/extracted/ncbi_dataset/data/<accession>/genomic.gff`
(the ephemeral cache `studies/bacteria/UHM_Koxytoca/bin/fetch_koxytoca_outgroup_dna.sh` already populated) —
copying them into `data_dir/gff3/<GCF_accession>.gff3` and setting each outgroup
row's `GFF3` column in `config.csv` is all that's needed.

## Running

```bash
sbatch studies/bacteria/UHM_Koxytoca/run_Koxytoca.sh
```

See `run_Koxytoca.sh`'s header comment for the config/data_dir rebuild caveat (do not let
`bin/run_study.sh`'s generic fallback rebuild this study — it calls the wrong
builder for this study's species.csv schema).
