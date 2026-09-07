# `Kpn78578` model organism provenance and regeneration

Model organism for `studies/bacteria/UHM_Koxytoca/modelorgs.yaml`: **Klebsiella
pneumoniae subsp. pneumoniae MGH 78578** (ATCC 700721), the closest
well-annotated reference for the study's Klebsiella oxytoca ingroup/outgroup.

## Sources (fetched 2026-09-07)

- **NCBI RefSeq assembly `GCF_000016305.1`** (paired GenBank `GCA_000016305.1`,
  BioProject PRJNA31, assembly `ASM1630v1`) -- Complete Genome, 6 replicons: the
  5.3 Mb chromosome `NC_009648.1` (5,049 CDS) plus 5 plasmids, one of which is
  `NC_009653.1` (a 3.5 kb, 5-CDS plasmid -- NOT the chromosome, despite an earlier
  assumption in this study's setup that it was; kept anyway as part of the whole
  assembly per the "use the whole assembly" decision below). Fetched as reference
  data only (genome + GFF3 + protein.faa, all replicons) -- lands in the
  ephemeral, gitignored `data/ncbi/GCF_000016305.1/` cache, never used directly by
  `modelorgs.yaml` (see next point for why).
- **UniProt reference proteome `UP000000265`** (taxid `272620`, mnemonic `KLEP7`)
  -- confirmed via `https://rest.uniprot.org/proteomes/search?query=taxonomy_id:573`
  to be the exact MGH 78578 strain proteome (not the species-level *K. pneumoniae*
  taxon 573), and its own `genomeAssembly.assemblyId` is `GCA_000016305.1` --
  same underlying genome as the NCBI pull above, confirming the two sources agree.
  Lands in `data/uniprot/UP000000265/`, gitignored.

Both fetched via `studies/bacteria/UHM_Koxytoca/bin/fetch_kpn78578_modelorg.py` (module `ncbi_datasets/18.30.1`
on PATH required for the NCBI half).

## Why the UniProt pull, not the NCBI one, drives `modelorgs.yaml`

`lib/model_organisms.py`'s `id_transform: diamond_fasta` (nf_NovInvenio) needs a
`gene_names_csv` and a `protein_fasta` in the **same ID space** -- the diamond hit's
subject ID must resolve, via a `gene=` field in `protein_fasta`'s header, straight
to a key in `gene_names_csv`. Using NCBI's protein.faa (`WP_*` accessions) as
`protein_fasta` but a UniProt-derived `gene_names_csv` (UniProt accessions) would
put those two IDs in different spaces and silently annotate nothing. So both the
lookup table and the diamond target are built from UniProt `UP000000265` alone:

- `config_support/modelorgs/Kpn78578_gene_names_UniProt.tsv` -- `bin/extract_dat_annotations.py`
  run against the UniProt `.dat.gz` (GN/DE/DR lines: gene_name, description,
  go_ids, pfam_ids, pfam_names, interpro_ids, ec_numbers, alphafold_id), keyed by
  UniProt accession.
- `config_support/modelorgs/Kpn78578_protein.faa` -- the UniProt `.fasta.gz`,
  header-rewritten from `>tr|ACCESSION|ENTRY_NAME description...` to
  `>ACCESSION gene=ACCESSION` (UniProt's stock header has no `gene=` field;
  this makes the post-diamond FASTA-header lookup step a same-ID round-trip).
- `config_support/modelorgs/Kpn78578_vs_UHM_ingroup.diamond.tsv` -- best diamond
  hit (`--max-target-seqs 1 --evalue 1e-5`) for every ingroup MAG protein
  (all 18 `koxytoca_ingroup_faa/*.faa` concatenated as query) against
  `Kpn78578_protein.faa`, via `studies/bacteria/UHM_Koxytoca/bin/build_kpn78578_diamond_hits.sh`. 85,053 query
  proteins -> 71,617 hits at that e-value. On a memory-constrained/interactive
  shell, diamond's default block size can get OOM-killed even with ample node
  RAM free -- use `--block-size 0.4 --index-chunks 4 --threads 2` (already set in
  the script) if that happens.

## Regeneration

```bash
source /etc/profile.d/modules.sh
module load ncbi_datasets/18.30.1
python3 studies/bacteria/UHM_Koxytoca/bin/fetch_kpn78578_modelorg.py          # NCBI GCF_000016305.1 (reference only)
                                                 # + UniProt UP000000265 (drives modelorgs.yaml)
module load diamond/2.1.12
studies/bacteria/UHM_Koxytoca/bin/build_kpn78578_diamond_hits.sh              # ingroup MAG proteins vs. Kpn78578_protein.faa
```

Re-run `build_kpn78578_diamond_hits.sh` whenever the ingroup `.faa` set changes;
re-run `fetch_kpn78578_modelorg.py` only if the UniProt release needs bumping
(this is a *new* pull with its own provenance/access_date, not a silent overwrite
of the existing one -- CLAUDE.md's provenance rule).
