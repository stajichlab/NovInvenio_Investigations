# `AkkMuc` model organism provenance and regeneration

Model organism for `studies/bacteria/UHM_Akkermansia/modelorgs.yaml`:
**Akkermansia muciniphila strain ATCC BAA-835** (the type strain of the genus
*Akkermansia*), the reference genome requested for this study's 6 candidate-genus
ingroup MAGs.

## Sources (fetched 2026-09-07)

- **NCBI RefSeq assembly `GCF_000020225.1`** (paired GenBank `GCA_000020225.1`,
  assembly `ASM2022v1`) -- Complete Genome, single replicon: chromosome
  `NC_010655.1` (GenBank `CP001071`). Fetched as reference data only (genome +
  GFF3 + protein.faa) -- lands in the ephemeral, gitignored
  `data/ncbi/GCF_000020225.1/` cache, never used directly by `modelorgs.yaml`
  (see next point for why).
- **UniProt reference proteome `UP000001031`** (taxid `349741`, mnemonic
  `AKKM8`) -- confirmed via
  `https://rest.uniprot.org/proteomes/search?query=taxonomy_id:349741` to be the
  exact ATCC BAA-835 strain proteome (2,137 proteins), whose sole component
  cross-references GenBank accession `CP001071` -- the same chromosome NCBI's
  `GCF_000020225.1` wraps (confirmed by `fetch_akkmuc_modelorg.py`'s own run:
  "GCA accession for genome/GFF3 pull: GCA_000020225.1"), confirming the two
  sources agree. Lands in `data/uniprot/UP000001031/`, gitignored.

Both fetched via `studies/bacteria/UHM_Akkermansia/bin/fetch_akkmuc_modelorg.py`
(module `ncbi_datasets` on PATH required for the NCBI half).

## Why the UniProt pull, not the NCBI one, drives `modelorgs.yaml`

Same reasoning as UHM_Koxytoca's Kpn78578 model organism
(`config_support/MODELORG_KPN78578_PROVENANCE.md`): `lib/model_organisms.py`'s
`id_transform: diamond_fasta` (nf_NovInvenio) needs a `gene_names_csv` and a
`protein_fasta` in the **same ID space** -- the diamond hit's subject ID must
resolve, via a `gene=` field in `protein_fasta`'s header, straight to a key in
`gene_names_csv`. Using NCBI's protein.faa (locus-tag/`WP_*` accessions) as
`protein_fasta` but a UniProt-derived `gene_names_csv` (UniProt accessions)
would put those two IDs in different spaces and silently annotate nothing. So
both the lookup table and the diamond target are built from UniProt
`UP000001031` alone:

- `config_support/modelorgs/AkkMuc_gene_names_UniProt.tsv` -- `bin/extract_dat_annotations.py`
  run against the UniProt `.dat.gz` (GN/DE/DR lines: gene_name, description,
  go_ids, pfam_ids, pfam_names, interpro_ids, ec_numbers, alphafold_id), keyed
  by UniProt accession. 2,137 proteins.
- `config_support/modelorgs/AkkMuc_protein.faa` -- the UniProt `.fasta.gz`,
  header-rewritten from `>tr|ACCESSION|ENTRY_NAME description...` to
  `>ACCESSION gene=ACCESSION` (UniProt's stock header has no `gene=` field;
  this makes the post-diamond FASTA-header lookup step a same-ID round-trip).
- `config_support/modelorgs/AkkMuc_vs_UHM_ingroup.diamond.tsv` -- best diamond
  hit (`--max-target-seqs 1 --evalue 1e-5`) for every candidate-genus ingroup
  protein (the 6 `data_dir/pep/{C286,C287,C288,C289,C294,C298}.pep.fa` files
  concatenated as query, already `<Short>__`-prefixed by
  `bin/build_akkermansia_config.py`) against `AkkMuc_protein.faa`, via
  `studies/bacteria/UHM_Akkermansia/bin/build_akkmuc_diamond_hits.sh`. 8,385
  hits.

Same `model_organisms.py` gotcha as Koxytoca's `modelorgs.yaml` (documented in
that file's header, repeated in this study's own `modelorgs.yaml`): a
`model_organisms` entry only fires for candidates whose OWN `source_proteome`
Short matches the entry's `short:` -- so `modelorgs.yaml` has one entry per
INGROUP Short (`C286`, `C287`, `C288`, `C289`, `C294`, `C298`), never `AkkMuc`
itself and never the 8 `clade0N` outgroup Shorts.

## Regeneration

```bash
source /etc/profile.d/modules.sh
module load ncbi_datasets
python3 studies/bacteria/UHM_Akkermansia/bin/fetch_akkmuc_modelorg.py         # NCBI GCF_000020225.1 (reference only)
                                                # + UniProt UP000001031 (drives modelorgs.yaml)
module load diamond
studies/bacteria/UHM_Akkermansia/bin/build_akkmuc_diamond_hits.sh             # 6 ingroup proteomes vs. AkkMuc_protein.faa
```

Re-run `build_akkmuc_diamond_hits.sh` whenever the ingroup Short set (or its
`.pep.fa` content) changes; re-run `fetch_akkmuc_modelorg.py` only if the
UniProt release needs bumping (this is a *new* pull with its own
provenance/access_date, not a silent overwrite of the existing one --
CLAUDE.md's provenance rule).
