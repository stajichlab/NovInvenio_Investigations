# `Ncra` model organism provenance and regeneration

Model organism for `studies/fungi/pezizo_set1/modelorgs.yaml`: **Neurospora
crassa OR74A** (UniProt reference proteome `UP000001805`, taxid `367110`).
Ncra is one of this study's own IN-group proteomes, not a separate comparison
reference (unlike `AkkMuc`/`Kpn78578` -- see `MODELORG_AKKMUC_PROVENANCE.md`).

## Sources

- **FungiDB gene-info export**, release **68**, downloaded **~May 2026**.
  `config_support/modelorgs/Neurospora_crassa_gene_names_FungiDB.csv`
  (9,759 rows: `Gene ID`, `Product Description`, `Gene Name or Symbol`, plus
  organism/location/pseudo/transcript-count columns).
  **Confidence: partial.** Release (68) and approximate month (May 2026)
  confirmed by the user (2026-09-08); the exact export URL/method (FungiDB
  Gene page batch download vs. GeneListDownload) and exact download date are
  not confirmed -- this file predates NII and was inherited from
  `/bigdata/stajichlab/jstajich/projects/NovInvenio/db/modelorgs/` without a
  sidecar record. If the exact date/URL is later found (e.g. in an old shell
  history or download log), update this record -- do not silently leave it
  approximate if a precise value turns up.
- **UniProt reference proteome `UP000001805`** (taxid `367110`) -- already
  fully provenance-tracked as part of this study's own data pull; see
  `studies/fungi/pezizo_set1/DATA_MANIFEST.yaml` (`UniProt reference
  proteomes release 2026_02 (10-Jun-2026)`, fetched 2026-09-07). Ncra.pep.fa
  and `studies/fungi/pezizo_set1/annotations/UP000001805.tsv` both derive from
  this pull.
- **`config_support/modelorgs/Ncra_self_id_crosswalk.tsv`** -- NOT fetched,
  derived. Built by `studies/fungi/pezizo_set1/bin/build_ncra_fungidb_crosswalk.py`
  from `Ncra.pep.fa`'s own headers (identity map, see that script's docstring
  for why an actual diamond search isn't needed here). 9,759 rows. Sidecar:
  `Ncra_self_id_crosswalk.tsv.provenance.yaml`.

## Not used by this study

`config_support/modelorgs/Ncra_vs_FungiDB_Ncra.diamond.tsv` -- a **real**
diamond blastp search output, but against a different, non-UniProt local
protein ID scheme (`FC69C3D3_000001-T1`-style query IDs). It does not match
`pezizo_set1`'s `Ncra.pep.fa` (UniProt `tr|ACCESSION|...` headers) at all --
using it here would silently annotate nothing. Left in place in case another
study's Ncra protein set matches its ID scheme; not referenced by
`pezizo_set1/modelorgs.yaml`.

## Why `id_transform: diamond_fasta` without an actual diamond search

`lib/model_organisms.py`'s `diamond_fasta` id_transform is a 2-hop lookup:
`protein_id -> (diamond_hits) -> ref_protein -> (protein_fasta header) ->
gene`. Normally the first hop needs a real similarity search (different ID
spaces, e.g. AkkMuc/Kpn78578's MAG/prodigal locus tags vs. UniProt
accessions). Here it doesn't: `Ncra.pep.fa`'s own UniProt headers already
carry the FungiDB `NCU` gene ID directly in a `GN=` field
(`>tr|A7UWL5|A7UWL5_NEUCR ... GN=NCU10683 ...`, confirmed present on 9,758 of
9,759 headers). So:

- `protein_fasta` = `Ncra.pep.fa` itself, `fasta_gene_field: "GN"` -- no
  synthesized FASTA needed, this is the study's actual proteome file.
- `diamond_hits` = `Ncra_self_id_crosswalk.tsv`, an **identity** map
  (`protein_id` -> itself) -- needed only because `ModelOrgAnnotator`'s code
  path always requires both fields set for `id_transform: diamond_fasta`;
  it is not a similarity search and is exact by construction, not an
  approximation.
- Final `gene_key` = the bare `NCU` accession, which is also what
  `gene_url_template` resolves against (`ModelOrgAnnotator.gene_url()`
  resolves the URL against `gene_key`, i.e. `gene_id_col`'s value -- never
  `gene_name_col`) -- so the FungiDB linkout
  (`https://fungidb.org/fungidb/app/record/gene/{gene}`) always fills in a
  real `NCU` gene ID, never a UniProt accession.

## Regeneration

```bash
python3 studies/fungi/pezizo_set1/bin/build_ncra_fungidb_crosswalk.py
```

Re-run only if `Ncra.pep.fa` changes (e.g. a newer UniProt release is pulled
for this study -- a new study run per `DESIGN.md` Sec 10, not an in-place
refresh). `Neurospora_crassa_gene_names_FungiDB.csv` only needs re-pulling if
FungiDB cuts a new release; that would be a new, separately provenance-dated
file, not a silent overwrite of the release-68 one above.
