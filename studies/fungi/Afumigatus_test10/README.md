# Afumigatus_test10: 10-genome A. fumigatus pangenome test set (ingroup only)

Purpose: a small, fast set to test the pangenome workflow end to end and to
iterate on the report figures. It is test set 1 of 2 (set 2:
`../Afumigatus_test_c1c7`).

- 10 A. fumigatus genomes, all `GROUP=IN`, no outgroup.
- `Asfu_Af293`, `Asfu_A1163`: reference strains.
- The best assembly of each Barber cluster 1-7, plus a second genome of
  cluster 1 (the most distinct population).
- "Best assembly": fewest contigs among genomes with BUSCO Complete >= 98.9%.
- `config.csv` is written by `bin/make_test_subsets.py` from
  `../Afumigatus_pangenome/config.csv`; provenance in `DATA_MANIFEST.yaml`.
- Data are read from `../Afumigatus_pangenome/data_dir` (see
  `pangenome_runs.yaml`); nothing is copied.
- No `publish.yaml`: not on the site.

Run: `bin/ni pangenome run --study-dir studies/fungi/Afumigatus_test10 --run in10_v1`
