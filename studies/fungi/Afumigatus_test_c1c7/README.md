# Afumigatus_test_c1c7: A. fumigatus Barber cluster 1 (IN) vs cluster 7 (OUT)

Test set 2 of 2 (set 1: `../Afumigatus_test10`). Tests the ingroup/outgroup
logic (per-group bins, `bin_out`, overlap figures) with two populations of one
species, where the expected pattern is easier to reason about than in the
two-species Coccidioides design.

- 5 genomes of Barber cluster 1 as `IN`, 5 of Barber cluster 7 as `OUT`.
- Clusters 1 and 7 are the most distant pair: mean mash distance 0.0075
  between them, 0.0009 and 0.0021 within.
- Genome choice: fewest contigs among genomes with BUSCO Complete >= 98.9%,
  by `../Afumigatus_test10/bin/make_test_subsets.py`; provenance in
  `DATA_MANIFEST.yaml`.
- Data are read from `../Afumigatus_pangenome/data_dir`; nothing is copied.
- No `publish.yaml`: not on the site.

Run: `bin/ni pangenome run --study-dir studies/fungi/Afumigatus_test_c1c7 --run c1_vs_c7_v1`
