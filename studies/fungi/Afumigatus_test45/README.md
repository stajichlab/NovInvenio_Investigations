# Afumigatus_test45: 45-genome A. fumigatus test set for island generation

Purpose: the smallest set that reliably gives accessory islands, for fast
iteration on island generation, the island/locus view and clinker pages.
`Afumigatus_test10` (9 representatives) gives 0 islands.

- 45 A. fumigatus genomes (`GROUP=IN`) + the 2 UniProt outgroup references
  (`Aslen_ref`, `Neofi_ref`, `GROUP=OUT`), as in `Afumigatus_pangenome`.
- The 45 are the N=45, seed=3 subset of the island genome-count sweep
  (`../Afumigatus_pangenome/analysis/island_genome_count/ISLAND_GENOME_COUNT.md`):
  representative strains of `full_v070`, Af293 forced in, 44 drawn at random.
  A1163 is included. No Barber cluster 6 genome.
- Sweep result for this subset (families reused from the 293-genome run, so a
  real run will differ): 75,194 co-occurring pairs, 3,152 islands, all 7
  high-confidence Af293 Starships hit, 10 of 13 scorable Starships matched by
  presence pattern.
- Sweep threshold: 0 islands at <= 20 representatives, > 1,500 at >= 30.
- `config.csv` is written by `bin/make_test45_config.py`; provenance in
  `DATA_MANIFEST.yaml`. Data are read from `../Afumigatus_pangenome/data_dir`.
- Published on the site (`publish.yaml`, study `Afumigatus_test45`), so island
  output can be compared across runs.

Run: `bin/ni pangenome run --study-dir studies/fungi/Afumigatus_test45 --run isl45_v1`
