# Locus summary figures

Quick graphics of the island loci from `island_loci.json`. They show, per locus, which strains have the locus content, in Mash-clade order. They also show where Af293 Starships overlap the loci.

Status: exploratory. Branch `locus-frequency-strata` (NovInvenio PR #237) adds the pipeline option used here.

## Files

- `run_stratified_loci.sh` (SLURM). Reruns ISLAND_DNA_TARGETS, the DNA check and ISLAND_LOCI offline with `--candidate_strata`. Output goes to `strat/`.
- `make_locus_figures.py`. Reads `strat/island_loci.json` and writes `fig/`. It works for any pangenome study. Run it with `/usr/bin/python3.12`. `zstd` must be on PATH.

## Outputs in `fig/` (not in git, about 120 MB)

- `locus_table.tsv`: one row per locus, with type, state counts and Starship overlap.
- `locus_overview.png`, `locus_strips/<key>.png`, `locus_strips.pdf`, `locus_strips_selected.png`.
- `starship_overlay.png`, only when `--starships` is given.
- `idiomorph_candidates.tsv`.

## Strip colours

Gene in place (codes 1, 3). DNA at the site, no gene (6). DNA checked, absent (7). Not DNA-checked (0). Elsewhere or contig break (2, 4, 5). The DNA check runs only for strains with intact flanks.

## Locus types (focal species only)

relocated, unresolved, monomorphic, absent_in_focal, rare_content, swap_candidate, presence_absence, content_variants. Rules are in `classify()`.

## Limits

- The `unresolved` type is the majority (137 of 248 Afumigatus loci, 146 of 250 Coccidioides loci).
- `swap_candidate` is verified on one locus only (MAT, L001 in the Afumigatus run). Treat the others as candidates.
- Starship overlap uses Af293 coordinates only (`../island_genome_count/inputs/af293_starships_table_s7.tsv`, provenance in `DATA_MANIFEST.yaml`). Coccidioides has no Starship list.
- Locus keys are renumbered per run.

## Run

    # Afumigatus (defaults in the script)
    sbatch run_stratified_loci.sh
    # then
    P=<results/full_v072/config/pangenome>
    /usr/bin/python3.12 make_locus_figures.py --loci strat/island_loci.json \
      --clades $P/clade_assignments.tsv --config $P/samplesheet.with_clades.csv \
      --species "Aspergillus fumigatus" --outdir fig --label "A. fumigatus full_v072" \
      --pick L001 --starships ../island_genome_count/inputs/af293_starships_table_s7.tsv \
      --ref_strain Asfu_Af293 --family_positions $P/family_positions.tsv.zst \
      --gene_positions $P/gene_positions.tsv.zst --rescue_positions $P/rescue_positions.tsv

For Coccidioides see `../../../coccidioides_pangenome/analysis/locus_summary/README.md`.
