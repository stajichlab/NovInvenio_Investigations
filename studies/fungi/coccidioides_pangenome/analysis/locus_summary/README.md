# Locus summary figures (C. immitis only)

Uses the scripts in `../../../Afumigatus_pangenome/analysis/locus_summary/` (see its README). Posadasii is not run here.

Strata: `STRATA=0.01,0.03,0.10,0.22` (a guess; the ingroup is about 169 of 529 strains). Run `run_stratified_loci.sh` from the Afumigatus folder with:
`T=<ISLAND_DNA_TARGETS task dir of immitis_in_posadasii_out_v5> CFG=config_immitis_in_posadasii_out.csv REFS=RS,WA_211 PROJECT=config_immitis_in_posadasii_out STRATA=0.01,0.03,0.10,0.22 OUT=<here>/strat`

    CP=<results/immitis_in_posadasii_out_v5/config_immitis_in_posadasii_out/pangenome>
    /usr/bin/python3.12 ../../../Afumigatus_pangenome/analysis/locus_summary/make_locus_figures.py \
      --loci strat/island_loci.json --clades $CP/clade_assignments.tsv \
      --config $CP/samplesheet.with_clades.csv --species "Coccidioides immitis" \
      --outdir fig --label "C. immitis v5 (immitis only)"

No Starship coordinate list exists for Coccidioides, so the overlay is absent.
