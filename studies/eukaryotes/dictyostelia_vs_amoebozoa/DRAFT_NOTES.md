Draft. Not yet run. No data pulled, so no DATA_MANIFEST.yaml yet (class 2 rule: add provenance when the pull runs). One proteome per species, best BUSCO, genera spread. Prune before running. BUSCO values come from the UniProt proteome record (UniProt does not state the lineage set in the TSV field used here).

Question: genes specific to Dictyostelia (aggregative multicellularity) vs other Amoebozoa. Literature (not re-verified): Glockner et al. 2016 Nat Commun (about 80 percent of needed proteins already in unicellular relatives; unique part mostly extracellular); Sucgang et al. 2011 Genome Biol. Weak design: only 3 non-Dictyostelia Amoebozoa reference proteomes were found (Acanthamoeba castellanii BUSCO 77.6, Planoprotostelium fungivorum, Entamoeba histolytica BUSCO 52.5). Physarum and Protostelium were not found in the Amoebozoa reference set. Opisthokonts are added as a distant outgroup tier. Non-reference proteomes could be added by hand after checking quality.

## Selected proteomes (UniProt reference proteomes, queried 2026-10-11)

| Short | Species | Group | TaxonGroup | UPID | BUSCO complete % | Proteins |
|---|---|---|---|---|---|---|
| Dictdisc | Dictyostelium discoideum | IN | Dictyostelia | UP000002195 | 93.7 | 12726 |
| Dictpurp | Dictyostelium purpureum | IN | Dictyostelia | UP000001064 | 89.0 | 12347 |
| Hetepall | Heterostelium pallidum | IN | Dictyostelia | UP000001396 | 86.7 | 12351 |
| Cavefasc | Cavenderia fasciculata | IN | Dictyostelia | UP000007797 | 90.6 | 12152 |
| Tieglact | Tieghemostelium lacteum | IN | Dictyostelia | UP000076078 | 91.0 | 10208 |
| Polyviol | Polysphondylium violaceum | IN | Dictyostelia | UP000695562 | 94.5 | 10615 |
| Dictfirm | Dictyostelium firmibasis | IN | Dictyostelia | UP001344447 | 93.3 | 10512 |
| Acancast | Acanthamoeba castellanii | OUT | Amoebozoa_other | UP000011083 | 77.6 | 14939 |
| Planfung | Planoprotostelium fungivorum | OUT | Amoebozoa_other | UP000241769 | 94.5 | 16856 |
| Entahist | Entamoeba histolytica | OUT | Amoebozoa_other | UP000001926 | 52.5 | 7967 |
| Homosapi | Homo sapiens | OUT | Opisthokonta | UP000005640 | 99.5 | 147520 |
| Sacccere | Saccharomyces cerevisiae | OUT | Opisthokonta | UP000002311 | 99.6 | 6067 |
| Monobrev | Monosiga brevicollis | OUT | Opisthokonta | UP000001357 | 78.8 | 9156 |
| Salprose | Salpingoeca rosetta | OUT | Opisthokonta | UP000007799 | 83.1 | 11698 |
