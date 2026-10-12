Draft. Not yet run. No data pulled, so no DATA_MANIFEST.yaml yet (class 2 rule: add provenance when the pull runs). One proteome per species, best BUSCO, genera spread. Prune before running. BUSCO values come from the UniProt proteome record (UniProt does not state the lineage set in the TSV field used here).

Question: genes gained in Malassezia (lipid-dependent skin yeasts) vs other Ustilaginomycotina. Positive controls from the literature (Wu et al. 2015 PLoS Genet; Ianiri et al. 2020 PNAS; not re-verified): secreted lipases and phospholipases, aspartyl proteases, flavohemoglobin (horizontal transfer). Expected loss signal: glycosyl hydrolases, fatty acid synthase. Run both directions. Risk: gain set is small (about 64 families reported).

## Selected proteomes (UniProt reference proteomes, queried 2026-10-11)

| Short | Species | Group | TaxonGroup | UPID | BUSCO complete % | Proteins |
|---|---|---|---|---|---|---|
| Malaglob | Malassezia globosa | IN | Malasseziales | UP000008837 | 73.8 | 4271 |
| Malasymp | Malassezia sympodialis | IN | Malasseziales | UP000186303 | 92.3 | 4501 |
| Malarest | Malassezia restricta | IN | Malasseziales | UP000269793 | 86.2 | 4097 |
| Malafurf | Malassezia furfur | IN | Malasseziales | UP000818624 | 81.0 | 4532 |
| Malapach | Malassezia pachydermatis | IN | Malasseziales | UP000037751 | 88.1 | 4201 |
| Malavesp | Malassezia vespertilionis | IN | Malasseziales | UP000232875 | 84.0 | 3791 |
| Malanana | Malassezia nana | IN | Malasseziales | UP001213623 | 81.5 | 3953 |
| Malabras | Malassezia brasiliensis | IN | Malasseziales | UP001216638 | 83.8 | 3911 |
| Malajapo | Malassezia japonica | IN | Malasseziales | UP001217754 | 88.8 | 4299 |
| Malaarun | Malassezia arunalokei | IN | Malasseziales | UP001217582 | 82.4 | 3785 |
| Mycomayd | Mycosarcoma maydis | OUT | Ustilaginales | UP000000561 | 99.4 | 6805 |
| Ustihord | Ustilago hordei | OUT | Ustilaginales | UP000006174 | 99.4 | 7111 |
| Sporreil | Sporisorium reilianum | OUT | Ustilaginales | UP000008867 | 99.8 | 6673 |
| Moesaphi | Moesziomyces aphidis | OUT | Ustilaginales | UP000019462 | 89.5 | 6011 |
| Kalmbras | Kalmanozyma brasiliensis | OUT | Ustilaginales | UP000019377 | 96.1 | 5650 |
| Pseufloc | Pseudozyma flocculosa | OUT | Ustilaginales | UP000323386 | 98.9 | 7767 |
| Sporgram | Sporisorium graminicola | OUT | Ustilaginales | UP000306050 | 97.3 | 6585 |
| Testcype | Testicularia cyperi | OUT | Ustilaginales | UP000246740 | 97.1 | 7036 |
| Tillindi | Tilletia indica | OUT | Tilletiales | UP000077521 | 95.2 | 9656 |
| Tillcari | Tilletia caries | OUT | Tilletiales | UP000836402 | 94.4 | 10036 |
| Meirmilt | Meira miltonrushii | OUT | Exobasidiomycetes | UP000245771 | 94.3 | 7452 |
| Acaringo | Acaromyces ingoldii | OUT | Exobasidiomycetes | UP000245768 | 95.8 | 8019 |
| Jamirose | Jaminaea rosea | OUT | Exobasidiomycetes | UP000245884 | 90.2 | 6858 |
| Micrgluc | Microstroma glucosiphilum | OUT | Exobasidiomycetes | UP000245942 | 91.6 | 6679 |
| Tillwash | Tilletiopsis washingtonensis | OUT | Exobasidiomycetes | UP000245946 | 90.7 | 7007 |
| Cerabomb | Ceraceosorus bombacis | OUT | Ceraceosorales_Violaceomycetales | UP000054845 | 86.8 | 8017 |
| Violpalu | Violaceomyces palustris | OUT | Ceraceosorales_Violaceomycetales | UP000245626 | 95.0 | 8058 |
| Tillanom | Tilletiaria anomala | OUT | Tilletiariales | UP000027361 | 95.7 | 6806 |
