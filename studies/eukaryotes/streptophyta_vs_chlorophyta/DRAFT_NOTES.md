Draft. Not yet run. No data pulled, so no DATA_MANIFEST.yaml yet (class 2 rule: add provenance when the pull runs). One proteome per species, best BUSCO, genera spread. Prune before running. BUSCO values come from the UniProt proteome record (UniProt does not state the lineage set in the TSV field used here).

Question: genes gained at the streptophyte node (terrestrialization: cell wall, hormone signalling, stress, plastid signalling) vs Chlorophyta. Literature (not re-verified): Hori et al. 2014 Nat Commun (Klebsormidium); Nishiyama et al. 2018 Cell (Chara); Kunz et al. 2025 Trends Genet warns many land-plant genes predate land plants. Only two reference proteomes exist for streptophyte algae (Klebsormidium nitens, Chara braunii). No Zygnematophyceae, Mesostigma, Chlorokybus or Coleochaete reference proteomes were found, so the node closest to land plants is not sampled. Chara has 35,576 proteins and C. richardii 64,680, so search cost is high. Group IN holds both algae and land plants; split by TaxonGroup if you want the algae-only node.

## Selected proteomes (UniProt reference proteomes, queried 2026-10-11)

| Short | Species | Group | TaxonGroup | UPID | BUSCO complete % | Proteins |
|---|---|---|---|---|---|---|
| Klebnite | Klebsormidium nitens | IN | Streptophyte_algae | UP000054558 | 96.9 | 16251 |
| Charbrau | Chara braunii | IN | Streptophyte_algae | UP000265515 | 72.2 | 35576 |
| Marcpoly | Marchantia polymorpha | IN | Embryophyta | UP000244005 | 87.4 | 22017 |
| Physpate | Physcomitrium patens | IN | Embryophyta | UP000006727 | 88.7 | 47782 |
| Selamoel | Selaginella moellendorffii | IN | Embryophyta | UP000001514 | 77.6 | 33150 |
| Cerarich | Ceratopteris richardii | IN | Embryophyta | UP000825935 | 84.6 | 64680 |
| Ambotric | Amborella trichopoda | IN | Embryophyta | UP000017836 | 90.3 | 27406 |
| Arabthal | Arabidopsis thaliana | IN | Embryophyta | UP000006548 | 100.0 | 39272 |
| Oryzsati | Oryza sativa | IN | Embryophyta | UP000059680 | 84.3 | 48896 |
| Chlarein | Chlamydomonas reinhardtii | OUT | Chlorophyta | UP000006906 | 98.9 | 18840 |
| Volvcart | Volvox carteri | OUT | Chlorophyta | UP000001058 | 96.6 | 14335 |
| Ostrtaur | Ostreococcus tauri | OUT | Chlorophyta | UP000009170 | 98.2 | 7751 |
| Micrpusi | Micromonas pusilla | OUT | Chlorophyta | UP000001876 | 97.0 | 10250 |
| Bathpras | Bathycoccus prasinos | OUT | Chlorophyta | UP000198341 | 96.6 | 7886 |
| Chlovari | Chlorella variabilis | OUT | Chlorophyta | UP000008141 | 92.9 | 9751 |
| Coccsube | Coccomyxa subellipsoidea | OUT | Chlorophyta | UP000007264 | 98.1 | 9799 |
| Auxeprot | Auxenochlorella protothecoides | OUT | Chlorophyta | UP000028924 | 93.9 | 7122 |
| Mononegl | Monoraphidium neglectum | OUT | Chlorophyta | UP000054498 | 74.1 | 16738 |
| Chloprim | Chloropicon primus | OUT | Chlorophyta | UP000316726 | 91.4 | 8712 |
| Tetrobli | Tetradesmus obliquus | OUT | Chlorophyta | UP000256970 | 92.8 | 18631 |
| Raphsubc | Raphidocelis subcapitata | OUT | Chlorophyta | UP000247498 | 96.1 | 13347 |
