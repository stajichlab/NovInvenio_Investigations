Draft. Not yet run. No data pulled, so no DATA_MANIFEST.yaml yet (class 2 rule: add provenance when the pull runs). One proteome per species, best BUSCO, genera spread. Prune before running. BUSCO values come from the UniProt proteome record (UniProt does not state the lineage set in the TSV field used here).

Question: genes specific to Serinales (CUG-Ser1 clade) vs Pichiales and Alaninales (sisters in Pichiomycetes) with Saccharomycetaceae as a more distant outgroup. Literature (Groenewald et al. 2023 Stud Mycol) defines the order; no paper found that lists gained families at its base. Shen et al. 2018 Cell and Feng et al. 2024 bioRxiv report loss-dominated yeast evolution, so expect many absent-in-sister calls to be sister loss. Exploratory. Alaninales has one reference proteome only. TaxonGroup separates the sister and the outgroup. The Metschnikowiaceae clade here is a subset of IN; see metschnikowiaceae_cah.

## Selected proteomes (UniProt reference proteomes, queried 2026-10-11)

| Short | Species | Group | TaxonGroup | UPID | BUSCO complete % | Proteins |
|---|---|---|---|---|---|---|
| Candalbi | Candida albicans | IN | Serinales_Debaryomycetaceae | UP000000559 | 96.4 | 6041 |
| Canddubl | Candida dubliniensis | IN | Serinales_Debaryomycetaceae | UP000002605 | 97.1 | 5829 |
| Candtrop | Candida tropicalis | IN | Serinales_Debaryomycetaceae | UP000002037 | 90.6 | 6226 |
| Loddelon | Lodderomyces elongisporus | IN | Serinales_Debaryomycetaceae | UP000001996 | 87.6 | 5779 |
| Meyeguil | Meyerozyma guilliermondii | IN | Serinales_Debaryomycetaceae | UP000001997 | 90.7 | 5919 |
| Debahans | Debaryomyces hansenii | IN | Serinales_Debaryomycetaceae | UP000000599 | 98.7 | 6284 |
| Schestip | Scheffersomyces stipitis | IN | Serinales_Debaryomycetaceae | UP000002258 | 96.4 | 5797 |
| Spatpass | Spathaspora passalidarum | IN | Serinales_Debaryomycetaceae | UP000000709 | 92.9 | 5973 |
| Hyphburt | Hyphopichia burtonii | IN | Serinales_Debaryomycetaceae | UP000095085 | 94.6 | 5994 |
| Suhotanz | Suhomyces tanzawaensis | IN | Serinales_Debaryomycetaceae | UP000094285 | 96.6 | 5882 |
| Candtenu | Candida tenuis | IN | Serinales_Debaryomycetaceae | UP000000707 | 91.1 | 6051 |
| Babjinos | Babjeviella inositovora | IN | Serinales_other | UP000094336 | 94.8 | 6394 |
| Diutrugo | Diutina rugosa | IN | Serinales_other | UP000449547 | 91.7 | 5811 |
| Yamajaro | Yamadazyma jaroonii | IN | Serinales_other | UP001152531 | 95.0 | 5820 |
| Candauri | Candidozyma auris | IN | Serinales_other | UP000230249 | 94.8 | 5433 |
| Clavlusi | Clavispora lusitaniae | IN | Serinales_other | UP000007703 | 88.0 | 5932 |
| Metsbicu | Metschnikowia bicuspidata | IN | Serinales_other | UP000092555 | 92.7 | 5827 |
| Lizaspar | Lizanozyma spartinae | IN | Serinales_other | UP000790833 | 89.6 | 5041 |
| Komaphaf | Komagataella phaffii | OUT | Pichiales | UP000000314 | 93.5 | 5075 |
| Ogatpara | Ogataea parapolymorpha | OUT | Pichiales | UP000008673 | 93.0 | 5319 |
| Kuracaps | Kuraishia capsulata | OUT | Pichiales | UP000019384 | 93.2 | 6027 |
| Pichmemb | Pichia membranifaciens | OUT | Pichiales | UP000094455 | 86.1 | 5536 |
| Pichkudr | Pichia kudriavzevii | OUT | Pichiales | UP000249293 | 86.5 | 5146 |
| Bretnaar | Brettanomyces naardenensis | OUT | Pichiales | UP000290900 | 87.3 | 5175 |
| Eeninana | Eeniella nana | OUT | Pichiales | UP000662931 | 85.6 | 4703 |
| Ogatpoly | Ogataea polymorpha | OUT | Pichiales | UP000788993 | 90.3 | 6801 |
| Pachtann | Pachysolen tannophilus | OUT | Alaninales | UP000094236 | 93.0 | 5645 |
| Sacccere | Saccharomyces cerevisiae | OUT | Saccharomycetaceae | UP000002311 | 99.6 | 6067 |
| Kluylact | Kluyveromyces lactis | OUT | Saccharomycetaceae | UP000000598 | 99.5 | 5071 |
| Lachther | Lachancea thermotolerans | OUT | Saccharomycetaceae | UP000002036 | 99.8 | 5093 |
| Candglab | Candida glabrata | OUT | Saccharomycetaceae | UP000002428 | 99.4 | 5200 |
| Torudelb | Torulaspora delbrueckii | OUT | Saccharomycetaceae | UP000005627 | 99.9 | 4974 |
| Eremgoss | Eremothecium gossypii | OUT | Saccharomycetaceae | UP000000591 | 99.5 | 4760 |

## Babjeviella inositovora (Babjinos) is a divergent outlier (2026-10-13)

Decision: keep Babjinos in the IN group. Flag it as an outlier when interpreting results. Source: `results/serinales_vs_pichiales/` (pairwise-diamond run).

Observations (counted from `novelties.*.tsv`, `presence_matrix.tsv`, `loss_presence_matrix.tsv`, `loss_tblastn_summary.tsv`):

- Novelty rows: Babjinos 23. The other 17 IN species: 143-205.
- The 2,995 IN candidates (`candidates.txt`): Babjinos carries 23.3%. Every other IN species carries 85-99%, including the other Serinales_other species.
- 2,298 candidates are absent from Babjinos. The other Serinales_other species carry about 90% of them, and the Debaryomycetaceae about 95%.
- Loss direction: 9,813 OUT proteins absent from the Babjinos proteome. This is the second lowest of 18 IN species (median about 12,300). Babjinos keeps more OUT-shared genes than most IN species.
- Completeness: UniProt BUSCO 94.8% (above the IN median) and 6,394 proteins (most of any IN species). BUSCO was not rerun here.

Interpretation: the pattern fits a lineage that diverged before the rest of the IN group. It does not fit an incomplete proteome. Open points: no phylogeny in this repo confirms the branching position. The tblastn rescue covered 169 of 76,006 loss rows, so it does not test genome-wide annotation completeness.

Use: IN-specific candidate counts across the other 17 species are not affected by Babjinos. Counts "shared by all IN" are lowered by Babjinos (about 2,300 candidates are absent from it).
