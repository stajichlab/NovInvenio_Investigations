Draft. Not yet run. No data pulled, so no DATA_MANIFEST.yaml yet (class 2 rule: add provenance when the pull runs). One proteome per species, best BUSCO, genera spread. Prune before running. BUSCO values come from the UniProt proteome record (UniProt does not state the lineage set in the TSV field used here).

Question: what is shared in Candidozyma (Candida) auris, C. haemuli complex, Clavispora lusitaniae and Metschnikowia that other Serinales lack. Motivation from the maintainer: C. auris is a major human pathogen and biofilm former, and these lineages evolve drug resistance. Literature (not re-verified): placement of the CAH clade in Metschnikowiaceae (Liu et al. 2024 Persoonia); shared small GC-poor centromeres (Narayanan et al. 2021 mBio); haploid state (Du et al. 2020 PLoS Pathog). Thermotolerance differs within the clade (Xiao et al. 2025), so it is not a clade-level target. No paper found tying drug resistance to gained gene families at the clade base. Metschnikowia pulcherrima UP000649328 (BUSCO 78.0) and Australozyma monospora (81.3) are left out for quality. Where drug resistance genes are tested, use curated lists (ERG11, CDR1, MDR1, FKS1, TAC1) as controls, to be sourced from a reference.

## Selected proteomes (UniProt reference proteomes, queried 2026-10-11)

| Short | Species | Group | TaxonGroup | UPID | BUSCO complete % | Proteins |
|---|---|---|---|---|---|---|
| Candauri | Candidozyma auris | IN | Candidozyma | UP000230249 | 94.8 | 5433 |
| Candhaem | Candidozyma haemuli | IN | Candidozyma | UP000244309 | 90.4 | 5240 |
| Candpseu | Candidozyma pseudohaemuli | IN | Candidozyma | UP000241107 | 89.9 | 5132 |
| Candduob | Candidozyma duobushaemuli | IN | Candidozyma | UP000244406 | 88.7 | 5170 |
| Clavlusi | Clavispora lusitaniae | IN | Metschnikowiaceae_other | UP000007703 | 88.0 | 5932 |
| Sunginte | Sungouiella intermedia | IN | Metschnikowiaceae_other | UP000182334 | 95.7 | 5921 |
| Metsbicu | Metschnikowia bicuspidata | IN | Metschnikowiaceae_other | UP000092555 | 92.7 | 5827 |
| Metsaff | Metschnikowia aff. | IN | Metschnikowiaceae_other | UP000292447 | 89.3 | 5795 |
| Austsaop | Australozyma saopauloensis | IN | Metschnikowiaceae_other | UP001338582 | 90.0 | 5012 |
| Candalbi | Candida albicans | OUT | Serinales_other | UP000000559 | 96.4 | 6041 |
| Candtrop | Candida tropicalis | OUT | Serinales_other | UP000002037 | 90.6 | 6226 |
| Loddelon | Lodderomyces elongisporus | OUT | Serinales_other | UP000001996 | 87.6 | 5779 |
| Meyeguil | Meyerozyma guilliermondii | OUT | Serinales_other | UP000001997 | 90.7 | 5919 |
| Debahans | Debaryomyces hansenii | OUT | Serinales_other | UP000000599 | 98.7 | 6284 |
| Schestip | Scheffersomyces stipitis | OUT | Serinales_other | UP000002258 | 96.4 | 5797 |
| Spatpass | Spathaspora passalidarum | OUT | Serinales_other | UP000000709 | 92.9 | 5973 |
| Diutrugo | Diutina rugosa | OUT | Serinales_other | UP000449547 | 91.7 | 5811 |
| Komaphaf | Komagataella phaffii | OUT | Pichiales | UP000000314 | 93.5 | 5075 |
| Ogatpara | Ogataea parapolymorpha | OUT | Pichiales | UP000008673 | 93.0 | 5319 |
| Kuracaps | Kuraishia capsulata | OUT | Pichiales | UP000019384 | 93.2 | 6027 |
| Pichkudr | Pichia kudriavzevii | OUT | Pichiales | UP000249293 | 86.5 | 5146 |
| Pachtann | Pachysolen tannophilus | OUT | Alaninales | UP000094236 | 93.0 | 5645 |
