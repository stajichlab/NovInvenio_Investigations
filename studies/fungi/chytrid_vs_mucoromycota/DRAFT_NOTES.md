Draft. Not yet run. No data pulled, so no DATA_MANIFEST.yaml yet (class 2 rule: add provenance when the pull runs). One proteome per species, best BUSCO, genera spread. Prune before running. BUSCO values come from the UniProt proteome record (UniProt does not state the lineage set in the TSV field used here).

Question: chytrids (zoosporic) vs zygomycetes (Mucoromycota, one Zoopagomycota). Literature controls (not re-verified): flagellar proteins (Galindo et al. 2020 Nat Commun; Amses et al. 2022 PNAS) as a loss test in Mucoromycota; CotH expansion in Mucoromycota (Wang et al. 2022 GBE); cell wall synthases (Ruiz-Herrera et al. 2010). Mostly a loss design: run both directions. Not included: Blastocladiomycota (Allomyces UP000054350, Catenaria UP000193411; zoosporic, not chytrid) and Dikarya outgroup. The existing study zoosporic_dikarya covers chytrid vs Dikarya. Zoopagomycota has 3 reference proteomes and 2 have BUSCO below 55 percent, so only one is kept.

## Selected proteomes (UniProt reference proteomes, queried 2026-10-11)

| Short | Species | Group | TaxonGroup | UPID | BUSCO complete % | Proteins |
|---|---|---|---|---|---|---|
| Batrdend | Batrachochytrium dendrobatidis | IN | Chytridiomycota | UP000007241 | 89.7 | 8610 |
| Batrsala | Batrachochytrium salamandrivorans | IN | Chytridiomycota | UP001648503 | 91.2 | 11223 |
| Spizpunc | Spizellomyces punctatus | IN | Chytridiomycota | UP000053201 | 95.6 | 9267 |
| Powehirt | Powellomyces hirtus | IN | Chytridiomycota | UP000318582 | 93.7 | 6533 |
| Geravari | Geranomyces variabilis | IN | Chytridiomycota | UP001212152 | 92.7 | 8696 |
| Syncendo | Synchytrium endobioticum | IN | Chytridiomycota | UP000317494 | 84.2 | 7895 |
| Syncmicr | Synchytrium microbalum | IN | Chytridiomycota | UP000319731 | 90.2 | 6304 |
| Chytconf | Chytriomyces confervae | IN | Chytridiomycota | UP000320333 | 83.6 | 10661 |
| Polystew | Polyrhizophydium stewartii | IN | Chytridiomycota | UP001527925 | 88.3 | 9229 |
| Rhizglob | Rhizoclosmatium globosum | IN | Chytridiomycota | UP000193642 | 82.2 | 16725 |
| Rhizrose | Rhizophlyctis rosea | IN | Chytridiomycota | UP001212841 | 76.0 | 11454 |
| Mucocirc | Mucor circinelloides | OUT | Mucoromycotina | UP000014254 | 99.7 | 12143 |
| Rhizdele | Rhizopus delemar | OUT | Mucoromycotina | UP000009138 | 80.2 | 16971 |
| Lichcory | Lichtheimia corymbifera | OUT | Mucoromycotina | UP000027586 | 91.0 | 13257 |
| Phycblak | Phycomyces blakesleeanus | OUT | Mucoromycotina | UP000077315 | 92.5 | 16155 |
| Absiglau | Absidia glauca | OUT | Mucoromycotina | UP000078561 | 96.7 | 14825 |
| Syncrace | Syncephalastrum racemosum | OUT | Mucoromycotina | UP000242180 | 97.2 | 11037 |
| Choacucu | Choanephora cucurbitarum | OUT | Mucoromycotina | UP000093000 | 95.6 | 11924 |
| Hessvesi | Hesseltinella vesiculosa | OUT | Mucoromycotina | UP000242146 | 96.7 | 11058 |
| Umbevina | Umbelopsis vinacea | OUT | Mortierellomycotina_Umbelopsidomycotina | UP000612746 | 95.2 | 9235 |
| Mortisab | Mortierella isabellina | OUT | Mortierellomycotina_Umbelopsidomycotina | UP000654370 | 94.7 | 9021 |
| Thamspha | Thamnocephalis sphaerospora | OUT | Zoopagomycota | UP000271241 | 73.2 | 6769 |
