Run complete (6 IN, 13 OUT; pairwise-diamond, SLURM 29876233, 2026-10-12). Published on the site. The 4 NCBI/JGI IN proteomes (Neocalif, Neolanat, Neoconst, PiroE2sp) were added after the first draft; see `species.csv`. See "Results" at the end. The proteome table below is the first draft (queried 2026-10-11); it lists 2 IN species only. BUSCO values come from the UniProt proteome record; lineage sets differ (n=255 eukaryota, n=954 metazoa, n=3285 etc.), so values are not comparable between clades.

Question: genes of anaerobic gut fungi (Neocallimastigomycota) absent from other early-diverging fungi. Literature (not re-verified): Murphy et al. 2019 Appl Environ Microbiol (277 horizontal transfer events in 27 transcriptomes; 91.7 percent specific to the group; donors mostly gut anaerobic bacteria); Pratt et al. 2024 Nat Commun. Expected recovery: CAZymes (cellulosome components, dockerin-containing proteins), fermentation and electron-disposal enzymes. Only 2 UniProt reference proteomes exist for the group (Piromyces finnis UP000193719, Anaeromyces robustus UP000193944). Neocallimastix californiae (UP000193920, BUSCO 85.4) has UniProt proteome type 'Excluded' and cannot be used as a reference. Two Neocallimastix proteomes (UP001798839, UP001897618) and Piromyces sp. E2 (UP000195826, BUSCO 47.4) have no usable quality values or are poor. With IN = 2 species, a novel family must appear in both; this is a very strict ingroup. [Superseded: the run used 6 IN species, the 2 UniProt proteomes plus 4 NCBI/JGI proteomes added by hand with provenance.] Options: accept this, add proteomes from JGI or NCBI by hand with provenance (leaves the UniProt-only rule), or use the 3 plus Neocallimastix californiae by non-reference pull (needs change to bin/fetch_uniprot_proteome.py, which only pulls the reference_proteomes FTP). Blastocladiomycota Catenaria has BUSCO 74.1.

## Selected proteomes (UniProt reference proteomes, queried 2026-10-11)

| Short | Species | Group | TaxonGroup | UPID | BUSCO complete % | Proteins |
|---|---|---|---|---|---|---|
| Pirofinn | Piromyces finnis | IN | Neocallimastigomycota | UP000193719 | 86.5 | 10787 |
| Anaerobu | Anaeromyces robustus | IN | Neocallimastigomycota | UP000193944 | 82.8 | 12666 |
| Spizpunc | Spizellomyces punctatus | OUT | Chytridiomycota | UP000053201 | 95.6 | 9267 |
| Batrdend | Batrachochytrium dendrobatidis | OUT | Chytridiomycota | UP000007241 | 89.7 | 8610 |
| Powehirt | Powellomyces hirtus | OUT | Chytridiomycota | UP000318582 | 93.7 | 6533 |
| Geravari | Geranomyces variabilis | OUT | Chytridiomycota | UP001212152 | 92.7 | 8696 |
| Chytconf | Chytriomyces confervae | OUT | Chytridiomycota | UP000320333 | 83.6 | 10661 |
| Syncendo | Synchytrium endobioticum | OUT | Chytridiomycota | UP000317494 | 84.2 | 7895 |
| Allomacr | Allomyces macrogynus | OUT | Blastocladiomycota | UP000054350 | 86.9 | 19092 |
| Cateangu | Catenaria anguillulae | OUT | Blastocladiomycota | UP000193411 | 74.1 | 12599 |
| Rhizdele | Rhizopus delemar | OUT | Mucoromycota | UP000009138 | 80.2 | 16971 |
| Mucocirc | Mucor circinelloides | OUT | Mucoromycota | UP000014254 | 99.7 | 12143 |
| Umbevina | Umbelopsis vinacea | OUT | Mucoromycota | UP000612746 | 95.2 | 9235 |
| Sacccere | Saccharomyces cerevisiae | OUT | Dikarya | UP000002311 | 99.6 | 6067 |
| Neurcras | Neurospora crassa | OUT | Dikarya | UP000001805 | 99.2 | 10255 |

## Results (added 2026-10-13)

Sources: `results/neocallimastigomycota_vs_early_fungi` (not tracked), `analysis/near_miss/outputs/`, `analysis/enrichment/outputs/`. Numbers are counts, not interpretation.

Run:
- IN (6): Pirofinn, Anaerobu (UniProt), Neocalif, Neolanat, Neoconst, PiroE2sp (NCBI/JGI, no UniProt annotation). OUT (13) as in the table above.
- Novelty rule: ingroup fraction 0.75, OUT hits = 0. 32,717 candidate proteins. Only about 9% carry a UniProt annotation, because 4 of 6 IN proteomes have none.
- Dockerin domains are enriched. Most expected CAZymes (cellulosome components) are not candidates.
- Near-miss table (at most 2 OUT species with a hit): the GH48 and PFL cases are explained by single-species OUT hits (Umbevina, Cateangu).

Loss direction (present in at least 0.75 of the 13 OUT species, absent from all 6 IN): 9,173 loss candidates (6.7 percent of the OUT proteome), 99.6 percent annotated. Primary enrichment, whole set, family level, FDR < 0.05: Pfam 288, InterPro 689, GO 731 (GO secondary).

Open items:
- Novelty enrichment is run for the 2 annotated IN species only. Pfam annotation for the 4 NCBI/JGI proteomes needs hmmscan or hmmsearch (pilot on PiroE2sp: 9,371 s on 16 cores for 14,612 proteins, output in `analysis/hmmscan/`). The other 3 are not run.
- Expected-gene control lists have no sourced references.
- The literature statements above are not re-verified.
