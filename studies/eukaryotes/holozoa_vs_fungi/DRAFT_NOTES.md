Draft. Not yet run. No data pulled, so no DATA_MANIFEST.yaml yet. One proteome per species. Prune before running. BUSCO values come from the UniProt proteome record; lineage sets differ (n=255 eukaryota, n=954 metazoa, n=3285 etc.), so values are not comparable between clades.

Question: genes gained in Holozoa (animals plus unicellular relatives) vs Fungi. Literature (not re-verified): Grau-Bove et al. 2017 eLife (burst of gene diversity in the holozoan ancestor; transcription factors, adhesion); Suga et al. 2013 Nat Commun (Capsaspora); Ocana-Pallares et al. 2022 Nature (gain at the animal root, loss at the fungal root); Sebe-Pedros et al. 2016 Cell argues regulation, not gene innovation. The literature is split, so a weak result is possible. Unicellular Holozoa are thin in UniProt: 2 choanoflagellates (Monosiga brevicollis, Salpingoeca rosetta), 1 ichthyosporean (Sphaeroforma arctica, BUSCO 63.1, low), and Capsaspora owczarzaki (Filasterea). Corallochytrium and Creolimax were not found. Metazoa are sampled lightly on purpose (7 of 500+ reference proteomes). Animal proteomes are large (the UniProt count for human is 147,520, far above the gene count; the cause, for example isoforms, was not checked). TaxonGroup splits the unicellular holozoans, early metazoans and bilaterians so a node-level call is possible (for example, unicellular Holozoa + animals vs fungi, or animals only).

## Selected proteomes (UniProt reference proteomes, queried 2026-10-11)

| Short | Species | Group | TaxonGroup | UPID | BUSCO complete % | Proteins |
|---|---|---|---|---|---|---|
| Monobrev | Monosiga brevicollis | IN | Choanoflagellata | UP000001357 | 78.8 | 9156 |
| Salprose | Salpingoeca rosetta | IN | Choanoflagellata | UP000007799 | 83.1 | 11698 |
| Capsowcz | Capsaspora owczarzaki | IN | Filasterea_Ichthyosporea | UP000008743 | 93.7 | 9794 |
| Sphaarct | Sphaeroforma arctica | IN | Filasterea_Ichthyosporea | UP000054560 | 63.1 | 18649 |
| Amphquee | Amphimedon queenslandica | IN | Metazoa_early | UP000007879 | 90.8 | 22831 |
| Tricadha | Trichoplax adhaerens | IN | Metazoa_early | UP000009022 | 89.3 | 11520 |
| Nemavect | Nematostella vectensis | IN | Metazoa_early | UP000001593 | 89.7 | 24435 |
| Mnemleid | Mnemiopsis leidyi | IN | Metazoa_early | UP001729289 | 80.9 | 25028 |
| Caeneleg | Caenorhabditis elegans | IN | Metazoa_bilateria | UP000001940 | 100.0 | 26629 |
| Drosmela | Drosophila melanogaster | IN | Metazoa_bilateria | UP000000803 | 100.0 | 21947 |
| Homosapi | Homo sapiens | IN | Metazoa_bilateria | UP000005640 | 99.5 | 147520 |
| Allomacr | Allomyces macrogynus | OUT | Fungi_zoosporic | UP000054350 | 86.9 | 19092 |
| Spizpunc | Spizellomyces punctatus | OUT | Fungi_zoosporic | UP000053201 | 95.6 | 9267 |
| Batrdend | Batrachochytrium dendrobatidis | OUT | Fungi_zoosporic | UP000007241 | 89.7 | 8610 |
| Pirofinn | Piromyces finnis | OUT | Fungi_zoosporic | UP000193719 | 86.5 | 10787 |
| Rhizdele | Rhizopus delemar | OUT | Fungi_Mucoromycota | UP000009138 | 80.2 | 16971 |
| Mucocirc | Mucor circinelloides | OUT | Fungi_Mucoromycota | UP000014254 | 99.7 | 12143 |
| Sacccere | Saccharomyces cerevisiae | OUT | Fungi_Ascomycota | UP000002311 | 99.6 | 6067 |
| Schipomb | Schizosaccharomyces pombe | OUT | Fungi_Ascomycota | UP000002485 | 81.8 | 5130 |
| Neurcras | Neurospora crassa | OUT | Fungi_Ascomycota | UP000001805 | 99.2 | 10255 |
| Aspenidu | Aspergillus nidulans | OUT | Fungi_Ascomycota | UP000000560 | 94.2 | 10561 |
| Mycomayd | Mycosarcoma maydis | OUT | Fungi_Basidiomycota | UP000000561 | 99.4 | 6805 |
| Crypneof | Cryptococcus neoformans | OUT | Fungi_Basidiomycota | UP000010091 | 99.6 | 7427 |
| Coprcine | Coprinopsis cinerea | OUT | Fungi_Basidiomycota | UP000001861 | 95.8 | 13335 |

## Actin-network literature (added 2026-10-12)

Source: Velle KB, Swafford AJM, Garner RM, Fritz-Laylin LK (2024) Actin network evolution as a key driver of eukaryotic diversification. J Cell Sci 137(15):jcs261660. https://journals.biologists.com/jcs/article/137/15/jcs261660/361509/Actin-network-evolution-as-a-key-driver-of

How this was read: the page text was fetched with WebFetch and summarised by a small model. Only the first 200,000 of 231,822 characters were covered. Quotes below were returned by that tool and were not checked against the PDF. Gene-level claims are not in the text I read unless stated.

What the paper says (as returned):
- "deeply branching fungi (which include chytrid fungi) also have similar phenotypes to animal cells, including the ability to crawl using branched actin networks."
- "The yeast lineage lost the ability to crawl concurrently with the loss of specific actin regulators." The text did not name the lost regulators.
- "Myosin 2 is phylogenetically limited to a single lineage, the Amorphea (which is the group that includes animals, fungi and Amoebozoa), with one exception."
- "Formin family proteins have undergone rampant domain swapping." No formin names were returned.
- Figure 3 shows presence and absence of core actin regulators across eukaryotes, with high- and low-confidence hits.
- It cites Prostak et al. 2021 ("The actin networks of chytrid fungi reveal evolutionary loss of cytoskeletal complexity in the fungal kingdom"), Kollmar and Muhlhausen 2017 (myosin repertoire expansion) and Chalkia et al. 2008 (formin evolution). Not read by me.
- Sampling advice in the paper: study neglected non-opisthokont lineages (SAR, Discoba, Amoebozoa) and early-diverging fungi.
- Not found in the text I read: statements about choanoflagellates, Capsaspora or Blastocladiomycota actin regulators.

What presence/absence can and cannot test here:
- Core regulators (Arp2/3 subunits, SCAR/WAVE, WASP, formins, profilin, cofilin, capping protein, myosins) are widespread in eukaryotes. They are called present in both groups, so they are not candidates.
- A candidate needs a family (or a domain-defined protein) that is in most IN species and in no OUT species. Candidates would be lineage-restricted regulators, for example ones in chytrids and animals but lost in yeasts and filamentous Dikarya.
- Domain swapping, duplication and sequence divergence are not seen by this method. Compare Pfam and UniProt annotations of the candidates and of the conserved families by hand.

Controls to check in the output (names come from general knowledge, not from the paper; verify before use): Arp2/3 subunits (ARPC1-5, ACTR2, ACTR3), WAVE complex (NCKAP1, CYFIP, ABI, BRK1), WASP family, formin FH2 domain proteins, Ena/VASP, coronin, cortactin, myosin 2 heavy chain. Expected: most are present in both groups. Record which are lost in Dikarya.

Sampling decision for this study (2026-10-12):
- Zoosporic fungi stay in OUT. Question here: genes gained in Holozoa and absent from all fungi. Animal and zoosporic-fungus shared genes are ancestral to Opisthokonta, so excluding them is correct for this question.
- For the actin question use the inverse design, `zoosporic_holozoa_vs_dikarya`. IN = Holozoa plus zoosporic fungi. OUT = Dikarya.
- Amoebozoa are not added as OUT here. They share myosin 2 and crawling with animals, so they would remove the actin signal.
