# zoosporic_dikarya

IN: 4 zoosporic fungi (Batrdend, Spizpunc, Syncendo, Cateangu). OUT: 5 Dikarya (Scer, Ncra, Anid, Spom, Ccin). Rerun on NovInvenio d03f5a2 on 2026-10-11; see results/zoosporic_dikarya_pre_d03f5a2 for the old run.

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

Sampling for the actin question: IN has only chytrids and Blastocladiomycota (4 species). It has no Holozoa. For the contrast animals + zoosporic vs Dikarya use `eukaryotes/zoosporic_holozoa_vs_dikarya`.
