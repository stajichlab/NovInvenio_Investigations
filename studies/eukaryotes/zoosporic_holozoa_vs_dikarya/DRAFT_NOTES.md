# zoosporic_holozoa_vs_dikarya

Run complete (2026-10-12, SLURM 29887139, pairwise-diamond, NovInvenio main d03f5a2). Species rows copied from `holozoa_vs_fungi`, `chytrid_vs_mucoromycota` and `zoosporic_dikarya` (UniProt reference proteomes release 2026_03; genomes from NCBI). `config.csv` and `DATA_MANIFEST.yaml` are committed. The report, alignments, enrichment (novelty and loss) are published on the site. See "Results" at the end.

Question: which protein families are shared by Holozoa and zoosporic (flagellated) fungi but are absent from Dikarya? This is the inverse of `holozoa_vs_fungi`. It targets traits kept in animals and early-diverging fungi and lost in yeasts and filamentous Dikarya, for example the actin regulators linked to branched-actin crawling and the flagellum.

Design:
- IN (25): 11 Holozoa (2 choanoflagellates, Capsaspora, Sphaeroforma, 4 early metazoans, 3 bilaterians), 11 Chytridiomycota, 2 Blastocladiomycota (Allomyces, Catenaria), 1 Neocallimastigomycota (Piromyces finnis).
- OUT (7): 4 Ascomycota, 3 Basidiomycota (Dikarya).
- Mucoromycota are left out. They lost the flagellum. Whether they kept crawling is not stated in the text I read. They could be added as context later.
- Proteomes: canonical FASTA only (`bin/fetch_uniprot_proteome.py` does not pull the `_additional` isoform file). Human has 20,652 entries.
- Default thresholds: ingroup fraction 0.75 means a family needs at least 19 of 25 IN species. Holozoa alone (11) cannot pass. Zoosporic alone (14) cannot pass. The rule does not force a candidate into both groups broadly: it needs at least 5 of 11 Holozoa and at least 8 of 14 zoosporic fungi. A per-clade minimum would need a separate filter.
- Weak points: Neocallimastigomycota and some chytrids are derived. Low-BUSCO proteomes (Sphaeroforma 63.1) are kept on purpose. One Dikarya homolog removes a family from the candidates. Only 7 Dikarya are sampled, so a gene lost in all 7 may still exist in other Dikarya.
- Not tested: Saccharomycotina other than S. cerevisiae, Ustilaginomycotina other than Mycosarcoma, Pucciniomycotina.

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
- Core regulators (Arp2/3 subunits, WASP, profilin, cofilin, capping protein, myosins) are widespread in eukaryotes. The run shows most of them present in both groups, so they are not candidates. This expectation did not hold for the WAVE complex and some formins; see "Results".
- A candidate needs a family (or a domain-defined protein) that is in most IN species and in no OUT species. Candidates would be lineage-restricted regulators, for example ones in chytrids and animals but lost in yeasts and filamentous Dikarya.
- Domain swapping, duplication and sequence divergence are not seen by this method. Compare Pfam and UniProt annotations of the candidates and of the conserved families by hand.

Controls to check in the output (names come from general knowledge, not from the paper; verify before use): Arp2/3 subunits (ARPC1-5, ACTR2, ACTR3), WAVE complex (NCKAP1, CYFIP, ABI, BRK1), WASP family, formin FH2 domain proteins, Ena/VASP, coronin, cortactin, myosin 2 heavy chain. Expected: most are present in both groups. Record which are lost in Dikarya.


## Results (added 2026-10-13)

Sources: `results/zoosporic_holozoa_vs_dikarya` (not tracked), `analysis/enrichment/outputs/`, `analysis/loss_near_miss/outputs/`. Numbers are counts, not interpretation.

Run:
- 1,272 tasks succeeded. One `BUILD_PRESENCE_MATRIX` task failed with out-of-memory (exit 137) and succeeded on retry.
- 32 proteomes. Human has 20,652 canonical entries.
- `config.csv` has no GFF3 for Amphimedon or human. The effect was not checked.
- `candidates.pfam.tblout` and `candidates.swissprot.tsv` are empty: no `--pfam_hmm` or SwissProt DB was passed. Pfam names come from UniProt annotation only.

Novelty direction (present in at least 19 of 25 IN, absent from all 7 OUT):
- 6,615 candidate proteins. The rule recomputed from `presence_matrix.tsv` gives the same set.
- 125 to 413 per IN species. 2,784 come from Holozoa proteomes and 3,831 from zoosporic fungi.
- The pipeline's 2,958 candidate clusters are not ortholog families: 2,161 are singletons.
- 670 candidates have `broad` outgroup evidence and 1,385 have a TBLASTN hit in at least one Dikarya genome. Some absences may be false.
- 723 candidates are in 6 or fewer Holozoa.

Actin controls (Pfam presence among the novelty candidates; `outputs/proteome_clustered/controls.tsv`):
- Candidates: NCKAP1 (PF09735) 22, CYFIP (PF05994) 26, SCAR_helical (PF27569) 15, WH2 motif 5, formin FH2 7, formin Drf domains 5, coronin 2, cortactin HS1_rep 1.
- Not candidates (0): ARPC2, ARPC3, ARPC4, ARPC5, WASP (WH1), ABI, VASP_tetra, Myosin_tail_1, Myosin head.
- So the WAVE complex core (NCKAP1, CYFIP) passes the novelty rule. The earlier expectation that core regulators are not candidates was wrong for these.
- FH2 and Drf formin hits have Dikarya near-misses (9 and 5 proteins in 1 OUT species). They may be domain-level artifacts. Not checked.
- The design cannot separate flagellum-linked loss from crawling-linked loss, because Dikarya lost both.

Enrichment, novelty direction (whole set, family level, FDR < 0.05 / terms tested): Pfam 297/6,027, InterPro 564/12,289, GO 532/13,705 (GO secondary, computational evidence). Candidate units with a Pfam annotation: 92%.

Loss direction (present in at least 6 of 7 Dikarya, absent from all 25 IN; the flipped question):
- 332 proteins, 248 families. The rebuild from the diamond search cache reproduces the pipeline's loss matrix and `loss_candidates.txt` exactly (symmetric difference 0).
- Looser sets (at most K IN species carry the gene): K=1 514, K=2 716, K=3 931 proteins.
- Enrichment, primary: Pfam 26/1,963, InterPro 49/3,927, GO 111/5,440. K<=2: Pfam 46/1,963, InterPro 86/3,930, GO 159/5,448.
- Small set, so low power. Some terms have k = n at 5 or 6 units.

Open items:
- The Velle 2024 quotes above are still unchecked against the PDF.
- Annotation provenance (UniProt and Pfam release behind `annotations/*.tsv.gz`) is not recorded.
- Per-species enrichment results were produced and not reviewed.
