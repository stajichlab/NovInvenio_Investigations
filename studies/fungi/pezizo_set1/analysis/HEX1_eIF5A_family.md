# HEX1 / eIF5A family — sequences in `HEX1_eIF5A_family.fa`

Source: pezizo_set1 diamond search results (`results/pezizo_set1/presence_matrix.targets.tsv`,
`results/pezizo_set1/self_hits/*.paralog_cutoffs.tsv`). Query protein: HEX1_NEUCR (P87252,
NCU08332 in FungiDB terms (corrected 2026-09-20: NCU08726 is fl, the C6 zinc-finger regulator fluffy -- NCU08332 is the 176aa woronin body major protein, exact length match to UniProt P87252, confirmed against db/modelorgs/Neurospora_crassa_gene_names_FungiDB.csv)) and its in-genome paralog eIF5A (P38672).

HEX1 (Woronin body major protein) is a lineage-specific duplicate of the universal translation
factor eIF5A, restricted to Pezizomycotina. Its presence-matrix row is 0 across every outgroup —
a genuine absence of signal (see A7UWR3 comparison below), not a paralog-competition artifact.

| Group | Species | HEX1-lineage copy | eIF5A copy (universal) |
|---|---|---|---|
| IN (Pezizomycotina) | Neurospora crassa (Ncra) | P87252 (HEX1_NEUCR) | P38672 (IF5A_NEUCR) |
| IN | Aspergillus fumigatus (Afum) | Q4WUL0 (HEXA_ASPFU) | Q4WK14 |
| IN | Arthrobotrys megalospora (Amega) | A0ACF5C2F1 | A0ACF5BXR3 |
| IN | Zymoseptoria tritici (Ztri) | F9X534 | F9XLE1 |
| IN | Coccidioides immitis (Cimm) | A0A0E1S226 | J3KEP1 |
| OUT | Coprinopsis cinerea (Ccin) | — (no ortholog) | A8NJF7 |
| OUT | Cryptococcus neoformans (CneoH99) | — | J9VQ24 |
| OUT | Mucor circinelloides (Mcir) | — | S2K2D9 + S2K710 (own in-genome paralog pair) |
| OUT | Neolecta irregularis (Nirr) | — | A0A1U7LWG5 |
| OUT | Saccharomyces cerevisiae (Scer) | — | P23301 (IF5A1) + P19211 (IF5A2) |
| OUT | Schizosaccharomyces pombe (Spom) | — | Q9UST4 (IF5A2) + P56289 (IF5A1) |

Mcir, Scer, and Spom each carry two in-genome eIF5A paralogs — likely WGD-derived pairs, not an
extraction artifact — so expect sister pairs for those three taxa in the resulting tree, not
singleton tips.
