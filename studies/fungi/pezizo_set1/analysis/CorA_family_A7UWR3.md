# CorA-family Mg2+/Mn2+ transporter — sequences in `CorA_family_A7UWR3.fa`

Source: pezizo_set1 diamond search results (`results/pezizo_set1/presence_matrix.targets.tsv`,
`results/pezizo_set1/self_hits/*.paralog_cutoffs.tsv`). Query protein: A7UWR3_NEUCR (NCU11312)
and its in-genome paralog Q7SEN2_NEUCR (NCU03312).

Unlike HEX1, A7UWR3 is *not* absent from the outgroups — diamond finds strong hits
(~1e-93..1e-97) in every one of the six outgroup proteomes. Those hits are disqualified from
the presence matrix only because the paralog Q7SEN2 wins head-to-head on the same outgroup
targets (paralog-competition filter 2). Each Pezizomycotina ingroup species independently
carries its own version of the same duplicate pair, which is the point of this table: the
ancestral CorA duplication predates the Pezizomycotina radiation, and each descendant lineage
kept both copies.

| Group | Species | A7UWR3-lineage copy | In-genome paralog (Q7SEN2-lineage) |
|---|---|---|---|
| IN (Pezizomycotina) | Neurospora crassa (Ncra) | A7UWR3 (NCU11312) | Q7SEN2 (NCU03312) |
| IN | Aspergillus fumigatus (Afum) | Q4WTR5 | Q4X211 |
| IN | Arthrobotrys megalospora (Amega) | A0ACF5C4E7 (annotated ALR1_2) | A0ACF5BNE9 (annotated MNR2) |
| IN | Zymoseptoria tritici (Ztri) | F9XHY1 | F9XDN8 |
| IN | Coccidioides immitis (Cimm) | J3K7R2 | J3KI06 |
| OUT | Coprinopsis cinerea (Ccin) | A8NYS6 (annotated MNR2) | A8N143 |
| OUT | Saccharomyces cerevisiae (Scer) | Q08269 (ALR1) / P43553 (ALR2) | P35724 (MNR2) |

Notes:
- Each IN-group pair (e.g. Afum Q4WTR5/Q4X211, mutual reciprocal-best e-value ~1e-99) is its
  own independent paralog pair, not a single tree-wide ortholog/paralog split — expect five
  separate sister-pairs among the ingroup taxa, not one ingroup clade beside one outgroup clade.
- Scer's family is three-way (ALR1/ALR2/MNR2): ALR1 and ALR2 are near-identical (e=0.0, almost
  certainly a *S. cerevisiae*-specific duplication), with MNR2 the more diverged of the three
  (~1e-70 from ALR1/ALR2). Which yeast copy pairs with "A7UWR3-lineage" vs "Q7SEN2-lineage" per
  species differs — see per-outgroup best-hit direction in `presence_matrix.targets.tsv` rather
  than assuming a fixed 1:1 correspondence across taxa.
- CneoH99, Mcir, Nirr, and Spom outgroup copies are not included here (this table only covers
  the two outgroups you asked for, Ccin and Scer); the full six-outgroup hit list is in the
  conversation record / `search_cache/Ncra_vs_*.diamond.tsv.gz` if you want to extend it.
