# serinales_vs_pichiales_v2 (2026-10-13)

Revision of `serinales_vs_pichiales`. Built from its data, not re-pulled (data_dir and annotations are hard links; the DATA_MANIFEST.yaml records are copied for the kept species, and the config.csv record is rewritten for this subset).

Changes against v1:
- Babjeviella inositovora (Babjinos): Group IN -> OUT, TaxonGroup `Serinales_Babjeviella_outgroup`. Basis: 23 novelty rows against 143-205 for the other IN species; 23.3% of the 2,995 candidates present against 85-99%; 789 OUT proteins found in no other IN species. See `../serinales_vs_pichiales/DRAFT_NOTES.md`.
- IN dropped: Candida dubliniensis (Jaccard 0.983 with C. albicans), Suhomyces tanzawaensis (0.972 with D. hansenii), Hyphopichia burtonii (0.960), Candida tropicalis (0.951).
- OUT dropped: Eremothecium gossypii, Ogataea parapolymorpha, Eeniella nana (first three in a greedy removal order; the next, Candida glabrata, was kept on request).
- Kept by request: Clavispora lusitaniae, Candidozyma auris, Candida glabrata.

Result: 13 IN, 13 OUT. v1 estimates: dropping the OUT species adds about +1.5% proteins with no OUT hit (a proxy, not the final candidate count); dropping the IN species and Babjinos from IN raises OUT proteins with no IN hit from 5,412 to 6,420.
