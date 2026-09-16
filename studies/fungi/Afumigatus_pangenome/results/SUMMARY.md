# Pangenome summary -- post-rescue, FULLY CORRECTED (2026-09-15: tblastn per-strain redesign + rescue-position fix for pair classification)

## Pangenome composition

| Bin | Families | % of total |
|---|---:|---:|
| core | 25,436 | 53.0% |
| soft_core | 643 | 1.3% |
| shell | 6,998 | 14.6% |
| cloud | 4,054 | 8.4% |
| singleton | 10,852 | 22.6% |
| **Total** | **47,983** | **100.0%** |

## Assembly completeness (BUSCO, genome mode)

| Metric | Value |
|---|---:|
| Strains | 295 |
| Complete BUSCO % (min / median / max) | 96.7 / 98.8 / 99.4 |
| Strains with >100 contigs (draft-assembly-consistent) | 279 |
| Strains with <=20 contigs (near-chromosome-consistent) | 13 |

## HAC / hacA targeted screen

| Locus | Present in | % of strains |
|---|---:|---:|
| hacA (Afu3g04070, UPR transcription factor) | 293/293 | 100.0% |
| hrmA ortholog (subtelomeric HAC, Starship-mobile) | 8/293 | 2.7% |

## Co-occurring pair classification

| Classification | Pairs | % of total |
|---|---:|---:|
| trans | 2,613,303 | 61.6% |
| trans_unconfirmed | 1,533,331 | 36.1% |
| unexplained_physical | 44,576 | 1.1% |
| ambiguous_linkage | 29,949 | 0.7% |
| insufficient_data | 13,643 | 0.3% |
| starship_explained | 8,890 | 0.2% |
| **Total FDR-significant pairs** | **4,243,692** | **100.0%** |
