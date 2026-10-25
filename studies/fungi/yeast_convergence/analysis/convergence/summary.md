# yeast_convergence: cross-contrast comparison

Scope `whole_set`, level `family`, significant = FDR <= 0.05, k >= 3, fold > 1. Recurrence: Poisson-binomial null, BH q <= 0.05, in >= 3 origins.

Contrasts: 9 origin, 5 pseudo, 2 reversal

Novelty lists not used (fewer than 2 ingroup species; losses still used): yc_rev_neolecta.

Primary call: pfam, pfam_clan over all contrasts. GO/InterPro (secondary) only over the 4 all-UniProt contrasts; excluded: yc_taphrinomycotina, yc_chaetothyriales, yc_dothideales, yc_capnodiales, yc_tremellomycetes, yc_pucciniomycotina, yc_mucoromycota, yc_pseudo_powdery_mildew, yc_pseudo_dermatophyte, yc_pseudo_arbuscular_mycorrhizal, yc_pseudo_lichen, yc_rev_neolecta.

Missing tables (left out): yc_pseudo_lichen/novelty, yc_pseudo_lichen/loss

## novelty

**yeast_convergent** (primary): 0 terms; secondary GO/InterPro cross-check: 0

**streamlining_generic** (primary): 3 terms; secondary GO/InterPro cross-check: 0
- pfam_clan CL0271 F-box: 3/8 origins, q=0.0248, controls 2/4 (deep 1/1)
- pfam PF00646 F-box: F-box domain: 3/8 origins, q=0.0478, controls 2/4 (deep 1/1)
- pfam PF12937 F-box-like: F-box-like: 3/9 origins, q=0.0478, controls 1/4 (deep 1/1)

Mean Pfam Jaccard of significant sets by role pair: origin-origin 0.031 (n=36), origin-pseudo 0.031 (n=36), origin-reversal 0.000 (n=9), pseudo-pseudo 0.060 (n=6), pseudo-reversal 0.000 (n=4)

## loss

**yeast_convergent** (primary): 0 terms; secondary GO/InterPro cross-check: 0

**streamlining_generic** (primary): 1 terms; secondary GO/InterPro cross-check: 0
- pfam_clan CL0219 RNase_H: 4/9 origins, q=0.00277, controls 1/4 (deep 0/1)

Mean Pfam Jaccard of significant sets by role pair: origin-origin 0.051 (n=36), origin-pseudo 0.034 (n=36), origin-reversal 0.002 (n=18), pseudo-pseudo 0.059 (n=6), pseudo-reversal 0.007 (n=8), reversal-reversal 0.000 (n=1)

## Reversals

- yc_rev_neolecta vs yc_taphrinomycotina (gained_in_yeast_origin__lost_in_reversal): pfam PF09271 LAG1-DNAbind: LAG1, DNA binding
- yc_rev_neolecta vs yc_taphrinomycotina (gained_in_yeast_origin__lost_in_reversal): pfam PF09816 EAF: RNA polymerase II transcription elongation factor

## Focus families (Nagy et al. 2014)

12 significant calls across contrasts; full grid in focus_terms.tsv.
- PF00172 Zn_clus (duplication expected): yc_tremellomycetes [origin] novelty, k=25, fdr=0.00542
- PF00172 Zn_clus (duplication expected): yc_pucciniomycotina [origin] novelty, k=20, fdr=0.0123
- PF04082 Fungal_trans (duplication expected): yc_pseudo_powdery_mildew [pseudo] loss, k=58, fdr=0.00191
- PF01185 Hydrophobin (loss expected): yc_tremellomycetes [origin] loss, k=6, fdr=0.0302
- PF06766 Hydrophobin_2 (loss expected): yc_capnodiales [origin] loss, k=3, fdr=0.0229
- PF00067 p450 (loss expected): yc_tremellomycetes [origin] loss, k=34, fdr=2.79e-11
- PF00734 CBM_1 (loss expected): yc_saccharomycotina [origin] loss, k=30, fdr=5.21e-06
- PF00734 CBM_1 (loss expected): yc_pseudo_powdery_mildew [pseudo] loss, k=15, fdr=2.82e-05
- PF00734 CBM_1 (loss expected): yc_pseudo_ectomycorrhizal [pseudo] loss, k=11, fdr=3.97e-09
- PF00734 CBM_1 (loss expected): yc_pseudo_dermatophyte [pseudo] loss, k=15, fdr=1.08e-10
- PF03443 AA9 (loss expected): yc_saccharomycotina [origin] loss, k=8, fdr=0.0246
- PF03443 AA9 (loss expected): yc_pseudo_ectomycorrhizal [pseudo] loss, k=5, fdr=0.000115
