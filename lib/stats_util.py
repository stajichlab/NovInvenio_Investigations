"""Shared statistics helpers for the ORA enrichment scripts (bin/go_enrichment.py,
bin/domain_enrichment.py). Kept dependency-free (no statsmodels) -- goatools' own
built-in fdr_bh method requires statsmodels' sandbox.stats.multicomp.multipletests,
which no longer exists in current statsmodels releases; using one BH implementation
here for both scripts avoids that broken dependency and keeps the correction
identical across GO/Pfam/InterPro results (DESIGN.md Sec 10: all three should be
comparable, tested with the same correction)."""


def bh_fdr(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR. Returns q-values in the same order as pvalues."""
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvalues[i])
    fdr = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=0):
        i = n - rank
        p = pvalues[idx]
        val = min(prev, p * n / i)
        fdr[idx] = val
        prev = val
    return fdr
