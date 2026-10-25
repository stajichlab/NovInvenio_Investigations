"""Name-based flags for transposable-element (TE) and repeat Pfam families.

A heuristic label for figures and tables. It is not a test and not a curated list.
It reads the Pfam name as written in the enrichment tables ("Plavaka: Plavaka transposase":
short name, then description) and returns one of:

    "TE"      transposable-element or mobile-element family (transposase, retrotransposon gag/pol,
              integrase, reverse transcriptase, hAT, MULE, Helitron, KDZ/Plavaka, IS605/TnpB,
              retrovirus-like protease, PEG10 and similar TE-derived families)
    "repeat"  generic repeat families that are not TEs (ankyrin, LRR, WD40, TPR, HEAT, Kelch, ...)
    "other"   everything else

HTH_psq and ZBED1 are flagged because, in the 106 yeast_convergence proteomes, 81% of HTH_psq
proteins and 56% of ZBED1_RNase-like proteins also carry a TE-flagged domain (Tc5/DDE_1 and
Dimer_Tnp_hAT). DUF659 is not flagged: 12% of its 115 proteins carry one.

Not flagged on purpose: RNase_H alone, generic "Aspartyl protease" families (PF13650,
PF09668: these also occur outside TEs), restriction endonucleases (Uma2, Mrr_cat),
homing endonucleases (LAGLIDADG), FAR1 (domesticated). Add a pattern here when a
reviewed term needs a flag, and add a test case for it in tests/test_te_flags.py.
"""
from __future__ import annotations

import re

TE_RE = re.compile(
    r"transpos|retrotranspos|retrotran|retrovir|retro-?element|\bTy[13]\b|Ty[13]_|copia|gypsy|\bLTR\b|"
    r"integrase|^(rve|rva)|\brve\b|reverse transcriptase|^(RVT|RT_)|zf-RVT|\bRVP\b|RVP_|"
    r"\bgag\b|gag[-_ ]|polyprotein|^Pol_|\bPol polyprotein|Helitron|(?-i:\bhAT\b)|Dimer_Tnp|\bTnp|_Tnp|Tnp_|"
    r"IS605|OrfB|TnpB|Cas12f1|MULE|Mutator|piggyBac|KDZ|Plavaka|\bCxC\d|^DDE_|DDE superfamily|"
    r"PEG10|Peptidase_A2B|Spuma|mobile element|^Tc[1-5]\b|HTH_Tnp|CENP-B|TNP-like|Plant_tran|"
    r"^HTH_psq|^ZBED\d|"
    r"\b(Hermes|Tc1|mariner|Sleeping Beauty|Harbinger|Pogo|Tigger)\b",
    re.I,
)
REPEAT_RE = re.compile(
    r"^(Ank|LRR|WD40|TPR|HEAT|Kelch|Sel1|PPR|Chromo|Exo_endo_phos$|Mito_fiss_reg|zf-CCHC|PDDEXK|DUF4219)",
    re.I,
)


def te_class(name: str) -> str:
    """Return "TE", "repeat" or "other" for a Pfam name ("ShortName: description")."""
    short = name.split(":")[0].strip()
    if TE_RE.search(name):
        return "TE"
    if REPEAT_RE.match(short):
        return "repeat"
    return "other"
