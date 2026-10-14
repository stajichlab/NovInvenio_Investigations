"""vocabulary_coverage() in the study enrichment script: which species have terms per vocabulary.

Script under test (study-specific, loaded by absolute path):
  studies/fungi/zoosporic_dikarya/analysis/enrichment/enrich_novelties_v2.py
"""
import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("scipy")

REPO = Path(__file__).resolve().parent.parent
ENRICH = REPO / "studies/fungi/zoosporic_dikarya/analysis/enrichment/enrich_novelties_v2.py"
spec = importlib.util.spec_from_file_location("enrich_novelties_v2_vc", ENRICH)
enrich = importlib.util.module_from_spec(spec)
spec.loader.exec_module(enrich)


def row(sp, subset, n, pfam, clan, ipr, go):
    return {"species": sp, "subset": subset, "n_proteins": n, "n_pfam": pfam, "n_pfam_clan": clan,
            "n_interpro": ipr, "n_go": go}


def test_pfam_only_species_listed_without_go_interpro():
    cov = [row("Uni", "universe", 100, 70, 40, 80, 60), row("Uni", "candidates", 10, 3, 2, 4, 2),
           row("Hmm", "universe", 50, 30, 20, 0, 0), row("Hmm", "candidates", 5, 1, 1, 0, 0),
           row("ALL_IN", "universe", 150, 100, 60, 80, 60)]
    vc = enrich.vocabulary_coverage(cov)
    assert vc["pfam"]["species_without_terms"] == [] and vc["pfam"]["n_with"] == 2
    assert vc["go"]["species_with_terms"] == ["Uni"] and vc["go"]["species_without_terms"] == ["Hmm"]
    assert vc["interpro"]["n_with"] == 1 and vc["interpro"]["n_species"] == 2


def test_all_annotated_has_no_missing_species():
    cov = [row("A", "universe", 10, 5, 4, 6, 3), row("B", "universe", 10, 5, 4, 6, 3)]
    assert all(not c["species_without_terms"] for c in enrich.vocabulary_coverage(cov).values())
