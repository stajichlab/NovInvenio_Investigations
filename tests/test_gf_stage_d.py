"""Tests for bin/gf_stage_d.py (stage D input and the crude class call used for calibration)."""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))


def _load():
    spec = importlib.util.spec_from_file_location("gf_stage_d", ROOT / "bin" / "gf_stage_d.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


sd = _load()


def test_compact_rows_orders_species_level_first_and_strips_tabs():
    recs = [
        {"sentence_id": "g", "match_level": "genus", "uid": "1", "text": "Genus only.", "year": 2001, "section": "A"},
        {"sentence_id": "a", "match_level": "article", "uid": "1", "text": "the\tfungus grows", "year": 2002, "section": "B"},
        {"sentence_id": "f", "match_level": "full", "uid": "2", "text": "Full\nname.", "year": 2003, "section": "C"},
    ]
    rows = sd.compact_rows(recs)
    assert [r["sentence_id"] for r in rows] == ["f", "a", "g"]
    assert "\t" not in rows[1]["text"] and "\n" not in rows[0]["text"]


def test_class_call_reduces_species_level_rows():
    r = lambda sup, lvl="species", ev="described": {"supports": sup, "organism_level": lvl, "evidence_type": ev}
    assert sd.class_call([r("filamentous")]) == "filamentous"
    assert sd.class_call([r("yeast")]) == "yeast"
    assert sd.class_call([r("yeast"), r("filamentous")]) == "dimorphic"            # 50% each
    assert sd.class_call([r("dimorphic"), r("yeast")]) == "dimorphic"
    assert sd.class_call([r("yeast", "genus")]) == "unsourced"                      # genus rows give no species call
    assert sd.class_call([r("other")]) == "unsourced"
    assert sd.class_call([]) == "unsourced"


def test_class_call_ignores_predicted_rows_and_stray_minorities():
    r = lambda sup, ev="described": {"supports": sup, "organism_level": "species", "evidence_type": ev}
    ten_fil = [r("filamentous") for _ in range(10)]
    assert sd.class_call(ten_fil + [r("dimorphic", "predicted")]) == "filamentous"   # a prediction is not an observation
    assert sd.class_call(ten_fil + [r("yeast")]) == "filamentous"                    # 1 of 11 is under 20%
    assert sd.class_call(ten_fil[:4] + [r("yeast"), r("yeast")]) == "dimorphic"      # 2 of 6 reach 20%
    assert sd.class_call([r("dimorphic", "predicted")]) == "unsourced"


def test_split_chunks_and_base_slug():
    rows = [{"i": str(i)} for i in range(250)]
    chunks = sd.split_chunks(rows, 120)
    assert [len(c) for c in chunks] == [120, 120, 10]
    assert sd.split_chunks(rows, 0) == [rows] and sd.split_chunks(rows[:5], 120) == [rows[:5]]
    assert sd.base_slug("mucor_circinelloides.part3") == "mucor_circinelloides"
    assert sd.base_slug("taphrina_deformans") == "taphrina_deformans"


def _row(sup, ev="described", lvl="species", phase="none", sid="s", quote="q"):
    return {"supports": sup, "evidence_type": ev, "organism_level": lvl, "phase_stated": phase,
            "sentence_id": sid, "quote": quote}


def test_phase_call():
    assert sd.phase_call([_row("dimorphic", phase="hyphal"), _row("dimorphic", phase="hyphal"),
                          _row("dimorphic", phase="yeast")]) == "hyphal"           # 2 of 3 = 0.67
    assert sd.phase_call([_row("dimorphic", phase="hyphal"), _row("dimorphic", phase="yeast")]) == "unclear"
    assert sd.phase_call([_row("yeast")]) == "none"
    assert sd.phase_call([_row("yeast", ev="predicted", phase="yeast")]) == "none"   # predicted rows do not count
    assert sd.phase_call([_row("yeast", lvl="genus", phase="yeast")]) == "none"


def test_compare_species_agreement_and_flags():
    a = [_row("filamentous", sid=f"a{i}") for i in range(5)]
    b = [_row("filamentous", sid=f"a{i}") for i in range(3)] + [_row("filamentous", sid="b9")]
    res = sd.compare_species({"p1": a, "p2": b})
    assert res["class_agree"] == "yes" and res["calls"] == {"p1": "filamentous", "p2": "filamentous"}
    assert res["needs_adjudication"] == "no"
    assert res["shared_sentences"] == 3 and res["only"] == {"p1": 2, "p2": 1}


def test_compare_species_flags_class_disagreement():
    p1 = [_row("filamentous", sid=f"f{i}") for i in range(4)]
    p2 = [_row("dimorphic", sid=f"d{i}") for i in range(4)]
    res = sd.compare_species({"p1": p1, "p2": p2})
    assert res["class_agree"] == "no" and "class_disagree" in res["needs_adjudication"]


def test_mixed_forms_is_information_not_a_flag():
    mixed = [_row("yeast", sid=f"y{i}") for i in range(3)] + [_row("filamentous", sid=f"f{i}") for i in range(3)]
    res = sd.compare_species({"p1": mixed, "p2": mixed})
    assert res["class_agree"] == "yes" and res["mixed"] == {"p1": True, "p2": True}
    assert res["needs_adjudication"] == "no"                  # a dimorphic species is not a disagreement


def test_thin_evidence_is_flagged_even_when_passes_agree():
    one = [_row("yeast", sid="y1")]
    res = sd.compare_species({"p1": one, "p2": one})
    assert res["class_agree"] == "yes" and "thin_evidence:p1" in res["needs_adjudication"]


def test_compare_species_flags_one_pass_without_species_evidence():
    p1 = [_row("filamentous", sid="f1")]
    p2 = [_row("filamentous", lvl="genus", sid="g1")]
    res = sd.compare_species({"p1": p1, "p2": p2})
    assert res["class_agree"] == "no" and "no_species_rows_in:p2" in res["needs_adjudication"]


def test_disagreement_rows_lists_minority_and_unique_rows():
    p1 = [_row("filamentous", sid=f"f{i}", quote=f"q{i}") for i in range(4)] + [_row("yeast", sid="y1", quote="loose")]
    p2 = [_row("filamentous", sid="f0", quote="q0"), _row("filamentous", sid="z", quote="only p2")]
    rows = sd.disagreement_rows({"p1": p1, "p2": p2})
    quotes = {(r["pass_tag"], r["quote"]) for r in rows}
    assert ("p1", "loose") in quotes                       # a minority form
    assert ("p2", "only p2") in quotes                     # unique to one pass
    assert ("p1", "q0") not in quotes                      # shared and in the majority form

