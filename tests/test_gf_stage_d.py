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
    r = lambda s, lvl="species": {"supports": s, "organism_level": lvl}
    assert sd.class_call([r("filamentous")]) == "filamentous"
    assert sd.class_call([r("yeast")]) == "yeast"
    assert sd.class_call([r("yeast"), r("filamentous")]) == "dimorphic"
    assert sd.class_call([r("dimorphic"), r("yeast")]) == "dimorphic"
    assert sd.class_call([r("yeast", "genus")]) == "unsourced"      # genus-level rows do not give a species call
    assert sd.class_call([r("other")]) == "unsourced"
    assert sd.class_call([]) == "unsourced"


def test_split_chunks_and_base_slug():
    rows = [{"i": str(i)} for i in range(250)]
    chunks = sd.split_chunks(rows, 120)
    assert [len(c) for c in chunks] == [120, 120, 10]
    assert sd.split_chunks(rows, 0) == [rows] and sd.split_chunks(rows[:5], 120) == [rows[:5]]
    assert sd.base_slug("mucor_circinelloides.part3") == "mucor_circinelloides"
    assert sd.base_slug("taphrina_deformans") == "taphrina_deformans"
