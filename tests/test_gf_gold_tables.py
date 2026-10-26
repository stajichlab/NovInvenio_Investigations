"""Tests for bin/gf_gold_tables.py (independent gold rows from tables of review and comparative papers)."""
import importlib.util
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))


def _load():
    spec = importlib.util.spec_from_file_location("gf_gold_tables", ROOT / "bin" / "gf_gold_tables.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gt = _load()

XML = b"""<?xml version="1.0"?>
<pmc-articleset><article><front><article-meta>
<article-id pub-id-type="pmc">123</article-id><article-id pub-id-type="pmid">456</article-id>
<title-group><article-title>Dimorphism in fungi</article-title></title-group></article-meta></front>
<body><sec>
<table-wrap id="T1"><label>Table 1</label><caption><title>Growth forms of selected fungi</title></caption>
<table><thead><tr><th>Species</th><th>Growth form</th></tr></thead>
<tbody>
<tr><td><italic>Yarrowia lipolytica</italic></td><td>dimorphic (yeast and hyphae)</td></tr>
<tr><td>Saccharomyces cerevisiae</td><td>yeast</td></tr>
<tr><td>Unknown species</td><td>filamentous</td></tr>
</tbody></table></table-wrap>
<table-wrap id="T2"><label>Table 2</label><caption><title>Drug susceptibility</title></caption>
<table><tbody><tr><td>Yarrowia lipolytica</td><td>MIC 0.5</td></tr></tbody></table></table-wrap>
</sec></body></article></pmc-articleset>"""


def test_table_rows_reads_label_caption_and_cells():
    meta, tables = gt.parse_tables(XML)
    assert meta["pmcid"] == "123" and meta["pmid"] == "456"
    assert [t["label"] for t in tables] == ["Table 1", "Table 2"]
    assert tables[0]["caption"] == "Growth forms of selected fungi"
    assert tables[0]["rows"][1] == ["Yarrowia lipolytica", "dimorphic (yeast and hyphae)"]


def test_match_rows_needs_species_and_a_growth_term():
    pats = {"Yarrowia lipolytica": [re.compile(r"Yarrowia lipolytica", re.I)],
            "Saccharomyces cerevisiae": [re.compile(r"Saccharomyces cerevisiae", re.I)]}
    _meta, tables = gt.parse_tables(XML)
    hits = gt.match_rows(tables, pats)
    got = {(h["species"], h["table"]) for h in hits}
    assert ("Yarrowia lipolytica", "Table 1") in got
    assert ("Saccharomyces cerevisiae", "Table 1") in got
    assert ("Yarrowia lipolytica", "Table 2") not in got          # no growth term in the row or caption
    assert not any(h["species"] == "Unknown species" for h in hits)
    assert all(len(h["row_text"].split()) <= 40 for h in hits)


def test_noise_captions_are_skipped_even_with_growth_words():
    xml = XML.replace(b"Growth forms of selected fungi", b"MIC results against yeast and filamentous fungi")
    pats = {"Yarrowia lipolytica": [re.compile(r"Yarrowia lipolytica", re.I)]}
    _meta, tables = gt.parse_tables(xml)
    assert gt.match_rows(tables, pats) == []


def test_sentence_hits_need_species_and_specific_term():
    pats = {"Taphrina deformans": [re.compile(r"Taphrina deformans", re.I)]}
    secs = [("Intro", "Taphrina deformans is a dimorphic fungus that infects peach. It grows on leaves. "
                      "Taphrina deformans was sequenced in 2013.")]
    hits = gt.sentence_hits(secs, pats)
    assert len(hits) == 1 and "dimorphic" in hits[0]["row_text"]
    assert gt.sentence_hits([("x", "Yeast and mold grow well.")], pats) == []


def test_short_phrase_is_a_row_substring():
    pats = {"Yarrowia lipolytica": [re.compile(r"Yarrowia lipolytica", re.I)]}
    _meta, tables = gt.parse_tables(XML)
    h = gt.match_rows(tables, pats)[0]
    assert h["row_text"].startswith("Yarrowia lipolytica | dimorphic")
