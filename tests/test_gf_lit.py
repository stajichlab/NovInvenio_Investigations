"""Tests for bin/gf_lit.py (growth-form literature pipeline)."""
import argparse
import importlib.util
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location("gf_lit", ROOT / "bin" / "gf_lit.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gf = _load()


NODES = (
    "1\t|\t1\t|\tno rank\t|\t.\t|\n"
    "2\t|\t1\t|\tdomain\t|\t.\t|\n"
    "4829\t|\t2\t|\tspecies\t|\t.\t|\n"
    "5341\t|\t2\t|\tspecies\t|\t.\t|\n"
    "1220926\t|\t4829\t|\tno rank\t|\t.\t|\n"
)

NAMES = (
    "4829\t|\tAbsidia glauca\t|\t\t|\tscientific name\t|\n"
    "4829\t|\tAbsidia glauca Hagem\t|\tauthority\t|\tscientific name\t|\n"
    "4829\t|\tTieghemella glauca\t|\t\t|\tsynonym\t|\n"
    "5341\t|\tAgaricus bisporus\t|\t\t|\tscientific name\t|\n"
    "5341\t|\tAgaricus brunnescens\t|\t\t|\tsynonym\t|\n"
    "5341\t|\tcultivated mushroom\t|\t\t|\tgenbank common name\t|\n"
)


def _taxonomy(tmp_path):
    d = tmp_path / "taxonomy"
    d.mkdir()
    (d / "nodes.dmp").write_text(NODES)
    (d / "names.dmp").write_text(NAMES)
    return d


def _species_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    header = "Species,Short,Taxon_ID,Group,TaxonGroup\n"
    path.write_text(header + "".join(",".join(r) + "\n" for r in rows))


def test_parse_nodes_climbs_and_stops_at_root(tmp_path):
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text(NODES)
    parents, ranks = gf.parse_nodes_for_taxids(nodes, {"1220926"})
    assert parents["1220926"] == "4829"
    assert ranks["4829"] == "species"
    assert gf.species_node_for("1220926", parents, ranks) == "4829"
    assert gf.species_node_for("4829", parents, ranks) == "4829"


def test_parse_nodes_missing_taxid_raises(tmp_path):
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text(NODES)
    with pytest.raises(SystemExit, match="not found"):
        gf.parse_nodes_for_taxids(nodes, {"999999"})


def test_read_names_skips_authority(tmp_path):
    names = tmp_path / "names.dmp"
    names.write_text(NAMES)
    got = gf.read_names_for(names, {"4829"})["4829"]
    assert got["scientific name"] == ["Absidia glauca"]
    assert "Absidia glauca Hagem" not in got["scientific name"]
    assert got["synonym"] == ["Tieghemella glauca"]


@pytest.mark.parametrize(
    "name,tier,usable",
    [
        ("Absidia glauca", "canonical", True),
        ("Blumeria graminis f. sp. hordei", "ncbi_synonym", True),
        ("brewer's yeast", "genbank_common", True),
        ("Absidia glauca Hagem, 1908", "ncbi_synonym", False),
        ("Nannizzia gypsea (Nannizzi) Stockdale", "ncbi_synonym", False),
        ("Pseudozyma aphidis DSM 70725", "ncbi_synonym", False),
        ("Saccharomyces boulardii (nom. inval.)", "ncbi_synonym", False),
        ("Absidia sp.", "ncbi_synonym", False),
        ("Absidia", "ncbi_synonym", False),
        ("Aureobasidium pullulans (de Bary & Lowenthal) G. Arnaud, 1918",
         "ncbi_synonym", False),
    ],
)
def test_is_usable_name(name, tier, usable):
    assert gf.is_usable_name(name, tier) is usable


def test_build_queries_uses_correct_fields():
    qs = gf.build_queries("Mucor circinelloides")
    assert len(qs) == 7
    by = {(q["db"], q["mode"]): q["term"] for q in qs}
    assert set(by) == {
        ("pubmed", "name"), ("pubmed", "growth"), ("pubmed", "dimorph"), ("pubmed", "review"),
        ("pmc", "name"), ("pmc", "growth"), ("pmc", "dimorph"),
    }
    assert "review[pt]" in by[("pubmed", "review")] and '"dimorphism"[Title/Abstract]' in by[("pubmed", "review")]
    assert by[("pubmed", "name")] == '"Mucor circinelloides"[Title/Abstract]'
    assert by[("pmc", "name")] == '"Mucor circinelloides"[Text Word]'
    assert "[All Fields]" in by[("pubmed", "growth")]
    assert "[Text Word]" in by[("pmc", "growth")]
    assert '"dimorphism"[Title/Abstract]' in by[("pubmed", "dimorph")]
    assert '"yeast-to-hypha"[Text Word]' in by[("pmc", "dimorph")]


def test_query_priority_puts_dimorphism_first():
    assert gf.QUERY_PRIORITY["dimorph"] < gf.QUERY_PRIORITY["growth"] < gf.QUERY_PRIORITY["name"]
    assert gf.QUERY_PRIORITY["review"] == gf.QUERY_PRIORITY["dimorph"]


def _rec(db, uid, priority, rank, n_queries=1, pmid="", pmcid=""):
    return {"db": db, "uid": uid, "priority": priority, "rank": rank, "n_queries": n_queries,
            "pmid": pmid, "pmcid": pmcid}


def test_merge_article_records_joins_pubmed_and_pmc():
    recs = [_rec("pubmed", "111", 0, 2, 2, pmid="111", pmcid="999"),
            _rec("pmc", "999", 1, 48, 1, pmcid="999")]
    merged = gf.merge_article_records(recs)
    assert len(merged) == 1
    m = merged[0]
    assert (m["db"], m["uid"]) == ("pmc", "999")          # full text wins
    assert m["priority"] == 0 and m["rank"] == 2           # best priority, its rank
    assert m["n_queries"] == 3 and m["pmid"] == "111"


def test_merge_normalizes_prefixed_pmcid():
    """The ID converter gives 'PMC999'; PMC rows use '999'. They are one article, and the merged
    uid is the bare number that the cache uses."""
    recs = [_rec("pubmed", "111", 0, 2, 2, pmid="111", pmcid="PMC999"),
            _rec("pmc", "999", 1, 48, 1, pmcid="999")]
    merged = gf.merge_article_records(recs)
    assert len(merged) == 1
    assert merged[0]["uid"] == "999" and merged[0]["pmcid"] == "999"
    assert merged[0]["n_queries"] == 3


def test_two_step_uid_is_cache_key_for_prefixed_pmcid():
    recs = [_rec("pubmed", "111", 0, 0, 1, pmid="111", pmcid="PMC999")]
    out = gf.select_two_step(recs, {}, ["Species name"], 5, 5, screen_articles=5)
    assert (out[0]["db"], out[0]["uid"]) == ("pmc", "999")


def test_score_selection_beats_rank_ties():
    """Many name queries return rank-0 rows; the paper that several queries agree on must win."""
    recs = [_rec("pubmed", str(100 + i), 0, 0, 1, pmid=str(100 + i)) for i in range(15)]
    recs.append(_rec("pubmed", "777", 0, 3, 4, pmid="777"))           # rank 3 but 4 queries agree
    by_rank = gf.select_for_species(recs, 12, 10, selection="rank")
    by_score = gf.select_for_species(recs, 12, 10, selection="score")
    assert "777" not in {r["uid"] for r in by_rank}                     # the old order drops it
    assert by_score[0]["uid"] == "777"


def test_selection_respects_caps_and_fulltext_limit():
    recs = [_rec("pmc", str(i), 0, i, 1, pmcid=str(i)) for i in range(30)]
    recs += [_rec("pubmed", str(1000 + i), 1, i, 1, pmid=str(1000 + i)) for i in range(30)]
    out = gf.select_for_species(recs, 12, 5, selection="score")
    assert len(out) == 12
    assert sum(1 for r in out if r["db"] == "pmc") == 5
    with pytest.raises(ValueError):
        gf.select_for_species(recs, 12, 5, selection="nope")


def test_forced_records_sort_first():
    recs = [_rec("pmc", "5", 0, 0, 3, pmcid="5"), _rec("pmc", "9", -1, -1, 1, pmcid="9")]
    assert gf.select_for_species(recs, 1, 1)[0]["uid"] == "9"


def test_gold_recall_rows_stages():
    gold = [
        {"species": "Sp a", "level": "species", "pmid": "1", "pmcid": "10", "phrase": "yeast-to-hypha transition"},
        {"species": "Sp a", "level": "species", "pmid": "2", "pmcid": "", "phrase": "never found"},
        {"species": "Sp a", "level": "species", "pmid": "3", "pmcid": "30", "phrase": "forced source"},
        {"species": "Sp a", "level": "species", "pmid": "4", "pmcid": "", "phrase": "x"},
    ]
    candidates = [{"species": "Sp a", "uid": "10"}, {"species": "Sp a", "uid": "2"}, {"species": "Sp a", "uid": "30"}]
    articles = [
        {"species": "Sp a", "pmid": "1", "pmcid": "10", "priority": "0"},
        {"species": "Sp a", "pmid": "2", "pmcid": "", "priority": "1"},
        {"species": "Sp a", "pmid": "3", "pmcid": "30", "priority": "-1"},
    ]
    sentences = {"sp_a": [
        {"pmid": "1", "pmcid": "10", "text": "Cells in the", "context_before": "", "context_after": "yeast-to-hypha transition."},
        {"pmid": "2", "pmcid": "", "text": "something else", "context_before": "", "context_after": ""},
    ]}
    rows = gf.gold_recall_rows(gold, candidates, articles, sentences)
    r = {x["pmid"]: x for x in rows}
    assert (r["1"]["retrieved"], r["1"]["selected"], r["1"]["sentence"], r["1"]["phrase"]) == ("yes",) * 4
    assert (r["2"]["selected"], r["2"]["sentence"], r["2"]["phrase"]) == ("yes", "yes", "no")
    assert r["3"]["selected"] == "forced_only"                       # forced inclusion does not count
    assert (r["4"]["retrieved"], r["4"]["selected"]) == ("no", "no")


def test_split_sentences_protects_abbreviations():
    parts = gf.split_sentences(
        "We used S. cerevisiae, e.g. strain X, and M. circinelloides. "
        "It grew as yeast. Fig. 1 shows this."
    )
    assert parts == [
        "We used S. cerevisiae, e.g. strain X, and M. circinelloides.",
        "It grew as yeast.",
        "Fig. 1 shows this.",
    ]


PMC = """<?xml version="1.0"?>
<pmc-articleset><article>
<front><article-meta>
<article-id pub-id-type="pmcid">PMC999</article-id>
<article-id pub-id-type="pmid">12345</article-id>
<article-id pub-id-type="doi">10.1/x</article-id>
<title-group><article-title>Test title</article-title></title-group>
<journal-title>Test journal</journal-title>
<pub-date><year>2026</year></pub-date>
<abstract><p>Abstract sentence about yeast.</p></abstract>
</article-meta></front>
<body>
<p>{body}</p>
<back><ref-list><ref><mixed-citation>M. circinelloides was dimorphic.</mixed-citation></ref></ref-list></back>
</body>
</article></pmc-articleset>
""".format(body="Body sentence about dimorphic growth. " * 20)


def test_parse_pmc_accepts_pmcid_and_excludes_refs():
    got = gf.parse_pmc_articles(PMC.encode())["999"]
    assert got["pmcid"] == "999"
    assert got["pmid"] == "12345"
    assert got["doi"] == "10.1/x"
    assert got["title"] == "Test title"
    assert got["has_fulltext"] is True
    body_text = " ".join(t for s, t in got["sections"] if s == "Body")
    assert "dimorphic growth" in body_text
    assert "was dimorphic" not in body_text


def test_split_pmc_elements_matches_numeric_uid():
    parts = gf.split_pmc_elements(PMC.encode())
    assert set(parts) == {"999"}
    assert gf.parse_pmc_articles(parts["999"])["999"]["pmcid"] == "999"


def test_parse_pmc_single_article_root():
    root = PMC.encode()
    article = root.split(b"<pmc-articleset>", 1)[1].rsplit(b"</pmc-articleset>", 1)[0]
    got = gf.parse_pmc_articles(article)["999"]
    assert got["pmcid"] == "999"


PUBMED = """<PubmedArticleSet><PubmedArticle>
<MedlineCitation><PMID>12345</PMID>
<Article><ArticleTitle>Pub title</ArticleTitle>
<Abstract><AbstractText Label="BACKGROUND">Abstract text.</AbstractText></Abstract>
<Journal><Title>Pub journal</Title>
<JournalIssue><PubDate><Year>2020</Year></PubDate></JournalIssue></Journal>
<ArticleIdList><ArticleId IdType="doi">10.2/y</ArticleId></ArticleIdList>
</Article></MedlineCitation>
</PubmedArticle></PubmedArticleSet>"""


def test_parse_pubmed_batch_and_single():
    got = gf.parse_pubmed_articles(PUBMED.encode())["12345"]
    assert got["pmid"] == "12345"
    assert got["title"] == "Pub title"
    assert got["journal"] == "Pub journal"
    assert got["year"] == "2020"
    assert got["doi"] == "10.2/y"
    assert got["sections"] == [("Abstract", "BACKGROUND: Abstract text.")]
    single = gf.split_pubmed_elements(PUBMED.encode())["12345"]
    assert gf.parse_pubmed_articles(single)["12345"]["title"] == "Pub title"


def test_quote_in_text_is_normalized():
    assert gf.quote_in_text(
        "Mucor circinelloides was DIMORPHIC.",
        "we found mucor circinelloides was dimorphic in culture",
    )
    assert not gf.quote_in_text(
        "Mucor circinelloides is filamentous",
        "Mucor circinelloides is dimorphic",
    )


def test_cmd_species_and_names(tmp_path):
    study = tmp_path / "study"
    _species_csv(
        study / "yc_a" / "species.csv",
        [("Strain one", "S1", "1220926", "IN", "Mucorales"),
         ("Species two", "S2", "5341", "OUT", "Agaricales")],
    )
    batches = study / "batches"
    batches.mkdir()
    (batches / "b.yaml").write_text(yaml.safe_dump({
        "batch": "b",
        "species": {
            "Strain one": {"growth_form": "dimorphic", "dominant_phase": "yeast"},
            "Species two": {"growth_form": "filamentous"},
        },
    }))
    taxonomy = _taxonomy(tmp_path)
    out = study / "growth_form"
    gf.cmd_species(argparse.Namespace(
        study_dir=str(study), out=str(out), taxonomy_dir=str(taxonomy)))
    species = {r["species"]: r for r in gf.read_tsv(out / "species.tsv")}
    assert species["Strain one"]["species_taxon_id"] == "4829"
    assert species["Strain one"]["taxon_rank"] == "no rank"
    assert species["Species two"]["species_taxon_id"] == "5341"

    (out / "aliases.tsv").write_text(
        "species\talias\tsource\tnote\treview_status\n"
        "Strain one\tMucor lusitanicus\tNCBI Taxonomy 20260128 taxid 29924\tt\tseed\n"
    )
    gf.cmd_names(argparse.Namespace(
        study_dir=str(study), out=str(out), taxonomy_dir=str(taxonomy)))
    rows = gf.read_tsv(out / "names.tsv")
    one = {r["name"]: r for r in rows if r["species"] == "Strain one"}
    two = {r["name"]: r for r in rows if r["species"] == "Species two"}
    assert one["Absidia glauca"]["name_tier"] == "canonical"
    assert one["Tieghemella glauca"]["name_tier"] == "ncbi_synonym"
    assert one["Mucor lusitanicus"]["name_tier"] == "literature_alias"
    assert "Absidia glauca Hagem" not in one
    assert two["cultivated mushroom"]["name_tier"] == "genbank_common"


def test_cmd_species_rejects_batch_without_csv_row(tmp_path):
    study = tmp_path / "study"
    _species_csv(study / "yc_a" / "species.csv",
                 [("Strain one", "S1", "1220926", "IN", "Mucorales")])
    batches = study / "batches"
    batches.mkdir()
    (batches / "b.yaml").write_text(yaml.safe_dump({
        "batch": "b",
        "species": {"Strain one": {"growth_form": "dimorphic"},
                    "Ghost": {"growth_form": "yeast"}},
    }))
    taxonomy = _taxonomy(tmp_path)
    with pytest.raises(SystemExit, match="without species.csv rows"):
        gf.cmd_species(argparse.Namespace(
            study_dir=str(study), out=str(study / "growth_form"),
            taxonomy_dir=str(taxonomy)))


def test_screen_score_rewards_dimorphism_and_name_in_title():
    names = ["Yarrowia lipolytica"]
    dim = gf.screen_score(
        "TORC1 represses yeast-to-hypha transition in Yarrowia lipolytica",
        "The yeast-to-hypha dimorphic transition is important. Hyphae form in glucose media.", names)
    plain = gf.screen_score("Lipid production by Yarrowia lipolytica", "We improved lipid yield.", names)
    other = gf.screen_score("A dimorphic fungus", "Dimorphism and hyphae in another species.", names)
    assert plain == 5.0                                          # name in title, no growth terms
    assert dim > other > plain
    assert gf.screen_score("", "", names) == 0.0
    assert gf.screen_score("T", "x", ["abc"]) == 0.0            # names shorter than 6 characters are ignored


def test_two_step_selection_finds_paper_beyond_the_cap():
    """A dimorphism paper that sits at rank 80 of a saturated species is found by the screen."""
    recs = [_rec("pubmed", str(100 + i), 0, i, 1, pmid=str(100 + i)) for i in range(100)]
    screen = {str(100 + i): ("Lipid yield of Yarrowia lipolytica", "Biotechnology.") for i in range(100)}
    screen["180"] = ("Yeast-to-hypha dimorphic transition in Yarrowia lipolytica", "Hyphae and yeast forms.")
    plain = gf.select_for_species(recs, 12, 10, selection="score")
    two = gf.select_two_step(recs, screen, ["Yarrowia lipolytica"], 12, 10, screen_articles=100)
    assert "180" not in {r["uid"] for r in plain}
    assert "180" in {r["uid"] for r in two}                 # found by the screen half
    assert {r["uid"] for r in two} >= {"100", "101", "102", "103", "104", "105"}   # base half kept
    assert len(two) == 12


def test_two_step_never_ranks_a_screened_out_paper_above_a_screened_in_one():
    recs = [_rec("pubmed", str(i), 0, i, 1, pmid=str(i)) for i in range(10)]
    screen = {"9": ("Dimorphic Yarrowia lipolytica", "yeast hyphae dimorphism")}      # best text, but rank 9
    out = gf.select_two_step(recs, screen, ["Yarrowia lipolytica"], 3, 3, screen_articles=5)
    assert "9" not in {r["uid"] for r in out}                                          # outside the pool of 5


def test_two_step_keeps_forced_records_and_handles_no_abstract():
    recs = [_rec("pmc", "5", 0, 0, 1, pmcid="5"), _rec("pmc", "9", -1, -1, 1, pmcid="9")]
    out = gf.select_two_step(recs, {}, ["Species name"], 2, 2, screen_articles=1)
    assert [r["uid"] for r in out][0] == "9"
    assert {r["uid"] for r in out} == {"5", "9"}


# ---------- article-level sentences and parenthetical genus names ----------
def test_parenthetical_genus_pattern_matches_wangiella_exophiala():
    names = [("Exophiala dermatitidis", "canonical"), ("Wangiella dermatitidis", "ncbi_synonym"),
             ("Hormiscium dermatitidis", "ncbi_synonym")]
    pats = gf.parenthetical_genus_patterns(names)
    assert pats, "same epithet under several genera should give a pattern"
    label, tier, pat = pats[0]
    assert label == "Exophiala dermatitidis" and tier == "canonical"
    assert pat.search("The pathogen Wangiella (Exophiala) dermatitidis grows as yeast.")
    assert pat.search("Exophiala (Wangiella) dermatitidis is a black yeast")
    assert not pat.search("Wangiella dermatitidis")          # plain forms are matched elsewhere
    assert not pat.search("Wangiella (Exophiala) spinifera")


def test_parenthetical_genus_pattern_needs_two_genera():
    assert gf.parenthetical_genus_patterns([("Mucor circinelloides", "canonical")]) == []


def test_is_about_species():
    pats = [("Taphrina deformans", "canonical", gf.name_pattern("Taphrina deformans"))]
    assert gf.is_about_species("Early events in the Taphrina deformans-Prunus interaction", "x", pats, "abstract")
    assert gf.is_about_species("Peach leaf curl", "Caused by Taphrina deformans.", pats, "abstract")
    assert not gf.is_about_species("Peach leaf curl", "Caused by Taphrina deformans.", pats, "title")
    assert not gf.is_about_species("A review of fungi", "Many fungi.", pats, "abstract")
    assert not gf.is_about_species("Taphrina deformans", "x", pats, "off")


def test_pick_article_level_prefers_dimorphism_sentences_and_caps():
    cands = [(1, 0, "a"), (0, 5, "b"), (1, 2, "c"), (0, 9, "d"), (1, 3, "e")]     # (priority, position, payload)
    assert [c[2] for c in gf.pick_article_level(cands, 3)] == ["b", "d", "a"]
    assert gf.pick_article_level(cands, 0) == []


def test_gold_recall_blank_phrase_is_not_applicable():
    gold = [{"species": "Sp a", "level": "species", "pmid": "1", "pmcid": "10", "phrase": ""}]
    candidates = [{"species": "Sp a", "uid": "10"}]
    articles = [{"species": "Sp a", "pmid": "1", "pmcid": "10", "priority": "0"}]
    r = gf.gold_recall_rows(gold, candidates, articles, {})[0]
    assert (r["retrieved"], r["selected"], r["sentence"], r["phrase"]) == ("yes", "yes", "n/a", "n/a")


def test_name_reasons_apply_only_to_species_level_sentences():
    full_ok = {"match_level": "full", "matched_name": "Taphrina deformans", "name_used_in_source": "yes"}
    assert gf.name_reasons(full_ok, "the dimorphic fungus Taphrina deformans") == []
    assert gf.name_reasons(full_ok, "the dimorphic fungus") == ["name_not_in_quote"]
    no_name = {"match_level": "full", "matched_name": "Taphrina deformans", "name_used_in_source": "no"}
    assert "name_not_verbatim_in_sentence" in gf.name_reasons(no_name, "Taphrina deformans is dimorphic")
    for level in ("article", "genus"):
        rec = {"match_level": level, "matched_name": "Taphrina deformans", "name_used_in_source": "no"}
        assert gf.name_reasons(rec, "the fungus switches between yeast and hyphae") == []


def test_read_evidence_accepts_jsonl(tmp_path):
    f = tmp_path / "e.jsonl"
    f.write_text('{"species": "Sp a", "sentence_id": "x1", "quote": "q", "n": 3, "none": null}\n\n'
                 '{"species": "Sp a", "sentence_id": "x2", "quote": "r"}\n')
    rows = gf.read_evidence(f)
    assert [r["sentence_id"] for r in rows] == ["x1", "x2"]
    assert rows[0]["n"] == "3" and rows[0]["none"] == ""


def test_name_reasons_accept_the_abbreviation_for_abbrev_sentences():
    rec = {"match_level": "abbrev", "matched_name": "Mucor circinelloides", "name_used_in_source": "yes"}
    assert gf.name_reasons(rec, "M. circinelloides is a dimorphic fungus") == []
    assert gf.name_reasons(rec, "Mucor circinelloides is a dimorphic fungus") == []
    assert gf.name_reasons(rec, "is a dimorphic fungus") == ["name_not_in_quote"]
    full = {"match_level": "full", "matched_name": "Mucor circinelloides", "name_used_in_source": "yes"}
    assert gf.name_reasons(full, "M. circinelloides is a dimorphic fungus") == ["name_not_in_quote"]

