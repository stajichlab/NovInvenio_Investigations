import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from starbase_crossvalidation import (
    parse_tblastn_spans,
    load_starbase_af293_starships,
    match_captains_to_starships,
)


def test_parse_tblastn_spans_collapses_hsps_per_query_contig():
    lines = [
        "q1\tctg1\t90\t100\t1e-50\t200\t1\t100\t500\t600\t150",
        "q1\tctg1\t80\t50\t1e-10\t80\t101\t150\t650\t700\t150",
        "q2\tctg1\t99\t80\t1e-40\t150\t1\t80\t1000\t1080\t80",
    ]
    spans = parse_tblastn_spans(lines)
    assert spans["q1"] == {"contig": "ctg1", "start": 500, "end": 700}
    assert spans["q2"] == {"contig": "ctg1", "start": 1000, "end": 1080}


def test_parse_tblastn_spans_handles_reverse_strand_coordinates():
    # sstart > send on the minus strand -- span must still be (min, max).
    lines = ["q1\tctg1\t90\t100\t1e-50\t200\t1\t100\t700\t600\t150"]
    spans = parse_tblastn_spans(lines)
    assert spans["q1"] == {"contig": "ctg1", "start": 600, "end": 700}


def test_parse_tblastn_spans_keeps_only_best_contig_per_query():
    # A query with HSPs on two different contigs: keep only the
    # contig containing the lowest-evalue hit, not a bogus merged span.
    lines = [
        "q1\tctgA\t99\t100\t1e-60\t200\t1\t100\t500\t600\t150",
        "q1\tctgB\t40\t50\t1e-05\t60\t1\t50\t9000\t9050\t150",
    ]
    spans = parse_tblastn_spans(lines)
    assert spans["q1"]["contig"] == "ctgA"


def test_parse_tblastn_spans_does_not_merge_distant_clusters_on_same_contig():
    # Real case found 2026-09-17: a query with a tight, high-identity HSP
    # cluster at one locus and a second, low-identity cluster ~445kb away
    # on the SAME contig -- two different loci, not one 445kb-wide span.
    lines = [
        "q1\tctg1\t85\t300\t1e-100\t500\t1\t300\t639573\t639873\t601",
        "q1\tctg1\t100\t173\t1e-60\t300\t429\t601\t641567\t641740\t601",
        "q1\tctg1\t42\t404\t1e-05\t100\t221\t574\t1087403\t1087807\t601",
    ]
    spans = parse_tblastn_spans(lines, max_gap=20_000)
    assert spans["q1"] == {"contig": "ctg1", "start": 639573, "end": 641740}


def test_load_starbase_af293_starships_dedupes_and_picks_widest_span(tmp_path):
    import sqlite3

    db_path = tmp_path / "starbase.sqlite"
    conn = sqlite3.connect(db_path)
    conn.executescript("""
        CREATE TABLE taxonomy (id INTEGER PRIMARY KEY, strain TEXT, taxID TEXT);
        CREATE TABLE joined_ships (id INTEGER PRIMARY KEY, tax_id INTEGER, ship_id INTEGER, ship_family_id INTEGER);
        CREATE TABLE family_names (id INTEGER PRIMARY KEY, familyName TEXT);
        CREATE TABLE starship_features (
            id INTEGER PRIMARY KEY, ship_id INTEGER, contigID TEXT,
            starshipID TEXT, captainID TEXT, elementBegin TEXT, elementEnd TEXT, boundaryType TEXT
        );
        INSERT INTO taxonomy VALUES (1, 'Af293', '330879');
        INSERT INTO joined_ships VALUES (1, 1, 100, 5);
        INSERT INTO family_names VALUES (5, 'Galactica');
        INSERT INTO starship_features VALUES
            (1, 100, 'aspfum5_CM000171.1', 's01', 'cap1', '639106', '688498', ''),
            (2, 100, 'aspfum5_CM000171.1', 's01', 'cap1', '639099', '688505', 'flank');
    """)
    conn.commit()
    conn.close()

    starships = load_starbase_af293_starships(str(db_path), af293_taxonomy_ids=[1])
    assert len(starships) == 1
    s = starships[0]
    assert s["contig"] == "CM000171.1"
    assert s["start"] == 639099
    assert s["end"] == 688505
    assert s["family"] == "Galactica"


def test_load_starbase_af293_starships_normalizes_chr_label_contigs(tmp_path):
    import sqlite3

    db_path = tmp_path / "starbase.sqlite"
    conn = sqlite3.connect(db_path)
    conn.executescript("""
        CREATE TABLE taxonomy (id INTEGER PRIMARY KEY, strain TEXT, taxID TEXT);
        CREATE TABLE joined_ships (id INTEGER PRIMARY KEY, tax_id INTEGER, ship_id INTEGER, ship_family_id INTEGER);
        CREATE TABLE family_names (id INTEGER PRIMARY KEY, familyName TEXT);
        CREATE TABLE starship_features (
            id INTEGER PRIMARY KEY, ship_id INTEGER, contigID TEXT,
            starshipID TEXT, captainID TEXT, elementBegin TEXT, elementEnd TEXT, boundaryType TEXT
        );
        INSERT INTO taxonomy VALUES (1, 'Af293', '330879');
        INSERT INTO joined_ships VALUES (1, 1, 100, NULL);
        INSERT INTO starship_features VALUES
            (1, 100, 'aspfum3_Chr3AfumigatusAf293', 's02', 'cap2', '1012235', '1089747', '');
    """)
    conn.commit()
    conn.close()

    starships = load_starbase_af293_starships(str(db_path), af293_taxonomy_ids=[1])
    assert starships[0]["contig"] == "CM000171.1"  # Chr3 -> CM000171.1 (chr1=CM000169.1)


def test_match_captains_to_starships_flags_span_containment():
    captain_spans = {
        "q1": {"contig": "CM000169.1", "start": 4397258, "end": 4399952},
        "q2": {"contig": "CM000169.1", "start": 325121, "end": 328217},
    }
    starships = [
        {"contig": "CM000169.1", "start": 4343061, "end": 4400403, "starshipID": "s01", "family": "Hephaestus"},
    ]
    matches = match_captains_to_starships(captain_spans, starships)
    assert matches["q1"]["matched_starship"] == "s01"
    assert matches["q2"]["matched_starship"] is None


def test_match_captains_to_starships_reports_unmatched_starships():
    captain_spans = {"q1": {"contig": "ctg1", "start": 100, "end": 200}}
    starships = [
        {"contig": "ctg1", "start": 90, "end": 210, "starshipID": "s01", "family": "F1"},
        {"contig": "ctg2", "start": 5000, "end": 6000, "starshipID": "s02", "family": "F2"},
    ]
    matches = match_captains_to_starships(captain_spans, starships)
    unmatched = [s["starshipID"] for s in starships if s["starshipID"] not in
                 {m["matched_starship"] for m in matches.values() if m["matched_starship"]}]
    assert unmatched == ["s02"]
