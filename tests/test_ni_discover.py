import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "lib"))
sys.path.insert(0, str(REPO / "bin"))
import ni_discover as nd  # noqa: E402

FIXTURES = REPO / "tests" / "fixtures"


def _load_fixture_lines(name: str) -> str:
    return (FIXTURES / name).read_text()


def test_query_species_genomes_parses_live_fixture():
    fake_run = MagicMock(stdout=_load_fixture_lines("fox_annotated.jsonl"), returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        genomes = nd.query_species_genomes("5507")
    assert len(genomes) > 0
    # Every record must carry the fields later tasks depend on.
    for g in genomes:
        assert "accession" in g
        assert "organism_name" in g
        assert "record_taxon_id" in g
        assert "has_annotation" in g
        assert g["has_annotation"] is True  # --annotated query -> every record has annotation_info
    # Spot-check a known record from the live fixture (adjust if the live
    # capture no longer contains this exact accession -- the point is that
    # SOME record's fields are checked against known-real values, not that
    # this specific accession is eternal).
    fo47 = [g for g in genomes if g["strain"] == "Fo47"]
    assert fo47, "expected at least one Fo47 record in the live fixture"

    # Live-verified 2026-09-12: annotation_info.stats.gene_counts.protein_coding
    # and annotation_info.busco.complete are real field paths -- at least one
    # record in the fixture has both populated.
    with_gene_count = [g for g in genomes if g["protein_coding_count"] is not None]
    assert with_gene_count, "expected at least one record with a protein_coding_count"
    with_busco = [g for g in genomes if g["busco_score"] is not None]
    assert with_busco, "expected at least one record with a busco_score"
    for g in with_busco:
        assert isinstance(g["busco_score"], float)

    # Live-verified: RefSeq (GCF_) records carry annotation_info.pipeline
    # ("NCBI Eukaryotic Annotation Propagation Pipeline"); GenBank-submitted
    # records do not (annotation_pipeline is None for those, which is
    # expected, not a parsing bug).
    refseq = [g for g in genomes if g["accession"].startswith("GCF_")]
    assert refseq
    assert any(g["annotation_pipeline"] for g in refseq)


def test_query_species_genomes_complex_fixture_is_superset():
    fake_run = MagicMock(stdout=_load_fixture_lines("fox_complex_annotated.jsonl"), returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        genomes = nd.query_species_genomes("171631")
    assert len(genomes) >= len(_load_fixture_lines("fox_annotated.jsonl").strip().splitlines())
    odoratissimum = [g for g in genomes if "odoratissimum" in g["organism_name"]]
    assert odoratissimum, "expected Fusarium odoratissimum records in the species-complex fixture"


def test_fetch_taxonomy_ranks_batched_call_shape():
    # Real shape verified live 2026-09-12 via:
    #   datasets summary taxonomy taxon 1229664 61366
    # `datasets summary taxonomy taxon` returns ONE JSON object (not
    # JSON-lines) shaped {"reports": [...], "total_count": N}. Each report is
    # {"query": ["<id>"], "taxonomy": {...}}. The taxonomy object's `rank` and
    # `tax_id` are direct fields; the taxon's own name is at
    # `current_scientific_name.name` (NOT a top-level "name" field); the
    # species-level classification is at `classification.species.name`
    # (absent/None for a taxon at or above species rank, e.g. a
    # FORMA_SPECIALIS or SPECIES_GROUP node has no "species" sub-key distinct
    # from itself, or the key is simply missing for high ranks); and
    # `parents` is the FULL ancestor lineage as a flat list of ints ordered
    # root-first, with the immediate parent LAST -- not a single-element
    # direct-parent field as an early sketch assumed.
    real_response = {
        "reports": [
            {
                "query": ["1229664"],
                "taxonomy": {
                    "tax_id": 1229664,
                    "rank": "STRAIN",
                    "current_scientific_name": {"name": "Fusarium oxysporum f. sp. cubense race 1"},
                    "classification": {"species": {"id": 5507, "name": "Fusarium oxysporum"}},
                    "parents": [171631, 5507, 61366],
                },
            },
            {
                "query": ["61366"],
                "taxonomy": {
                    "tax_id": 61366,
                    "rank": "FORMA_SPECIALIS",
                    "current_scientific_name": {"name": "Fusarium oxysporum f. sp. cubense"},
                    "classification": {"species": {"id": 5507, "name": "Fusarium oxysporum"}},
                    "parents": [5506, 171631, 5507],
                },
            },
        ],
        "total_count": 2,
    }
    fake_run = MagicMock(stdout=json.dumps(real_response), returncode=0)
    with patch("subprocess.run", return_value=fake_run) as mock_run:
        result = nd.fetch_taxonomy_ranks(["1229664", "61366"])

    # One batched call, both ids in a single command.
    assert mock_run.call_count == 1
    called_cmd = mock_run.call_args[0][0]
    assert called_cmd[:4] == ["datasets", "summary", "taxonomy", "taxon"]
    assert "1229664" in called_cmd and "61366" in called_cmd

    assert result == {
        "1229664": {
            "rank": "STRAIN",
            "name": "Fusarium oxysporum f. sp. cubense race 1",
            "parents": ["171631", "5507", "61366"],
            "species_name": "Fusarium oxysporum",
        },
        "61366": {
            "rank": "FORMA_SPECIALIS",
            "name": "Fusarium oxysporum f. sp. cubense",
            "parents": ["5506", "171631", "5507"],
            "species_name": "Fusarium oxysporum",
        },
    }


def test_fetch_taxonomy_ranks_handles_missing_species_classification():
    # Live-verified: a SPECIES_GROUP-rank taxon (171631, "Fusarium oxysporum
    # species complex") has no "species" key under classification at all.
    real_response = {
        "reports": [
            {
                "query": ["171631"],
                "taxonomy": {
                    "tax_id": 171631,
                    "rank": "SPECIES_GROUP",
                    "current_scientific_name": {"name": "Fusarium oxysporum species complex"},
                    "classification": {"genus": {"id": 5506, "name": "Fusarium"}},
                    "parents": [5125, 110618, 5506],
                },
            }
        ],
        "total_count": 1,
    }
    fake_run = MagicMock(stdout=json.dumps(real_response), returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        result = nd.fetch_taxonomy_ranks(["171631"])
    assert result["171631"]["species_name"] == ""


def test_fetch_taxonomy_ranks_empty_input_makes_no_call():
    with patch("subprocess.run") as mock_run:
        result = nd.fetch_taxonomy_ranks([])
    assert result == {}
    mock_run.assert_not_called()


def test_find_forma_specialis_group_walks_rank_lookup():
    # cubense race 1 (taxid 1229664) -> parent 61366 "Fusarium oxysporum f.
    # sp. cubense" (rank FORMA_SPECIALIS) -- live-verified in the design
    # spec's review. Build rank_lookup with exactly the real shape confirmed
    # in fetch_taxonomy_ranks's own test above.
    rank_lookup = {
        "1229664": {"rank": "STRAIN", "name": "Fusarium oxysporum f. sp. cubense race 1", "parents": ["61366"], "species_name": "Fusarium oxysporum"},
        "61366": {"rank": "FORMA_SPECIALIS", "name": "Fusarium oxysporum f. sp. cubense", "parents": ["5507"], "species_name": "Fusarium oxysporum"},
    }
    label, used_fallback = nd.find_forma_specialis_group("1229664", rank_lookup, "Fusarium oxysporum f. sp. cubense race 1")
    assert label == "cubense"
    assert used_fallback is False


def test_find_forma_specialis_group_walks_realistic_multi_element_parents():
    # Regression test for the parents[-1] vs parents[0] walk direction.
    # Every other rank_lookup fixture in this file uses single-element (or
    # empty) `parents` lists, for which [-1] and [0] are identical -- this
    # test uses the REAL, full-length root-first `parents` arrays live-
    # captured for taxid 1229664 (Fusarium oxysporum f. sp. cubense race 1)
    # and its immediate parent 61366 (see task-1-report.md for the full
    # `datasets summary taxonomy taxon 1229664 61366` transcript this was
    # transcribed from). Swapping parents[-1] back to parents[0] in
    # find_forma_specialis_group would walk toward taxid "1" (the root) on
    # the very first hop instead of to 61366, and "1" is not in rank_lookup,
    # so the walk would break and fall through to the organism_name regex.
    # organism_name is deliberately given WITHOUT any "f. sp." text so that
    # fallback path can't coincidentally produce the right label -- a
    # parents[0] regression must show up as a wrong label/used_fallback, not
    # be masked by the fallback also finding "cubense".
    rank_lookup = {
        "1229664": {
            "rank": "STRAIN",
            "name": "Fusarium oxysporum f. sp. cubense race 1",
            "species_name": "Fusarium oxysporum",
            "parents": [
                "1", "131567", "2759", "33154", "4751", "451864", "4890",
                "716545", "147538", "716546", "715989", "147550", "222543",
                "5125", "110618", "5506", "171631", "5507", "61366",
            ],
        },
        "61366": {
            "rank": "FORMA_SPECIALIS",
            "name": "Fusarium oxysporum f. sp. cubense",
            "species_name": "Fusarium oxysporum",
            "parents": [
                "1", "131567", "2759", "33154", "4751", "451864", "4890",
                "716545", "147538", "716546", "715989", "147550", "222543",
                "5125", "110618", "5506", "171631", "5507",
            ],
        },
    }
    label, used_fallback = nd.find_forma_specialis_group(
        "1229664", rank_lookup, "Fusarium isolate FOC-race1"
    )
    assert label == "cubense"
    assert used_fallback is False


def test_find_forma_specialis_group_falls_back_to_regex():
    # A taxid with no FORMA_SPECIALIS ancestor in rank_lookup at all.
    rank_lookup = {
        "999999": {"rank": "STRAIN", "name": "Some organism f. sp. madeup", "parents": ["1"], "species_name": "Some organism"},
        "1": {"rank": "SPECIES", "name": "Some organism", "parents": [], "species_name": "Some organism"},
    }
    label, used_fallback = nd.find_forma_specialis_group("999999", rank_lookup, "Some organism f. sp. madeup")
    assert label == "madeup"
    assert used_fallback is True


def test_find_forma_specialis_group_no_fsp_anywhere():
    rank_lookup = {
        "660027": {"rank": "STRAIN", "name": "Fusarium oxysporum Fo47", "parents": ["5507"], "species_name": "Fusarium oxysporum"},
        "5507": {"rank": "SPECIES", "name": "Fusarium oxysporum", "parents": ["171631"], "species_name": "Fusarium oxysporum"},
    }
    label, used_fallback = nd.find_forma_specialis_group("660027", rank_lookup, "Fusarium oxysporum Fo47")
    assert label == "no-fsp-in-name"


def test_find_species_complex_detects_species_group_ancestor():
    rank_lookup = {
        "5507": {"rank": "SPECIES", "name": "Fusarium oxysporum", "parents": ["171631"], "species_name": "Fusarium oxysporum"},
        "171631": {"rank": "SPECIES_GROUP", "name": "Fusarium oxysporum species complex", "parents": [], "species_name": ""},
    }
    complex_info = nd.find_species_complex("5507", rank_lookup)
    assert complex_info == {"taxon_id": "171631", "name": "Fusarium oxysporum species complex"}


def test_find_species_complex_walks_realistic_multi_element_parents_multi_hop():
    # Regression test for the parents[-1] vs parents[0] walk direction,
    # requiring TWO hops (61366 -> 5507 -> 171631) before reaching the
    # SPECIES_GROUP ancestor, using the real full-length root-first
    # `parents` arrays live-captured for taxids 61366, 5507, and 171631
    # (see task-1-report.md). Swapping parents[-1] back to parents[0] would
    # walk toward taxid "1" on the first hop, which is not in rank_lookup,
    # so the walk would break and return None instead of finding 171631.
    rank_lookup = {
        "61366": {
            "rank": "FORMA_SPECIALIS",
            "name": "Fusarium oxysporum f. sp. cubense",
            "species_name": "Fusarium oxysporum",
            "parents": [
                "1", "131567", "2759", "33154", "4751", "451864", "4890",
                "716545", "147538", "716546", "715989", "147550", "222543",
                "5125", "110618", "5506", "171631", "5507",
            ],
        },
        "5507": {
            "rank": "SPECIES",
            "name": "Fusarium oxysporum",
            "species_name": "Fusarium oxysporum",
            "parents": [
                "1", "131567", "2759", "33154", "4751", "451864", "4890",
                "716545", "147538", "716546", "715989", "147550", "222543",
                "5125", "110618", "5506", "171631",
            ],
        },
        "171631": {
            "rank": "SPECIES_GROUP",
            "name": "Fusarium oxysporum species complex",
            "species_name": "",
            "parents": [
                "1", "131567", "2759", "33154", "4751", "451864", "4890",
                "716545", "147538", "716546", "715989", "147550", "222543",
                "5125", "110618", "5506",
            ],
        },
    }
    complex_info = nd.find_species_complex("61366", rank_lookup)
    assert complex_info == {"taxon_id": "171631", "name": "Fusarium oxysporum species complex"}


def test_find_species_complex_none_when_no_species_group_ancestor():
    rank_lookup = {
        "5507": {"rank": "SPECIES", "name": "Fusarium oxysporum", "parents": ["4751"], "species_name": "Fusarium oxysporum"},
        "4751": {"rank": "KINGDOM", "name": "Fungi", "parents": [], "species_name": ""},
    }
    assert nd.find_species_complex("5507", rank_lookup) is None


def test_derive_species_name_uses_taxonomy_not_string_surgery():
    genome = {"record_taxon_id": "660027", "organism_name": "Fusarium oxysporum Fo47"}
    rank_lookup = {"660027": {"rank": "STRAIN", "name": "Fusarium oxysporum Fo47", "parents": ["5507"], "species_name": "Fusarium oxysporum"}}
    assert nd.derive_species_name(genome, rank_lookup) == "Fusarium oxysporum"


def test_pick_short_prefers_isolate_over_race_token_strain():
    genome = {"strain": "TR4", "isolate": "UK0001"}
    short = nd.pick_short(genome, "Foxy", set())
    assert "TR4" not in short
    assert "UK0001" in short


def test_pick_short_uses_strain_when_not_a_race_token():
    genome = {"strain": "Fo47", "isolate": None}
    short = nd.pick_short(genome, "Foxy", set())
    assert "Fo47" in short


def test_pick_short_deduplicates_on_collision():
    used = set()
    first = nd.pick_short({"strain": "X1", "isolate": None}, "Foxy", used)
    used.add(first)
    second = nd.pick_short({"strain": "X1", "isolate": None}, "Foxy", used)
    assert first != second


def test_pick_short_never_bare_numeric():
    short = nd.pick_short({"strain": "9", "isolate": None}, "Foxy", set())
    assert not short.isdigit()


def test_detect_duplicate_groups_unions_biosample_and_strain_match():
    # Fo5176 registered 3 times with 3 DIFFERENT biosample accessions
    # (live-verified in the design spec's review) -- must still collapse via
    # the strain-string match arm of the union, not the biosample arm.
    genomes = [
        {"accession": "GCA_1", "biosample_accession": "SAMN_A", "strain": "Fo5176", "isolate": None},
        {"accession": "GCA_2", "biosample_accession": "SAMN_B", "strain": "Fo5176", "isolate": None},
        {"accession": "GCA_3", "biosample_accession": "SAMN_C", "strain": "Fo5176", "isolate": None},
        {"accession": "GCA_4", "biosample_accession": "SAMN_D", "strain": "SomethingElse", "isolate": None},
    ]
    groups = nd.detect_duplicate_groups(genomes)
    sizes = sorted(len(g) for g in groups)
    assert sizes == [1, 3]


def test_detect_duplicate_groups_unions_via_biosample_when_strain_differs():
    genomes = [
        {"accession": "GCA_1", "biosample_accession": "SAMN_SHARED", "strain": "NameA", "isolate": None},
        {"accession": "GCA_2", "biosample_accession": "SAMN_SHARED", "strain": "NameB", "isolate": None},
    ]
    groups = nd.detect_duplicate_groups(genomes)
    assert len(groups) == 1
    assert len(groups[0]) == 2


def test_build_report_includes_group_counts():
    genomes = [
        {"accession": "GCA_1", "record_taxon_id": "1", "organism_name": "X f. sp. a", "strain": "s1", "isolate": None, "assembly_level": "Scaffold", "annotation_pipeline": "p1", "protein_coding_count": 100, "busco_score": 90.0},
        {"accession": "GCA_2", "record_taxon_id": "2", "organism_name": "X f. sp. a", "strain": "s2", "isolate": None, "assembly_level": "Scaffold", "annotation_pipeline": "p1", "protein_coding_count": 200, "busco_score": 95.0},
        {"accession": "GCA_3", "record_taxon_id": "3", "organism_name": "X", "strain": "s3", "isolate": None, "assembly_level": "Chromosome", "annotation_pipeline": "p2", "protein_coding_count": 150, "busco_score": None},
    ]
    rank_lookup = {}  # forma_group precomputed and stashed on each genome dict for this test
    for g, grp in zip(genomes, ["a", "a", "no-fsp-in-name"]):
        g["_forma_group"] = grp
    report = nd.build_report(genomes, rank_lookup, duplicate_groups=[[g] for g in genomes], species_complex=None, complex_genome_count=0)
    assert "a" in report
    assert "2" in report  # count for group "a"
    assert "no-fsp-in-name" in report


def test_build_report_notes_duplicate_groups():
    genomes = [
        {"accession": "GCA_1", "strain": "Fo5176", "isolate": None, "_forma_group": "no-fsp-in-name"},
        {"accession": "GCA_2", "strain": "Fo5176", "isolate": None, "_forma_group": "no-fsp-in-name"},
    ]
    report = nd.build_report(genomes, {}, duplicate_groups=[genomes], species_complex=None, complex_genome_count=0)
    assert "GCA_1" in report and "GCA_2" in report
    assert "duplicate" in report.lower() or "merged" in report.lower()


def test_build_report_notes_species_complex():
    report = nd.build_report(
        [], {}, duplicate_groups=[],
        species_complex={"taxon_id": "171631", "name": "Fusarium oxysporum species complex"},
        complex_genome_count=4,
    )
    assert "Fusarium oxysporum species complex" in report
    assert "4" in report
