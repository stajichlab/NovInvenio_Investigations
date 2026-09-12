import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "lib"))
sys.path.insert(0, str(REPO / "bin"))
import ni_resolve as nr  # noqa: E402


def _mock_response(payload, link_header=None):
    """A context-manager mock matching urllib.request.urlopen's return value.

    NOTE: the design brief's original version of this helper assigned
    `resp.headers.get = lambda ...` onto a plain dict, which raises
    "AttributeError: 'dict' object attribute 'get' is read-only" (dict.get is
    a bound method on the type, not an instance attribute you can overwrite).
    Fixed here by giving `headers` its own MagicMock with `.get` wired via
    side_effect, which behaves like urllib's real (case-insensitive)
    HTTPMessage headers object closely enough for this test's purposes.
    """
    resp = MagicMock()
    resp.__enter__.return_value = resp
    resp.read.return_value = json.dumps(payload).encode("utf-8")
    resp.headers = MagicMock()
    resp.headers.get = MagicMock(
        side_effect=lambda k, default=None: (link_header if k == "Link" else default)
    )
    return resp


def test_resolve_taxon_id_single_match():
    payload = {"results": [{"taxonId": 5334, "scientificName": "Schizophyllum commune"}]}
    with patch("urllib.request.urlopen", return_value=_mock_response(payload)):
        hits = nr.resolve_taxon_id("Schizophyllum commune")
    assert hits == [{"taxon_id": "5334", "scientific_name": "Schizophyllum commune"}]


def test_resolve_taxon_id_no_match():
    payload = {"results": []}
    with patch("urllib.request.urlopen", return_value=_mock_response(payload)):
        hits = nr.resolve_taxon_id("Xylaria sp.")
    assert hits == []


def test_search_uniprot_proteomes_excludes_mycovirus_via_taxon_scoping():
    # Taxon-ID-scoped query should never see the mycovirus in the first place --
    # this test asserts the QUERY URL built uses taxonomy_id, not organism_name,
    # which is the actual fix (the mycovirus contamination happens at the API
    # level for organism_name queries; scoping by taxon_id is what avoids it).
    #
    # NOTE on shape (verified live against rest.uniprot.org 2026-09-11): UniProt
    # puts "strain" as a TOP-LEVEL field on the proteome record, not nested under
    # "taxonomy" as the design sketch assumed -- confirmed via
    # `curl 'https://rest.uniprot.org/proteomes/search?query=taxonomy_id:5334+AND+proteome_type:REFERENCE&format=json'`,
    # which for UP000001805 (Neurospora crassa reference proteome) returns
    # "strain": "ATCC 24698 / 74-OR23-1A / CBS 708.71 / DSM 1257 / FGSC 987" at
    # the top level, alongside (not inside) "taxonomy": {"scientificName": ...,
    # "taxonId": ..., "mnemonic": ...}. superkingdom is also lowercase
    # ("eukaryota") in real responses, not "Eukaryota" -- harmless here since
    # check_ftp_available/DOMAIN_FOLDERS already normalizes via .capitalize().
    payload = {
        "results": [
            {
                "id": "UP000001805",
                "taxonomy": {"scientificName": "Neurospora crassa", "taxonId": 5334, "mnemonic": "NEUCR"},
                "strain": "ATCC 24698 / 74-OR23-1A / CBS 708.71 / DSM 1257 / FGSC 987",
                "superkingdom": "eukaryota",
                "genomeAssembly": {"assemblyId": "GCA_000182925.2"},
                "proteomeCompletenessReport": {"buscoReport": {"score": 98.5}},
            }
        ]
    }
    captured_url = {}

    def fake_urlopen(req, *a, **kw):
        url = req.full_url if hasattr(req, "full_url") else req
        captured_url["url"] = url
        return _mock_response(payload)

    with patch("urllib.request.urlopen", side_effect=fake_urlopen):
        candidates = nr.search_uniprot_proteomes("5334", reference_only=True)

    assert "taxonomy_id:5334" in captured_url["url"]
    assert "organism_name" not in captured_url["url"]
    assert len(candidates) == 1
    assert candidates[0]["proteome_id"] == "UP000001805"
    assert candidates[0]["genome_assembly_id"] == "GCA_000182925.2"
    assert candidates[0]["busco_score"] == 98.5
    assert candidates[0]["strain"] == "ATCC 24698 / 74-OR23-1A / CBS 708.71 / DSM 1257 / FGSC 987"


def test_search_uniprot_proteomes_paginates():
    page1 = {"results": [{"id": "UP1", "strain": None, "superkingdom": "eukaryota"}]}
    page2 = {"results": [{"id": "UP2", "strain": None, "superkingdom": "eukaryota"}]}
    responses = [
        _mock_response(page1, link_header='<https://rest.uniprot.org/proteomes/search?cursor=abc>; rel="next"'),
        _mock_response(page2),
    ]

    def fake_urlopen(req, *a, **kw):
        return responses.pop(0)

    with patch("urllib.request.urlopen", side_effect=fake_urlopen):
        candidates = nr.search_uniprot_proteomes("999", reference_only=False)

    assert [c["proteome_id"] for c in candidates] == ["UP1", "UP2"]


def test_check_ftp_available_true_on_200():
    resp = MagicMock()
    resp.__enter__.return_value = resp
    resp.status = 200
    with patch("urllib.request.urlopen", return_value=resp):
        assert nr.check_ftp_available("UP000001805", "eukaryota") is True


def test_check_ftp_available_false_on_404():
    import urllib.error
    with patch("urllib.request.urlopen", side_effect=urllib.error.HTTPError(
        "url", 404, "not found", {}, None
    )):
        assert nr.check_ftp_available("UP001497681", "eukaryota") is False


def test_strain_matches_substring_case_insensitive():
    assert nr.strain_matches("H4-8 / ATCC 38548", "h4-8") is True
    assert nr.strain_matches("Tattone D", "H4-8") is False
    assert nr.strain_matches(None, "H4-8") is False
    assert nr.strain_matches("", "H4-8") is False
