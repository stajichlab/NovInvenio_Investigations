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


# NCBI Datasets fixtures below mirror the REAL live shape observed 2026-09-11
# via `pixi run datasets summary genome taxon 5334 --as-json-lines` against a
# GCA/GCF reference-genome pair (Schizophyllum commune, GCA_000143185.2 /
# GCF_000143185.2). Two field paths differ from this task's brief sketch:
#
# - `paired_accession` is a TOP-LEVEL field on the record (`r["paired_accession"]`),
#   not nested under `assembly_info` (`assembly_info.paired_accession` does not
#   exist; there's a differently-shaped `assembly_info.paired_assembly.accession`
#   object instead, which we don't use).
# - There is no `assembly_info.submission_date` field at all. The closest
#   analog present on every record is `assembly_info.release_date` (the date the
#   assembly was released to the public database); we use that for the
#   "submission_date" key in our returned candidate dicts. (A `submission_date`
#   does exist, but only nested three levels down under
#   `assembly_info.biosample.submission_date`, which describes the biological
#   sample, not the assembly record -- not what this field is meant to capture.)
NCRASSA_PAIR_JSONL = "\n".join([
    json.dumps({
        "accession": "GCA_000182925.2",
        "paired_accession": "GCF_000182925.2",
        "assembly_info": {
            "refseq_category": "reference genome",
            "assembly_level": "Chromosome",
            "assembly_status": "current",
            "submitter": "Broad Institute",
            "release_date": "2014-05-01",
        },
        "organism": {"infraspecific_names": {"strain": "OR74A"}},
        "annotation_info": {"name": "GCF_000182925.2-RS_2014_05"},
    }),
    json.dumps({
        "accession": "GCF_000182925.2",
        "paired_accession": "GCA_000182925.2",
        "assembly_info": {
            "refseq_category": "reference genome",
            "assembly_level": "Chromosome",
            "assembly_status": "current",
            "submitter": "Broad Institute",
            "release_date": "2014-05-01",
        },
        "organism": {"infraspecific_names": {"strain": "OR74A"}},
        "annotation_info": {"name": "GCF_000182925.2-RS_2014_05"},
    }),
])


def test_query_ncbi_assemblies_parses_jsonl():
    fake_run = MagicMock(stdout=NCRASSA_PAIR_JSONL, returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        candidates = nr.query_ncbi_assemblies("5334", require_reference=False)
    assert len(candidates) == 2
    accessions = {c["accession"] for c in candidates}
    assert accessions == {"GCA_000182925.2", "GCF_000182925.2"}


def test_rank_ncbi_candidates_collapses_gca_gcf_pair_to_one():
    fake_run = MagicMock(stdout=NCRASSA_PAIR_JSONL, returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        candidates = nr.query_ncbi_assemblies("5334", require_reference=False)
    top = nr.rank_ncbi_candidates(candidates, require_annotation=False)
    assert len(top) == 1  # NOT 2 -- this is the pair-collapse fix
    assert top[0]["accession"] == "GCA_000182925.2"  # GCA preferred for recording


def test_rank_ncbi_candidates_requires_annotation_when_asked():
    unannotated = json.dumps({
        "accession": "GCA_999999999.1",
        "paired_accession": None,
        "assembly_info": {
            "refseq_category": None,
            "assembly_level": "Complete Genome",
            "assembly_status": "current",
            "submitter": "Some Lab",
            "release_date": "2025-01-01",
        },
        "organism": {"infraspecific_names": {}},
    })  # no annotation_info key at all
    fake_run = MagicMock(stdout=unannotated, returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        candidates = nr.query_ncbi_assemblies("1", require_reference=False)
    assert nr.rank_ncbi_candidates(candidates, require_annotation=True) == []
    assert len(nr.rank_ncbi_candidates(candidates, require_annotation=False)) == 1


def test_rank_ncbi_candidates_ties_when_no_reference_designation():
    two_unranked = "\n".join([
        json.dumps({
            "accession": f"GCA_{i}00000000.1",
            "paired_accession": None,
            "assembly_info": {
                "refseq_category": None, "assembly_level": "Chromosome",
                "assembly_status": "current",
                "submitter": f"Lab {i}", "release_date": "2025-01-01",
            },
            "organism": {"infraspecific_names": {"strain": f"strain{i}"}},
            "annotation_info": {"name": "x"},
        }) for i in (1, 2)
    ])
    fake_run = MagicMock(stdout=two_unranked, returncode=0)
    with patch("subprocess.run", return_value=fake_run):
        candidates = nr.query_ncbi_assemblies("1", require_reference=False)
    tied = nr.rank_ncbi_candidates(candidates, require_annotation=False)
    assert len(tied) == 2  # genuine tie, not collapsed -- different accessions, no ranking signal


def test_resolve_row_uniprot_reference_match(monkeypatch):
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [{"taxon_id": "5334", "scientific_name": s}])
    monkeypatch.setattr(nr, "search_uniprot_proteomes", lambda tid, reference_only: (
        [{"proteome_id": "UP000001805", "strain": "OR74A / FGSC 987", "busco_score": 98.5,
          "genome_assembly_id": "GCA_000182925.2", "superkingdom": "Eukaryota"}]
        if reference_only else []
    ))
    monkeypatch.setattr(nr, "check_ftp_available", lambda pid, sk: True)
    monkeypatch.setattr(nr, "query_ncbi_assemblies", lambda tid, require_reference: [])

    row = nr.resolve_row("Ncra", "Neurospora crassa", "OR74A")

    assert row.resolved is True
    assert row.protein_source == "uniprot"
    assert row.protein_accession == "UP000001805"
    assert row.taxon_id == "5334"
    assert row.genome_source == "ncbi"
    assert row.genome_accession == "GCA_000182925.2"


def test_resolve_row_uniprot_ftp_not_yet_available_falls_back_to_ncbi(monkeypatch):
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [{"taxon_id": "1", "scientific_name": s}])
    monkeypatch.setattr(nr, "search_uniprot_proteomes", lambda tid, reference_only: (
        [{"proteome_id": "UP001497681", "strain": "Tattone D", "busco_score": 97.0,
          "genome_assembly_id": "GCA_023508785.1", "superkingdom": "Eukaryota"}]
        if reference_only else []
    ))
    monkeypatch.setattr(nr, "check_ftp_available", lambda pid, sk: False)  # not in FTP snapshot yet
    monkeypatch.setattr(nr, "query_ncbi_assemblies", lambda tid, require_reference: [
        {"accession": "GCA_023508785.1", "paired_accession": "GCF_023508785.1",
         "refseq_category": "reference genome", "assembly_level": "Chromosome",
         "assembly_status": "current", "strain": "H4-8", "has_annotation": True,
         "submitter": "x", "submission_date": "2022-01-01"},
        {"accession": "GCF_023508785.1", "paired_accession": "GCA_023508785.1",
         "refseq_category": "reference genome", "assembly_level": "Chromosome",
         "assembly_status": "current", "strain": "H4-8", "has_annotation": True,
         "submitter": "x", "submission_date": "2022-01-01"},
    ] if not require_reference else [])

    row = nr.resolve_row("Scom", "Schizophyllum commune", "H4-8")

    assert row.resolved is True
    assert row.protein_source == "ncbi"  # NOT uniprot -- FTP wasn't available
    assert row.protein_accession == "GCA_023508785.1"
    assert row.genome_source == "ncbi"
    assert row.genome_accession == "GCA_023508785.1"
    assert any("FTP" in line for line in row.report_lines)


def test_resolve_row_no_match_anywhere_leaves_blank(monkeypatch):
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [{"taxon_id": "1", "scientific_name": s}])
    monkeypatch.setattr(nr, "search_uniprot_proteomes", lambda tid, reference_only: [])
    monkeypatch.setattr(nr, "query_ncbi_assemblies", lambda tid, require_reference: [])

    row = nr.resolve_row("Xsp", "Xylaria sp.", "")

    assert row.resolved is False
    assert row.protein_source == ""
    assert row.genome_source == ""


def test_resolve_row_unresolvable_taxon_name(monkeypatch):
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [])

    row = nr.resolve_row("Xsp", "Xylaria sp.", "")

    assert row.resolved is False
    assert any("taxon" in line.lower() for line in row.report_lines)


def test_resolve_row_strain_targeted_widening(monkeypatch):
    # Default reference-only UniProt search finds nothing; strain is given;
    # widened (non-reference) search matches by strain name.
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [{"taxon_id": "1", "scientific_name": s}])

    def fake_search(tid, reference_only):
        if reference_only:
            return []
        return [
            {"proteome_id": "UP1", "strain": "Fo47", "busco_score": 95.0,
             "genome_assembly_id": None, "superkingdom": "Eukaryota"},
            {"proteome_id": "UP2", "strain": "race 4", "busco_score": 93.0,
             "genome_assembly_id": None, "superkingdom": "Eukaryota"},
        ]

    monkeypatch.setattr(nr, "search_uniprot_proteomes", fake_search)
    monkeypatch.setattr(nr, "check_ftp_available", lambda pid, sk: True)
    monkeypatch.setattr(nr, "query_ncbi_assemblies", lambda tid, require_reference: [])

    row = nr.resolve_row("Foxy1", "Fusarium oxysporum", "Fo47")

    assert row.resolved is True
    assert row.protein_source == "uniprot"
    assert row.protein_accession == "UP1"
    assert row.genome_source == ""
    assert any("strain-targeted" in line.lower() for line in row.report_lines)
    assert any("no genome found" in line.lower() for line in row.report_lines)


def test_resolve_row_ncbi_strain_widening_requires_annotation_for_protein(monkeypatch):
    # Step 1 finds nothing; Step 2's default (require_reference) search finds
    # nothing; Step 2b widens by strain and finds exactly one match -- but it
    # has no annotation, so it cannot supply Protein_Source. Must NOT be
    # silently accepted.
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [{"taxon_id": "1", "scientific_name": s}])
    monkeypatch.setattr(nr, "search_uniprot_proteomes", lambda tid, reference_only: [])

    def fake_query(tid, require_reference):
        if require_reference:
            return []
        return [
            {"accession": "GCA_111111111.1", "paired_accession": None,
             "refseq_category": None, "assembly_level": "Contig",
             "assembly_status": "current", "strain": "XYZ123", "has_annotation": False,
             "submitter": "y", "submission_date": "2020-01-01"},
        ]

    monkeypatch.setattr(nr, "query_ncbi_assemblies", fake_query)

    row = nr.resolve_row("Tst1", "Test species", "XYZ123")

    assert row.resolved is False
    assert row.protein_source == ""
    assert any("annotation" in line.lower() for line in row.report_lines)


def test_resolve_row_uniprot_genome_flagged_superseded(monkeypatch):
    # UniProt's proteome record links a genome_assembly_id directly (Step 1's
    # one-query-resolves-both-columns path). If that accession is no longer
    # NCBI's current assembly for the taxon, this must be flagged, not
    # silently accepted (concrete live example from the spec: UP000007431 ->
    # GCA_000143185.1 while NCBI's current is .2).
    monkeypatch.setattr(nr, "resolve_taxon_id", lambda s: [{"taxon_id": "5", "scientific_name": s}])
    monkeypatch.setattr(nr, "search_uniprot_proteomes", lambda tid, reference_only: (
        [{"proteome_id": "UP000007431", "strain": "X", "busco_score": 90.0,
          "genome_assembly_id": "GCA_000143185.1", "superkingdom": "Eukaryota"}]
        if reference_only else []
    ))
    monkeypatch.setattr(nr, "check_ftp_available", lambda pid, sk: True)
    monkeypatch.setattr(nr, "query_ncbi_assemblies", lambda tid, require_reference: [
        {"accession": "GCA_000143185.1", "paired_accession": "GCF_000143185.1",
         "refseq_category": "reference genome", "assembly_level": "Chromosome",
         "assembly_status": "replaced", "strain": "X", "has_annotation": True,
         "submitter": "y", "submission_date": "2010-01-01"},
    ])

    row = nr.resolve_row("Scom2", "Schizophyllum commune", "X")

    assert row.resolved is True
    assert row.genome_accession == "GCA_000143185.1"
    assert any("superseded" in line.lower() for line in row.report_lines)


import csv


HEADER = [
    "Short", "Species", "Strain", "Group", "TaxonGroup",
    "Protein_Source", "Protein_Accession", "Taxon_ID",
    "Genome_Source", "Genome_Accession", "GFF3_Source", "GFF3_Accession",
]


def _write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER)
        w.writeheader()
        w.writerows(rows)


def test_resolve_species_csv_fills_blank_rows_only(tmp_path, monkeypatch, capsys):
    study_dir = tmp_path / "studies" / "fungi" / "toy"
    study_dir.mkdir(parents=True)
    already_row = {  # already filled -- must NOT be touched
        "Short": "Already", "Species": "X", "Strain": "", "Group": "IN", "TaxonGroup": "X",
        "Protein_Source": "local_faa", "Protein_Accession": "/some/path.faa", "Taxon_ID": "",
        "Genome_Source": "local_genome", "Genome_Accession": "/some/path.fna",
        "GFF3_Source": "", "GFF3_Accession": "",
    }
    _write_csv(study_dir / "species.csv", [
        already_row,
        {  # blank -- should resolve
            "Short": "Ncra", "Species": "Neurospora crassa", "Strain": "OR74A", "Group": "OUT",
            "TaxonGroup": "Pezizomycotina",
            "Protein_Source": "", "Protein_Accession": "", "Taxon_ID": "",
            "Genome_Source": "", "Genome_Accession": "", "GFF3_Source": "", "GFF3_Accession": "",
        },
    ])

    def fake_resolve_row(short, species, strain):
        assert short == "Ncra"
        return nr.ResolvedRow(
            short=short, resolved=True, protein_source="uniprot",
            protein_accession="UP000001805", taxon_id="5334",
            genome_source="ncbi", genome_accession="GCA_000182925.2",
            report_lines=[],
        )

    monkeypatch.setattr(nr, "resolve_row", fake_resolve_row)
    nr.resolve_species_csv(study_dir)

    with open(study_dir / "species.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))

    already = next(r for r in rows if r["Short"] == "Already")
    assert already == already_row  # every field on this row is untouched, byte for byte

    ncra = next(r for r in rows if r["Short"] == "Ncra")
    assert ncra["Protein_Source"] == "uniprot"
    assert ncra["Protein_Accession"] == "UP000001805"
    assert ncra["Genome_Source"] == "ncbi"
    assert ncra["Genome_Accession"] == "GCA_000182925.2"

    out = capsys.readouterr().out
    assert "Resolved 1 of 1" in out


def test_resolve_species_csv_leaves_unresolved_blank_and_reports(tmp_path, monkeypatch, capsys):
    study_dir = tmp_path / "studies" / "fungi" / "toy2"
    study_dir.mkdir(parents=True)
    _write_csv(study_dir / "species.csv", [
        {"Short": "Xsp", "Species": "Xylaria sp.", "Strain": "", "Group": "OUT", "TaxonGroup": "X",
         "Protein_Source": "", "Protein_Accession": "", "Taxon_ID": "",
         "Genome_Source": "", "Genome_Accession": "", "GFF3_Source": "", "GFF3_Accession": ""},
    ])

    monkeypatch.setattr(nr, "resolve_row", lambda short, species, strain: nr.ResolvedRow(
        short=short, resolved=False, report_lines=["not a resolvable NCBI taxon name"]
    ))
    nr.resolve_species_csv(study_dir)

    with open(study_dir / "species.csv", newline="") as fh:
        row = next(csv.DictReader(fh))
    assert row["Protein_Source"] == ""

    out = capsys.readouterr().out
    assert "need a decision" in out
    assert "not a resolvable NCBI taxon name" in out


def test_resolve_species_csv_missing_species_csv_exits_cleanly(tmp_path):
    study_dir = tmp_path / "studies" / "fungi" / "nocsv"
    study_dir.mkdir(parents=True)
    # No species.csv written here.
    try:
        nr.resolve_species_csv(study_dir)
        assert False, "expected SystemExit"
    except SystemExit as e:
        assert str(study_dir / "species.csv") in str(e)
        assert "not found" in str(e)
