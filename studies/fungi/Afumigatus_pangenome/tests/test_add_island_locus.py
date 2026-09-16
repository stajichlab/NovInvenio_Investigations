import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from add_island_locus import (
    load_needed_families,
    load_family_members_by_strain,
    load_gene_positions_for_strains,
    compute_locus,
    annotate_islands,
)


def test_load_needed_families_collects_all_member_families(tmp_path):
    path = tmp_path / "islands.tsv"
    path.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\n"
        "1\tAsfu_A\t2\tfamA,famB\n"
        "1\tAsfu_B\t2\tfamB,famC\n"
    )
    assert load_needed_families(str(path)) == {"famA", "famB", "famC"}


def test_load_family_members_by_strain_filters_to_needed_families(tmp_path):
    path = tmp_path / "cluster.tsv"
    path.write_text(
        "famA\tAsfu_X|p1\n"
        "famA\tAsfu_Y|p2\n"
        "famZ\tAsfu_X|p3\n"  # famZ not needed -- must be skipped
    )
    result = load_family_members_by_strain(str(path), needed_families={"famA"})
    assert result == {"famA": {"Asfu_X": ["p1"], "Asfu_Y": ["p2"]}}


def test_load_gene_positions_for_strains_filters_to_needed_strains(tmp_path):
    path = tmp_path / "gene_positions.tsv"
    path.write_text(
        "Short\tprotein_id\tcontig\tstart\tend\n"
        "Asfu_X\tp1\tctg1\t100\t500\n"
        "Asfu_Y\tp2\tctg1\t600\t900\n"
        "Asfu_Z\tp9\tctg2\t1\t50\n"  # Asfu_Z not needed
    )
    result = load_gene_positions_for_strains(str(path), needed_strains={"Asfu_X", "Asfu_Y"})
    assert result == {
        ("Asfu_X", "p1"): ("ctg1", 100, 500),
        ("Asfu_Y", "p2"): ("ctg1", 600, 900),
    }


def test_compute_locus_spans_all_resolved_members():
    family_members_by_strain = {
        "famA": {"Asfu_X": ["p1"]},
        "famB": {"Asfu_X": ["p2"]},
    }
    gene_positions = {
        ("Asfu_X", "p1"): ("ctg1", 100, 500),
        ("Asfu_X", "p2"): ("ctg1", 600, 900),
    }
    locus = compute_locus("Asfu_X", ["famA", "famB"], family_members_by_strain, gene_positions)
    assert locus["contig"] == "ctg1"
    assert locus["start"] == 100
    assert locus["end"] == 900
    assert locus["n_resolved"] == 2
    assert locus["n_total"] == 2
    assert locus["locus_id"] == "Asfu_X:ctg1:100-900"


def test_compute_locus_handles_unresolved_members_without_crashing():
    # famC has no entry for Asfu_X at all (e.g. a rescue-pass-only call for
    # this strain, with no annotated protein_id to resolve a position from).
    family_members_by_strain = {"famA": {"Asfu_X": ["p1"]}}
    gene_positions = {("Asfu_X", "p1"): ("ctg1", 100, 500)}
    locus = compute_locus("Asfu_X", ["famA", "famC"], family_members_by_strain, gene_positions)
    assert locus["n_resolved"] == 1
    assert locus["n_total"] == 2
    assert locus["contig"] == "ctg1"
    assert locus["start"] == 100
    assert locus["end"] == 500


def test_compute_locus_returns_none_locus_when_nothing_resolves():
    locus = compute_locus("Asfu_X", ["famA"], {}, {})
    assert locus["n_resolved"] == 0
    assert locus["locus_id"] == ""


def test_compute_locus_uses_widest_span_across_multiple_contigs_flagged():
    # A pathological case: members resolve to two different contigs (should
    # not normally happen for a real accessory island, since it's built from
    # one strain's single-contig consecutive run) -- must not silently
    # average across contigs; report the example (first-seen) contig and
    # flag it via n_contigs > 1 rather than crash or produce a bogus span.
    family_members_by_strain = {
        "famA": {"Asfu_X": ["p1"]},
        "famB": {"Asfu_X": ["p2"]},
    }
    gene_positions = {
        ("Asfu_X", "p1"): ("ctg1", 100, 500),
        ("Asfu_X", "p2"): ("ctg2", 600, 900),
    }
    locus = compute_locus("Asfu_X", ["famA", "famB"], family_members_by_strain, gene_positions)
    assert locus["n_contigs"] == 2


def test_annotate_islands_adds_locus_columns():
    rows = [{"example_strain": "Asfu_X", "member_families": "famA,famB"}]
    family_members_by_strain = {
        "famA": {"Asfu_X": ["p1"]},
        "famB": {"Asfu_X": ["p2"]},
    }
    gene_positions = {
        ("Asfu_X", "p1"): ("ctg1", 100, 500),
        ("Asfu_X", "p2"): ("ctg1", 600, 900),
    }
    annotated = annotate_islands(rows, family_members_by_strain, gene_positions)
    assert annotated[0]["locus_id"] == "Asfu_X:ctg1:100-900"
    assert annotated[0]["n_members_with_coordinates"] == "2"
