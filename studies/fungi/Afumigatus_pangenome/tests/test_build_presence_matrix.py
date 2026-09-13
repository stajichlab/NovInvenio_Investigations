import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PRESENT, ABSENT
from build_presence_matrix import (
    split_member_id,
    count_family_members,
    build_matrix,
)


def test_split_member_id_splits_on_first_separator_only():
    assert split_member_id("s1|prot_1") == ("s1", "prot_1")
    # A protein ID that itself contains the separator survives intact.
    assert split_member_id("s1|prot|with|pipes") == ("s1", "prot|with|pipes")
    # No separator at all -> no strain.
    assert split_member_id("bare_protein") == ("", "bare_protein")


def test_count_family_members_counts_copies_and_flags_unknown_strains():
    families = {
        "s1|repA": ["s1|repA", "s1|paralogA", "s2|orthoA"],
        "s2|repB": ["s2|repB", "nosuchstrain|x", "bare_id"],
    }
    counts, unrecognized = count_family_members(families, ["s1", "s2"])
    assert counts[("s1|repA", "s1")] == 2      # two s1 proteins in the family
    assert counts[("s1|repA", "s2")] == 1
    assert counts[("s2|repB", "s2")] == 1
    assert unrecognized == 2                    # nosuchstrain|x and bare_id


def test_build_matrix_sets_presence_copy_number_and_absence():
    families = {
        "s1|famCore": ["s1|famCore", "s2|g1", "s3|g1"],
        "s1|famDup": ["s1|famDup", "s1|famDup_paralog"],
    }
    matrix, unrecognized = build_matrix(families, ["s1", "s2", "s3"])
    assert unrecognized == 0
    assert matrix.strains == ["s1", "s2", "s3"]
    assert matrix.families == ["s1|famCore", "s1|famDup"]

    assert matrix.call("s1|famCore", "s2") == PRESENT
    assert matrix.copy_number[("s1|famCore", "s1")] == 1
    assert matrix.strain_count("s1|famCore") == 3

    # famDup is a two-copy, single-strain family: present only in s1.
    assert matrix.copy_number[("s1|famDup", "s1")] == 2
    assert matrix.call("s1|famDup", "s2") == ABSENT
    assert matrix.strain_count("s1|famDup") == 1


def test_build_matrix_round_trips_through_tsv_with_copy_numbers(tmp_path):
    families = {"s1|famA": ["s1|famA", "s1|famA_b", "s2|x"]}
    matrix, _ = build_matrix(families, ["s1", "s2"])
    out = tmp_path / "presence_matrix.tsv"
    matrix.to_tsv(out)

    from pangenome_matrix import PresenceMatrix
    loaded = PresenceMatrix.from_tsv(out)
    assert loaded.call("s1|famA", "s1") == PRESENT
    assert loaded.copy_number[("s1|famA", "s1")] == 2
    assert loaded.copy_number[("s1|famA", "s2")] == 1


def test_build_matrix_honours_a_custom_separator():
    families = {"s1::famA": ["s1::famA", "s2::g9"]}
    matrix, unrecognized = build_matrix(families, ["s1", "s2"], id_sep="::")
    assert unrecognized == 0
    assert matrix.strain_count("s1::famA") == 2
