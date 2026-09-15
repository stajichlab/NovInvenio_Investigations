import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import (
    PRESENT, GENOME_ONLY, ABSENT,
    read_cluster_tsv, build_families, PresenceMatrix, copy_number_path,
)


def test_read_cluster_tsv_parses_rep_member_pairs(tmp_path):
    p = tmp_path / "clusters.tsv"
    p.write_text("repA\trepA\nrepA\tmemberA2\nrepB\trepB\n")
    result = read_cluster_tsv(p)
    assert result == {"repA": "repA", "memberA2": "repA", "repB": "repB"}


def test_build_families_groups_by_rep_including_singletons():
    member_to_rep = {"repA": "repA", "memberA2": "repA", "repB": "repB"}
    families = build_families(member_to_rep)
    assert families == {"repA": ["memberA2", "repA"], "repB": ["repB"]}


def test_presence_matrix_set_and_query():
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2", "s3"])
    pm.set_call("famA", "s1", PRESENT)
    pm.set_call("famA", "s2", GENOME_ONLY)
    # famA absent in s3 (never set) -- default ABSENT
    assert pm.call("famA", "s1") == PRESENT
    assert pm.call("famA", "s3") == ABSENT
    assert pm.is_present("famA", "s1") is True
    assert pm.is_present("famA", "s2") is True   # genome_only still counts as present
    assert pm.is_present("famA", "s3") is False
    assert pm.presence_vector("famA") == [True, True, False]
    assert pm.strain_count("famA") == 2
    assert pm.frequency("famA") == 2 / 3


def test_presence_matrix_tsv_roundtrip(tmp_path):
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", PRESENT)
    pm.set_call("famB", "s2", GENOME_ONLY)
    out = tmp_path / "matrix.tsv"
    pm.to_tsv(out)
    loaded = PresenceMatrix.from_tsv(out)
    assert loaded.strains == ["s1", "s2"]
    assert loaded.families == ["famA", "famB"]
    assert loaded.call("famA", "s1") == PRESENT
    assert loaded.call("famB", "s2") == GENOME_ONLY
    assert loaded.call("famA", "s2") == ABSENT


def test_tsv_roundtrip_preserves_copy_number(tmp_path):
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", PRESENT, copies=3)
    pm.set_call("famA", "s2", PRESENT, copies=1)
    pm.set_call("famB", "s2", GENOME_ONLY)
    out = tmp_path / "matrix.tsv"
    pm.to_tsv(out)

    assert copy_number_path(out).exists()
    loaded = PresenceMatrix.from_tsv(out)
    assert loaded.copy_number[("famA", "s1")] == 3
    assert loaded.copy_number[("famA", "s2")] == 1
    # A cell with no copies recorded stays 0, and the states round-trip unchanged.
    assert loaded.copy_number.get(("famB", "s2"), 0) == 0
    assert loaded.call("famA", "s1") == PRESENT
    assert loaded.call("famB", "s2") == GENOME_ONLY


@pytest.mark.parametrize("suffix", [".gz", ".zst"])
def test_presence_matrix_tsv_roundtrip_compressed(tmp_path, suffix):
    """to_tsv actually compresses when given a .gz/.zst path (not just
    accepts one on read) -- the write-side half of this study's general
    storage-compression convention, matching from_tsv's existing transparent
    read support."""
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", PRESENT)
    pm.set_call("famB", "s2", GENOME_ONLY)
    out = tmp_path / f"matrix.tsv{suffix}"
    pm.to_tsv(out)

    assert out.exists()
    # A real compressed file is not readable as plain UTF-8 text -- decoding
    # its raw bytes as text must fail (or trivially "succeed" into garbage
    # that isn't the TSV header), otherwise to_tsv silently wrote plain text
    # under a misleading .gz/.zst name.
    raw = out.read_bytes()
    assert not raw.startswith(b"family\t")

    loaded = PresenceMatrix.from_tsv(out)
    assert loaded.strains == ["s1", "s2"]
    assert loaded.families == ["famA", "famB"]
    assert loaded.call("famA", "s1") == PRESENT
    assert loaded.call("famB", "s2") == GENOME_ONLY
    assert loaded.call("famA", "s2") == ABSENT


@pytest.mark.parametrize("suffix", [".gz", ".zst"])
def test_compressed_matrix_gets_compressed_copy_number_sidecar(tmp_path, suffix):
    pm = PresenceMatrix(families=["famA"], strains=["s1"])
    pm.set_call("famA", "s1", PRESENT, copies=4)
    out = tmp_path / f"matrix.tsv{suffix}"
    pm.to_tsv(out)

    sidecar = copy_number_path(out)
    assert sidecar.exists()
    assert str(sidecar).endswith(suffix), (
        f"sidecar {sidecar} should carry the same compression suffix as {out}"
    )
    assert not sidecar.read_bytes().startswith(b"family\t")

    loaded = PresenceMatrix.from_tsv(out)
    assert loaded.copy_number[("famA", "s1")] == 4


def test_state_only_matrix_writes_no_sidecar_and_loads_fine(tmp_path):
    """Backward compatibility: a matrix with no copy numbers produces exactly
    one file, and loading a matrix whose sidecar is absent is not an error."""
    pm = PresenceMatrix(families=["famA"], strains=["s1"])
    pm.set_call("famA", "s1", PRESENT)
    out = tmp_path / "matrix.tsv"
    pm.to_tsv(out)
    assert not copy_number_path(out).exists()
    loaded = PresenceMatrix.from_tsv(out)
    assert loaded.copy_number == {}
    assert loaded.call("famA", "s1") == PRESENT


def test_from_tsv_rejects_unknown_state(tmp_path):
    p = tmp_path / "corrupt.tsv"
    p.write_text("family\ts1\ts2\nfamA\tpresent\tPRESENT\n")
    try:
        PresenceMatrix.from_tsv(p)
        assert False, "expected ValueError"
    except ValueError as exc:
        assert "invalid presence state" in str(exc)


def test_set_call_rejects_invalid_state():
    pm = PresenceMatrix(families=["famA"], strains=["s1"])
    try:
        pm.set_call("famA", "s1", "not_a_real_state")
        assert False, "expected ValueError"
    except ValueError:
        pass
