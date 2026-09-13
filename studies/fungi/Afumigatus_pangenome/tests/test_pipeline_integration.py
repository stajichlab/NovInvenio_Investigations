"""End-to-end synthetic chain test.

Every other test file exercises one function in isolation, which structurally
cannot catch a SEAM mismatch between two steps -- and two real ones existed:
nothing produced the `presence_matrix.tsv` every downstream step read, and
copy numbers silently vanished across `to_tsv`/`from_tsv`. This test walks a
tiny hand-built cluster TSV through the real chain:

    cluster TSV -> build_presence_matrix -> to_tsv/from_tsv
                -> rescue_pass.apply_rescue
                -> frequency_bins.compute_frequency_table
                -> cooccurrence.find_cooccurring_pairs

3 strains, 4 families, entirely in tmp_path. No external tools.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY, ABSENT, build_families, read_cluster_tsv
from build_presence_matrix import build_matrix
from rescue_pass import apply_rescue
from frequency_bins import compute_frequency_table
from cooccurrence import find_cooccurring_pairs

# s1, s2, s3 are ingroup; out1 is the outgroup reference used for polarization.
STRAINS = ["s1", "s2", "s3"]
ALL_COLUMNS = STRAINS + ["out1"]

# rep<TAB>member, Short-prefixed exactly as build_presence_matrix.py requires.
CLUSTER_TSV = "\n".join([
    # famCore: every strain plus the outgroup; s1 carries two copies.
    "s1|core\ts1|core",
    "s1|core\ts1|core_paralog",
    "s1|core\ts2|core",
    "s1|core\ts3|core",
    "s1|core\tout1|core",
    # famShellA / famShellB: the same two strains -> a perfectly co-occurring pair.
    "s1|shellA\ts1|shellA",
    "s1|shellA\ts2|shellA",
    "s1|shellB\ts1|shellB",
    "s1|shellB\ts2|shellB",
    # famSingleton: s3 only.
    "s3|only\ts3|only",
]) + "\n"


def _matrix(tmp_path) -> PresenceMatrix:
    cluster = tmp_path / "tier1_cluster.tsv"
    cluster.write_text(CLUSTER_TSV)
    families = build_families(read_cluster_tsv(cluster))
    matrix, unrecognized = build_matrix(families, ALL_COLUMNS)
    assert unrecognized == 0
    return matrix


def test_cluster_tsv_to_presence_matrix_file(tmp_path):
    matrix = _matrix(tmp_path)
    out = tmp_path / "presence_matrix.tsv"
    matrix.to_tsv(out)

    loaded = PresenceMatrix.from_tsv(out)
    assert loaded.strains == ALL_COLUMNS
    assert loaded.families == ["s1|core", "s1|shellA", "s1|shellB", "s3|only"]
    assert loaded.call("s1|core", "s2") == PRESENT
    assert loaded.call("s1|shellA", "s3") == ABSENT
    # The seam that used to drop copy numbers on the floor.
    assert loaded.copy_number[("s1|core", "s1")] == 2


def test_rescue_then_binning_then_cooccurrence(tmp_path):
    matrix = _matrix(tmp_path)
    path = tmp_path / "presence_matrix.tsv"
    matrix.to_tsv(path)
    matrix = PresenceMatrix.from_tsv(path)

    # Rescue pass: a genomic tblastn hit puts famShellA into s3 as genome_only.
    applied, skipped = apply_rescue(matrix, {("s1|shellA", "s3")})
    assert (applied, skipped) == (1, 0)
    assert matrix.call("s1|shellA", "s3") == GENOME_ONLY

    rescued = tmp_path / "presence_matrix.rescued.tsv"
    matrix.to_tsv(rescued)
    matrix = PresenceMatrix.from_tsv(rescued)
    assert matrix.copy_number[("s1|core", "s1")] == 2   # survives the rescue write

    # Binning over the INGROUP strains only: the outgroup column must not be
    # part of the denominator (3 strains, not 4).
    table = compute_frequency_table(matrix, strains=STRAINS)
    by_family = {row["family"]: row for row in table}
    assert by_family["s1|core"]["frequency"] == 1.0
    assert by_family["s1|core"]["bin"] == "core"
    # shellA was rescued into all three strains; shellB is still 2/3.
    assert by_family["s1|shellA"]["strain_count"] == 3
    assert by_family["s1|shellB"]["strain_count"] == 2
    assert by_family["s1|shellB"]["bin"] == "shell"
    assert by_family["s3|only"]["bin"] == "singleton"

    # Co-occurrence over the same ingroup strain list. Force shellA back to its
    # pre-rescue state so the pair is a clean 2/3-vs-2/3 co-occurrence.
    matrix.set_call("s1|shellA", "s3", ABSENT)
    table = compute_frequency_table(matrix, strains=STRAINS)
    outgroup_presence = {
        fam: (int(matrix.is_present(fam, "out1")), 1) for fam in matrix.families
    }
    pairs = find_cooccurring_pairs(
        matrix,
        table,
        clade_of_strain={s: "cladeX" for s in STRAINS},
        outgroup_presence=outgroup_presence,
        min_strain_count=2,
        fdr_alpha=1.0,       # tiny fixture: keep every tested pair
        n_perms=20,
        strains=STRAINS,
    )
    tested = {(p["family_a"], p["family_b"]) for p in pairs}
    assert ("s1|shellA", "s1|shellB") in tested
    # Core and singleton families are not eligible (wrong bin / below the floor).
    assert all("s1|core" not in pair and "s3|only" not in pair for pair in tested)

    shell_pair = next(p for p in pairs if p["family_a"] == "s1|shellA")
    assert shell_pair["jaccard"] == 1.0
    # Absent from the outgroup reference -> a gain in the strains that carry it.
    assert shell_pair["direction_a"] == "gain"
    assert shell_pair["clade_composition"] == {"cladeX": 2}
