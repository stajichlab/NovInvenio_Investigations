import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from strain_inventory import read_representative_shorts

HEADER = "Short\tn_contigs\tn50\ttotal_length\tdedup_group\tis_representative\n"


def test_reads_only_representative_strains(tmp_path):
    p = tmp_path / "strain_inventory.tsv"
    p.write_text(
        HEADER
        + "s1\t10\t100\t1000\t0\t1\n"
        + "s1_dup\t50\t20\t990\t0\t0\n"
        + "s2\t8\t300\t1100\t1\t1\n"
    )
    assert read_representative_shorts(p) == ["s1", "s2"]


def test_rejects_a_file_that_is_not_an_inventory(tmp_path):
    p = tmp_path / "wrong.tsv"
    p.write_text("family\ts1\ts2\nfamA\tpresent\tabsent\n")
    with pytest.raises(ValueError):
        read_representative_shorts(p)


def test_rejects_an_inventory_with_no_representatives(tmp_path):
    p = tmp_path / "empty.tsv"
    p.write_text(HEADER + "s1\t10\t100\t1000\t0\t0\n")
    with pytest.raises(ValueError):
        read_representative_shorts(p)
