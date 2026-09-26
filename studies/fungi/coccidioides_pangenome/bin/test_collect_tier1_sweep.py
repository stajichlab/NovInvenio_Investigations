import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from collect_tier1_sweep import (
    NONE_LABEL,
    aligned_label_arrays,
    ami,
    ari,
    parse_point,
    protein_module_labels,
    read_family_to_module,
    read_protein_to_family,
)


def _write(path, text):
    path.write_text(text)
    return path


def _labels(tmp_path, name, cluster_rows, module_rows):
    cl = _write(tmp_path / f"{name}.cluster.tsv", "".join(f"{r}\t{m}\n" for r, m in cluster_rows))
    fm = _write(tmp_path / f"{name}.family_modules.tsv",
                "family\tmodule_id\tmodule_size\n" + "".join(f"{f}\t{m}\t0\n" for f, m in module_rows))
    return protein_module_labels(read_protein_to_family(cl), read_family_to_module(fm))


def test_same_partition_with_different_family_ids_scores_one(tmp_path):
    # Run A names families by p1/p3; run B picks different representatives
    # (p2/p4) and different module IDs. The protein grouping is identical,
    # so a protein-level comparison must give AMI = ARI = 1.
    a = _labels(tmp_path, "a", [("p1", "p1"), ("p1", "p2"), ("p3", "p3"), ("p3", "p4")],
                [("p1", 0), ("p3", 1)])
    b = _labels(tmp_path, "b", [("p2", "p1"), ("p2", "p2"), ("p4", "p3"), ("p4", "p4")],
                [("p2", 7), ("p4", 3)])
    xa, xb = aligned_label_arrays(a, b, drop_none=False)
    assert ami(xa, xb) == pytest.approx(1.0)
    assert ari(xa, xb) == pytest.approx(1.0)


def test_family_without_module_maps_to_none_and_can_be_dropped(tmp_path):
    a = _labels(tmp_path, "a", [("p1", "p1"), ("p1", "p2"), ("p3", "p3"), ("p4", "p4")],
                [("p1", 0), ("p3", 1)])
    assert a == {"p1": "0", "p2": "0", "p3": "1", "p4": NONE_LABEL}
    b = _labels(tmp_path, "b", [("p1", "p1"), ("p1", "p2"), ("p3", "p3"), ("p3", "p4")],
                [("p1", 0), ("p3", 1)])
    xa, xb = aligned_label_arrays(a, b, drop_none=False)
    assert len(xa) == 4
    xa, xb = aligned_label_arrays(a, b, drop_none=True)
    assert len(xa) == 3  # p4 is 'none' in run A
    assert ami(xa, xb) == pytest.approx(1.0)


def test_protein_sets_must_match(tmp_path):
    a = _labels(tmp_path, "a", [("p1", "p1"), ("p1", "p2")], [("p1", 0)])
    b = _labels(tmp_path, "b", [("p1", "p1")], [("p1", 0)])
    with pytest.raises(ValueError):
        aligned_label_arrays(a, b, drop_none=False)


def test_scores_match_sklearn_reference_values():
    # Reference values from sklearn 1.9.1 adjusted_mutual_info_score /
    # adjusted_rand_score on the same inputs.
    a = np.array([0, 0, 1, 2])
    b = np.array([0, 0, 1, 1])
    assert ami(a, b) == pytest.approx(0.5714285714285714)
    assert ari(a, b) == pytest.approx(0.5714285714285714)
    assert ami(np.array([0, 0, 0, 0]), np.array([0, 1, 2, 3])) == pytest.approx(0.0)


def test_parse_point():
    assert parse_point("coccidioides_sweep_id0.95_cov0.5") == (0.95, 0.5)
