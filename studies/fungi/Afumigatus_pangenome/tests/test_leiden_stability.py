import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from leiden_stability import align_labels, pairwise_ami, singleton_fraction


def test_align_labels_orders_by_shared_family_list():
    a = {"f1": 0, "f2": 0, "f3": 1}
    b = {"f1": 5, "f2": 9, "f3": 9}
    families = ["f1", "f2", "f3"]
    la, lb = align_labels(a, b, families)
    assert la == [0, 0, 1]
    assert lb == [5, 9, 9]


def test_pairwise_ami_identical_partitions_is_one():
    families = ["f1", "f2", "f3", "f4"]
    partitions = [
        {"f1": 0, "f2": 0, "f3": 1, "f4": 1},
        {"f1": 7, "f2": 7, "f3": 3, "f4": 3},  # same grouping, different labels
    ]
    mean_ami, std_ami = pairwise_ami(partitions, families)
    assert mean_ami == 1.0
    assert std_ami == 0.0


def test_pairwise_ami_unrelated_partitions_is_low():
    families = [f"f{i}" for i in range(20)]
    # partition 1: two even blocks; partition 2: interleaved (unrelated structure)
    p1 = {f: (0 if i < 10 else 1) for i, f in enumerate(families)}
    p2 = {f: (i % 2) for i, f in enumerate(families)}
    mean_ami, std_ami = pairwise_ami([p1, p2], families)
    assert mean_ami < 0.2


def test_singleton_fraction():
    partition = {"f1": 0, "f2": 0, "f3": 1, "f4": 2}
    # modules {f3} and {f4} are singletons, out of 3 modules total
    assert singleton_fraction(partition) == 2 / 3
