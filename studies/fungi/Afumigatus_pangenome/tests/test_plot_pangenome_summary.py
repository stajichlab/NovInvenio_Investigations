import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT
from plot_pangenome_summary import (
    band_counts, matrix_to_binary_array, accumulation_curve, read_frequency_table,
)


def test_read_frequency_table_parses_rows(tmp_path):
    f = tmp_path / "freq.tsv"
    f.write_text("family\tfrequency\tstrain_count\tbin\nfamA\t1.0000\t10\tcore\n")
    rows = read_frequency_table(str(f))
    assert rows == [{"family": "famA", "frequency": "1.0000", "strain_count": "10", "bin": "core"}]


def test_band_counts_tallies_by_bin():
    rows = [
        {"bin": "core"}, {"bin": "core"}, {"bin": "shell"},
        {"bin": "singleton"}, {"bin": "singleton"}, {"bin": "singleton"},
    ]
    counts = band_counts(rows)
    assert counts["core"] == 2
    assert counts["shell"] == 1
    assert counts["singleton"] == 3
    assert counts["cloud"] == 0
    assert counts["soft_core"] == 0


def test_matrix_to_binary_array_matches_presence_matrix():
    strains = ["s1", "s2", "s3"]
    pm = PresenceMatrix(families=["famA", "famB"], strains=strains)
    pm.set_call("famA", "s1", PRESENT)
    pm.set_call("famA", "s2", PRESENT)
    pm.set_call("famB", "s3", PRESENT)
    arr = matrix_to_binary_array(pm, ["famA", "famB"], strains)
    expected = np.array([[True, True, False], [False, False, True]])
    assert np.array_equal(arr, expected)


def test_accumulation_curve_core_never_exceeds_pangenome_size():
    # famA present in all 4 strains (core); famB present in only strain 0.
    arr = np.array([
        [True, True, True, True],
        [True, False, False, False],
    ])
    pan_mean, pan_std, core_mean, core_std = accumulation_curve(arr, n_permutations=10, seed=0)
    assert len(pan_mean) == 4
    # At n=4 strains sampled, pangenome must include both families (union=2)
    assert pan_mean[-1] == 2
    # Core (intersection over ALL 4) can only ever be famA -> core size = 1 at n=4
    assert core_mean[-1] == 1
    # Core count is always <= pangenome count at every sample size
    assert all(c <= p for c, p in zip(core_mean, pan_mean))


def test_accumulation_curve_is_deterministic_given_seed():
    arr = np.array([[True, False, True], [False, True, True], [True, True, False]])
    result_a = accumulation_curve(arr, n_permutations=5, seed=42)
    result_b = accumulation_curve(arr, n_permutations=5, seed=42)
    for a, b in zip(result_a, result_b):
        assert np.array_equal(a, b)


def test_accumulation_curve_pangenome_never_decreases():
    arr = np.random.default_rng(0).random((20, 10)) > 0.5
    pan_mean, _, _, _ = accumulation_curve(arr, n_permutations=15, seed=1)
    assert all(pan_mean[i] <= pan_mean[i + 1] for i in range(len(pan_mean) - 1))
