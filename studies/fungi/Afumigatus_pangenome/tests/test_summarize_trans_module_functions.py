import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from summarize_trans_module_functions import load_modules


def test_load_modules_excludes_singletons_by_default(tmp_path):
    f = tmp_path / "family_modules.tsv"
    f.write_text(
        "family\tmodule_id\tmodule_size\n"
        "famA\t0\t1\n"
        "famB\t1\t2\n"
        "famC\t1\t2\n"
        "famD\t2\t3\n"
        "famE\t2\t3\n"
        "famF\t2\t3\n"
    )
    modules = load_modules(str(f), min_module_size=2)
    assert modules == {
        "1": ["famB", "famC"],
        "2": ["famD", "famE", "famF"],
    }


def test_load_modules_respects_custom_min_size(tmp_path):
    f = tmp_path / "family_modules.tsv"
    f.write_text(
        "family\tmodule_id\tmodule_size\n"
        "famA\t0\t2\n"
        "famB\t0\t2\n"
        "famC\t1\t5\n"
        "famD\t1\t5\n"
    )
    modules = load_modules(str(f), min_module_size=5)
    assert set(modules.keys()) == {"1"}
    assert sorted(modules["1"]) == ["famC", "famD"]


def test_load_modules_empty_file_returns_empty(tmp_path):
    f = tmp_path / "family_modules.tsv"
    f.write_text("family\tmodule_id\tmodule_size\n")
    assert load_modules(str(f), min_module_size=2) == {}
