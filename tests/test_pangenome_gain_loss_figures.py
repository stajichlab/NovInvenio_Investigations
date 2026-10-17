"""Tests for the fig4 class grouping in bin/pangenome_gain_loss_figures.py."""
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "bin" / "pangenome_gain_loss_figures.py"
spec = importlib.util.spec_from_file_location("pangenome_gain_loss_figures", SCRIPT)
glf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(glf)


def _row(gain, loss, amb):
    return {"direction_gain": str(gain), "direction_loss": str(loss),
            "direction_ambiguous": str(amb), "direction_conserved": "0"}


BYCLS = {
    "core": _row(0, 0, 0),
    "soft_core": _row(1, 2, 3),
    "shell": _row(10, 20, 30),
    "cloud": _row(100, 0, 5),
    "singleton": _row(1000, 1, 0),
    "nonrep_only": _row(7, 7, 7),
}


def test_class_group_counts_sums_only_member_classes():
    g = glf.class_group_counts(BYCLS)
    assert list(g) == ["shell+soft_core", "cloud+singleton"]
    assert g["shell+soft_core"] == {"gain": 11, "loss": 22, "ambiguous": 33}
    assert g["cloud+singleton"] == {"gain": 1100, "loss": 1, "ambiguous": 5}
    # core and nonrep_only never enter either group
    total = sum(sum(v.values()) for v in g.values())
    assert total == 66 + 1106


def test_class_group_rows_fractions():
    rows = glf.class_group_rows([("S", {"bycls": BYCLS})])
    assert len(rows) == 6
    d = {(r[1], r[2]): (r[3], float(r[4])) for r in rows}
    assert d[("shell+soft_core", "loss")] == (22, round(22 / 66, 6))
    assert d[("cloud+singleton", "gain")][0] == 1100
    for grp in ("shell+soft_core", "cloud+singleton"):
        assert abs(sum(f for (g, _), (_, f) in d.items() if g == grp) - 1.0) < 1e-5


def test_empty_group_fraction_is_zero():
    empty = dict(BYCLS, cloud=_row(0, 0, 0), singleton=_row(0, 0, 0))
    rows = glf.class_group_rows([("S", {"bycls": empty})])
    assert [r[4] for r in rows if r[1] == "cloud+singleton"] == ["0.000000"] * 3


def test_limitation_text_only_for_small_outgroups():
    def st(n):
        return {"prov": {"strain_counts": {"out_all": n}}}
    txt = glf.limitation_text([("A", st(2)), ("B", st(360)), ("C", st(2))])
    assert "A and C: outgroup = 2 genomes" in txt and "0/2, 1/2, 2/2" in txt
    assert "B" not in txt
    assert glf.limitation_text([("B", st(360))]) == ""
