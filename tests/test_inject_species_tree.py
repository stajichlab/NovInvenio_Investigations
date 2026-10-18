"""bin/inject_species_tree.py: card insertion, replacement and run.json hash upkeep."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "bin" / "inject_species_tree.py"
spec = importlib.util.spec_from_file_location("inject_species_tree", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

PAGE = (
    "<html><body>\n  <section class=\"card\">\n    <h2 class=\"card-title\">Reports</h2>\n  </section>\n\n"
    "  <section class=\"card\">\n    <h2 class=\"card-title\">Run summary</h2>\n  </section>\n</body></html>\n"
)
INFO = {"method": "outgroup", "n_out": 5, "n_out_outside": 0}


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def test_inject_adds_card_before_run_summary():
    out = mod.inject(PAGE, mod.build_card(b"png", 10, INFO))
    assert out.count('id="species-tree"') == 1
    assert out.index('id="species-tree"') < out.index("Run summary")
    assert "data:image/png;base64," in out
    assert "Rooted on the OUT group (monophyletic)." in out


def test_inject_twice_replaces_the_card():
    once = mod.inject(PAGE, mod.build_card(b"one", 10, INFO))
    twice = mod.inject(once, mod.build_card(b"two", 7, INFO))
    assert twice.count('id="species-tree"') == 1
    assert "top 7 markers" in twice and "top 10 markers" not in twice


def test_inject_without_anchor_leaves_page_alone():
    with pytest.raises(SystemExit):
        mod.inject("<html><body>no cards</body></html>", "<section>x</section>")


def test_run_json_hash_updated_when_it_matched(tmp_path):
    old, new = b"old report", b"new report"
    meta = tmp_path / "run.json"
    meta.write_text(json.dumps({"report_sha256": _sha(old), "n_species": 9}))
    note = mod.update_run_json(meta, old, new)
    data = json.loads(meta.read_text())
    assert data["report_sha256"] == _sha(new)
    assert data["n_species"] == 9
    assert "updated" in note


def test_run_json_hash_left_alone_when_it_did_not_match(tmp_path):
    meta = tmp_path / "run.json"
    meta.write_text(json.dumps({"report_sha256": "deadbeef"}))
    note = mod.update_run_json(meta, b"old", b"new")
    assert json.loads(meta.read_text())["report_sha256"] == "deadbeef"
    assert "WARNING" in note


def test_run_json_missing_is_not_an_error(tmp_path):
    assert mod.update_run_json(tmp_path / "run.json", b"a", b"b") == ""


def test_main_writes_card_and_updates_run_json(tmp_path, monkeypatch):
    run = tmp_path / "run"
    run.mkdir()
    report = run / "report.html"
    report.write_text(PAGE)
    (run / "run.json").write_text(json.dumps({"report_sha256": _sha(PAGE.encode())}))
    tree = tmp_path / "tree"
    (tree / "filter").mkdir(parents=True)
    (tree / "filter" / "a.aa.mfa").write_text(">x\nAA\n")
    (tree / "filter" / "treeness.tsv").write_text("x\n")
    (tree / "species_tree.svg").write_text("<svg/>")
    (tree / "species_tree.rooting.json").write_text(json.dumps(INFO))
    monkeypatch.setattr(mod, "render_png", lambda svg, density=150: b"\x89PNG")
    assert mod.main(["--report", str(report), "--tree-dir", str(tree)]) == 0
    html = report.read_text()
    assert html.count('id="species-tree"') == 1
    assert "top 1 markers" in html
    assert json.loads((run / "run.json").read_text())["report_sha256"] == _sha(report.read_bytes())
