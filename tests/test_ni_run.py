import importlib.machinery
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def _load_ni_module():
    """bin/ni has no .py extension, so it needs an explicit file-based import
    (importlib.util) rather than a plain `import ni` -- same reasoning as any
    extensionless CLI entry point in this repo."""
    loader = importlib.machinery.SourceFileLoader("ni_cli", str(REPO / "bin" / "ni"))
    spec = importlib.util.spec_from_file_location("ni_cli", REPO / "bin" / "ni", loader=loader)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_run_subcommand_invokes_bin_run_study_sh(tmp_path, monkeypatch):
    """Regression test: `ni run` previously built cmd[0] as
    <repo_root>/run_study.sh (which does not exist -- the script lives at
    bin/run_study.sh) and crashed with FileNotFoundError on every real
    invocation. `run` was never actually exercised end-to-end before this,
    so the bug shipped silently."""
    ni = _load_ni_module()

    study_dir = REPO / "studies" / "bacteria" / "cyanobacteria"
    calls = []

    def fake_run(cmd, *args, **kwargs):
        calls.append(cmd)

        class FakeResult:
            returncode = 0
        return FakeResult()

    monkeypatch.setattr(ni.subprocess, "run", fake_run)
    monkeypatch.setattr(sys, "argv", [
        "ni", "run", "--study-dir", str(study_dir), "--", "-profile", "slurm",
    ])

    rc = ni.main()

    assert rc == 0
    assert len(calls) == 1
    cmd = calls[0]
    assert cmd[0] == str(ni.BIN / "run_study.sh")
    assert Path(cmd[0]).exists(), f"{cmd[0]} does not exist -- ni run would crash for real"
    assert cmd[1] == "bacteria/cyanobacteria"
    assert "-profile" in cmd and "slurm" in cmd
    assert "--" not in cmd  # the literal "--" separator is stripped, not forwarded


def test_run_subcommand_rejects_study_dir_outside_studies_root(tmp_path, monkeypatch):
    ni = _load_ni_module()
    outside = tmp_path / "not_under_studies"
    outside.mkdir()

    monkeypatch.setattr(sys, "argv", ["ni", "run", "--study-dir", str(outside)])

    rc = ni.main()
    assert rc == 1
