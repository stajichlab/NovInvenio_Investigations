import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "bin" / "pangenome_outgroup_polarity.py"


def _run(tmp_path, freq_text):
    (tmp_path / "m.tsv").write_text("family\ti1\to1\nfA\tpresent\tpresent\nfB\tabsent\tpresent\n")
    (tmp_path / "f.tsv").write_text(freq_text)
    (tmp_path / "p.tsv").write_text("family_a\tfamily_b\tdirection_a\nfA\tfB\tgain\n")
    return subprocess.run([sys.executable, str(SCRIPT), "--matrix", str(tmp_path / "m.tsv"),
                           "--freq", str(tmp_path / "f.tsv"), "--pairs", str(tmp_path / "p.tsv"),
                           "--outgroup", "o1"], capture_output=True, text=True, check=True).stdout


def test_reads_old_four_column_table(tmp_path):
    out = _run(tmp_path, "family\tfrequency\tstrain_count\tbin\nfA\t1.0\t1\tcore\nfB\t0.0\t0\tsingleton\n")
    assert "\ncore\t1\t" in out


def test_reads_per_group_table_and_reports_new_classes(tmp_path):
    # NovInvenio #212: 7-column frequency_table with outgroup_only / nonrep_only.
    out = _run(tmp_path,
               "family\tfrequency\tstrain_count\tbin\tfrequency_out\tstrain_count_out\tbin_out\n"
               "fA\t1.0\t1\tcore\t-\t-\t-\nfB\t0.0\t0\toutgroup_only\t-\t-\t-\n")
    assert "\ncore\t1\t" in out
    assert "\noutgroup_only\t1\t" in out
