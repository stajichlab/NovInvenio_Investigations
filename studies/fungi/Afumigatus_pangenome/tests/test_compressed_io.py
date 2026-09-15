import gzip
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from compressed_io import open_maybe_compressed


def test_open_maybe_compressed_reads_plain_file(tmp_path):
    p = tmp_path / "plain.txt"
    p.write_text("hello\nworld\n")
    with open_maybe_compressed(p) as fh:
        assert fh.read() == "hello\nworld\n"


def test_open_maybe_compressed_reads_gz_file(tmp_path):
    p = tmp_path / "data.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("hello\nworld\n")
    with open_maybe_compressed(p) as fh:
        assert fh.read() == "hello\nworld\n"


def test_open_maybe_compressed_reads_zst_file(tmp_path):
    p = tmp_path / "data.txt.zst"
    plain = tmp_path / "data.txt"
    plain.write_text("hello\nworld\n")
    subprocess.run(["zstd", "-q", "-f", str(plain), "-o", str(p)], check=True)
    with open_maybe_compressed(p) as fh:
        assert fh.read() == "hello\nworld\n"
