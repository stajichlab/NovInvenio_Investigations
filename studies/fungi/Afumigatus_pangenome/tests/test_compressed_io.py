import gzip
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from compressed_io import open_maybe_compressed, open_maybe_compressed_write


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


def test_open_maybe_compressed_write_plain_file(tmp_path):
    p = tmp_path / "plain.txt"
    with open_maybe_compressed_write(p) as fh:
        fh.write("hello\nworld\n")
    assert p.read_text() == "hello\nworld\n"


def test_open_maybe_compressed_write_gz_roundtrip(tmp_path):
    p = tmp_path / "data.txt.gz"
    with open_maybe_compressed_write(p) as fh:
        fh.write("hello\nworld\n")
    with gzip.open(p, "rt") as fh:
        assert fh.read() == "hello\nworld\n"
    # Round-trips through the read-side helper too.
    with open_maybe_compressed(p) as fh:
        assert fh.read() == "hello\nworld\n"


def test_open_maybe_compressed_write_zst_roundtrip(tmp_path):
    p = tmp_path / "data.txt.zst"
    with open_maybe_compressed_write(p) as fh:
        fh.write("hello\nworld\n")
    decompressed = subprocess.run(
        ["zstd", "-dc", str(p)], capture_output=True, text=True, check=True,
    )
    assert decompressed.stdout == "hello\nworld\n"
    with open_maybe_compressed(p) as fh:
        assert fh.read() == "hello\nworld\n"


def test_open_maybe_compressed_write_zst_raises_on_subprocess_failure(tmp_path):
    """A write to a path whose parent directory doesn't exist must raise,
    not silently produce a truncated/missing .zst file -- see
    _ZstdWriter.__exit__'s returncode check."""
    bad_path = tmp_path / "no_such_dir" / "data.txt.zst"
    with pytest.raises(subprocess.CalledProcessError):
        with open_maybe_compressed_write(bad_path) as fh:
            fh.write("hello\n")
