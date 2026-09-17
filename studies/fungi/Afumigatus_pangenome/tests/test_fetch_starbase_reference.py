import hashlib
import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

import fetch_starbase_reference as fsr


def _fake_urlretrieve_factory(content_by_name: dict[str, bytes]):
    def _fake(url, dest, *a, **kw):
        for name, content in content_by_name.items():
            if name in url:
                Path(dest).write_bytes(content)
                return
        raise AssertionError(f"unexpected url {url}")
    return _fake


def test_fetch_skips_existing_file_with_correct_checksum(tmp_path):
    name = "starship-scan.py"
    content = b"already here"
    md5 = hashlib.md5(content).hexdigest()
    with patch.dict(fsr.FILES, {name: {"url": "http://x/" + name, "md5": md5, "size": len(content)}}, clear=True):
        (tmp_path / name).write_bytes(content)
        with patch("urllib.request.urlretrieve") as mock_dl:
            fetched = fsr.fetch_starbase_files(tmp_path)
            mock_dl.assert_not_called()
        assert fetched == [tmp_path / name]


def test_fetch_downloads_and_verifies_checksum(tmp_path):
    name = "starship-scan.py"
    content = b"downloaded content"
    md5 = hashlib.md5(content).hexdigest()
    with patch.dict(fsr.FILES, {name: {"url": "http://x/" + name, "md5": md5, "size": len(content)}}, clear=True):
        with patch("urllib.request.urlretrieve", side_effect=_fake_urlretrieve_factory({name: content})):
            fetched = fsr.fetch_starbase_files(tmp_path)
        assert fetched == [tmp_path / name]
        assert (tmp_path / name).read_bytes() == content


def test_fetch_raises_on_checksum_mismatch(tmp_path):
    name = "starship-scan.py"
    content = b"corrupted"
    with patch.dict(fsr.FILES, {name: {"url": "http://x/" + name, "md5": "0" * 32, "size": len(content)}}, clear=True):
        with patch("urllib.request.urlretrieve", side_effect=_fake_urlretrieve_factory({name: content})):
            try:
                fsr.fetch_starbase_files(tmp_path)
                assert False, "expected RuntimeError"
            except RuntimeError as e:
                assert "checksum mismatch" in str(e)


def test_fetch_redownloads_with_force_even_if_present(tmp_path):
    name = "starship-scan.py"
    old_content = b"stale"
    new_content = b"fresh"
    md5 = hashlib.md5(new_content).hexdigest()
    with patch.dict(fsr.FILES, {name: {"url": "http://x/" + name, "md5": md5, "size": len(new_content)}}, clear=True):
        (tmp_path / name).write_bytes(old_content)
        with patch("urllib.request.urlretrieve", side_effect=_fake_urlretrieve_factory({name: new_content})):
            fsr.fetch_starbase_files(tmp_path, force=True)
        assert (tmp_path / name).read_bytes() == new_content
