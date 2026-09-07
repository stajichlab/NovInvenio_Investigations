"""Shared provenance-record helpers for NII's data-pull scripts.

Every tracked/curated data file NII commits must carry a provenance record (see
DESIGN.md Sec 4): source URL, source version/release, access date, checksum, license,
and (if derived) the script+parameters that produced it. Raw ephemeral pulls under
data/ are gitignored and never need this -- this module is for the sidecar records that
accompany anything that gets promoted out of data/ into a tracked location
(config_support/, studies/<domain>/<set>/).
"""
from __future__ import annotations

import hashlib
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


def sha256_of(path: Path, chunk_size: int = 1 << 20) -> str:
    """Return the sha256 hex digest of a file, streamed (safe for multi-GB inputs)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_record(
    *,
    source_url: str,
    source_release: str,
    license: str,
    local_path: Path,
    checksum: str | None = None,
    derived_by: str | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build one provenance record dict, ready to write via write_record()/append_manifest().

    checksum: pass explicitly if already computed (avoids re-hashing large files);
        otherwise computed here from local_path.
    derived_by: free-text description of the script+parameters that produced local_path,
        when it's not a raw pull (e.g. "bin/build_study_config.py --study pezizo_set1").
    """
    return {
        "local_path": str(local_path),
        "source_url": source_url,
        "source_release": source_release,
        "access_date": now_utc_iso(),
        "checksum_sha256": checksum or sha256_of(local_path),
        "license": license,
        "derived_by": derived_by,
        "pulled_with": f"{Path(sys.argv[0]).name} on {platform.node()}",
        **(extra or {}),
    }


def write_record(record: dict[str, Any], sidecar_path: Path) -> None:
    """Write a single provenance record as a YAML sidecar next to the data file it describes."""
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    with open(sidecar_path, "w") as fh:
        yaml.safe_dump(record, fh, sort_keys=False)


def append_manifest(records: list[dict[str, Any]], manifest_path: Path) -> None:
    """Append/merge records into a study-level DATA_MANIFEST.yaml (list of records, keyed
    by local_path so re-running a pull updates its own entry rather than duplicating it)."""
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    existing: dict[str, Any] = {}
    if manifest_path.exists():
        with open(manifest_path) as fh:
            loaded = yaml.safe_load(fh) or []
        existing = {r["local_path"]: r for r in loaded}
    for r in records:
        existing[r["local_path"]] = r
    with open(manifest_path, "w") as fh:
        yaml.safe_dump(list(existing.values()), fh, sort_keys=False)
