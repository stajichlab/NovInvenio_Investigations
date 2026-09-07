#!/usr/bin/env python3
"""Pull go-basic.obo (GO-DAG structure, needed for goatools ORA -- DESIGN.md Sec 5/7).

Distinct source from UniProt: geneontology.org, not the UniProt FTP. Same ephemeral/
recipe-driven treatment (DESIGN.md Sec 4): lands under --outdir (default data/go/,
gitignored), never archived, provenance sidecar written alongside it.
"""
import argparse
import re
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from provenance import build_record, write_record  # noqa: E402

GO_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"
LICENSE = "CC-BY-4.0 (Gene Ontology Consortium, http://geneontology.org/docs/go-citation-policy/)"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data/go", help="Output directory (default: data/go)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    dest = outdir / "go-basic.obo"

    print(f"downloading {GO_OBO_URL} -> {dest}", file=sys.stderr)
    # purl.obolibrary.org's WAF 403s on urllib's default User-Agent -- set one explicitly.
    req = urllib.request.Request(GO_OBO_URL, headers={"User-Agent": "NovInvenio_Investigations/1.0"})
    with urllib.request.urlopen(req) as resp, open(dest, "wb") as fh:
        fh.write(resp.read())

    # go-basic.obo's own header carries "data-version: releases/YYYY-MM-DD"
    header = dest.read_text(encoding="utf-8", errors="replace")[:2000]
    m = re.search(r"data-version:\s*(\S+)", header)
    release = m.group(1) if m else "unknown"

    record = build_record(
        source_url=GO_OBO_URL,
        source_release=f"go-basic.obo {release}",
        license=LICENSE,
        local_path=dest,
    )
    write_record(record, dest.with_suffix(dest.suffix + ".provenance.yaml"))
    print(f"OK  release={release}  sha256={record['checksum_sha256'][:12]}...", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
