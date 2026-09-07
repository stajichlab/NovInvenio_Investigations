#!/usr/bin/env python3
"""Pull one UniProt reference proteome's protein FASTA + annotation flat file.

Ephemeral/recipe-driven pull (DESIGN.md Sec 4/5): downloads land under --outdir
(default data/uniprot/, gitignored) and are never archived. This script is the
checked-in artifact; its output is disposable and re-fetchable on demand.

Per proteome, downloads exactly two files from the UniProt reference_proteomes FTP
archive:
  - {Proteome_ID}_{taxid}.fasta.gz  -- canonical protein sequences
  - {Proteome_ID}_{taxid}.dat.gz    -- SwissProt/TrEMBL flat file: carries taxonomy
        (OC/OX lines), gene names (GN), and GO/Pfam/InterPro cross-references (DR
        lines) all in one file -- see DESIGN.md Sec 5 for why this replaces a
        separate GOA pull for v1.

A provenance sidecar (<output>.provenance.yaml) is written next to each downloaded
file recording the exact source URL, the UniProt release, access date, and a sha256
checksum, per DESIGN.md's provenance rule.
"""
import argparse
import re
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from provenance import build_record, now_utc_iso, sha256_of, write_record  # noqa: E402

UNIPROT_REST = "https://rest.uniprot.org"
UNIPROT_FTP_BASE = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/reference_proteomes"
)
LICENSE = "CC-BY-4.0 (UniProt, https://www.uniprot.org/help/license)"

# UniProt's proteome JSON 'superkingdom' field values -> reference_proteomes FTP folder
DOMAIN_FOLDERS = {
    "Eukaryota": "Eukaryota",
    "Bacteria": "Bacteria",
    "Archaea": "Archaea",
    "Viruses": "Viruses",
}


def fetch_json(proteome_id: str) -> dict:
    url = f"{UNIPROT_REST}/proteomes/{proteome_id}.json"
    with urllib.request.urlopen(url) as resp:
        import json

        return json.load(resp)


def fetch_release_string() -> str:
    """Parse 'Release 2026_02, 10-Jun-2026' out of the reference_proteomes README."""
    url = f"{UNIPROT_FTP_BASE}/README"
    with urllib.request.urlopen(url) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    m = re.search(r"Release\s+(\S+),\s+([\d-]+-\w+-\d+)", text)
    if not m:
        return "unknown (README format changed, check manually)"
    return f"{m.group(1)} ({m.group(2)})"


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, dest)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--proteome-id", required=True, help="UniProt proteome ID, e.g. UP000001805")
    ap.add_argument("--taxid", help="NCBI taxon ID (looked up from UniProt if omitted)")
    ap.add_argument("--outdir", default="data/uniprot", help="Root output directory (default: data/uniprot)")
    ap.add_argument("--short", help="Short species code, used only for a human-readable log line")
    args = ap.parse_args()

    upid = args.proteome_id
    print(f"[{upid}] fetching proteome metadata...", file=sys.stderr)
    meta = fetch_json(upid)
    taxid = args.taxid or str(meta["taxonomy"]["taxonId"])
    superkingdom = meta.get("superkingdom", "Eukaryota")
    domain_folder = DOMAIN_FOLDERS.get(superkingdom, "Eukaryota")
    genome_assembly = meta.get("genomeAssembly", {})
    release = fetch_release_string()

    outdir = Path(args.outdir) / upid
    stem = f"{upid}_{taxid}"
    records = []
    for ext in ("fasta.gz", "dat.gz"):
        remote = f"{UNIPROT_FTP_BASE}/{domain_folder}/{upid}/{stem}.{ext}"
        local = outdir / f"{stem}.{ext}"
        label = args.short or upid
        print(f"[{label}] downloading {remote} -> {local}", file=sys.stderr)
        download(remote, local)
        checksum = sha256_of(local)
        record = build_record(
            source_url=remote,
            source_release=f"UniProt reference proteomes release {release}",
            license=LICENSE,
            local_path=local,
            checksum=checksum,
            extra={
                "proteome_id": upid,
                "taxon_id": taxid,
                "organism": meta.get("taxonomy", {}).get("scientificName"),
                "genome_assembly": genome_assembly,
            },
        )
        write_record(record, local.with_suffix(local.suffix + ".provenance.yaml"))
        records.append(record)
        print(f"[{label}] OK  sha256={checksum[:12]}...  ({local.stat().st_size:,} bytes)", file=sys.stderr)

    print(
        f"[{upid}] done. GCA accession for genome/GFF3 pull: "
        f"{genome_assembly.get('assemblyId', 'NOT FOUND')}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
