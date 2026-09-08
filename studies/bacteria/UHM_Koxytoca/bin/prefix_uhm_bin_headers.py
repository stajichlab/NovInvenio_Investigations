#!/usr/bin/env python3
"""Prefix UHM MAG protein/DNA FASTA headers with their bin's Short code.

Why: the 18 UHM ingroup MAGs (Group=IN in config.csv) are each a separately
binned+prodigal-called assembly (build_koxytoca_config.py's docstring), so
megahit's contig ids (k141_<N>) restart per bin. Across the 18 bins' pep FASTAs
that collides for real: 341 duplicate protein ids (e.g. k141_4769_1 appearing
in two unrelated bins) out of ~85k total, which is silently ambiguous
everywhere downstream that keys on protein_id alone (presence matrix,
novelties/core/losses reports, annotation merge) -- two different proteins
from two different bins can get merged under one id.

Rewrites data_dir/pep/<Short>.pep.fa and data_dir/dna/<Short>.dna.fa in place
(config.csv's own Protein/DNA columns, filtered to GROUP=IN), prepending
"<Short>__" to the id token of every header line. Idempotent: a header already
carrying its own Short prefix is left alone, so re-running after a partial
edit or a fresh build_koxytoca_config.py copy is safe.

This is a study-specific script (lives under this study's own bin/, not NII's
shared bin/ -- see CLAUDE.md's "Where new code goes"), so it hardcodes NII_ROOT
below rather than deriving paths from __file__/cwd.
"""
import csv
import re
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))
from provenance import append_manifest, build_record, sha256_of  # noqa: E402

STUDY_DIR = NII_ROOT / "studies/bacteria/UHM_Koxytoca"

HEADER_RE = re.compile(r"^>(\S+)(.*)$")


def prefix_file(path: Path, short: str) -> bool:
    """Rewrite path's headers in place with '<short>__' prepended to the id.
    Returns True if the file was modified, False if already prefixed (no-op)."""
    lines = path.read_text().splitlines(keepends=True)
    changed = False
    out = []
    for line in lines:
        if line.startswith(">"):
            m = HEADER_RE.match(line.rstrip("\n"))
            ident, rest = m.group(1), m.group(2)
            if ident.startswith(f"{short}__"):
                out.append(line)
                continue
            newline = "\n" if line.endswith("\n") else ""
            out.append(f">{short}__{ident}{rest}{newline}")
            changed = True
        else:
            out.append(line)
    if changed:
        path.write_text("".join(out))
    return changed


def main() -> int:
    config_csv = STUDY_DIR / "config.csv"
    if not config_csv.exists():
        sys.exit(f"ERROR: {config_csv} not found -- run build_koxytoca_config.py first")

    pep_dir, dna_dir = STUDY_DIR / "data_dir" / "pep", STUDY_DIR / "data_dir" / "dna"
    manifest_records = []
    n_modified = 0

    with open(config_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["GROUP"] != "IN":
                continue
            short = row["Short"]

            for fname, subdir in ((row["Protein"], pep_dir), (row["DNA"], dna_dir)):
                if not fname:
                    continue
                fpath = subdir / fname
                if not fpath.exists():
                    sys.exit(f"[{short}] ERROR: expected {fpath} not found")
                modified = prefix_file(fpath, short)
                status = "prefixed" if modified else "already prefixed, skipped"
                print(f"[{short}] {fpath.name}: {status}", file=sys.stderr)
                if modified:
                    n_modified += 1
                    manifest_records.append(build_record(
                        source_url="(derived, no external source)",
                        source_release="n/a",
                        license="n/a",
                        local_path=fpath,
                        checksum=sha256_of(fpath),
                        derived_by=(
                            f"headers prefixed with '{short}__' to disambiguate megahit "
                            "contig ids (k141_N) that collide across separately-assembled "
                            "UHM MAG bins; studies/bacteria/UHM_Koxytoca/bin/prefix_uhm_bin_headers.py"
                        ),
                    ))

    if manifest_records:
        append_manifest(manifest_records, STUDY_DIR / "DATA_MANIFEST.yaml")
        print(f"\nUpdated {STUDY_DIR / 'DATA_MANIFEST.yaml'} ({len(manifest_records)} records)", file=sys.stderr)

    print(f"\n{n_modified} file(s) modified", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
