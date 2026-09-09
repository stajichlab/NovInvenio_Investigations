#!/usr/bin/env python3
"""Build the Ncra model-organism crosswalk for pezizo_set1's modelorgs.yaml.

Neurospora crassa (Ncra) is one of this study's own IN-group proteomes
(studies/fungi/pezizo_set1/data_dir/pep/Ncra.pep.fa), fetched from UniProt
reference proteome UP000001805 -- so its FASTA headers are plain UniProt
format: ">tr|A7UWL5|A7UWL5_NEUCR ... GN=NCU10683 ...". lib/model_organisms.py
(nf_NovInvenio) needs a gene_names_csv keyed by FungiDB NCU gene ID to attach
FungiDB's own Product Description/Gene Name text -- config_support/modelorgs/
Neurospora_crassa_gene_names_FungiDB.csv (FungiDB release 68, see
config_support/MODELORG_NCRA_PROVENANCE.md) -- but that means resolving each
Ncra protein's raw ID down to its bare NCU accession first.

Two options existed:
  1. A real diamond blastp search against a FungiDB Ncra proteome (the
     approach config_support/modelorgs/Ncra_vs_FungiDB_Ncra.diamond.tsv used
     for some *other*, non-UniProt local protein ID scheme -- its query IDs
     are "FC69C3D3_000001-T1"-style, which don't match this study's
     Ncra.pep.fa UniProt headers at all; that file is NOT used here).
  2. Since Ncra.pep.fa's own headers already carry the NCU ID directly in a
     "GN=" field, no alignment is needed -- id_transform: diamond_fasta's
     query->subject->gene 2-hop lookup can run entirely self-referentially:
     protein_fasta = Ncra.pep.fa itself (fasta_gene_field: "GN", matching its
     own "GN=NCU10683" tag -- no synthesized FASTA needed), and diamond_hits
     only needs to be an IDENTITY map (protein_id -> itself) to satisfy
     ModelOrgAnnotator's 2-hop code path, which always requires both fields
     set for id_transform: diamond_fasta (lib/model_organisms.py's
     _validate()). This script produces exactly that identity map -- it is
     NOT a diamond search output despite the field name, and is exact (a
     literal ID match, not a similarity search).

Output: config_support/modelorgs/Ncra_self_id_crosswalk.tsv (3 cols, no
header: protein_id, protein_id, "100" -- one row per Ncra.pep.fa sequence).
The 3rd column is a harmless placeholder, not a real percent-identity value
-- lib/model_organisms.py's _load_diamond_hits() splits each line on '\t' and
keeps parts[1] AS-IS (no rstrip), so with only 2 columns the subject-ID field
would retain the line's trailing '\n' and silently fail every lookup (real
diamond outfmt6 output never hits this because its subject-ID column always
has more columns after it, e.g. evalue/bitscore, absorbing the newline
instead) -- confirmed by hand: a 2-column version produced 0/9759 resolved
gene_keys end-to-end despite every intermediate map loading correctly.

This is a study-specific script (hardcodes this study's Ncra.pep.fa; see
CLAUDE.md's "Where new code goes" -- study-specific vs. shared scripts), so
NII_ROOT is hardcoded rather than derived from __file__/cwd.
"""
import re
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))
from provenance import build_record, sha256_of, write_record  # noqa: E402

NCRA_PEP_FASTA = NII_ROOT / "studies" / "fungi" / "pezizo_set1" / "data_dir" / "pep" / "Ncra.pep.fa"
OUT_TSV = NII_ROOT / "config_support" / "modelorgs" / "Ncra_self_id_crosswalk.tsv"

HEADER_RE = re.compile(r"^>(\S+)")


def build_identity_map(fasta_path: Path, out_path: Path) -> int:
    n = 0
    with open(fasta_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if not line.startswith(">"):
                continue
            m = HEADER_RE.match(line)
            if not m:
                sys.exit(f"ERROR: unrecognized FASTA header: {line!r}")
            pid = m.group(1)
            fout.write(f"{pid}\t{pid}\t100\n")
            n += 1
    return n


def main() -> int:
    if not NCRA_PEP_FASTA.exists():
        sys.exit(f"ERROR: {NCRA_PEP_FASTA} not found")

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    n = build_identity_map(NCRA_PEP_FASTA, OUT_TSV)
    print(f"[Ncra crosswalk] wrote {OUT_TSV} ({n} proteins, identity map)", file=sys.stderr)

    checksum = sha256_of(OUT_TSV)
    record = build_record(
        source_url="",
        source_release="derived, not fetched -- see module docstring",
        license="N/A (derived from this study's own UniProt-sourced Ncra.pep.fa)",
        local_path=OUT_TSV,
        checksum=checksum,
        derived_by=(
            "studies/fungi/pezizo_set1/bin/build_ncra_fungidb_crosswalk.py, "
            "identity map from studies/fungi/pezizo_set1/data_dir/pep/Ncra.pep.fa "
            "headers (UniProt proteome UP000001805, see this study's own "
            "DATA_MANIFEST.yaml for that pull's provenance)"
        ),
    )
    write_record(record, OUT_TSV.with_suffix(OUT_TSV.suffix + ".provenance.yaml"))
    print(f"[Ncra crosswalk] OK  sha256={checksum[:12]}...", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
