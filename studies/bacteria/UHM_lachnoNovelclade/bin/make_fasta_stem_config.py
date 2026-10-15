#!/usr/bin/env python3
"""Write a derived copy of config.csv whose Short equals the protein FASTA stem.

Why: in studies/bacteria/UHM_lachnoNovelclade/config.csv, 6 of 9 Short values
(UHM1088b56, UHM1210b85, UHM896b22, UHM904b119, UHM207b71, Lachimisp1) differ from the
protein FASTA stems (Ub56, Ub85, Ub22, Ub119, Ub71, Lsp1). The pipeline results
(results/UHM_lachnoNovelclade: presence-matrix columns, candidates.txt, novelties.<Short>.tsv)
and species.csv use the FASTA stems. enrich_novelties_v2.py reads Short from --config, so it
needs a config whose Short is the stem. The tracked config.csv is not changed.

What: reads config.csv, sets Short = Protein without ".pep.fa" (all other columns
unchanged), writes OUT (default analysis/enrichment/config.fasta_stems.csv) and
OUT.provenance.yaml (source path + sha256, the Short changes, command, date).
Stops if a stem is not a species.csv Short or not a presence_matrix.tsv column.

Study-specific script (CLAUDE.md "Where new code goes"): absolute NII_ROOT, no __file__.
"""
import argparse
import csv
import datetime
import hashlib
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
STUDY = NII_ROOT / "studies/bacteria/UHM_lachnoNovelclade"
RESULTS = NII_ROOT / "results/UHM_lachnoNovelclade"
sys.path.insert(0, str(NII_ROOT / "lib"))
import yaml  # noqa: E402


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=STUDY / "analysis/enrichment/config.fasta_stems.csv")
    args = ap.parse_args()
    src = STUDY / "config.csv"
    rows = list(csv.DictReader(open(src)))
    fields = list(rows[0].keys())
    sp_short = {r["Short"] for r in csv.DictReader(open(STUDY / "species.csv"))}
    with open(RESULTS / "presence_matrix.tsv") as fh:
        cols = set(fh.readline().rstrip("\n").split("\t"))
    changes = []
    for r in rows:
        stem = r["Protein"][: -len(".pep.fa")] if r["Protein"].endswith(".pep.fa") else None
        if not stem:
            sys.exit(f"Protein does not end in .pep.fa: {r['Protein']}")
        if stem not in sp_short or stem not in cols:
            sys.exit(f"stem {stem} not in species.csv Short or presence_matrix.tsv columns")
        if stem != r["Short"]:
            changes.append({"config_short": r["Short"], "fasta_stem": stem})
        r["Short"] = stem
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    rec = {
        "local_path": str(args.out),
        "checksum_sha256": sha256(args.out),
        "derived_from": {"path": str(src), "sha256": sha256(src)},
        "derived_by": " ".join([sys.executable] + sys.argv),
        "script_sha256": sha256(Path(sys.argv[0]).resolve()),
        "date": datetime.datetime.now().isoformat(timespec="seconds"),
        "rule": "Short = Protein without '.pep.fa'; all other columns unchanged",
        "why": ("config.csv Short differs from the protein FASTA stem for some species; the "
                "pipeline results and species.csv use the stem. Used as enrich_novelties_v2.py "
                "--config. The tracked config.csv is unchanged."),
        "short_changes": changes,
    }
    Path(str(args.out) + ".provenance.yaml").write_text(yaml.safe_dump(rec, sort_keys=False))
    print(f"wrote {args.out}; {len(changes)} Short values changed: {changes}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
