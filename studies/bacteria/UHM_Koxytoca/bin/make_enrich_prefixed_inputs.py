#!/usr/bin/env python3
"""Make enrichment inputs for UHM_Koxytoca whose IN protein ids match the pipeline results.

Why: the pipeline run in results/UHM_Koxytoca (2026-09-10) read IN (UHM MAG) proteomes
whose headers had the "<Short>__" prefix (bin/prefix_uhm_bin_headers.py; megahit contig
ids k141_<N> repeat across bins). The 2026-09-11 rebuild by bin/build_study_config.py
wrote data_dir/pep/<Short>.pep.fa again WITHOUT the prefix, and the Pfam annotation files
(annotations/<name>.tsv.gz, 2026-10-14) were made from those unprefixed FASTAs. So the
candidate ids (UHM125bin22::UHM125bin22__k141_154965_41) are not in data_dir/pep and not
in annotations/. enrich_novelties_v2.py stops with "candidates not in universe".

What: writes, into OUT_DIR (default studies/bacteria/UHM_Koxytoca/data_dir_enrich_prefixed,
gitignored by studies/*/*/data_dir_*):
  pep/<Short>.pep.fa             IN: header id -> "<Short>__<id>" (same rule as
                                 bin/prefix_uhm_bin_headers.py; already-prefixed ids kept).
                                 OUT: symlink to data_dir/pep/<Short>.pep.fa (ids unchanged).
  annotations/<name>.tsv.gz      IN: accession -> "<Short>__<accession>"; other columns
                                 unchanged. OUT: copy of annotations/<name>.tsv.gz (a symlink is
                                 refused by lib/annotation_names.annotation_path).
  annotations/<name>.tsv.gz.provenance.yaml
                                 IN: the original sidecar with local_path, checksum_sha256,
                                 derived_from and derived_by updated. OUT: copy.
  PROVENANCE.yaml                inputs and outputs with sha256, command, date.
Nothing in data_dir/ or annotations/ is changed. The prefix is checked: every IN id must
be unique within its bin, and the script stops if a prefixed id would collide.

Study-specific script (CLAUDE.md "Where new code goes"): absolute NII_ROOT, no __file__.
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import os
import shutil
import sys
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
sys.path.insert(0, str(NII_ROOT / "lib"))
import yaml  # noqa: E402
from annotation_names import annotation_names, annotation_path  # noqa: E402

STUDY = NII_ROOT / "studies/bacteria/UHM_Koxytoca"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def link(src: Path, dst: Path) -> None:
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    os.symlink(src, dst)


def prefix_fasta(src: Path, dst: Path, short: str) -> int:
    seen = set()
    n = 0
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                ident, _, rest = line[1:].rstrip("\n").partition(" ")
                new = ident if ident.startswith(f"{short}__") else f"{short}__{ident}"
                if new in seen:
                    sys.exit(f"duplicate id {new} in {src}")
                seen.add(new)
                n += 1
                fo.write(f">{new} {rest}\n" if rest else f">{new}\n")
            else:
                fo.write(line)
    return n


def prefix_annotation(src: Path, dst: Path, short: str) -> int:
    n = 0
    with gzip.open(src, "rt", newline="") as fi, gzip.open(dst, "wt", newline="") as fo:
        rdr = csv.DictReader(fi, delimiter="\t")
        w = csv.DictWriter(fo, fieldnames=rdr.fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rdr:
            a = row["accession"]
            row["accession"] = a if a.startswith(f"{short}__") else f"{short}__{a}"
            w.writerow(row)
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=STUDY / "data_dir_enrich_prefixed")
    args = ap.parse_args()
    out = args.out_dir
    (out / "pep").mkdir(parents=True, exist_ok=True)
    (out / "annotations").mkdir(parents=True, exist_ok=True)
    cfg = list(csv.DictReader(open(STUDY / "config.csv")))
    names = annotation_names(list(csv.DictReader(open(STUDY / "species.csv"))))
    cmd = " ".join([sys.executable] + sys.argv)
    rec = {"generated": datetime.datetime.now().isoformat(timespec="seconds"),
           "script": {"path": str(Path(sys.argv[0]).resolve()), "sha256": sha256(Path(sys.argv[0]).resolve())},
           "command": cmd,
           "class": "ephemeral derived input (gitignored; regenerate with the command above)",
           "why": "IN protein ids in results/UHM_Koxytoca carry a '<Short>__' prefix; "
                  "data_dir/pep and annotations/ do not (see script docstring).",
           "files": []}
    for r in cfg:
        sp = r["Short"]
        src_fa = STUDY / "data_dir/pep" / r["Protein"]
        dst_fa = out / "pep" / r["Protein"]
        src_an = annotation_path(STUDY / "annotations", names[sp])
        dst_an = out / "annotations" / src_an.name
        side_src = Path(str(src_an) + ".provenance.yaml")
        side_dst = Path(str(dst_an) + ".provenance.yaml")
        if not src_an.exists():
            sys.exit(f"missing annotation file for {sp}: {src_an}")
        if r["GROUP"] == "IN":
            if dst_fa.is_symlink():
                dst_fa.unlink()
            if dst_an.is_symlink():
                dst_an.unlink()
            nf = prefix_fasta(src_fa, dst_fa, sp)
            na = prefix_annotation(src_an, dst_an, sp)
            side = yaml.safe_load(side_src.read_text()) if side_src.exists() else {}
            side["local_path"] = str(dst_an)
            side["checksum_sha256"] = sha256(dst_an)
            side["derived_from"] = {"path": str(src_an), "sha256": sha256(src_an)}
            side["derived_by"] = f"{cmd} (accession prefixed with '{sp}__'); original: {side.get('derived_by', '')}"
            if side_dst.is_symlink():
                side_dst.unlink()
            side_dst.write_text(yaml.safe_dump(side, sort_keys=False))
            action = "prefixed"
        else:
            link(src_fa, dst_fa)
            # Annotation files are copied, not linked: lib/annotation_names.annotation_path
            # resolves symlinks and refuses a file outside the --annotations dir.
            for s, d in ((src_an, dst_an), (side_src, side_dst)):
                if d.is_symlink() or d.exists():
                    d.unlink()
                if s.exists():
                    shutil.copy2(s, d)
            nf = na = None
            action = "fasta symlink, annotation copy"
        rec["files"].append({"short": sp, "group": r["GROUP"], "action": action,
                             "fasta_src": str(src_fa), "fasta_src_sha256": sha256(src_fa),
                             "fasta_out": str(dst_fa), "fasta_records": nf,
                             "annotation_src": str(src_an), "annotation_src_sha256": sha256(src_an),
                             "annotation_out": str(dst_an), "annotation_rows": na})
        print(f"{sp}\t{r['GROUP']}\t{action}\tfasta={nf}\tannotation_rows={na}")
    (out / "PROVENANCE.yaml").write_text(yaml.safe_dump(rec, sort_keys=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
