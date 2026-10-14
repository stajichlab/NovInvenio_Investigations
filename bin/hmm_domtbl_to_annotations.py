#!/usr/bin/env python3
"""Convert a HMMER Pfam domain table (--domtblout) into a per-proteome annotation TSV.

Why: proteomes from NCBI or JGI have no UniProt record, so bin/extract_dat_annotations.py
cannot annotate them. A Pfam-A hmmsearch (or hmmscan) run gives Pfam domains for these
proteins. This script writes them in the same file format as the UniProt annotation
files (<study>/annotations/<name>.tsv.gz), so the enrichment scripts read both sources
the same way. <name> comes from lib/annotation_names.py: the Protein_Accession for
UniProt/NCBI accessions; for a local FASTA path, its basename without the
proteome/compression suffix. With --study-dir and --short, --out defaults to that file,
and an --out outside <study-dir>/annotations/ is refused.

Input
  --domtbl  one or more HMMER 3 domain tables: plain text, .gz or .zst (read with
            `zstd -dc`). Several files (for example per-chunk outputs of one proteome)
            are read as one table. The layout is detected per file:
              hmmsearch  col 1 = protein, col 4 = Pfam name, col 5 = Pfam accession
              hmmscan    col 1 = Pfam name, col 2 = Pfam accession, col 4 = protein
            Detection uses the "# Program:" trailer line. Without a trailer, the column
            that holds a PFnnnnn accession decides. --tool forces a layout; a conflict
            with the trailer stops the script.
  --fasta   optional protein FASTA of the proteome. With it, the output has one row per
            protein in FASTA order (proteins without a Pfam hit get empty Pfam columns,
            as in the UniProt files), the description comes from the FASTA header, and
            a domain-table protein that is not in the FASTA stops the script.

Hits are used as given. The script does not apply a new threshold: the thresholds are
the ones of the HMMER run (normally --cut_ga). The Pfam version suffix is removed
(PF00244.27 -> PF00244).

Output columns (tab-separated, gzip when the name ends in .gz). The first 11 columns are
the columns of bin/extract_dat_annotations.py output, in the same order:
  accession      first token of the FASTA header (the key the enrichment script uses)
  taxon_id       --taxon-id, or empty
  gene_name      empty
  description    FASTA header text after the first token, without a trailing
                 " [Organism]" tag; empty if not known
  go_ids         empty (not available from a Pfam search)
  pfam_ids       Pfam accessions, "|"-separated, one per family, ordered by the
                 envelope start of the family's first domain
  pfam_names     Pfam short names (Pfam ID, e.g. "14-3-3"), same order as pfam_ids
  interpro_ids   empty (not available from a Pfam search)
  ec_numbers, alphafold_id, xrefs   empty
Two more columns follow (extra; readers that select columns by name ignore them):
  pfam_best_ievalue       best (lowest) domain i-Evalue per family, same order
  pfam_best_domain_score  best (highest) domain bit score per family, same order

A Pfam clan column is not written: the enrichment script
(studies/fungi/zoosporic_dikarya/analysis/enrichment/enrich_novelties_v2.py) maps Pfam
families to clans itself from Pfam-A.hmm.dat (--pfam-dat).

A provenance sidecar <out>.provenance.yaml is written with lib/provenance.py
(build_record/write_record): source domain tables with sha256, HMMER program, version and
recorded command, Pfam release and path, the protein FASTA sha256, counts, the command of
this script, and the statement that GO and InterPro are not available for this source.

Example
  bin/hmm_domtbl_to_annotations.py \\
      --domtbl analysis/hmmscan/Neocalif.pfam.hmmsearch.domtbl.zst \\
      --fasta data_dir/pep/Neocalif.pep.fa --taxon-id 1754190 \\
      --out annotations/GCA_002104975.1.tsv.gz
  # same, with the default name:  --study-dir studies/fungi/<study> --short Neocalif
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from annotation_names import annotation_names, annotation_path, ensure_within  # noqa: E402
from provenance import build_record, sha256_of, write_record  # noqa: E402

UNIPROT_COLUMNS = ["accession", "taxon_id", "gene_name", "description", "go_ids", "pfam_ids",
                   "pfam_names", "interpro_ids", "ec_numbers", "alphafold_id", "xrefs"]
EXTRA_COLUMNS = ["pfam_best_ievalue", "pfam_best_domain_score"]
COLUMNS = UNIPROT_COLUMNS + EXTRA_COLUMNS
PFAM_ACC_RE = re.compile(r"^PF\d{5}(\.\d+)?$")
ORGANISM_TAG_RE = re.compile(r"\s+\[[^\[\]]*\]$")
NOT_AVAILABLE = ("GO and InterPro are not available for this source: the annotation comes "
                 "from a Pfam-A HMM search only. The go_ids and interpro_ids columns are empty "
                 "for every protein.")
PFAM_LICENSE = "Pfam data: CC0 1.0 (https://www.ebi.ac.uk/interpro/about/pfam/)"


# --- reading -------------------------------------------------------------------


def read_text(path: Path) -> list[str]:
    """Return the lines of a plain, .gz or .zst text file."""
    path = Path(path)
    if path.suffix == ".zst":
        zstd = shutil.which("zstd")
        if zstd is None:
            raise SystemExit(f"ERROR: {path} is .zst but no zstd binary is in PATH")
        res = subprocess.run([zstd, "-dc", str(path)], capture_output=True, check=False)
        if res.returncode != 0:
            raise SystemExit(f"ERROR: zstd -dc {path} failed: {res.stderr.decode()[-500:]}")
        return res.stdout.decode().splitlines()
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            return fh.read().splitlines()
    with open(path) as fh:
        return fh.read().splitlines()


def parse_trailer(lines: list[str]) -> dict:
    """Read the '# Program:', '# Version:', '# Option settings:', '# Query file:',
    '# Target file:' and '# Date:' lines that HMMER writes at the end of a table."""
    keys = {"Program": "program", "Version": "version", "Option settings": "option_settings",
            "Query file": "query_file", "Target file": "target_file", "Date": "date"}
    out: dict = {}
    for line in lines:
        if not line.startswith("# "):
            continue
        key, sep, val = line[2:].partition(":")
        if sep and key.strip() in keys and keys[key.strip()] not in out:
            out[keys[key.strip()]] = val.strip()
    return out


def detect_layout(lines: list[str], trailer: dict, tool: str) -> str:
    """Return 'hmmsearch' or 'hmmscan'."""
    program = trailer.get("program", "")
    if program not in ("", "hmmsearch", "hmmscan"):
        raise SystemExit(f"ERROR: unsupported HMMER program {program!r} (need hmmsearch or hmmscan)")
    if tool != "auto":
        if program and program != tool:
            raise SystemExit(f"ERROR: --tool {tool} but the table trailer says {program}")
        return tool
    if program:
        return program
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) >= 5 and PFAM_ACC_RE.match(f[1]):
            return "hmmscan"
        if len(f) >= 5 and PFAM_ACC_RE.match(f[4]):
            return "hmmsearch"
        raise SystemExit(f"ERROR: cannot detect the table layout from line: {line[:120]!r}; use --tool")
    return "hmmsearch"  # empty table: the layout does not change the (empty) result


def parse_domtbl(lines: list[str], layout: str) -> list[dict]:
    """Return one dict per domain row."""
    rows = []
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split(None, 22)
        if len(f) < 22:
            raise SystemExit(f"ERROR: domain table row has {len(f)} fields (need >= 22): {line[:120]!r}")
        if layout == "hmmsearch":
            protein, name, acc = f[0], f[3], f[4]
            desc = f[22] if len(f) > 22 else ""
        else:
            protein, name, acc = f[3], f[0], f[1]
            desc = ""  # hmmscan: col 23 describes the Pfam family, not the protein
        if not PFAM_ACC_RE.match(acc):
            raise SystemExit(f"ERROR: {acc!r} is not a Pfam accession ({layout} layout); line: {line[:120]!r}")
        rows.append({"protein": protein, "pfam": acc.split(".")[0], "name": name,
                     "ievalue": float(f[12]), "score": float(f[13]), "env_from": int(f[19]),
                     "desc": desc.strip()})
    return rows


def read_fasta_headers(path: Path) -> list[tuple[str, str]]:
    """Return [(first token, rest of header)] in file order."""
    out = []
    for line in read_text(path):
        if line.startswith(">"):
            head = line[1:].strip()
            tok, _, rest = head.partition(" ")
            out.append((tok, rest.strip()))
    return out


def clean_description(text: str) -> str:
    return ORGANISM_TAG_RE.sub("", text.strip())


# --- conversion ----------------------------------------------------------------


def fmt_num(x: float) -> str:
    return f"{x:.3g}"


def build_rows(domains: list[dict], fasta: list[tuple[str, str]] | None, taxon_id: str = "") -> list[dict]:
    """One output row per protein. Families per protein are ordered by the envelope start
    of their first domain (ties: accession). Best i-Evalue and best domain score per family."""
    per_prot: dict[str, dict[str, dict]] = {}
    desc_of: dict[str, str] = {}
    for d in domains:
        fam = per_prot.setdefault(d["protein"], {}).get(d["pfam"])
        if fam is None:
            per_prot[d["protein"]][d["pfam"]] = {"name": d["name"], "ievalue": d["ievalue"],
                                                 "score": d["score"], "start": d["env_from"]}
        else:
            fam["ievalue"] = min(fam["ievalue"], d["ievalue"])
            fam["score"] = max(fam["score"], d["score"])
            fam["start"] = min(fam["start"], d["env_from"])
        if d["desc"]:
            desc_of.setdefault(d["protein"], d["desc"])
    if fasta is not None:
        known = {tok for tok, _ in fasta}
        missing = sorted(set(per_prot) - known)
        if missing:
            raise SystemExit(f"ERROR: {len(missing)} domain-table proteins are not in the FASTA, "
                             f"e.g. {missing[:3]}")
        order = [(tok, rest) for tok, rest in fasta]
    else:
        order = [(p, desc_of.get(p, "")) for p in sorted(per_prot)]
    rows = []
    for prot, desc in order:
        fams = sorted(per_prot.get(prot, {}).items(), key=lambda kv: (kv[1]["start"], kv[0]))
        rows.append({
            "accession": prot, "taxon_id": taxon_id, "gene_name": "",
            "description": clean_description(desc), "go_ids": "",
            "pfam_ids": "|".join(acc for acc, _ in fams),
            "pfam_names": "|".join(v["name"] for _, v in fams),
            "interpro_ids": "", "ec_numbers": "", "alphafold_id": "", "xrefs": "",
            "pfam_best_ievalue": "|".join(fmt_num(v["ievalue"]) for _, v in fams),
            "pfam_best_domain_score": "|".join(f"{v['score']:.1f}" for _, v in fams),
        })
    return rows


def write_tsv(rows: list[dict], out: Path) -> None:
    """Write the table. A .gz name gives gzip with mtime 0, so the same input gives the same bytes."""
    out.parent.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    data = buf.getvalue().encode()
    if out.suffix == ".gz":
        with open(out, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz:
            gz.write(data)
    else:
        out.write_bytes(data)


def pfam_release_from_path(path: str) -> str:
    """'/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm' -> 'Pfam38.2 (from directory name
    2026-01-27-Pfam38.2)'. The HMM file itself has no release line."""
    if not path:
        return "not recorded"
    parent = Path(path).parent.name
    m = re.search(r"Pfam[-_]?(\d+(?:\.\d+)?)", parent, re.I)
    return f"Pfam {m.group(1)} (from directory name {parent})" if m else f"not recorded (directory {parent})"


# --- main ----------------------------------------------------------------------


def convert(domtbls: list[Path], out: Path, fasta: Path | None = None, tool: str = "auto",
            taxon_id: str = "", pfam_hmm: str = "", pfam_release: str = "",
            argv: list[str] | None = None) -> dict:
    """Convert, write the TSV and the sidecar. Return summary counts."""
    domains: list[dict] = []
    sources = []
    hmm_paths = []
    for p in domtbls:
        lines = read_text(p)
        trailer = parse_trailer(lines)
        layout = detect_layout(lines, trailer, tool)
        rows = parse_domtbl(lines, layout)
        domains.extend(rows)
        opts = trailer.get("option_settings", "")
        hmm_paths.append(trailer.get("query_file" if layout == "hmmsearch" else "target_file", ""))
        sources.append({"path": str(Path(p).resolve()), "sha256": sha256_of(Path(p)), "layout": layout,
                        "program": trailer.get("program", "not recorded"),
                        "hmmer_version": trailer.get("version", "not recorded"),
                        "hmmer_command": opts or "not recorded",
                        "cut_ga_in_command": ("--cut_ga" in shlex.split(opts)) if opts else None,
                        "run_date": trailer.get("date", "not recorded"),
                        "n_domain_rows": len(rows)})
    if not pfam_hmm:
        pfam_hmm = next((h for h in hmm_paths if h), "")
    fasta_rows = read_fasta_headers(fasta) if fasta else None
    rows = build_rows(domains, fasta_rows, taxon_id)
    write_tsv(rows, out)

    n_with = sum(1 for r in rows if r["pfam_ids"])
    versions = sorted({s["hmmer_version"] for s in sources})
    programs = sorted({s["program"] for s in sources})
    release = pfam_release or pfam_release_from_path(pfam_hmm)
    summary = {"n_rows": len(rows), "n_proteins_with_pfam": n_with,
               "n_domain_rows": len(domains), "n_pfam_families": len({d["pfam"] for d in domains}),
               "n_fasta_proteins": len(fasta_rows) if fasta_rows is not None else None}
    cmd = " ".join(shlex.quote(x) for x in ([sys.executable] + (argv if argv is not None else sys.argv)))
    record = build_record(
        source_url=("file://" + sources[0]["path"]) if len(sources) == 1
        else "; ".join("file://" + s["path"] for s in sources),
        source_release=f"{release}; HMMER {', '.join(versions)} ({', '.join(programs)})",
        license=PFAM_LICENSE,
        local_path=out.resolve(),
        derived_by=f"bin/hmm_domtbl_to_annotations.py: {cmd}",
        extra={
            "annotation_source": f"Pfam-A HMM search ({', '.join(programs)}); Pfam only",
            "vocabularies_available": ["pfam", "pfam_clan (derived from Pfam-A.hmm.dat by the reader)"],
            "vocabularies_not_available": ["go", "interpro"],
            "statement": NOT_AVAILABLE,
            "thresholds": "used as given in the domain table (no new threshold applied); see "
                          "hmmer_command (--cut_ga = Pfam gathering thresholds)",
            "pfam_hmm": pfam_hmm or "not recorded",
            "pfam_release": release,
            "source_domtbl": sources,
            "protein_fasta": ({"path": str(Path(fasta).resolve()), "sha256": sha256_of(Path(fasta))}
                              if fasta else None),
            "taxon_id": taxon_id or None,
            "columns": COLUMNS,
            "counts": summary,
        },
    )
    write_record(record, Path(str(out) + ".provenance.yaml"))
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domtbl", required=True, nargs="+", type=Path,
                    help="HMMER --domtblout file(s): plain, .gz or .zst")
    ap.add_argument("--out", type=Path,
                    help="output TSV (.tsv.gz writes gzip). Default (needs --study-dir and --short): "
                         "<study-dir>/annotations/<name>.tsv.gz, <name> from lib/annotation_names.py")
    ap.add_argument("--study-dir", type=Path,
                    help="study dir with species.csv; gives the default --out and the guard "
                         "(--out must then be inside <study-dir>/annotations/)")
    ap.add_argument("--short", help="species Short name in <study-dir>/species.csv (with --study-dir)")
    ap.add_argument("--fasta", type=Path, help="protein FASTA of the proteome (plain, .gz or .zst)")
    ap.add_argument("--tool", choices=["auto", "hmmsearch", "hmmscan"], default="auto",
                    help="table layout (default: detect from the '# Program:' trailer or the columns)")
    ap.add_argument("--taxon-id", default="", help="NCBI taxon id for the taxon_id column")
    ap.add_argument("--pfam-hmm", default="",
                    help="Pfam-A.hmm path for the provenance record (default: from the table trailer)")
    ap.add_argument("--pfam-release", default="",
                    help="Pfam release text (default: parsed from the Pfam-A.hmm directory name)")
    args = ap.parse_args(argv)
    if args.study_dir is not None:
        ann_dir = (args.study_dir / "annotations").resolve()
        if args.out is None:
            if not args.short:
                ap.error("--study-dir without --out needs --short")
            with open(args.study_dir / "species.csv", newline="") as fh:
                rows = list(csv.DictReader(fh))
            names = annotation_names(rows)
            if args.short not in names:
                ap.error(f"--short {args.short} is not in {args.study_dir / 'species.csv'}")
            args.out = annotation_path(ann_dir, names[args.short])
        try:
            ensure_within(args.out, [ann_dir])
        except ValueError as e:
            raise SystemExit(f"ERROR: {e}")
    elif args.out is None:
        ap.error("give --out, or --study-dir and --short")
    s = convert(args.domtbl, args.out, fasta=args.fasta, tool=args.tool, taxon_id=args.taxon_id,
                pfam_hmm=args.pfam_hmm, pfam_release=args.pfam_release, argv=argv)
    pct = f"{100 * s['n_proteins_with_pfam'] / s['n_rows']:.1f}%" if s["n_rows"] else "NA"
    print(f"wrote {args.out}: {s['n_rows']} rows, {s['n_proteins_with_pfam']} with >= 1 Pfam ({pct}), "
          f"{s['n_domain_rows']} domain rows, {s['n_pfam_families']} Pfam families", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
