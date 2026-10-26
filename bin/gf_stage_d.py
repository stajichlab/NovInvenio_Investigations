#!/usr/bin/env python3
"""Stage D helpers: write the compact input for an extractor and score extractor passes.

  gf_stage_d.py input --out <growth_form dir> --species "Taphrina deformans" ...
      -> <out>/stage_d_input/<slug>.tsv  (sentence_id, match_level, year, section, text), one line per
         candidate sentence, for a model to read. Context sentences and article metadata are left out
         to keep the input short; `verify` still checks quotes against the full source.
  gf_stage_d.py score --out <growth_form dir> --pass-tag haiku-v1 [--reference <file>] [--gold <file>]
      -> reads <out>/evidence/<pass_tag>/<slug>.jsonl, runs `gf_lit.py verify` on each, and writes
         <out>/evidence/<pass_tag>/calibration_rows.tsv (one row per species) and appends one summary
         row to <out>/calibration.tsv.

The class call used for scoring is deliberately crude: from the verified rows that the model labels
`organism_level = species`, collect the `supports` values and reduce them (section 14.5 of the plan).
It is a calibration device, not the adjudicator (stage F).
Stdlib only (plus gf_lit.py).
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gf_lit as gf  # noqa: E402


def compact_rows(records: list[dict]) -> list[dict[str, str]]:
    """One input line per candidate sentence; species-level matches first, then article, then genus."""
    order = {"full": 0, "abbrev": 0, "article": 1, "genus": 2}
    rows = []
    for r in sorted(records, key=lambda r: (order.get(r["match_level"], 3), r["uid"], r["sentence_id"])):
        rows.append(
            {
                "sentence_id": r["sentence_id"],
                "match_level": r["match_level"],
                "year": str(r.get("year", "")),
                "section": (r.get("section", "") or "").replace("\t", " ")[:60],
                "text": (r.get("text", "") or "").replace("\t", " ").replace("\n", " "),
            }
        )
    return rows


def class_call(rows: list[dict[str, str]]) -> str:
    """Reduce verified species-level rows to yeast | filamentous | dimorphic | unsourced."""
    sup = {
        r.get("supports", "").strip().lower()
        for r in rows
        if r.get("organism_level", "").strip().lower() == "species"
    }
    sup &= {"yeast", "filamentous", "dimorphic"}
    if not sup:
        return "unsourced"
    if "dimorphic" in sup or {"yeast", "filamentous"} <= sup:
        return "dimorphic"
    return next(iter(sup))


PART_RE = re.compile(r"\.part\d+$")


def split_chunks(rows: list[dict[str, str]], size: int) -> list[list[dict[str, str]]]:
    """Cut the rows into chunks of at most `size` lines. A model that reads a long file in pages
    may stop early (a pilot run read 350 of 775 lines), so each call gets one short chunk."""
    if size <= 0 or len(rows) <= size:
        return [rows]
    return [rows[i : i + size] for i in range(0, len(rows), size)]


def base_slug(stem: str) -> str:
    """'mucor_circinelloides.part3' -> 'mucor_circinelloides'."""
    return PART_RE.sub("", stem)


def cmd_input(args: argparse.Namespace) -> None:
    out = Path(args.out)
    dest = out / "stage_d_input"
    dest.mkdir(parents=True, exist_ok=True)
    fields = ["sentence_id", "match_level", "year", "section", "text"]
    for sp in args.species:
        slug = gf.slugify(sp)
        path = out / "sentences" / f"{slug}.jsonl"
        recs = [json.loads(x) for x in open(path) if x.strip()]
        rows = compact_rows(recs)
        chunks = split_chunks(rows, args.chunk_size)
        names = []
        for i, chunk in enumerate(chunks, 1):
            name = slug if len(chunks) == 1 else f"{slug}.part{i}"
            gf.write_tsv(dest / f"{name}.tsv", chunk, fields)
            names.append(name)
        tokens = sum(len(r["text"]) for r in rows) // 4
        print(f"{sp}: {len(rows)} sentences, about {tokens} tokens, {len(chunks)} chunk(s): {', '.join(names)}")


def cmd_score(args: argparse.Namespace) -> None:
    out = Path(args.out)
    tag = args.pass_tag
    pass_dir = out / "evidence" / tag
    reference = {r["species"]: r for r in gf.read_tsv(Path(args.reference))} if args.reference else {}
    gold = gf.read_tsv(Path(args.gold)) if args.gold else []
    rows_out = []
    by_slug: dict[str, list[Path]] = {}
    for path in sorted(pass_dir.glob("*.jsonl")):
        if path.name.endswith(".merged.jsonl"):
            continue
        by_slug.setdefault(base_slug(path.stem), []).append(path)
    for slug, paths in sorted(by_slug.items()):
        evidence = [r for path in paths for r in gf.read_evidence(path)]
        species = evidence[0]["species"] if evidence else slug
        prefix = slug
        merged = pass_dir / f"{slug}.merged.jsonl"
        with open(merged, "w") as fh:
            for r in evidence:
                fh.write(json.dumps(r, ensure_ascii=False) + "\n")
        # verify through the same code path as the CLI
        ns = argparse.Namespace(study_dir=args.study_dir, out=str(out), evidence=str(merged),
                                out_prefix=prefix, verified_dir=str(pass_dir))
        gf.cmd_verify(ns)
        ver = gf.read_tsv(pass_dir / f"{prefix}.verified.tsv", required=False)
        rej = gf.read_tsv(pass_dir / f"{prefix}.rejected.tsv", required=False)
        reasons: dict[str, int] = {}
        for r in rej:
            for reason in r["rejection_reason"].split(";"):
                reasons[reason.split("(")[0]] = reasons.get(reason.split("(")[0], 0) + 1
        call = class_call(ver)
        ref = reference.get(species, {})
        accepted = {x.strip() for x in ref.get("accept", ref.get("class", "")).split("|") if x.strip()}
        found = {(g["pmid"], g["pmcid"]) for g in gold if g["species"] == species}
        got = {(r.get("pmid", ""), r.get("pmcid", "")) for r in ver}
        gold_found = sum(1 for (pm, pc) in found if any((pm and pm == g[0]) or (pc and pc == g[1]) for g in got))
        rows_out.append(
            {
                "pass_tag": tag, "species": species, "n_rows": len(evidence), "n_verified": len(ver),
                "n_rejected": len(rej), "rejection_reasons": ";".join(f"{k}:{v}" for k, v in sorted(reasons.items())),
                "class_call": call, "reference": ref.get("class", ""), "match": ("yes" if call in accepted else "no") if ref else "",
                "gold_sources": len(found), "gold_sources_in_verified": gold_found,
            }
        )
    gf.write_tsv(pass_dir / "calibration_rows.tsv", rows_out, list(rows_out[0].keys()))
    n = len(rows_out)
    tot = sum(r["n_rows"] for r in rows_out)
    ver = sum(r["n_verified"] for r in rows_out)
    summ = {
        "pass_tag": tag, "species": n, "rows": tot, "verified_rate": f"{ver / tot:.3f}" if tot else "",
        "class_match": f"{sum(1 for r in rows_out if r['match'] == 'yes')}/{sum(1 for r in rows_out if r['match'])}",
        "gold_sources_found": f"{sum(r['gold_sources_in_verified'] for r in rows_out)}/{sum(r['gold_sources'] for r in rows_out)}",
        "date": gf.utc_now(),
    }
    cal = out / "calibration.tsv"
    existing = gf.read_tsv(cal, required=False) if cal.exists() else []
    existing = [r for r in existing if r["pass_tag"] != tag] + [summ]
    gf.write_tsv(cal, existing, list(summ.keys()))
    print(json.dumps(summ))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("input")
    p.add_argument("--out", required=True)
    p.add_argument("--species", nargs="+", required=True)
    p.add_argument("--chunk-size", type=int, default=120,
                   help="sentences per chunk file (0 = one file); default 120 (about 6,000 tokens)")
    p.set_defaults(func=cmd_input)
    p = sub.add_parser("score")
    p.add_argument("--out", required=True)
    p.add_argument("--study-dir", default=str(gf.DEFAULT_STUDY))
    p.add_argument("--pass-tag", required=True)
    p.add_argument("--reference", default="")
    p.add_argument("--gold", default="")
    p.set_defaults(func=cmd_score)
    args = ap.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
