#!/usr/bin/env python3
"""Stage D helpers: write the compact input for an extractor and score extractor passes.

  gf_stage_d.py input --out <growth_form dir> --species "Taphrina deformans" ...
      -> <out>/stage_d_input/<slug>.tsv  (sentence_id, match_level, year, section, text), one line per
         candidate sentence, for a model to read. Context sentences and article metadata are left out
         to keep the input short; `verify` still checks quotes against the full source.
  gf_stage_d.py compare --out <growth_form dir> --pass-tags haiku-v2 sonnet-v2
      -> <out>/compare/<tags>.tsv (per species: each pass's class and phase call, agreement, shared and
         one-pass-only sentences, needs_adjudication reasons) and <tags>.disagreements.tsv (the rows an
         adjudicator should read). Run `score` for each pass first; compare reads <slug>.verified.tsv.
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
    """Reduce verified species-level rows to yeast | filamentous | dimorphic | unsourced.

    Only rows the model labelled organism_level `species` and evidence_type `tested` or `described` count:
    a `predicted` row (an inference by sequence or analogy) is not an observation of growth form. Each
    row votes by `supports`. A `dimorphic` vote share of at least 20% gives dimorphic; so does a yeast
    share and a filamentous share that are each at least 20%. Otherwise the larger share wins. One stray
    row against a dozen no longer flips the call. A calibration device, not the stage F adjudicator.
    """
    votes = {"yeast": 0, "filamentous": 0, "dimorphic": 0}
    for r in rows:
        if r.get("organism_level", "").strip().lower() != "species":
            continue
        if r.get("evidence_type", "").strip().lower() not in {"tested", "described"}:
            continue
        s = r.get("supports", "").strip().lower()
        if s in votes:
            votes[s] += 1
    total = sum(votes.values())
    if total == 0:
        return "unsourced"
    share = {k: v / total for k, v in votes.items()}
    if share["dimorphic"] >= 0.2 or (share["yeast"] >= 0.2 and share["filamentous"] >= 0.2):
        return "dimorphic"
    return max(votes, key=votes.get)


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


COUNTED_EVIDENCE = {"tested", "described"}
THIN_ROWS = 2   # a pass with this many counted rows or fewer is thin evidence


def counted_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Rows that bear on the species' own growth form: organism_level species, evidence tested or described."""
    return [
        r for r in rows
        if r.get("organism_level", "").strip().lower() == "species"
        and r.get("evidence_type", "").strip().lower() in COUNTED_EVIDENCE
        and r.get("supports", "").strip().lower() in {"yeast", "filamentous", "dimorphic"}
    ]


def phase_call(rows: list[dict[str, str]]) -> str:
    """The normal phase the counted rows state: yeast, hyphal, unclear (no 60% majority) or none."""
    votes = {"yeast": 0, "hyphal": 0}
    for r in counted_rows(rows):
        ph = r.get("phase_stated", "").strip().lower()
        if ph in votes:
            votes[ph] += 1
    total = sum(votes.values())
    if total == 0:
        return "none"
    top = max(votes, key=votes.get)
    return top if votes[top] / total >= 0.6 else "unclear"


def _shares(rows: list[dict[str, str]]) -> dict[str, float]:
    c = counted_rows(rows)
    n = len(c)
    return {k: (sum(1 for r in c if r["supports"].strip().lower() == k) / n if n else 0.0)
            for k in ("yeast", "filamentous", "dimorphic")}


def compare_species(rows_by_tag: dict[str, list[dict[str, str]]]) -> dict:
    """Compare the verified rows of several passes for one species.

    Reports each pass's class and phase call, whether the passes agree, the sentences they share and
    those found by one pass only, and `needs_adjudication`, a list of reasons: class_disagree,
    phase_disagree, thin_evidence:<tag> (a pass has one or two counted rows), no_species_rows_in:<tag>
    (a pass has none while another has some). `mixed` says per pass whether it has both forms.
    """
    calls = {tag: class_call(rows) for tag, rows in rows_by_tag.items()}
    phases = {tag: phase_call(rows) for tag, rows in rows_by_tag.items()}
    sents = {tag: {r["sentence_id"] for r in counted_rows(rows)} for tag, rows in rows_by_tag.items()}
    shared = set.intersection(*sents.values()) if sents else set()
    reasons: list[str] = []
    class_agree = len(set(calls.values())) == 1
    if not class_agree:
        reasons.append("class_disagree")
    known = {p for p in phases.values() if p in {"yeast", "hyphal"}}
    if len(known) > 1:
        reasons.append("phase_disagree")
    # mixed_forms is information, not a disagreement: a dimorphic species has yeast and hyphal rows.
    mixed = {
        tag: bool(
            (sh := _shares(rows))["dimorphic"] >= 0.2
            or (sh["yeast"] >= 0.2 and sh["filamentous"] >= 0.2)
        )
        for tag, rows in rows_by_tag.items()
    }
    # One or two rows are too few for a stable call: runs of the same model flipped between none and
    # one row on such species.
    for tag, rows in rows_by_tag.items():
        if 1 <= len(counted_rows(rows)) <= THIN_ROWS:
            reasons.append(f"thin_evidence:{tag}")
    have = {tag for tag, rows in rows_by_tag.items() if counted_rows(rows)}
    if have and len(have) < len(rows_by_tag):
        reasons += [f"no_species_rows_in:{tag}" for tag in rows_by_tag if tag not in have]
    return {
        "calls": calls, "phases": phases, "class_agree": "yes" if class_agree else "no",
        "phase_agree": "yes" if len(known) <= 1 else "no",
        "shared_sentences": len(shared), "only": {tag: len(s - shared) for tag, s in sents.items()},
        "n_counted": {tag: len(counted_rows(rows)) for tag, rows in rows_by_tag.items()},
        "mixed": mixed,
        "needs_adjudication": ";".join(reasons) if reasons else "no",
    }


def disagreement_rows(rows_by_tag: dict[str, list[dict[str, str]]], cap: int = 40) -> list[dict[str, str]]:
    """The rows an adjudicator should read: counted rows in a minority form of their pass, and counted
    rows whose sentence only one pass found. At most `cap` rows."""
    sents = {tag: {r["sentence_id"] for r in counted_rows(rows)} for tag, rows in rows_by_tag.items()}
    out: list[dict[str, str]] = []
    for tag, rows in rows_by_tag.items():
        c = counted_rows(rows)
        if not c:
            continue
        sh = _shares(rows)
        major = max(sh, key=sh.get)
        others = set.union(*(s for t, s in sents.items() if t != tag)) if len(sents) > 1 else set()
        for r in c:
            minority = r["supports"].strip().lower() != major
            unique = r["sentence_id"] not in others if len(sents) > 1 else False
            if minority or unique:
                out.append({"pass_tag": tag, "sentence_id": r["sentence_id"], "supports": r["supports"],
                            "phase_stated": r.get("phase_stated", ""), "evidence_type": r["evidence_type"],
                            "why": ("minority_form" if minority else "") + (";" if minority and unique else "")
                                   + ("only_this_pass" if unique else ""),
                            "quote": r["quote"]})
    return out[:cap]


def cmd_compare(args: argparse.Namespace) -> None:
    out = Path(args.out)
    tags = args.pass_tags
    by_species: dict[str, dict[str, list[dict[str, str]]]] = {}
    names: dict[str, str] = {}
    for tag in tags:
        for path in sorted((out / "evidence" / tag).glob("*.verified.tsv")):
            slug = path.name[: -len(".verified.tsv")]
            rows = gf.read_tsv(path, required=False)
            by_species.setdefault(slug, {})[tag] = rows
            if rows:
                names[slug] = rows[0]["species"]
    dest = out / "compare"
    dest.mkdir(parents=True, exist_ok=True)
    label = "+".join(tags)
    table, detail = [], []
    for slug in sorted(by_species):
        rows_by_tag = {tag: by_species[slug].get(tag, []) for tag in tags}
        res = compare_species(rows_by_tag)
        sp = names.get(slug, slug)
        table.append({"species": sp, **{f"call:{t}": res["calls"][t] for t in tags},
                      **{f"phase:{t}": res["phases"][t] for t in tags},
                      **{f"rows:{t}": res["n_counted"][t] for t in tags},
                      **{f"mixed:{t}": "yes" if res["mixed"][t] else "no" for t in tags},
                      "class_agree": res["class_agree"], "phase_agree": res["phase_agree"],
                      "shared_sentences": res["shared_sentences"],
                      **{f"only:{t}": res["only"][t] for t in tags},
                      "needs_adjudication": res["needs_adjudication"]})
        if res["needs_adjudication"] != "no":
            for r in disagreement_rows(rows_by_tag):
                detail.append({"species": sp, **r})
    gf.write_tsv(dest / f"{label}.tsv", table, list(table[0].keys()))
    if detail:
        gf.write_tsv(dest / f"{label}.disagreements.tsv", detail, list(detail[0].keys()))
    n_adj = sum(1 for r in table if r["needs_adjudication"] != "no")
    print(f"compare {label}: {len(table)} species; class agree {sum(1 for r in table if r['class_agree'] == 'yes')}; "
          f"need adjudication {n_adj}; disagreement rows {len(detail)} -> {dest}")


def cmd_phase_input(args: argparse.Namespace) -> None:
    """Input for the phase pass: the counted rows of an extraction pass, quote only (short)."""
    out = Path(args.out)
    dest = out / "stage_d_phase_input"
    dest.mkdir(parents=True, exist_ok=True)
    for sp in args.species:
        slug = gf.slugify(sp)
        rows = counted_rows(gf.read_tsv(out / "evidence" / args.pass_tag / f"{slug}.verified.tsv", required=False))
        seen, keep = set(), []
        for r in rows:
            if r["sentence_id"] in seen:
                continue
            seen.add(r["sentence_id"])
            keep.append({"sentence_id": r["sentence_id"], "evidence_type": r["evidence_type"],
                         "supports": r["supports"], "quote": r["quote"].replace("\t", " ")})
        gf.write_tsv(dest / f"{slug}.tsv", keep, ["sentence_id", "evidence_type", "supports", "quote"])
        print(f"{sp}: {len(keep)} quotes, about {sum(len(r['quote']) for r in keep) // 4} tokens")


def phase_compare(results: dict[str, dict], valid_ids: set[str]) -> dict:
    """Compare phase calls of several passes for one species. `results` maps pass tag to the JSON object."""
    calls = {tag: r.get("dominant_phase_in_culture", "") for tag, r in results.items()}
    bad_basis = {tag: [b for b in r.get("basis", []) if b not in valid_ids] for tag, r in results.items()}
    filled = {tag: sum(1 for x in r.get("rows", []) if x.get("phase") in {"yeast", "hyphal", "both"})
              for tag, r in results.items()}
    return {"calls": calls, "agree": "yes" if len(set(calls.values())) == 1 else "no",
            "invalid_basis": {t: b for t, b in bad_basis.items() if b}, "phase_stated_rows": filled}


def cmd_phase_compare(args: argparse.Namespace) -> None:
    out = Path(args.out)
    table = []
    for sp in args.species:
        slug = gf.slugify(sp)
        valid = {r["sentence_id"]: r for r in gf.read_tsv(out / "stage_d_phase_input" / f"{slug}.tsv", required=False)}
        results = {}
        for tag in args.pass_tags:
            path = out / "evidence_phase" / tag / f"{slug}.json"
            if path.exists():
                results[tag] = json.loads(path.read_text())
        if not results:
            continue
        res = phase_compare(results, set(valid))
        table.append({"species": sp, **{f"call:{t}": res["calls"].get(t, "") for t in args.pass_tags},
                      **{f"conf:{t}": results[t].get("confidence", "") if t in results else "" for t in args.pass_tags},
                      **{f"phase_rows:{t}": res["phase_stated_rows"].get(t, "") for t in args.pass_tags},
                      "agree": res["agree"], "invalid_basis": json.dumps(res["invalid_basis"])})
        for tag, r in results.items():
            for sid in r.get("basis", []):
                if sid in valid:
                    print(f"  [{sp[:14]}|{tag}|{r.get('dominant_phase_in_culture')}] {valid[sid]['quote'][:150]}")
    gf.write_tsv(out / "phase_compare.tsv", table, list(table[0].keys()))
    print(f"phase_compare: {len(table)} species, agree {sum(1 for r in table if r['agree'] == 'yes')}")


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
    p = sub.add_parser("phase-input", help="write the quotes of an extraction pass as input for the phase pass")
    p.add_argument("--out", required=True)
    p.add_argument("--pass-tag", required=True)
    p.add_argument("--species", nargs="+", required=True)
    p.set_defaults(func=cmd_phase_input)
    p = sub.add_parser("phase-compare", help="compare phase calls of several phase passes")
    p.add_argument("--out", required=True)
    p.add_argument("--pass-tags", nargs="+", required=True)
    p.add_argument("--species", nargs="+", required=True)
    p.set_defaults(func=cmd_phase_compare)
    p = sub.add_parser("compare", help="compare verified rows of two or more passes per species")
    p.add_argument("--out", required=True)
    p.add_argument("--pass-tags", nargs="+", required=True)
    p.set_defaults(func=cmd_compare)
    args = ap.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
