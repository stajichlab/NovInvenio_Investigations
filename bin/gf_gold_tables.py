#!/usr/bin/env python3
"""Build independent gold rows from tables of review and comparative papers.

The 19-row gold set in growth_form/gold_sources.tsv came from the sources a pilot search found, so a
method that copies the pilot scores well on it. This script builds a second gold set by a different
route: it searches PMC for review and comparative papers on dimorphism and growth form, reads the
tables of their full text, and keeps every table row that names one of the study species together
with a growth term (in the row or in the table caption). Each (species, paper) pair is a gold row at
the article level: the pipeline should retrieve the paper and select it. These rows have no phrase: they measure
whether the pipeline retrieves and selects the paper. The recall command reports sentence and phrase
as n/a for them.

Usage: pixi run python bin/gf_gold_tables.py --study-dir studies/fungi/yeast_convergence \\
           [--out <growth_form dir>] [--retmax 100]
Each gold row has a `strength`: `table` (a table row names the species and a growth-form term), `about`
(the species' full name is in the paper's title or abstract and a body sentence has a growth-form term),
or `passing` (only a body sentence mentions the species). A recall below 100% on `passing` rows is not a
failure: a review that names 30 species once is not a source the pipeline should select for each.
Output (in <out>): gold_tables_candidates.tsv (one row per matching table row or body sentence),
gold_sources_independent.tsv (one row per species and paper; same columns as gold_sources.tsv),
gold_tables_summary.tsv (queries, papers, rows). Papers already in gold_sources.tsv are excluded.
Stdlib only (plus gf_lit.py).
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gf_lit as gf  # noqa: E402

QUERIES = {
    "dimorphism_review": '("dimorphic fungi"[Title/Abstract] OR "fungal dimorphism"[Title/Abstract] '
    'OR "yeast-to-hypha"[Title/Abstract] OR "dimorphic transition"[Title/Abstract]) AND review[pt]',
    "dimorphic_comparative": '("dimorphic"[Title/Abstract] AND ("comparative genomic*"[Title/Abstract] '
    'OR "phylogen*"[Title/Abstract] OR "evolution of"[Title/Abstract]) AND (fungi[Title/Abstract] '
    'OR fungal[Title/Abstract]))',
    "yeast_origins": '("yeast"[Title/Abstract] AND ("convergent"[Title/Abstract] OR "multiple origins"'
    '[Title/Abstract] OR "independent origin*"[Title/Abstract] OR "evolution of yeast"[Title/Abstract]) '
    'AND fungi[Title/Abstract])',
    "growth_form_evolution": '("growth form"[Title/Abstract] OR "growth forms"[Title/Abstract] OR '
    '"cell morphology"[Title/Abstract]) AND fung*[Title/Abstract] AND ("evolution"[Title/Abstract] '
    'OR "traits"[Title/Abstract])',
    "mucorales_dimorphism": '(Mucorales[Title/Abstract] OR Mucor[Title/Abstract]) AND '
    '(dimorphism[Title/Abstract] OR dimorphic[Title/Abstract])',
    "black_yeasts": '("black yeasts"[Title/Abstract] OR Chaetothyriales[Title/Abstract]) AND '
    '(morpholog*[Title/Abstract] OR dimorph*[Title/Abstract])',
    "basidiomycete_yeasts": '(Basidiomycota[Title/Abstract] AND yeast[Title/Abstract] AND '
    '(dimorphic[Title/Abstract] OR "life cycle"[Title/Abstract] OR morphology[Title/Abstract]))',
}
# A table is about growth form only if its caption says so ...
CAPTION_TERM_RE = re.compile(
    r"dimorph|growth form|morpholog|yeast[- ]form|yeast[- ]like|hyphal form|life cycle|lifestyle|pleomorph", re.I
)
# ... and a row counts if the caption does, or the row itself names a growth form.
ROW_TERM_RE = re.compile(
    r"dimorph|yeast[- ]like|yeast form|yeast phase|hyphal|pseudohyph|filamentous|mycelial|pleomorph", re.I
)
# Tables of drug tests, strains, sequences and clinical series are not growth-form tables, whatever
# else their captions say.
NOISE_CAPTION_RE = re.compile(
    r"\bMICs?\b|\bMFCs?\b|susceptib|antifungal|resistan|primer|accession|genbank|sequence|"
    r"strain|isolate|patient|clinical|inhibit|toxic|virulence|plasmid|oligonucleotide|\bdrug|"
    r"treatment|\bdose\b|worm|mice|mouse",
    re.I,
)
# Body sentences: species and a specific growth-form term in one sentence.
SENT_TERM_RE = re.compile(
    gf.DIMORPH_RE.pattern + r"|yeast[- ]like|yeast form|yeast phase|hyphal form|pseudohyph\w*|"
    r"filamentous form|yeast and hyph\w*|hyph\w* and yeast",
    re.I,
)
MAX_ROW_WORDS = 40


def local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def text_of(el: ET.Element) -> str:
    return re.sub(r"\s+", " ", " ".join(el.itertext())).strip()


def parse_tables(xml_bytes: bytes) -> tuple[dict[str, str], list[dict]]:
    """Return (article ids and title, tables). A table is {label, caption, rows: [[cell, ...]]}."""
    root = ET.fromstring(xml_bytes)
    meta = {"pmcid": "", "pmid": "", "title": ""}
    for el in root.iter():
        name = local(el.tag)
        if name == "article-id":
            kind = el.get("pub-id-type", "")
            if kind in {"pmc", "pmcid"} and not meta["pmcid"]:
                meta["pmcid"] = gf.normalize_pmcid((el.text or "").strip())
            elif kind == "pmid" and not meta["pmid"]:
                meta["pmid"] = (el.text or "").strip()
        elif name == "article-title" and not meta["title"]:
            meta["title"] = text_of(el)
    tables = []
    for tw in root.iter():
        if local(tw.tag) != "table-wrap":
            continue
        label = caption = ""
        rows: list[list[str]] = []
        for el in tw.iter():
            n = local(el.tag)
            if n == "label" and not label:
                label = text_of(el)
            elif n == "caption" and not caption:
                caption = text_of(el)
            elif n == "tr":
                cells = [text_of(c) for c in el if local(c.tag) in {"td", "th"}]
                if cells:
                    rows.append(cells)
        tables.append({"label": label, "caption": caption, "rows": rows})
    return meta, tables


def match_rows(tables: list[dict], patterns: dict[str, list[re.Pattern[str]]]) -> list[dict]:
    """Table rows that name a species and carry a growth term in the row or the table caption."""
    hits = []
    for t in tables:
        cap = f"{t['label']} {t['caption']}"
        if NOISE_CAPTION_RE.search(cap):
            continue
        caption_hit = bool(CAPTION_TERM_RE.search(cap))
        for cells in t["rows"]:
            row = " | ".join(cells)
            for species, pats in patterns.items():
                if not any(p.search(row) for p in pats):
                    continue
                if not (caption_hit or ROW_TERM_RE.search(row)):
                    continue
                words = row.split()
                hits.append(
                    {
                        "species": species,
                        "table": t["label"],
                        "caption": t["caption"][:160],
                        "row_text": " ".join(words[:MAX_ROW_WORDS]),
                    }
                )
    return hits


def sentence_hits(
    sections: list[tuple[str, str]], patterns: dict[str, list[re.Pattern[str]]]
) -> list[dict]:
    """Body sentences that name a species and a specific growth-form term (dimorphic, yeast-like,
    hyphal form, ...). Gold at the article level: they show that the paper discusses the species'
    growth form. They are found by a rule that is independent of the pipeline's search."""
    hits = []
    for section, text in sections:
        for sentence in gf.split_sentences(text):
            if not SENT_TERM_RE.search(sentence):
                continue
            for species, pats in patterns.items():
                if any(p.search(sentence) for p in pats):
                    hits.append(
                        {
                            "species": species,
                            "table": f"text: {section[:40]}",
                            "caption": "",
                            "row_text": " ".join(sentence.split()[:MAX_ROW_WORDS]),
                        }
                    )
    return hits


def species_patterns(names_rows: list[dict[str, str]]) -> dict[str, list[re.Pattern[str]]]:
    out: dict[str, list[re.Pattern[str]]] = {}
    for r in names_rows:
        if r["name_tier"] == "genus" or len(r["name"]) < 6:
            continue
        # Full names only. An abbreviation ("C. cinerea") can belong to another genus in the same
        # paper, and a gold row must be right more often than the pipeline it tests.
        out.setdefault(r["species"], []).append(gf.name_pattern(r["name"]))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study-dir", default=str(gf.DEFAULT_STUDY))
    ap.add_argument("--out", default="", help="growth_form dir (default <study-dir>/growth_form)")
    ap.add_argument("--retmax", type=int, default=100)
    ap.add_argument("--max-per-species", type=int, default=8,
                    help="keep at most this many papers per species (papers with a table row first), "
                         "so that a heavily studied species does not dominate")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""))
    args = ap.parse_args(argv)
    out = Path(args.out) if args.out else Path(args.study_dir) / "growth_form"
    cache = out / "gold_cache"
    cache.mkdir(parents=True, exist_ok=True)
    client = gf.HttpClient(
        api_key=os.environ.get("NCBI_API_KEY") or None, email=args.email or None
    )
    pilot = set()
    gold_path = out / "gold_sources.tsv"
    if gold_path.exists():
        for r in gf.read_tsv(gold_path):
            if r.get("pmcid"):
                pilot.add(gf.normalize_pmcid(r["pmcid"]))

    theme_of: dict[str, list[str]] = {}
    summary: list[dict] = []
    for theme, term in QUERIES.items():
        res = client.esearch("pmc", term, args.retmax, sort="relevance")
        for uid in res["ids"]:
            theme_of.setdefault(uid, []).append(theme)
        summary.append({"item": f"query:{theme}", "value": f"count={res['count']} returned={len(res['ids'])}"})
    ids = sorted(u for u in theme_of if u not in pilot)
    gf.fetch_and_cache(client, cache, "pmc", ids)

    patterns = species_patterns(gf.read_tsv(out / "names.tsv"))
    cand: list[dict] = []
    about: set[tuple[str, str]] = set()          # (species, pmcid): full name in the title or abstract
    n_parsed = n_tables = 0
    for uid in ids:
        xml_path, _meta = gf.cache_paths(cache, "pmc", uid)
        if not xml_path.exists():
            continue
        try:
            meta, tables = parse_tables(gf.read_zst(xml_path))
        except ET.ParseError:
            continue
        n_parsed += 1
        n_tables += len(tables)
        parsed = gf.parse_pmc_articles(gf.read_zst(xml_path)).get(uid) or {}
        head = f"{parsed.get('title', '')} {parsed.get('abstract', '')}"
        for sp, pats in patterns.items():
            if any(p.search(head) for p in pats):
                about.add((sp, meta["pmcid"] or uid))
        rows = [("table", h) for h in match_rows(tables, patterns)] + [
            ("text", h) for h in sentence_hits(parsed.get("sections") or [], patterns)
        ]
        for kind, h in rows:
            cand.append(
                {**h, "kind": kind, "pmcid": meta["pmcid"] or uid, "pmid": meta["pmid"],
                 "title": meta["title"], "themes": ";".join(theme_of[uid])}
            )
    cand.sort(key=lambda r: (r["species"], r["pmcid"], r["kind"] != "table", r["table"], r["row_text"]))
    gf.write_tsv(out / "gold_tables_candidates.tsv", cand,
                 ["species", "kind", "pmcid", "pmid", "title", "table", "caption", "row_text", "themes"])

    gold: dict[tuple[str, str], dict] = {}
    table_keys0 = {(r["species"], r["pmcid"]) for r in cand if r["kind"] == "table"}
    for r in cand:
        key = (r["species"], r["pmcid"])
        if key not in gold:
            gold[key] = {"species": r["species"], "level": "species", "pmid": r["pmid"], "pmcid": r["pmcid"],
                         "phrase": "", "source": r["title"][:90],
                         "strength": "table" if key in table_keys0 else
                         ("about" if key in about else "passing"),
                         "note": f"{r['kind']} ({r['table']}): {r['row_text'][:140]}"}
    table_keys = {(r["species"], r["pmcid"]) for r in cand if r["kind"] == "table"}
    n_hits = {}
    for r in cand:
        n_hits[(r["species"], r["pmcid"])] = n_hits.get((r["species"], r["pmcid"]), 0) + 1
    per_species: dict[str, list[dict]] = {}
    for key, row in gold.items():
        per_species.setdefault(key[0], []).append((key, row))
    gold_rows = []
    for sp, items in per_species.items():
        order = {"table": 0, "about": 1, "passing": 2}
        items.sort(key=lambda kr: (order[kr[1]["strength"]], -n_hits[kr[0]], kr[0][1]))
        gold_rows += [row for _k, row in items[: args.max_per_species]]
    gold_rows.sort(key=lambda r: (r["species"], r["pmcid"]))
    gf.write_tsv(out / "gold_sources_independent.tsv", gold_rows,
                 ["species", "level", "pmid", "pmcid", "phrase", "source", "strength", "note"])
    species_all = {r["species"] for r in gf.read_tsv(out / "species.tsv")}
    covered = {r["species"] for r in gold_rows}
    summary += [
        {"item": "papers_searched", "value": len(theme_of)},
        {"item": "papers_excluded_already_in_pilot_gold", "value": len(set(theme_of) & pilot)},
        {"item": "papers_parsed", "value": n_parsed},
        {"item": "tables", "value": n_tables},
        {"item": "matching_table_rows", "value": sum(1 for r in cand if r["kind"] == "table")},
        {"item": "matching_body_sentences", "value": sum(1 for r in cand if r["kind"] == "text")},
        {"item": "gold_rows_species_paper", "value": len(gold_rows)},
        {"item": "species_with_gold", "value": f"{len(covered)} of {len(species_all)}"},
    ]
    gf.write_tsv(out / "gold_tables_summary.tsv", summary, ["item", "value"])
    for s in summary:
        print(f"{s['item']}: {s['value']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
