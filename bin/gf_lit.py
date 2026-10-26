#!/usr/bin/env python3
"""Literature retrieval and quote verification for growth-form coding.

This script implements stages A-C of GROWTH_FORM_LITERATURE_PLAN.md:

  species   build growth_form/species.tsv from the study's species.csv files
            and batch YAML codings (one row per species; resolves strain-level
            Taxon_IDs to species nodes with the local NCBI taxonomy dump)
  names     build growth_form/names.tsv (canonical, synonym, equivalent,
            GenBank common names plus curated literature aliases)
  search    run PubMed and PMC E-utilities queries per name; write
            growth_form/queries.tsv and growth_form/candidates.tsv
  fetch     resolve PMID<->PMCID, fetch abstracts and PMC full text with a
            per-article cache; write growth_form/articles.tsv
  sentences extract growth-form sentences from fetched text, excluding
            reference lists and back matter; write
            growth_form/sentences/<species>.jsonl and sentences_summary.tsv
  verify    check evidence rows against cached source text and write
            evidence.verified.tsv / evidence.rejected.tsv
  gold      check gold.tsv against retrieved sources and sentences

The script uses only the Python standard library plus PyYAML.  Compressed
caches use the zstd command-line tool (subprocess), not a Python module.

Examples:
  bin/gf_lit.py species --study-dir studies/fungi/yeast_convergence
  bin/gf_lit.py names   --study-dir studies/fungi/yeast_convergence
  bin/gf_lit.py search  --study-dir studies/fungi/yeast_convergence --max-per-query 200
  bin/gf_lit.py fetch   --study-dir studies/fungi/yeast_convergence \\
      --max-articles 60 --max-fulltext 30
  bin/gf_lit.py sentences --study-dir studies/fungi/yeast_convergence
  bin/gf_lit.py gold    --study-dir studies/fungi/yeast_convergence
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from provenance import build_record, sha256_of  # noqa: E402

NII_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_STUDY = NII_ROOT / "studies" / "fungi" / "yeast_convergence"
DEFAULT_TAXONOMY = Path("/srv/projects/db/ncbi/preformatted/20260128")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
IDCONV = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"
TOOL = "nii_gf_lit"
SCRIPT_VERSION = "gf_lit 2026-10-25"

NAME_TIER_ORDER = {
    "canonical": 0,
    "canonical_input": 1,
    "ncbi_synonym": 2,
    "genbank_common": 3,
    "literature_alias": 4,
    "genus": 5,
}
NAME_CLASSES = (
    "scientific name",
    "synonym",
    "equivalent name",
    "genbank common name",
)

GROWTH_TERMS = (
    "yeast",
    "yeasts",
    "yeast-like",
    "hypha",
    "hyphae",
    "hyphal",
    "filamentous",
    "filament",
    "dimorphic",
    "dimorphism",
    "pseudohypha",
    "pseudohyphae",
    "mycelium",
    "mycelial",
    "mold",
    "mould",
    "blastoconidia",
    "arthroconidia",
    "chlamydospore",
    "unicellular",
    "single cell",
    "single-cell",
    "budding",
    "germ tube",
    "germ-tube",
)

# Terms for the question itself (a switch between yeast and hyphal growth). Searched in a separate,
# top-priority query because the generic growth terms match thousands of papers for well-studied
# species and the relevance-sorted top results then miss the dimorphism papers.
DIMORPH_TERMS = (
    "dimorphic",
    "dimorphism",
    "yeast-to-hypha",
    "yeast to hypha",
    "yeast-hypha",
    "hypha-to-yeast",
    "yeast-to-mycelium",
    "yeast-mycelium",
    "hyphal transition",
    "polymorphic",
    "pseudohyphal growth",
    "filamentous growth",
    "morphological transition",
)
# Lower value = used first when articles are selected.
QUERY_PRIORITY = {"dimorph": 0, "growth": 1, "name": 2}

GROWTH_RE = re.compile(
    r"\b("
    r"yeasts?|yeast[- ]like|hyph\w*|filament\w*|dimorph\w*|pseudohyph\w*|"
    r"mycel\w*|molds?|moulds?|blastoconidi\w*|arthroconidi\w*|chlamydospor\w*|"
    r"unicellular|single[- ]cell\w*|budding|germ[- ]tube\w*|sporulat\w*"
    r")\b",
    re.I,
)
DIMORPH_RE = re.compile(
    "|".join(re.escape(term) for term in DIMORPH_TERMS), re.I
)

ABBREV_RE = re.compile(
    r"\b(?:"
    r"et al|Fig|Figs|Tab|Ref|Refs|No|vs|cf|approx|ca|e\.g|i\.e|"
    r"sp|spp|var|subsp|f|pv|ser|strain|isolate|"
    r"ATCC|CBS|NRRL|DSM|JCM|IFO|IMI|UAMH|MB|GenBank|NCBI|"
    r"DNA|RNA|rRNA|ITS|LSU|SSU|ORF|kDa|bp|kb|mya|Ma|min|max"
    r")\.",
    re.I,
)

SENT_SPLIT_RE = re.compile(r"(?<=[.!?])\s+(?=[A-Z0-9(\[\"'])")
SENT_PROTECT = "\u00a7"


# ---------------------------------------------------------------------------
# small utilities
# ---------------------------------------------------------------------------


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def norm_space(text: str | None) -> str:
    return re.sub(r"\s+", " ", text or "").strip()


def norm_quote(text: str | None) -> str:
    t = unicodedata.normalize("NFKC", text or "")
    for a, b in (
        ("\u2010", "-"),
        ("\u2011", "-"),
        ("\u2012", "-"),
        ("\u2013", "-"),
        ("\u2014", "-"),
        ("\u2018", "'"),
        ("\u2019", "'"),
        ("\u201c", '"'),
        ("\u201d", '"'),
    ):
        t = t.replace(a, b)
    return norm_space(t).lower()


def slugify(text: str) -> str:
    s = unicodedata.normalize("NFKD", text)
    s = s.encode("ascii", "ignore").decode("ascii").lower()
    s = re.sub(r"[^a-z0-9]+", "_", s).strip("_")
    return s or "x"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_text(text: str) -> str:
    return sha256_bytes(text.encode("utf-8"))


def chunks(items: list[Any], size: int) -> Iterable[list[Any]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_tsv(path: Path, required: bool = True) -> list[dict[str, str]]:
    if not path.exists():
        if required:
            raise SystemExit(f"missing file: {path}")
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_zst(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        ["zstd", "-q", "-f", "-o", str(path), "-"],
        input=data,
        check=False,
        capture_output=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"zstd failed for {path}: {proc.stderr.decode()}")


def read_zst(path: Path) -> bytes:
    proc = subprocess.run(
        ["zstd", "-dc", str(path)],
        check=False,
        capture_output=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"zstd failed for {path}: {proc.stderr.decode()}")
    return proc.stdout


def elem_text(el: ET.Element | None) -> str:
    if el is None:
        return ""
    return norm_space("".join(el.itertext()))


def localname(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def normalize_pmcid(value: str) -> str:
    value = norm_space(value)
    if value.upper().startswith("PMC"):
        value = value[3:]
    return value


def word_count(text: str) -> int:
    return len(re.findall(r"\b[\w'-]+\b", text or ""))


# ---------------------------------------------------------------------------
# HTTP with rate limiting, retries, and on-disk cache
# ---------------------------------------------------------------------------


class HttpClient:
    def __init__(
        self,
        api_key: str | None = None,
        email: str | None = None,
        tool: str = TOOL,
    ) -> None:
        self.api_key = api_key
        self.email = email
        self.tool = tool
        self.min_interval = 0.11 if api_key else 0.35
        self._last = 0.0
        self.calls = 0

    def _wait(self) -> None:
        delta = time.monotonic() - self._last
        if delta < self.min_interval:
            time.sleep(self.min_interval - delta)
        self._last = time.monotonic()

    def _params(self) -> dict[str, str]:
        p = {"tool": self.tool}
        if self.api_key:
            p["api_key"] = self.api_key
        if self.email:
            p["email"] = self.email
        return p

    def get(
        self,
        url: str,
        params: dict[str, str] | list[tuple[str, str]] | None = None,
        *,
        cache_path: Path | None = None,
        retries: int = 4,
    ) -> bytes:
        if cache_path is not None and cache_path.exists():
            return cache_path.read_bytes()
        if isinstance(params, list):
            query = urllib.parse.urlencode(params)
        else:
            query = urllib.parse.urlencode(params or {})
        full = f"{url}?{query}" if query else url
        last: Exception | None = None
        for attempt in range(retries):
            self._wait()
            try:
                req = urllib.request.Request(
                    full,
                    headers={"User-Agent": f"{self.tool}/1.0 (growth-form curation)"},
                )
                with urllib.request.urlopen(req, timeout=90) as resp:
                    data = resp.read()
                self.calls += 1
                if cache_path is not None:
                    cache_path.parent.mkdir(parents=True, exist_ok=True)
                    cache_path.write_bytes(data)
                return data
            except urllib.error.HTTPError as exc:
                last = exc
                if exc.code in (429, 500, 502, 503, 504):
                    time.sleep(min(60, 2**attempt))
                    continue
                raise
            except (urllib.error.URLError, TimeoutError) as exc:
                last = exc
                time.sleep(min(60, 2**attempt))
        raise RuntimeError(f"request failed after {retries}: {full}: {last}")

    def post(
        self,
        url: str,
        params: dict[str, str],
        *,
        retries: int = 4,
    ) -> bytes:
        body = urllib.parse.urlencode(params).encode("utf-8")
        last: Exception | None = None
        for attempt in range(retries):
            self._wait()
            try:
                req = urllib.request.Request(
                    url,
                    data=body,
                    headers={
                        "User-Agent": f"{self.tool}/1.0 (growth-form curation)",
                        "Content-Type": "application/x-www-form-urlencoded",
                    },
                    method="POST",
                )
                with urllib.request.urlopen(req, timeout=180) as resp:
                    data = resp.read()
                self.calls += 1
                return data
            except urllib.error.HTTPError as exc:
                last = exc
                if exc.code in (429, 500, 502, 503, 504):
                    time.sleep(min(60, 2**attempt))
                    continue
                raise
            except (urllib.error.URLError, TimeoutError) as exc:
                last = exc
                time.sleep(min(60, 2**attempt))
        raise RuntimeError(f"request failed after {retries}: {url}: {last}")

    def esearch(
        self,
        db: str,
        term: str,
        retmax: int,
        *,
        sort: str | None = "relevance",
        cache_path: Path | None = None,
    ) -> dict[str, Any]:
        params = {
            "db": db,
            "term": term,
            "retmax": str(retmax),
            "retmode": "json",
            **self._params(),
        }
        if sort:
            params["sort"] = sort
        data = self.get(f"{EUTILS}/esearch.fcgi", params, cache_path=cache_path)
        payload = json.loads(data)
        result = payload.get("esearchresult", {})
        return {
            "count": int(result.get("count", 0)),
            "ids": [str(x) for x in result.get("idlist", [])],
            "querytranslation": result.get("querytranslation", ""),
        }

    def efetch(
        self,
        db: str,
        ids: list[str],
        *,
        retmode: str = "xml",
    ) -> bytes:
        params = {
            "db": db,
            "id": ",".join(ids),
            "retmode": retmode,
            **self._params(),
        }
        if db == "pmc":
            params["rettype"] = "xml"
        return self.post(f"{EUTILS}/efetch.fcgi", params)

    def idconv(self, ids: list[str], cache_dir: Path) -> dict[str, dict[str, Any]]:
        out: dict[str, dict[str, Any]] = {}
        for chunk in chunks(ids, 100):
            key = sha256_text("|".join(chunk))[:16]
            path = cache_dir / "idconv" / f"{key}.json"
            params = [
                ("ids", ",".join(chunk)),
                ("format", "json"),
                ("tool", TOOL),
            ]
            if self.email:
                params.append(("email", self.email))
            data = self.get(IDCONV, params, cache_path=path)
            payload = json.loads(data)
            for rec in payload.get("records", []):
                pmid = str(rec.get("pmid", "")).strip()
                if pmid:
                    out[pmid] = rec
        return out


# ---------------------------------------------------------------------------
# taxonomy
# ---------------------------------------------------------------------------


def parse_nodes_for_taxids(
    nodes_path: Path, taxids: set[str]
) -> tuple[dict[str, str], dict[str, str]]:
    """Return (parent, rank) for the requested taxids and all their ancestors."""
    parents: dict[str, str] = {}
    ranks: dict[str, str] = {}
    wanted = {str(t) for t in taxids}
    while wanted:
        found = False
        with open(nodes_path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = [p.strip() for p in line.split("\t|\t")]
                if len(parts) < 3:
                    continue
                taxid, parent, rank = parts[0], parts[1], parts[2]
                if taxid in wanted:
                    parents[taxid] = parent
                    ranks[taxid] = rank
                    wanted.discard(taxid)
                    if parent not in parents:
                        wanted.add(parent)
                    found = True
        if not found:
            missing = ", ".join(sorted(wanted))
            raise SystemExit(f"taxids not found in {nodes_path}: {missing}")
    return parents, ranks


def species_node_for(
    taxid: str, parents: dict[str, str], ranks: dict[str, str]
) -> str:
    cur = str(taxid)
    seen = set()
    while cur and cur not in seen:
        seen.add(cur)
        if ranks.get(cur) == "species":
            return cur
        cur = parents.get(cur, "1")
    return str(taxid)


def read_names_for(
    names_path: Path, taxids: set[str]
) -> dict[str, dict[str, list[str]]]:
    out: dict[str, dict[str, list[str]]] = {str(t): {} for t in taxids}
    wanted = set(out)
    with open(names_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip() for p in line.rstrip("\n").split("\t|\t")]
            if not parts:
                continue
            parts[-1] = parts[-1].rstrip("\t|")
            if len(parts) < 4:
                continue
            taxid, name, unique, name_class = parts[0], parts[1], parts[2], parts[3]
            if taxid not in wanted or name_class not in NAME_CLASSES:
                continue
            if name_class == "scientific name" and unique == "authority":
                continue
            out[taxid].setdefault(name_class, []).append(name)
    return out


def is_usable_name(name: str, tier: str) -> bool:
    """Reject authority strings and other names that are poor search keys."""
    name = norm_space(name)
    if not name or len(name) < 4 or len(name) > 80:
        return False
    low = name.lower()
    if any(
        token in low
        for token in (
            "nom.",
            "sensu",
            "auct.",
            "complex",
            "unclassified",
        )
    ):
        return False
    if re.search(r"\b(1[6-9]\d{2}|20\d{2})\b", name):
        return False
    if "&" in name or "[" in name or "]" in name or '"' in name:
        return False
    if "," in name:
        return False
    if re.search(r"\)\s+[A-Z][a-z]", name):
        return False
    if re.search(
        r"\b(?:ATCC|CBS|DSM|JCM|NRRL|IFO|IMI|UAMH|CGMCC|PYCC|MUCL)\s*\d+\b",
        name,
    ):
        return False
    if re.search(r"\b[A-Z]{2,6}\s?\d{3,}\b", name):
        return False
    low_without_fsp = low.replace("f. sp.", "f_sp").replace("f. spp.", "f_spp")
    if re.search(r"\b(sp|spp|cf|aff)\.", low_without_fsp):
        return False
    if len(name.split()) == 1 and tier != "genbank_common":
        return False
    return True


def load_aliases(path: Path) -> dict[str, list[dict[str, str]]]:
    rows = read_tsv(path, required=False)
    out: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        species = norm_space(row.get("species"))
        alias = norm_space(row.get("alias"))
        if not species or not alias:
            continue
        out.setdefault(species, []).append(row)
    return out


# ---------------------------------------------------------------------------
# species table
# ---------------------------------------------------------------------------


def collect_species_rows(study_dir: Path) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for csv_path in sorted(study_dir.glob("*/species.csv")):
        contrast = csv_path.parent.name
        with open(csv_path, newline="") as fh:
            for row in csv.DictReader(fh):
                species = norm_space(row.get("Species"))
                if not species:
                    continue
                rec = rows.setdefault(
                    species,
                    {
                        "species": species,
                        "taxon_ids": set(),
                        "shorts": set(),
                        "contrasts": set(),
                        "roles": set(),
                        "taxon_groups": set(),
                    },
                )
                rec["taxon_ids"].add(norm_space(row.get("Taxon_ID")))
                rec["shorts"].add(norm_space(row.get("Short")))
                rec["contrasts"].add(contrast)
                rec["roles"].add(norm_space(row.get("Group")))
                rec["taxon_groups"].add(norm_space(row.get("TaxonGroup")))
    return rows


def load_batch_codings(study_dir: Path) -> dict[str, dict[str, str]]:
    codings: dict[str, dict[str, str]] = {}
    for path in sorted((study_dir / "batches").glob("*.yaml")):
        spec = yaml.safe_load(path.read_text()) or {}
        for species, values in (spec.get("species") or {}).items():
            values = values or {}
            values = {
                k: v
                for k, v in values.items()
                if k in ("growth_form", "dominant_phase")
            }
            if species in codings and codings[species] != values:
                raise SystemExit(
                    f"conflicting batch coding for {species}: "
                    f"{codings[species]} vs {values}"
                )
            codings[species] = values
    return codings


def cmd_species(args: argparse.Namespace) -> None:
    study = Path(args.study_dir)
    out = Path(args.out) if args.out else study / "growth_form"
    taxonomy = Path(args.taxonomy_dir)
    out.mkdir(parents=True, exist_ok=True)

    rows = collect_species_rows(study)
    codings = load_batch_codings(study)
    missing = set(codings) - set(rows)
    extra = set(rows) - set(codings)
    if missing:
        raise SystemExit(f"batch species without species.csv rows: {sorted(missing)}")
    if extra:
        raise SystemExit(f"species.csv rows without batch coding: {sorted(extra)}")

    taxids: set[str] = set()
    for species, rec in rows.items():
        if len(rec["taxon_ids"]) != 1:
            raise SystemExit(
                f"{species}: expected one Taxon_ID, got {sorted(rec['taxon_ids'])}"
            )
        rec["taxon_id"] = next(iter(rec["taxon_ids"]))
        taxids.add(rec["taxon_id"])

    parents, ranks = parse_nodes_for_taxids(
        taxonomy / "nodes.dmp", taxids
    )

    table: list[dict[str, Any]] = []
    climbed = 0
    for species in sorted(rows):
        rec = rows[species]
        species_taxid = species_node_for(rec["taxon_id"], parents, ranks)
        if species_taxid != rec["taxon_id"]:
            climbed += 1
        coding = codings[species]
        table.append(
            {
                "species": species,
                "taxon_id": rec["taxon_id"],
                "species_taxon_id": species_taxid,
                "taxon_rank": ranks.get(rec["taxon_id"], ""),
                "species_rank": ranks.get(species_taxid, ""),
                "current_growth_form": coding.get("growth_form", ""),
                "current_dominant_phase": coding.get("dominant_phase", ""),
                "taxon_group": "; ".join(
                    sorted(x for x in rec["taxon_groups"] if x)
                ),
                "short_names": "; ".join(sorted(x for x in rec["shorts"] if x)),
                "contrasts": "; ".join(sorted(rec["contrasts"])),
                "n_contrasts": len(rec["contrasts"]),
                "roles": "; ".join(sorted(x for x in rec["roles"] if x)),
            }
        )

    write_tsv(
        out / "species.tsv",
        table,
        [
            "species",
            "taxon_id",
            "species_taxon_id",
            "taxon_rank",
            "species_rank",
            "current_growth_form",
            "current_dominant_phase",
            "taxon_group",
            "short_names",
            "contrasts",
            "n_contrasts",
            "roles",
        ],
    )
    print(
        f"species.tsv: {len(table)} species; "
        f"{climbed} strain-level Taxon_IDs climbed to species nodes; "
        f"{len(table) - climbed} already species rank"
    )


# ---------------------------------------------------------------------------
# names table
# ---------------------------------------------------------------------------


def cmd_names(args: argparse.Namespace) -> None:
    study = Path(args.study_dir)
    out = Path(args.out) if args.out else study / "growth_form"
    taxonomy = Path(args.taxonomy_dir)
    species_rows = read_tsv(out / "species.tsv")
    aliases = load_aliases(out / "aliases.tsv")
    species_nodes = {r["species_taxon_id"] for r in species_rows}
    names_by_node = read_names_for(taxonomy / "names.dmp", species_nodes)

    table: list[dict[str, Any]] = []
    only_canonical = 0
    for row in species_rows:
        species = row["species"]
        node = row["species_taxon_id"]
        by_class = names_by_node.get(node, {})
        canonical_list = by_class.get("scientific name", [])
        if not canonical_list:
            raise SystemExit(f"no scientific name for taxid {node} ({species})")
        canonical = canonical_list[0]
        candidates: list[tuple[str, str, str]] = []
        candidates.append((canonical, "canonical", f"names.dmp:{node}"))
        candidates.append((species, "canonical_input", "species.csv"))
        for cls in ("synonym", "equivalent name", "genbank common name"):
            for name in by_class.get(cls, []):
                tier = "genbank_common" if cls == "genbank common name" else "ncbi_synonym"
                candidates.append((name, tier, f"names.dmp:{cls}"))
        for alias_row in aliases.get(species, []):
            candidates.append(
                (
                    alias_row["alias"],
                    "literature_alias",
                    f"aliases.tsv:{alias_row.get('source', 'curated')}",
                )
            )

        seen: dict[str, dict[str, str]] = {}
        for name, tier, source in candidates:
            name = norm_space(name)
            if not is_usable_name(name, tier):
                continue
            key = name.lower()
            if key in seen:
                if NAME_TIER_ORDER[tier] < NAME_TIER_ORDER[seen[key]["name_tier"]]:
                    seen[key] = {
                        "name": name,
                        "name_tier": tier,
                        "name_source": source,
                    }
                continue
            seen[key] = {
                "name": name,
                "name_tier": tier,
                "name_source": source,
            }
        ordered = sorted(
            seen.values(),
            key=lambda r: (NAME_TIER_ORDER[r["name_tier"]], r["name"].lower()),
        )
        if len(ordered) == 1:
            only_canonical += 1
        for rec in ordered:
            table.append(
                {
                    "species": species,
                    "species_taxon_id": node,
                    **rec,
                }
            )

    write_tsv(
        out / "names.tsv",
        table,
        ["species", "species_taxon_id", "name", "name_tier", "name_source"],
    )
    print(
        f"names.tsv: {len(table)} names for {len(species_rows)} species; "
        f"{only_canonical} species have only the canonical name; "
        f"{sum(len(v) for v in aliases.values())} curated aliases"
    )


# ---------------------------------------------------------------------------
# queries and candidates
# ---------------------------------------------------------------------------


def build_queries(name: str) -> list[dict[str, str]]:
    growth_pubmed = " OR ".join(f'"{t}"[All Fields]' for t in GROWTH_TERMS)
    growth_pmc = " OR ".join(f'"{t}"[Text Word]' for t in GROWTH_TERMS)
    dimorph_pubmed = " OR ".join(f'"{t}"[Title/Abstract]' for t in DIMORPH_TERMS)
    dimorph_pmc = " OR ".join(f'"{t}"[Text Word]' for t in DIMORPH_TERMS)
    return [
        {
            "db": "pubmed",
            "mode": "dimorph",
            "term": f'"{name}"[Title/Abstract] AND ({dimorph_pubmed})',
        },
        {
            "db": "pmc",
            "mode": "dimorph",
            "term": f'"{name}"[Text Word] AND ({dimorph_pmc})',
        },
        {
            "db": "pubmed",
            "mode": "name",
            "term": f'"{name}"[Title/Abstract]',
        },
        {
            "db": "pubmed",
            "mode": "growth",
            "term": f'"{name}"[Title/Abstract] AND ({growth_pubmed})',
        },
        {
            "db": "pmc",
            "mode": "name",
            "term": f'"{name}"[Text Word]',
        },
        {
            "db": "pmc",
            "mode": "growth",
            "term": f'"{name}"[Text Word] AND ({growth_pmc})',
        },
    ]


def cmd_search(args: argparse.Namespace) -> None:
    study = Path(args.study_dir)
    out = Path(args.out) if args.out else study / "growth_form"
    names = read_tsv(out / "names.tsv")
    if args.species:
        keep = set(args.species)
        names = [r for r in names if r["species"] in keep]
    cache = out / "cache" / "search"
    cache.mkdir(parents=True, exist_ok=True)
    client = HttpClient(
        api_key=os.environ.get("NCBI_API_KEY") or None,
        email=args.email or os.environ.get("NCBI_EMAIL") or None,
    )

    query_rows: list[dict[str, Any]] = []
    by_species: dict[str, list[dict[str, Any]]] = {}
    seen_queries: set[str] = set()
    for row in names:
        species = row["species"]
        name = row["name"]
        for spec in build_queries(name):
            qid = sha256_text(
                f"{spec['db']}|{spec['mode']}|{spec['term']}|"
                f"{args.max_per_query}|{args.sort or ''}"
            )[:16]
            if qid in seen_queries:
                result = {
                    "count": int(
                        next(
                            q["count"]
                            for q in query_rows
                            if q["query_id"] == qid
                        )
                    ),
                    "ids": [],
                }
                cached = cache / f"{qid}.json"
                if cached.exists():
                    payload = json.loads(cached.read_text())
                    result["ids"] = [
                        str(x)
                        for x in payload.get("esearchresult", {}).get("idlist", [])
                    ]
            else:
                seen_queries.add(qid)
                result = client.esearch(
                    spec["db"],
                    spec["term"],
                    args.max_per_query,
                    sort=args.sort or None,
                    cache_path=cache / f"{qid}.json",
                )
                query_rows.append(
                    {
                        "query_id": qid,
                        "db": spec["db"],
                        "mode": spec["mode"],
                        "term": spec["term"],
                        "sort": args.sort or "",
                        "retmax": args.max_per_query,
                        "count": result["count"],
                        "n_returned": len(result["ids"]),
                        "querytranslation": result["querytranslation"],
                    }
                )
            by_species.setdefault(species, []).append(
                {
                    "query_id": qid,
                    "db": spec["db"],
                    "mode": spec["mode"],
                    "ids": result["ids"],
                }
            )

    candidate_rows: list[dict[str, Any]] = []
    for species in sorted(by_species):
        best: dict[tuple[str, str], dict[str, Any]] = {}
        for spec in by_species[species]:
            priority = QUERY_PRIORITY[spec["mode"]]
            for pos, uid in enumerate(spec["ids"]):
                key = (spec["db"], str(uid))
                rec = best.get(key)
                if rec is None:
                    best[key] = {
                        "species": species,
                        "db": spec["db"],
                        "uid": str(uid),
                        "best_mode": spec["mode"],
                        "priority": priority,
                        "rank": pos,
                        "query_ids": [spec["query_id"]],
                    }
                else:
                    if priority < rec["priority"]:
                        rec["priority"] = priority
                        rec["best_mode"] = spec["mode"]
                        rec["rank"] = pos
                    elif priority == rec["priority"]:
                        rec["rank"] = min(rec["rank"], pos)
                    if spec["query_id"] not in rec["query_ids"]:
                        rec["query_ids"].append(spec["query_id"])
        for rec in best.values():
            rec["query_ids"] = ";".join(rec["query_ids"])
            rec["n_queries"] = len(rec["query_ids"].split(";"))
            candidate_rows.append(rec)

    candidate_rows.sort(
        key=lambda r: (r["species"], r["priority"], r["rank"], r["db"], r["uid"])
    )
    write_tsv(
        out / "queries.tsv",
        query_rows,
        [
            "query_id",
            "db",
            "mode",
            "term",
            "sort",
            "retmax",
            "count",
            "n_returned",
            "querytranslation",
        ],
    )
    write_tsv(
        out / "candidates.tsv",
        candidate_rows,
        [
            "species",
            "db",
            "uid",
            "best_mode",
            "priority",
            "rank",
            "query_ids",
            "n_queries",
        ],
    )
    total_hits = sum(int(r["count"]) for r in query_rows)
    print(
        f"queries.tsv: {len(query_rows)} unique queries; "
        f"{total_hits} total hits; candidates.tsv: {len(candidate_rows)} "
        f"species-article candidates; API calls: {client.calls}"
    )


# ---------------------------------------------------------------------------
# article fetching and parsing
# ---------------------------------------------------------------------------


def parse_pubmed_articles(xml_bytes: bytes) -> dict[str, dict[str, Any]]:
    root = ET.fromstring(xml_bytes)
    articles = root.findall(".//PubmedArticle")
    if not articles and localname(root.tag) == "PubmedArticle":
        articles = [root]
    out: dict[str, dict[str, Any]] = {}
    for art in articles:
        pmid = elem_text(art.find(".//PMID"))
        if not pmid:
            continue
        title = elem_text(art.find(".//ArticleTitle"))
        abstract_parts = []
        for ab in art.findall(".//AbstractText"):
            label = ab.get("Label") or ""
            text = elem_text(ab)
            if text:
                abstract_parts.append(f"{label}: {text}" if label else text)
        journal = elem_text(art.find(".//Journal/Title")) or elem_text(
            art.find(".//Journal/ISOAbbreviation")
        )
        year = elem_text(art.find(".//JournalIssue/PubDate/Year")) or elem_text(
            art.find(".//JournalIssue/PubDate/MedlineDate")
        )
        doi = ""
        for aid in art.findall(".//ArticleId"):
            if aid.get("IdType") == "doi":
                doi = elem_text(aid)
                break
        out[pmid] = {
            "pmid": pmid,
            "pmcid": "",
            "doi": doi,
            "title": title,
            "journal": journal,
            "year": year,
            "abstract": norm_space(" ".join(abstract_parts)),
            "sections": (
                [("Abstract", norm_space(" ".join(abstract_parts)))]
                if abstract_parts
                else []
            ),
            "has_fulltext": False,
        }
    return out


def block_text(el: ET.Element) -> str:
    tag = localname(el.tag)
    if tag in ("table", "table-wrap"):
        rows = []
        for tr in el.findall(".//tr"):
            cells = [
                norm_space("".join(c.itertext()))
                for c in tr
                if localname(c.tag) in ("td", "th")
            ]
            if cells:
                rows.append(" | ".join(cells))
        return norm_space(" ".join(rows))
    return norm_space("".join(el.itertext()))


def collect_sections(
    node: ET.Element, section: str = "Body", out: list[tuple[str, str]] | None = None
) -> list[tuple[str, str]]:
    if out is None:
        out = []
    for child in node:
        tag = localname(child.tag)
        if tag == "sec":
            title = elem_text(child.find("title"))
            name = f"{section} / {title}" if title else section
            collect_sections(child, name, out)
        elif tag in ("p", "disp-quote"):
            text = block_text(child)
            if text:
                out.append((section, text))
        elif tag in ("table-wrap", "fig", "boxed-text"):
            collect_sections(child, section, out)
        elif tag == "caption":
            text = block_text(child)
            if text:
                out.append((section, text))
        elif tag == "title":
            continue
        else:
            collect_sections(child, section, out)
    return out


def parse_pmc_articles(xml_bytes: bytes) -> dict[str, dict[str, Any]]:
    root = ET.fromstring(xml_bytes)
    articles = root.findall(".//article")
    if not articles and localname(root.tag) == "article":
        articles = [root]
    out: dict[str, dict[str, Any]] = {}
    for art in articles:
        pmcid = ""
        pmid = ""
        doi = ""
        for aid in art.findall(".//article-id"):
            kind = aid.get("pub-id-type")
            value = elem_text(aid)
            if kind in {"pmc", "pmcid"}:
                pmcid = normalize_pmcid(value)
            elif kind == "pmid":
                pmid = value
            elif kind == "doi":
                doi = value
        if not pmcid:
            continue
        title = elem_text(art.find(".//article-title"))
        journal = elem_text(art.find(".//journal-title"))
        year = elem_text(art.find(".//pub-date/year"))
        sections: list[tuple[str, str]] = []
        abstract = art.find(".//article-meta/abstract")
        if abstract is not None:
            sections.extend(collect_sections(abstract, "Abstract"))
        body = art.find(".//body")
        if body is not None:
            sections.extend(collect_sections(body, "Body"))
        body_len = sum(len(text) for _, text in sections)
        out[pmcid] = {
            "pmid": pmid,
            "pmcid": pmcid,
            "doi": doi,
            "title": title,
            "journal": journal,
            "year": year,
            "abstract": norm_space(
                " ".join(text for sec, text in sections if sec == "Abstract")
            ),
            "sections": sections,
            "has_fulltext": body is not None and body_len > 200,
        }
    return out


def cache_paths(cache: Path, db: str, uid: str) -> tuple[Path, Path]:
    return cache / db / f"{uid}.xml.zst", cache / db / f"{uid}.json"


def load_cached_article(
    cache: Path, db: str, uid: str
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    xml_path, meta_path = cache_paths(cache, db, uid)
    if not xml_path.exists() or not meta_path.exists():
        return None, {}
    meta = json.loads(meta_path.read_text())
    raw = read_zst(xml_path)
    if db == "pmc":
        parsed = parse_pmc_articles(raw).get(uid)
    else:
        parsed = parse_pubmed_articles(raw).get(uid)
    return parsed, meta


def split_pubmed_elements(xml_bytes: bytes) -> dict[str, bytes]:
    root = ET.fromstring(xml_bytes)
    articles = root.findall(".//PubmedArticle")
    if not articles and localname(root.tag) == "PubmedArticle":
        articles = [root]
    out: dict[str, bytes] = {}
    for art in articles:
        pmid = elem_text(art.find(".//PMID"))
        if pmid:
            out[pmid] = ET.tostring(art, encoding="utf-8")
    return out


def split_pmc_elements(xml_bytes: bytes) -> dict[str, bytes]:
    root = ET.fromstring(xml_bytes)
    articles = root.findall(".//article")
    if not articles and localname(root.tag) == "article":
        articles = [root]
    out: dict[str, bytes] = {}
    for art in articles:
        pmcid = ""
        for aid in art.findall(".//article-id"):
            if aid.get("pub-id-type") in {"pmc", "pmcid"}:
                pmcid = normalize_pmcid(elem_text(aid))
                break
        if pmcid:
            out[pmcid] = ET.tostring(art, encoding="utf-8")
    return out


def fetch_and_cache(
    client: HttpClient, cache: Path, db: str, uids: list[str]
) -> dict[str, str]:
    """Fetch uncached articles, write per-article cache, return uid->status."""
    status: dict[str, str] = {}
    if not uids:
        return status
    batch_size = 20 if db == "pmc" else 100
    for group in chunks(uids, batch_size):
        missing = [
            u
            for u in group
            if not cache_paths(cache, db, u)[0].exists()
            or not cache_paths(cache, db, u)[1].exists()
        ]
        if not missing:
            continue
        try:
            raw = client.efetch(db, missing)
            elements = (
                split_pmc_elements(raw)
                if db == "pmc"
                else split_pubmed_elements(raw)
            )
        except Exception as exc:  # noqa: BLE001
            if len(missing) == 1:
                status[missing[0]] = f"error: {exc}"
                continue
            for uid in missing:
                status.update(fetch_and_cache(client, cache, db, [uid]))
            continue
        for uid in missing:
            art_raw = elements.get(uid)
            if art_raw is None:
                status[uid] = "not_in_response"
                continue
            parsed = (
                parse_pmc_articles(art_raw)
                if db == "pmc"
                else parse_pubmed_articles(art_raw)
            ).get(uid)
            if parsed is None:
                status[uid] = "parse_error"
                continue
            xml_path, meta_path = cache_paths(cache, db, uid)
            write_zst(xml_path, art_raw)
            meta = {
                "uid": uid,
                "db": db,
                "retrieved_utc": utc_now(),
                "xml_sha256": sha256_bytes(art_raw),
                "bytes": len(art_raw),
                "pmid": parsed.get("pmid", ""),
                "pmcid": parsed.get("pmcid", ""),
                "title": parsed.get("title", ""),
                "journal": parsed.get("journal", ""),
                "year": parsed.get("year", ""),
                "has_fulltext": bool(parsed.get("has_fulltext")),
            }
            meta_path.parent.mkdir(parents=True, exist_ok=True)
            meta_path.write_text(json.dumps(meta, indent=2) + "\n")
            status[uid] = "ok"
    return status


def merge_article_records(
    records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """One record per article. A paper found in PubMed (uid = PMID) and in PMC (uid = PMCID)
    is the same article: merge them under the PMCID when it is known, so that agreement between
    the two databases counts and one article does not take two selection slots.

    The merged record keeps the best priority, the best rank within that priority, and the sum of
    n_queries. It is a PMC record (full text) when a PMCID is known, else a PubMed record.
    """
    merged: dict[str, dict[str, Any]] = {}
    for rec in records:
        # The ID converter returns "PMC123"; PMC records and the cache use the bare number "123".
        pmcid = normalize_pmcid(str(rec.get("pmcid", "") or ""))
        key = f"pmc:{pmcid}" if pmcid else f"{rec['db']}:{rec['uid']}"
        cur = merged.get(key)
        if cur is None:
            merged[key] = {
                **rec,
                "db": "pmc" if pmcid else rec["db"],
                "uid": pmcid if pmcid else rec["uid"],
                "pmcid": pmcid,
                "n_queries": int(rec.get("n_queries", 1) or 1),
            }
            continue
        cur["n_queries"] += int(rec.get("n_queries", 1) or 1)
        if not cur.get("pmid") and rec.get("pmid"):
            cur["pmid"] = rec["pmid"]
        if rec["priority"] < cur["priority"]:
            cur["priority"], cur["rank"] = rec["priority"], rec["rank"]
        elif rec["priority"] == cur["priority"]:
            cur["rank"] = min(cur["rank"], rec["rank"])
    return list(merged.values())


def select_for_species(
    records: list[dict[str, Any]],
    max_articles: int,
    max_fulltext: int,
    selection: str = "score",
) -> list[dict[str, Any]]:
    """Pick the articles to fetch for one species.

    selection "score" (default): merge duplicate articles, then order by
    (priority, -n_queries, rank): a paper that several queries returned, and that a dimorphism or
    growth query returned, comes first. This replaces the earlier order below, in which rank ties
    among the many name queries of one species filled the small cap with arbitrary papers.
    selection "rank": the earlier order (priority, rank, db, uid), no merging.
    Forced records (priority -1) always sort first.
    """
    if selection == "score":
        ordered = sorted(
            merge_article_records(records),
            key=lambda r: (r["priority"], -r["n_queries"], r["rank"], r["db"], r["uid"]),
        )
    elif selection == "rank":
        ordered = []
        seen: set[str] = set()
        for rec in sorted(
            records, key=lambda r: (r["priority"], r["rank"], r["db"], r["uid"])
        ):
            key = f"{rec['db']}:{rec['uid']}"
            if key in seen:
                continue
            seen.add(key)
            ordered.append(rec)
    else:
        raise ValueError(f"selection must be 'score' or 'rank', got {selection!r}")
    return apply_caps(ordered, max_articles, max_fulltext)


def apply_caps(
    ordered: list[dict[str, Any]], max_articles: int, max_fulltext: int
) -> list[dict[str, Any]]:
    """Take records in order until max_articles; at most max_fulltext of them PMC (full text)."""
    out: list[dict[str, Any]] = []
    n_ft = 0
    for rec in ordered:
        if len(out) >= max_articles:
            break
        if rec["db"] == "pmc":
            if n_ft >= max_fulltext:
                continue
            n_ft += 1
        out.append(rec)
    return out


def screen_score(title: str, abstract: str, names: list[str]) -> float:
    """How likely an article answers the growth-form question, from its title and abstract.

    +3 for each distinct dimorphism term, +1 for each distinct growth term (at most 5 growth terms),
    +5 if a usable name of the species is in the title, +2 if it is in the abstract. A heuristic for
    ordering only; it never decides that a coding is true.
    """
    low_t, low_a = (title or "").lower(), (abstract or "").lower()
    both = f"{low_t} {low_a}"
    dim = sum(1 for term in DIMORPH_TERMS if term.lower() in both)
    growth = min(5, sum(1 for term in GROWTH_TERMS if term.lower() in both))
    score = 3.0 * dim + 1.0 * growth
    lower_names = [n.lower() for n in names if len(n) >= 6]
    if any(n in low_t for n in lower_names):
        score += 5.0
    elif any(n in low_a for n in lower_names):
        score += 2.0
    return score


def select_two_step(
    records: list[dict[str, Any]],
    screen: dict[str, tuple[str, str]],
    names: list[str],
    max_articles: int,
    max_fulltext: int,
    screen_articles: int,
    base_share: float = 0.5,
) -> list[dict[str, Any]]:
    """Two-step selection, as a hybrid of two orders.

    Base order: merged records ordered by (priority, -n_queries, rank), as in `score` selection.
    Screen order: the first `screen_articles` of the base order, re-ordered by the title/abstract
    screen score (`screen` maps PMID to (title, abstract); no abstract scores 0), ties by base order.
    The first `base_share` of the slots go to the base order and the rest to the screen order. Forced
    records (priority < 0) always come first. Why a hybrid: in the 19-row gold test the base order missed
    a dimorphism paper that the screen found (Liang 2017), and the screen missed a paper whose abstract
    carries no dimorphism term (Tahiri 2026) that the base order found.
    """
    base = sorted(
        merge_article_records(records),
        key=lambda r: (r["priority"], -r["n_queries"], r["rank"], r["db"], r["uid"]),
    )
    pool = [r for i, r in enumerate(base) if i < screen_articles or r["priority"] < 0]
    order_idx = {id(r): i for i, r in enumerate(base)}
    for r in pool:
        title, abstract = screen.get(str(r.get("pmid", "") or ""), ("", ""))
        r["screen_score"] = screen_score(title, abstract, names)
    by_screen = sorted(
        pool, key=lambda r: (r["priority"] >= 0, -r["screen_score"], order_idx[id(r)])
    )
    n_base = int(max_articles * base_share)
    ordered: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    forced = [r for r in base if r["priority"] < 0]
    for r in forced + [r for r in base if r["priority"] >= 0][:n_base] + by_screen:
        key = (r["db"], r["uid"])
        if key in seen:
            continue
        seen.add(key)
        ordered.append(r)
    return apply_caps(ordered, max_articles, max_fulltext)


def cmd_fetch(args: argparse.Namespace) -> None:
    study = Path(args.study_dir)
    out = Path(args.out) if args.out else study / "growth_form"
    cache = out / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    candidates = read_tsv(out / "candidates.tsv")
    if args.species:
        keep = set(args.species)
        candidates = [r for r in candidates if r["species"] in keep]
    client = HttpClient(
        api_key=os.environ.get("NCBI_API_KEY") or None,
        email=args.email or os.environ.get("NCBI_EMAIL") or None,
    )

    all_pmids = sorted({r["uid"] for r in candidates if r["db"] == "pubmed"})
    idmap = client.idconv(all_pmids, cache) if all_pmids else {}

    by_species: dict[str, list[dict[str, Any]]] = {}
    for row in candidates:
        species = row["species"]
        db = row["db"]
        uid = row["uid"]
        rec = {
            "db": db,
            "uid": uid,
            "priority": int(row["priority"]),
            "rank": int(row.get("rank", 0) or 0),
            "n_queries": int(row.get("n_queries", 1) or 1),
            "pmid": uid if db == "pubmed" else "",
            "pmcid": uid if db == "pmc" else "",
        }
        if db == "pubmed":
            conv = idmap.get(uid, {})
            rec["pmcid"] = str(conv.get("pmcid", "") or "")
        by_species.setdefault(species, []).append(rec)

    include_all: list[str] = []
    include_by_species: dict[str, list[str]] = {}
    for item in args.include_pmcid or []:
        if "=" in item:
            sp, pmcid = item.split("=", 1)
            include_by_species.setdefault(sp.strip(), []).append(
                normalize_pmcid(pmcid)
            )
        else:
            include_all.append(normalize_pmcid(item))
    for species in by_species:
        for pmcid in include_all + include_by_species.get(species, []):
            if pmcid:
                by_species[species].append(
                    {
                        "db": "pmc",
                        "uid": pmcid,
                        "priority": -1,
                        "rank": -1,
                        "pmid": "",
                        "pmcid": pmcid,
                    }
                )

    selected: list[dict[str, Any]] = []
    if args.selection == "two-step":
        usable: dict[str, list[str]] = {}
        for r in read_tsv(out / "names.tsv"):
            if r.get("name_tier") != "genus":
                usable.setdefault(r["species"], []).append(r["name"])
        pools: dict[str, list[dict[str, Any]]] = {}
        screen_pmids: set[str] = set()
        for species in sorted(by_species):
            base = sorted(
                merge_article_records(by_species[species]),
                key=lambda r: (r["priority"], -r["n_queries"], r["rank"], r["db"], r["uid"]),
            )
            pools[species] = base[: args.screen_articles]
            screen_pmids.update(str(r["pmid"]) for r in pools[species] if r.get("pmid"))
        fetch_and_cache(client, cache, "pubmed", sorted(screen_pmids))
        screen: dict[str, tuple[str, str]] = {}
        for pmid in sorted(screen_pmids):
            parsed, _meta = load_cached_article(cache, "pubmed", pmid)
            if parsed:
                screen[pmid] = (parsed.get("title", ""), parsed.get("abstract", ""))
        for species in sorted(by_species):
            for rec in select_two_step(
                by_species[species],
                screen,
                usable.get(species, [species]),
                args.max_articles,
                args.max_fulltext,
                args.screen_articles,
                args.base_share,
            ):
                rec["species"] = species
                selected.append(rec)
    else:
        for species in sorted(by_species):
            for rec in select_for_species(
                by_species[species],
                args.max_articles,
                args.max_fulltext,
                args.selection,
            ):
                rec["species"] = species
                selected.append(rec)

    selected_pmc = {r["uid"] for r in selected if r["db"] == "pmc"}
    selected = [
        r
        for r in selected
        if not (r["db"] == "pubmed" and r["pmcid"] in selected_pmc)
    ]

    pmc_uids = sorted(
        {r["uid"] for r in selected if r["db"] == "pmc"}
    )
    pubmed_uids = sorted({r["uid"] for r in selected if r["db"] == "pubmed"})
    pmc_status = fetch_and_cache(client, cache, "pmc", pmc_uids)
    pubmed_status = fetch_and_cache(client, cache, "pubmed", pubmed_uids)

    article_rows: list[dict[str, Any]] = []
    failures: list[str] = []
    for rec in selected:
        parsed, meta = load_cached_article(cache, rec["db"], rec["uid"])
        if parsed is None or not meta:
            failures.append(
                f"{rec['species']}\t{rec['db']}\t{rec['uid']}\t"
                f"{pmc_status.get(rec['uid']) or pubmed_status.get(rec['uid']) or 'missing'}"
            )
            continue
        article_rows.append(
            {
                "species": rec["species"],
                "source_db": rec["db"],
                "uid": rec["uid"],
                "pmid": parsed.get("pmid", "") or rec["pmid"],
                "pmcid": parsed.get("pmcid", "") or rec["pmcid"],
                "doi": parsed.get("doi", ""),
                "title": parsed.get("title", ""),
                "journal": parsed.get("journal", ""),
                "year": parsed.get("year", ""),
                "has_fulltext": "yes" if parsed.get("has_fulltext") else "no",
                "has_abstract": "yes" if parsed.get("abstract") else "no",
                "xml_sha256": meta.get("xml_sha256", ""),
                "bytes": meta.get("bytes", ""),
                "retrieved_utc": meta.get("retrieved_utc", ""),
                "priority": rec["priority"],
                "read_level": "full_text"
                if rec["db"] == "pmc" and parsed.get("has_fulltext")
                else "abstract",
            }
        )

    article_rows.sort(key=lambda r: (r["species"], r["source_db"], r["uid"]))
    write_tsv(
        out / "articles.tsv",
        article_rows,
        [
            "species",
            "source_db",
            "uid",
            "pmid",
            "pmcid",
            "doi",
            "title",
            "journal",
            "year",
            "has_fulltext",
            "has_abstract",
            "xml_sha256",
            "bytes",
            "retrieved_utc",
            "priority",
            "read_level",
        ],
    )
    failure_path = out / "fetch_failures.tsv"
    if failures:
        write_tsv(
            failure_path,
            [{"failure": f} for f in failures],
            ["failure"],
        )
    elif failure_path.exists():
        failure_path.unlink()
    print(
        f"articles.tsv: {len(article_rows)} rows for {len(by_species)} species; "
        f"PMC fetched: {len(pmc_uids)}; PubMed abstracts fetched: "
        f"{len(pubmed_uids)}; failures: {len(failures)}; API calls: "
        f"{client.calls}"
    )


# ---------------------------------------------------------------------------
# sentence extraction
# ---------------------------------------------------------------------------


def split_sentences(text: str) -> list[str]:
    text = norm_space(text)
    if not text:
        return []
    protected = ABBREV_RE.sub(
        lambda m: m.group(0).replace(".", SENT_PROTECT), text
    )
    parts = SENT_SPLIT_RE.split(protected)
    return [
        part.replace(SENT_PROTECT, ".").strip()
        for part in parts
        if part.strip()
    ]


def name_pattern(name: str) -> re.Pattern[str]:
    return re.compile(
        rf"(?<![A-Za-z]){re.escape(name)}(?![A-Za-z])", re.I
    )


def abbreviation_for(name: str) -> str | None:
    parts = name.split()
    if len(parts) != 2:
        return None
    if not parts[0][:1].isupper() or not parts[1][:1].islower():
        return None
    if parts[1].endswith(".") or parts[0].lower() in {"sp", "spp"}:
        return None
    return f"{parts[0][0]}. {parts[1]}"


def article_sections(source_db: str, parsed: dict[str, Any]) -> list[tuple[str, str]]:
    sections = parsed.get("sections") or []
    if source_db == "pubmed" and not sections and parsed.get("abstract"):
        return [("Abstract", parsed["abstract"])]
    return sections


def parenthetical_genus_patterns(
    names: list[tuple[str, str]],
) -> list[tuple[str, str, re.Pattern[str]]]:
    """Patterns for names written with both genera, for example "Wangiella (Exophiala) dermatitidis".

    `names` holds (name, tier) pairs of one species. Names with the same epithet under two or more
    genera (a current name and its synonyms) give one pattern: "G1 (G2) epithet" with G1 and G2 any of
    those genera. The label and tier are those of the best-tier name of the group.
    """
    groups: dict[str, list[tuple[str, str]]] = {}
    for name, tier in names:
        parts = name.split()
        if len(parts) == 2 and parts[0][:1].isupper() and parts[1][:1].islower():
            groups.setdefault(parts[1], []).append((parts[0], tier, name))
    out: list[tuple[str, str, re.Pattern[str]]] = []
    for epithet, members in sorted(groups.items()):
        genera = sorted({g for g, _t, _n in members})
        if len(genera) < 2:
            continue
        best = min(members, key=lambda m: (NAME_TIER_ORDER.get(m[1], 9), m[2]))
        alt = "|".join(re.escape(g) for g in genera)
        pat = re.compile(
            rf"(?<![A-Za-z])(?:{alt})\s+\((?:{alt})\)\s+{re.escape(epithet)}(?![A-Za-z])", re.I
        )
        out.append((best[2], best[1], pat))
    return out


def is_about_species(
    title: str,
    abstract: str,
    patterns: list[tuple[str, str, re.Pattern[str]]],
    mode: str,
) -> bool:
    """True if the article is about the species: a full name in the title (mode "title"), or in the
    title or abstract (mode "abstract"). Mode "off" is always False. Genus-only matches do not count."""
    if mode == "off":
        return False
    text = title or ""
    if mode == "abstract":
        text = f"{title or ''} {abstract or ''}"
    return any(pat.search(text) for _name, _tier, pat in patterns)


def pick_article_level(
    candidates: list[tuple[int, int, Any]], cap: int
) -> list[tuple[int, int, Any]]:
    """Keep at most `cap` article-level candidates, (priority, position, payload) each. Priority 0
    (the sentence has a dimorphism term) comes before priority 1; document order breaks ties."""
    if cap <= 0:
        return []
    return sorted(candidates, key=lambda c: (c[0], c[1]))[:cap]


def cmd_sentences(args: argparse.Namespace) -> None:
    study = Path(args.study_dir)
    out = Path(args.out) if args.out else study / "growth_form"
    cache = out / "cache"
    names_rows = read_tsv(out / "names.tsv")
    articles = read_tsv(out / "articles.tsv")
    if args.species:
        keep = set(args.species)
        names_rows = [r for r in names_rows if r["species"] in keep]
        articles = [r for r in articles if r["species"] in keep]

    names_by_species: dict[str, list[dict[str, str]]] = {}
    for row in names_rows:
        names_by_species.setdefault(row["species"], []).append(row)

    articles_by_species: dict[str, list[dict[str, str]]] = {}
    for row in articles:
        articles_by_species.setdefault(row["species"], []).append(row)

    summary: list[dict[str, Any]] = []
    sent_dir = out / "sentences"
    sent_dir.mkdir(parents=True, exist_ok=True)

    for species in sorted(names_by_species):
        name_recs = names_by_species[species]
        patterns = [
            (r["name"], r["name_tier"], name_pattern(r["name"]))
            for r in name_recs
            if len(r["name"]) >= 4
        ]
        abbrev_patterns = []
        for name, tier, _pat in patterns:
            ab = abbreviation_for(name)
            if ab:
                abbrev_patterns.append((name, tier, name_pattern(ab)))
        # "Wangiella (Exophiala) dermatitidis": one name written with both genera.
        patterns += parenthetical_genus_patterns(
            [(r["name"], r["name_tier"]) for r in name_recs if len(r["name"]) >= 4]
        )
        genus_to_names: dict[str, list[str]] = {}
        for name, _tier, _pat in patterns:
            genus = name.split()[0]
            genus_to_names.setdefault(genus, []).append(name)

        records: list[dict[str, Any]] = []
        seen: dict[tuple[str, str, str], dict[str, Any]] = {}
        n_sent = 0
        n_fulltext_articles = 0
        n_abstract_articles = 0
        for article in articles_by_species.get(species, []):
            source_db = article["source_db"]
            uid = article["uid"]
            parsed, meta = load_cached_article(cache, source_db, uid)
            if parsed is None:
                continue
            if parsed.get("has_fulltext"):
                n_fulltext_articles += 1
            else:
                n_abstract_articles += 1
            sections = article_sections(source_db, parsed)
            article_text = norm_space(
                " ".join(text for _sec, text in sections)
            )
            full_names_found = {
                name
                for name, _tier, pat in patterns
                if pat.search(article_text)
            }
            # Article-level: for a paper about this species (name in the title or abstract), growth
            # sentences that do not name it ("the fungus ...") are kept too, labelled match_level
            # "article", at most --article-cap per article, dimorphism sentences first.
            about = is_about_species(
                parsed.get("title", ""), parsed.get("abstract", ""), patterns, args.article_level
            )
            article_cands: list[tuple[int, int, dict[str, Any]]] = []
            position = 0
            for section, section_text in sections:
                sentences = split_sentences(section_text)
                for idx, sentence in enumerate(sentences):
                    n_sent += 1
                    if not GROWTH_RE.search(sentence):
                        continue
                    hits: list[tuple[str, str, str]] = []
                    for name, tier, pat in patterns:
                        if pat.search(sentence):
                            hits.append((name, tier, "full"))
                    for name, tier, pat in abbrev_patterns:
                        if name not in full_names_found:
                            continue
                        if any(h[0] == name for h in hits):
                            continue
                        if pat.search(sentence):
                            hits.append((name, tier, "abbrev"))
                    position += 1
                    if not hits:
                        genus = species.split()[0]
                        if not re.search(
                            rf"(?<![A-Za-z]){re.escape(genus)}(?![A-Za-z])",
                            sentence,
                            re.I,
                        ):
                            if about:
                                best = min(
                                    (p for p in patterns if p[0] in full_names_found),
                                    key=lambda p: (NAME_TIER_ORDER.get(p[1], 9), -len(p[0])),
                                    default=(species, "canonical", None),
                                )
                                article_cands.append(
                                    (
                                        0 if DIMORPH_RE.search(sentence) else 1,
                                        position,
                                        {
                                            "section": section,
                                            "sentence": sentence,
                                            "idx": idx,
                                            "sentences": sentences,
                                            "name": best[0],
                                            "tier": best[1],
                                        },
                                    )
                                )
                            continue
                        other = [
                            n
                            for n in genus_to_names.get(genus, [])
                            if n != species
                            and any(
                                p.search(sentence)
                                for nm, _t, p in patterns
                                if nm == n
                            )
                        ]
                        if other:
                            continue
                        hits.append((genus, "genus", "genus"))

                    # Keep the best hit for this species (full > abbrev > genus,
                    # then canonical/nearer tiers).
                    level_rank = {"full": 0, "abbrev": 1, "genus": 2}
                    hits.sort(
                        key=lambda h: (
                            level_rank[h[2]],
                            NAME_TIER_ORDER.get(h[1], 9),
                            -len(h[0]),
                        )
                    )
                    matched_name, tier, level = hits[0]
                    text = sentence
                    key = (source_db, uid, norm_quote(text))
                    if key in seen:
                        old = seen[key]
                        if (
                            level_rank[level] < level_rank[old["match_level"]]
                            or (
                                level_rank[level] == level_rank[old["match_level"]]
                                and NAME_TIER_ORDER.get(tier, 9)
                                < NAME_TIER_ORDER.get(old["match_tier"], 9)
                            )
                        ):
                            seen[key] = {
                                **old,
                                "matched_name": matched_name,
                                "match_tier": tier,
                                "match_level": level,
                                "name_used_in_source": "yes"
                                if level in ("full", "abbrev")
                                else "no",
                            }
                        continue
                    rec = {
                        "sentence_id": sha256_text(
                            f"{source_db}|{uid}|{norm_quote(text)}"
                        )[:16],
                        "species": species,
                        "source_type": "pmc" if source_db == "pmc" else "pubmed",
                        "uid": uid,
                        "pmid": parsed.get("pmid", ""),
                        "pmcid": parsed.get("pmcid", ""),
                        "doi": parsed.get("doi", ""),
                        "title": parsed.get("title", ""),
                        "journal": parsed.get("journal", ""),
                        "year": parsed.get("year", ""),
                        "section": section,
                        "read_level": "full_text"
                        if parsed.get("has_fulltext")
                        else "abstract",
                        "matched_name": matched_name,
                        "match_tier": tier,
                        "match_level": level,
                        "name_used_in_source": "yes"
                        if level in ("full", "abbrev")
                        else "no",
                        "text": text,
                        "context_before": sentences[idx - 1]
                        if idx > 0
                        else "",
                        "context_after": sentences[idx + 1]
                        if idx + 1 < len(sentences)
                        else "",
                    }
                    seen[key] = rec
                    records.append(rec)

            for _prio, _pos, c in pick_article_level(article_cands, args.article_cap):
                text = c["sentence"]
                key = (source_db, uid, norm_quote(text))
                if key in seen:
                    continue
                sents = c["sentences"]
                rec = {
                    "sentence_id": sha256_text(
                        f"{source_db}|{uid}|{norm_quote(text)}"
                    )[:16],
                    "species": species,
                    "source_type": "pmc" if source_db == "pmc" else "pubmed",
                    "uid": uid,
                    "pmid": parsed.get("pmid", ""),
                    "pmcid": parsed.get("pmcid", ""),
                    "doi": parsed.get("doi", ""),
                    "title": parsed.get("title", ""),
                    "journal": parsed.get("journal", ""),
                    "year": parsed.get("year", ""),
                    "section": c["section"],
                    "read_level": "full_text" if parsed.get("has_fulltext") else "abstract",
                    "matched_name": c["name"],
                    "match_tier": c["tier"],
                    "match_level": "article",
                    "name_used_in_source": "no",
                    "text": text,
                    "context_before": sents[c["idx"] - 1] if c["idx"] > 0 else "",
                    "context_after": sents[c["idx"] + 1] if c["idx"] + 1 < len(sents) else "",
                }
                seen[key] = rec
                records.append(rec)

        records.sort(
            key=lambda r: (r["source_type"], r["uid"], r["section"], r["sentence_id"])
        )
        path = sent_dir / f"{slugify(species)}.jsonl"
        with open(path, "w") as fh:
            for rec in records:
                fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
        approx_tokens = sum(len(r["text"]) for r in records) // 4
        summary.append(
            {
                "species": species,
                "n_names": len(patterns),
                "n_articles_fulltext": n_fulltext_articles,
                "n_articles_abstract": n_abstract_articles,
                "n_sentences_scanned": n_sent,
                "n_records": len(records),
                "n_species_level": sum(
                    1 for r in records if r["match_level"] in ("full", "abbrev")
                ),
                "n_genus_level": sum(
                    1 for r in records if r["match_level"] == "genus"
                ),
                "n_article_level": sum(
                    1 for r in records if r["match_level"] == "article"
                ),
                "approx_text_tokens": approx_tokens,
            }
        )

    write_tsv(
        out / "sentences_summary.tsv",
        summary,
        [
            "species",
            "n_names",
            "n_articles_fulltext",
            "n_articles_abstract",
            "n_sentences_scanned",
            "n_records",
            "n_species_level",
            "n_genus_level",
            "n_article_level",
            "approx_text_tokens",
        ],
    )
    print(
        f"sentences_summary.tsv: {len(summary)} species; "
        f"{sum(r['n_records'] for r in summary)} candidate sentences; "
        f"{sum(r['approx_text_tokens'] for r in summary)} approx tokens"
    )


# ---------------------------------------------------------------------------
# quote verification
# ---------------------------------------------------------------------------


def source_text_for(species: str, article: dict[str, str], cache: Path) -> str:
    parsed, _meta = load_cached_article(
        cache, article["source_db"], article["uid"]
    )
    if parsed is None:
        return ""
    return norm_space(
        " ".join(text for _sec, text in article_sections(article["source_db"], parsed))
    )


def quote_in_text(quote: str, text: str) -> bool:
    q = norm_quote(quote).strip(" .,;:!?")
    t = norm_quote(text)
    return bool(q) and q in t


def cmd_verify(args: argparse.Namespace) -> None:
    study = Path(args.study_dir)
    out = Path(args.out) if args.out else study / "growth_form"
    cache = out / "cache"
    evidence = read_tsv(Path(args.evidence))
    articles = read_tsv(out / "articles.tsv", required=False)
    article_index: dict[tuple[str, str, str], dict[str, str]] = {
        (r["species"], r["source_db"], r["uid"]): r for r in articles
    }
    sentence_index: dict[str, dict[str, dict[str, Any]]] = {}
    sentence_dir = out / "sentences"
    if sentence_dir.exists():
        for path in sentence_dir.glob("*.jsonl"):
            by_id: dict[str, dict[str, Any]] = {}
            with open(path) as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    rec = json.loads(line)
                    by_id[rec["sentence_id"]] = rec
            sentence_index[path.stem] = by_id

    verified: list[dict[str, Any]] = []
    rejected: list[dict[str, Any]] = []
    for row in evidence:
        reasons: list[str] = []
        species = row.get("species", "").strip()
        sentence_id = row.get("sentence_id", "").strip()
        quote = norm_space(row.get("quote", ""))
        if not species:
            reasons.append("missing_species")
        if not sentence_id:
            reasons.append("missing_sentence_id")
        if not quote:
            reasons.append("missing_quote")
        if word_count(quote) > 40:
            reasons.append(f"quote_too_long({word_count(quote)})")

        rec = None
        if species and sentence_id:
            rec = sentence_index.get(slugify(species), {}).get(sentence_id)
            if rec is None:
                reasons.append("sentence_id_not_found")
        if rec is not None:
            key = (rec["species"], rec["source_type"], rec["uid"])
            article = article_index.get(key)
            if article is None:
                reasons.append("article_not_fetched")
            else:
                text = source_text_for(rec["species"], article, cache)
                if not text:
                    reasons.append("source_text_empty")
                elif not quote_in_text(quote, text):
                    reasons.append("quote_not_in_source")
                if quote and rec.get("matched_name"):
                    if not quote_in_text(rec["matched_name"], quote):
                        reasons.append("name_not_in_quote")
                if not rec.get("name_used_in_source") == "yes":
                    reasons.append("name_not_verbatim_in_sentence")

        out_row = {
            **row,
            "quote_word_count": word_count(quote),
            "verified": "no" if reasons else "yes",
            "rejection_reason": ";".join(reasons),
            "sentence_text": rec.get("text", "") if rec else "",
            "source_type": rec.get("source_type", "") if rec else "",
            "pmid": rec.get("pmid", "") if rec else "",
            "pmcid": rec.get("pmcid", "") if rec else "",
            "doi": rec.get("doi", "") if rec else "",
            "read_level": rec.get("read_level", "") if rec else "",
            "source_title": rec.get("title", "") if rec else "",
            "name_used_in_source": rec.get("name_used_in_source", "") if rec else "",
            "verified_utc": utc_now(),
            "verifier": SCRIPT_VERSION,
        }
        (verified if not reasons else rejected).append(out_row)

    fields = list(evidence[0].keys()) if evidence else ["species"]
    extra_fields = [
        "quote_word_count",
        "verified",
        "rejection_reason",
        "sentence_text",
        "source_type",
        "pmid",
        "pmcid",
        "doi",
        "read_level",
        "source_title",
        "name_used_in_source",
        "verified_utc",
        "verifier",
    ]
    write_tsv(out / "evidence.verified.tsv", verified, fields + extra_fields)
    write_tsv(out / "evidence.rejected.tsv", rejected, fields + extra_fields)
    print(
        f"verify: {len(verified)} verified, {len(rejected)} rejected "
        f"(of {len(evidence)})"
    )


# ---------------------------------------------------------------------------
# gold check
# ---------------------------------------------------------------------------


def cmd_gold(args: argparse.Namespace) -> None:
    study = Path(args.study_dir)
    out = Path(args.out) if args.out else study / "growth_form"
    cache = out / "cache"
    gold = read_tsv(out / "gold.tsv")
    articles = read_tsv(out / "articles.tsv", required=False)
    article_index: dict[tuple[str, str, str], dict[str, str]] = {
        (r["species"], r["source_db"], r["uid"]): r for r in articles
    }
    sentence_index: dict[str, list[dict[str, Any]]] = {}
    sent_dir = out / "sentences"
    if sent_dir.exists():
        for path in sent_dir.glob("*.jsonl"):
            with open(path) as fh:
                sentence_index[path.stem] = [
                    json.loads(line) for line in fh if line.strip()
                ]

    rows: list[dict[str, Any]] = []
    for gold_row in gold:
        species = gold_row["species"].strip()
        pmcid = gold_row.get("pmcid", "").strip()
        pmid = gold_row.get("pmid", "").strip()
        quote = norm_space(gold_row.get("quote", ""))
        found_source = False
        found_sentence = False
        quote_verified = False
        source_key: tuple[str, str, str] | None = None
        for (sp, db, uid), article in article_index.items():
            if sp != species:
                continue
            if (pmcid and article.get("pmcid") == pmcid) or (
                pmid and article.get("pmid") == pmid
            ):
                found_source = True
                source_key = (sp, db, uid)
                text = source_text_for(sp, article, cache)
                quote_verified = quote_in_text(quote, text)
                break
        for rec in sentence_index.get(slugify(species), []):
            if source_key and (rec["source_type"], rec["uid"]) != (
                source_key[1],
                source_key[2],
            ):
                continue
            if quote_in_text(quote, rec.get("text", "")):
                found_sentence = True
                break
        rows.append(
            {
                **gold_row,
                "found_source": "yes" if found_source else "no",
                "found_sentence": "yes" if found_sentence else "no",
                "quote_verified": "yes" if quote_verified else "no",
                "checked_utc": utc_now(),
            }
        )
    write_tsv(
        out / "gold_check.tsv",
        rows,
        list(gold[0].keys())
        + ["found_source", "found_sentence", "quote_verified", "checked_utc"],
    )
    print(
        f"gold_check.tsv: {sum(1 for r in rows if r['quote_verified'] == 'yes')}"
        f"/{len(rows)} quotes verified; "
        f"{sum(1 for r in rows if r['found_sentence'] == 'yes')}/{len(rows)} "
        "found as candidate sentences"
    )


def gold_recall_rows(
    gold: list[dict[str, str]],
    candidates: list[dict[str, str]],
    articles: list[dict[str, str]],
    sentences: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    """Stage-by-stage recall of known sources. For each gold row (species, level, pmid, pmcid,
    phrase) report whether the pipeline

      retrieved   the article among the search candidates,
      selected    it for fetching by the search ranking alone (an article forced in with
                  --include-pmcid, priority -1, does not count),
      sentence    gave at least one candidate sentence from it, and
      phrase      kept the gold phrase in a candidate sentence or its context.

    `gold` is not circular: a gold source must be found without being forced in.
    """
    out: list[dict[str, Any]] = []
    for g in gold:
        sp = g["species"].strip()
        pmid = g.get("pmid", "").strip()
        pmcid = normalize_pmcid(g.get("pmcid", "")) if g.get("pmcid", "").strip() else ""
        ids = {x for x in (pmid, pmcid) if x}
        retrieved = any(c["species"] == sp and c["uid"] in ids for c in candidates)
        arts = [
            a
            for a in articles
            if a["species"] == sp
            and ((pmid and a.get("pmid") == pmid) or (pmcid and a.get("pmcid") == pmcid))
        ]
        selected = any(str(a.get("priority")) != "-1" for a in arts)
        forced_only = bool(arts) and not selected
        recs = [
            s
            for s in sentences.get(slugify(sp), [])
            if (pmid and s.get("pmid") == pmid) or (pmcid and s.get("pmcid") == pmcid)
        ]
        phrase = g.get("phrase", "").strip()
        found_phrase = bool(phrase) and any(
            quote_in_text(
                phrase,
                " ".join(
                    [s.get("context_before", ""), s.get("text", ""), s.get("context_after", "")]
                ),
            )
            for s in recs
        )
        out.append(
            {
                **g,
                "retrieved": "yes" if retrieved else "no",
                "selected": "yes" if selected else ("forced_only" if forced_only else "no"),
                "sentence": "n/a" if not phrase else ("yes" if recs else "no"),
                "phrase": "n/a" if not phrase else ("yes" if found_phrase else "no"),
            }
        )
    return out


def cmd_recall(args: argparse.Namespace) -> None:
    study = Path(args.study_dir)
    out = Path(args.out) if args.out else study / "growth_form"
    gold = read_tsv(Path(args.gold) if args.gold else out / "gold_sources.tsv")
    sentences: dict[str, list[dict[str, Any]]] = {}
    sent_dir = out / "sentences"
    if sent_dir.exists():
        for path in sent_dir.glob("*.jsonl"):
            with open(path) as fh:
                sentences[path.stem] = [json.loads(x) for x in fh if x.strip()]
    rows = gold_recall_rows(
        gold,
        read_tsv(out / "candidates.tsv"),
        read_tsv(out / "articles.tsv"),
        sentences,
    )
    write_tsv(
        Path(args.recall_out) if args.recall_out else out / "gold_recall.tsv",
        rows,
        list(rows[0].keys()),
    )
    n = len(rows)
    for stage in ("retrieved", "selected", "sentence", "phrase"):
        k = sum(1 for r in rows if r[stage] == "yes")
        n_stage = sum(1 for r in rows if r[stage] != "n/a")
        print(f"{stage:10s} {k}/{n_stage}")
    forced = sum(1 for r in rows if r["selected"] == "forced_only")
    if forced:
        print(f"{forced} gold sources were only present because they were forced in")


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------


def cmd_provenance(args: argparse.Namespace) -> None:
    study = Path(args.study_dir)
    out = Path(args.out) if args.out else study / "growth_form"
    taxonomy = Path(args.taxonomy_dir)
    generated = utc_now()
    records: list[dict[str, Any]] = []

    common_extra = {
        "script": "bin/gf_lit.py",
        "script_sha256": sha256_of(Path(__file__)),
        "script_version": SCRIPT_VERSION,
        "generated_utc": generated,
        "taxonomy_dir": str(taxonomy),
        "eutils_base": EUTILS,
        "ncbi_tool": TOOL,
    }
    for name, source_url, release, license_text, derived_by in (
        (
            "species.tsv",
            str(study / "*/species.csv"),
            "study state at generation",
            "Internal / project data",
            "bin/gf_lit.py species",
        ),
        (
            "names.tsv",
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
            "NCBI Taxonomy 2026-01-28",
            "NCBI Taxonomy data (public)",
            "bin/gf_lit.py names",
        ),
        (
            "queries.tsv",
            EUTILS,
            "NCBI E-utilities (live queries)",
            "NCBI data use policy",
            "bin/gf_lit.py search",
        ),
        (
            "candidates.tsv",
            EUTILS,
            "NCBI E-utilities (live queries)",
            "NCBI data use policy",
            "bin/gf_lit.py search",
        ),
        (
            "articles.tsv",
            EUTILS,
            "NCBI E-utilities / PMC OA subset",
            "Mixed; see each article",
            "bin/gf_lit.py fetch",
        ),
        (
            "sentences_summary.tsv",
            EUTILS,
            "NCBI E-utilities / PMC OA subset",
            "Mixed; see each article",
            "bin/gf_lit.py sentences",
        ),
    ):
        path = out / name
        if not path.exists():
            continue
        records.append(
            build_record(
                source_url=source_url,
                source_release=release,
                license=license_text,
                local_path=path,
                derived_by=derived_by,
                extra=common_extra,
            )
        )
    for name, source_url, release, license_text in (
        ("aliases.tsv", "curated://growth_form/aliases.tsv", "curated seed", "Internal"),
        ("gold.tsv", "curated://growth_form/gold.tsv", "curated seed", "Internal"),
    ):
        path = out / name
        if path.exists():
            records.append(
                build_record(
                    source_url=source_url,
                    source_release=release,
                    license=license_text,
                    local_path=path,
                    derived_by="hand-curated",
                    extra={"script_version": SCRIPT_VERSION},
                )
            )
    path = out / "PROVENANCE.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "generated_utc": generated,
                "script": "bin/gf_lit.py",
                "script_sha256": sha256_of(Path(__file__)),
                "script_version": SCRIPT_VERSION,
                "taxonomy_dir": str(taxonomy),
                "taxonomy_nodes_sha256": sha256_of(taxonomy / "nodes.dmp")
                if (taxonomy / "nodes.dmp").exists()
                else "",
                "taxonomy_names_sha256": sha256_of(taxonomy / "names.dmp")
                if (taxonomy / "names.dmp").exists()
                else "",
                "records": records,
            },
            sort_keys=False,
        )
    )
    print(f"PROVENANCE.yaml: {len(records)} records; {path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--study-dir",
        default=str(DEFAULT_STUDY),
        help="study directory (default: %(default)s)",
    )
    parser.add_argument(
        "--out",
        default="",
        help="growth_form output directory (default: <study-dir>/growth_form)",
    )
    parser.add_argument(
        "--taxonomy-dir",
        default=str(DEFAULT_TAXONOMY),
        help="directory with nodes.dmp and names.dmp (default: %(default)s)",
    )
    parser.add_argument(
        "--species",
        nargs="+",
        default=None,
        help="limit to these species names",
    )
    parser.add_argument(
        "--email",
        default=os.environ.get("NCBI_EMAIL", ""),
        help="contact email for NCBI E-utilities (default: $NCBI_EMAIL)",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Literature retrieval and quote verification for growth-form coding."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("species", help="build species.tsv")
    add_common(p)
    p.set_defaults(func=cmd_species)

    p = sub.add_parser("names", help="build names.tsv")
    add_common(p)
    p.set_defaults(func=cmd_names)

    p = sub.add_parser("search", help="run PubMed/PMC searches")
    add_common(p)
    p.add_argument("--max-per-query", type=int, default=200)
    p.add_argument(
        "--sort",
        default="relevance",
        help="E-utilities sort order (default: relevance; empty for default)",
    )
    p.set_defaults(func=cmd_search)

    p = sub.add_parser("fetch", help="fetch abstracts and PMC full text")
    add_common(p)
    p.add_argument(
        "--selection",
        choices=["score", "rank", "two-step"],
        default="score",
        help="score: merge duplicates, order by (priority, -n_queries, rank); "
        "rank: earlier order (priority, rank, db, uid); "
        "two-step: fetch abstracts of the top --screen-articles candidates, order them by a "
        "title/abstract screen score, then fetch full text for the best",
    )
    p.add_argument("--screen-articles", type=int, default=300)
    p.add_argument(
        "--base-share",
        type=float,
        default=0.5,
        help="two-step: share of the article slots filled from the score order; the rest follow the "
        "title/abstract screen order",
    )
    p.add_argument("--max-articles", type=int, default=200)
    p.add_argument("--max-fulltext", type=int, default=100)
    p.add_argument(
        "--include-pmcid",
        nargs="+",
        default=None,
        help=(
            "force PMC IDs into the selection; use PMCID for all requested "
            "species, or Species=PMCID for one species"
        ),
    )
    p.set_defaults(func=cmd_fetch)

    p = sub.add_parser("sentences", help="extract growth-form sentences")
    add_common(p)
    p.add_argument(
        "--article-level",
        choices=["off", "title", "abstract"],
        default="abstract",
        help="for an article about the species (full name in the title, or in the title or abstract), "
        "also keep growth sentences that do not name it (match_level 'article')",
    )
    p.add_argument(
        "--article-cap",
        type=int,
        default=12,
        help="at most this many article-level sentences per article, dimorphism sentences first",
    )
    p.set_defaults(func=cmd_sentences)

    p = sub.add_parser("verify", help="verify evidence quotes")
    add_common(p)
    p.add_argument(
        "--evidence",
        default="",
        help="evidence TSV (default: <out>/evidence.tsv)",
    )
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser("gold", help="check gold.tsv")
    add_common(p)
    p.set_defaults(func=cmd_gold)

    p = sub.add_parser(
        "recall",
        help="stage recall of gold sources (retrieved, selected, sentence, phrase)",
    )
    add_common(p)
    p.add_argument("--gold", help="gold_sources.tsv (default <out>/gold_sources.tsv)")
    p.add_argument("--recall-out", default="", help="output TSV (default <out>/gold_recall.tsv)")
    p.set_defaults(func=cmd_recall)

    p = sub.add_parser("provenance", help="write PROVENANCE.yaml")
    add_common(p)
    p.set_defaults(func=cmd_provenance)

    args = parser.parse_args(argv)
    if args.command == "verify" and not args.evidence:
        args.evidence = str(
            (Path(args.out) if args.out else Path(args.study_dir) / "growth_form")
            / "evidence.tsv"
        )
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
