#!/usr/bin/env python3
"""Check the yeast-convergence batch specs and render one NII sub-study per contrast.

The batch YAMLs in ../batches/ are the single source of truth for which species go
into which pairwise run. They stay compatible with NovInvenio's
bin/build_targeted_configs.py (batch / outgroup_pools / studies[focal,
ingroup_extra, outgroup_pool]); the extra keys this script reads (study `id`,
`origin`, `contrast`, `notes`, and the top-level `species` map) are ignored there.

Subcommands
-----------
check   Validate every batch: unique study ids, focal not in its own outgroup,
        explicit members disjoint from the outgroup, every species has a
        taxon_group and a legal growth_form, and the growth forms fit the
        study's role (origin / reversal / pseudo-origin). Optional:
        --master-pool  report species missing from NovInvenio's master pool
        --traits       report species with no or a conflicting cell_morphology row
render  Write <folder>/<id>/species.csv + run_params.txt for every study, plus
        <folder>/study_index.tsv (read by run_all.sh, run_enrichment_all.sh and
        convergence_compare.py). Source columns are prefilled only where another
        tracked studies/*/*/species.csv already uses exactly one UniProt
        proteome for that species; every other row is left blank for
        `bin/ni resolve`. An existing species.csv is never overwritten
        unless --force is given (it may hold resolve output or hand edits).

Paths default to the cluster checkout (absolute, per CLAUDE.md "Study-specific
vs. shared scripts"); pass --nii-root to run elsewhere.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import yaml

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
FOLDER_REL = Path("studies/fungi/yeast_convergence")

SPECIES_FIELDS = [
    "Short", "Species", "Strain", "Group", "TaxonGroup", "Protein_Source",
    "Protein_Accession", "Taxon_ID", "Genome_Source", "Genome_Accession",
    "GFF3_Source", "GFF3_Accession",
]
INDEX_FIELDS = ["id", "batch", "role", "origin", "contrast", "compare_with", "n_in", "n_out",
                "n_prefilled", "n_to_resolve"]
RUN_PARAMS = """\
# Extra nextflow params for this study, applied by bin/run_study.sh before any
# args passed on the command line. One flag/value per line.
# Written by studies/fungi/yeast_convergence/bin/yc_batches.py render.

# diamond for pairwise search and self-search, same as every other
# yeast_convergence contrast so the runs stay comparable.
--run_tool diamond
"""
DOMINANT_PHASES = {"yeast", "hyphal"}
HYPHAL_STATES = {"filamentous", "dimorphic_hyphal"}


class SpecError(Exception):
    pass


def life_state(growth_form: str, dominant_phase: str | None) -> str:
    """Collapse growth_form + dominant_phase into the 3-state coding:
    filamentous / dimorphic_hyphal (yeast capacity, hypha dominant) /
    yeast_dominant (yeast, or dimorphic with the yeast phase dominant)."""
    if growth_form == "filamentous":
        return "filamentous"
    if growth_form == "yeast" or dominant_phase == "yeast":
        return "yeast_dominant"
    return "dimorphic_hyphal"


def role_of(contrast: str) -> str:
    if contrast.startswith("reversal"):
        return "reversal"
    if contrast.startswith("pseudo_origin"):
        return "pseudo"
    return "origin"


def load_growth_forms(trait_definitions: Path) -> set[str]:
    with open(trait_definitions) as fh:
        defs = yaml.safe_load(fh)
    return set(defs["traits"]["cell_morphology"]["values"])


def load_batches(batch_dir: Path) -> list[dict]:
    paths = sorted(batch_dir.glob("*.yaml"))
    if not paths:
        raise SpecError(f"no batch YAMLs in {batch_dir}")
    batches = []
    for p in paths:
        with open(p) as fh:
            spec = yaml.safe_load(fh)
        spec["_path"] = p
        batches.append(spec)
    return batches


def iter_studies(batches):
    """Yield (batch, study, in_species, out_species) for every study."""
    for b in batches:
        pools = {name: list(cfg["members"]) for name, cfg in b.get("outgroup_pools", {}).items()}
        for s in b["studies"]:
            extra = s.get("ingroup_extra") or {}
            if extra.get("mode") != "explicit":
                raise SpecError(f"{b['_path'].name}:{s.get('id')}: only mode: explicit is supported here")
            if s["outgroup_pool"] not in pools:
                raise SpecError(f"{b['_path'].name}:{s.get('id')}: outgroup_pool {s['outgroup_pool']!r} not defined")
            ins = [s["focal"]] + list(extra.get("members") or [])
            yield b, s, ins, pools[s["outgroup_pool"]]


def check(batches, growth_forms: set[str]) -> list[str]:
    """Return a list of problems (empty list = all good)."""
    problems = []
    seen_ids: dict[str, str] = {}
    for b in batches:
        for key in ("batch", "outgroup_pools", "studies", "species"):
            if key not in b:
                problems.append(f"{b['_path'].name}: missing top-level key {key!r}")
    if problems:
        return problems

    for b, s, ins, outs in iter_studies(batches):
        where = f"{b['_path'].name}:{s.get('id', s['focal'])}"
        sid = s.get("id")
        if not sid or not re.fullmatch(r"yc_[a-z0-9_]+", sid):
            problems.append(f"{where}: study id must match yc_[a-z0-9_]+ (got {sid!r})")
        elif sid in seen_ids:
            problems.append(f"{where}: duplicate study id (also in {seen_ids[sid]})")
        else:
            seen_ids[sid] = b["_path"].name
        for key in ("origin", "contrast"):
            if not s.get(key):
                problems.append(f"{where}: missing {key!r}")
        if role_of(s.get("contrast", "")) == "reversal" and not s.get("compare_with"):
            problems.append(f"{where}: reversal study needs compare_with (an origin study id)")
        if not (s.get("ingroup_extra") or {}).get("reason"):
            problems.append(f"{where}: explicit mode needs a reason")
        if len(set(ins)) != len(ins):
            problems.append(f"{where}: duplicate species in ingroup")
        if len(set(outs)) != len(outs):
            problems.append(f"{where}: duplicate species in outgroup pool")
        overlap = set(ins) & set(outs)
        if overlap:
            problems.append(f"{where}: species in both IN and OUT: {sorted(overlap)}")

        species_map = b["species"]
        states = {}
        for sp in ins + outs:
            meta = species_map.get(sp)
            if not meta:
                problems.append(f"{where}: {sp!r} has no entry in the batch's species map")
                continue
            if not meta.get("taxon_group"):
                problems.append(f"{where}: {sp!r} has no taxon_group")
            gf = meta.get("growth_form")
            if gf not in growth_forms:
                problems.append(f"{where}: {sp!r} growth_form {gf!r} not in {sorted(growth_forms)}")
                continue
            phase = meta.get("dominant_phase")
            if gf == "dimorphic" and phase not in DOMINANT_PHASES:
                problems.append(f"{where}: dimorphic {sp!r} needs dominant_phase in {sorted(DOMINANT_PHASES)}")
                continue
            if gf != "dimorphic" and phase is not None:
                problems.append(f"{where}: {sp!r} has dominant_phase but is not dimorphic")
            states[sp] = life_state(gf, phase)
        if any(sp not in states for sp in ins + outs):
            continue

        role = role_of(s.get("contrast", ""))
        contrast = s.get("contrast", "")
        in_states = {states[sp] for sp in ins}
        out_states = {states[sp] for sp in outs}
        if role == "origin":
            # A capacity origin (contrast dimorphic_*) only needs a yeast phase in IN;
            # every other origin needs at least one yeast-dominant IN species.
            wanted = {"yeast_dominant", "dimorphic_hyphal"} if contrast.startswith("dimorphic") else {"yeast_dominant"}
            if not in_states & wanted:
                problems.append(f"{where}: origin IN has no {' or '.join(sorted(wanted))} species")
            if not out_states <= HYPHAL_STATES:
                problems.append(f"{where}: origin OUT must be filamentous or hypha-dominant "
                                f"(got {sorted(out_states)})")
        elif role == "reversal":
            if in_states != {"filamentous"}:
                problems.append(f"{where}: reversal IN must be filamentous (got {sorted(in_states)})")
            if "yeast_dominant" not in out_states:
                problems.append(f"{where}: reversal OUT has no yeast-dominant species")
        elif role == "pseudo":
            if in_states | out_states != {"filamentous"}:
                problems.append(f"{where}: pseudo-origin control must be filamentous on both sides")
    reversal_targets = [(b["_path"].name, s["id"], s["compare_with"])
                        for b, s, _, _ in iter_studies(batches) if s.get("compare_with")]
    for fname, sid, target in reversal_targets:
        if target not in seen_ids:
            problems.append(f"{fname}:{sid}: compare_with {target!r} is not a study id")
    return problems


def check_master_pool(batches, master_pool: Path) -> list[str]:
    with open(master_pool, newline="") as fh:
        have = {row["Species"] for row in csv.DictReader(fh)}
    missing = set()
    for _, _, ins, outs in iter_studies(batches):
        missing |= {sp for sp in ins + outs if sp not in have}
    return sorted(missing)


def check_traits(batches, traits_csv: Path) -> list[str]:
    curated: dict[str, set[str]] = {}
    with open(traits_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["trait"] == "cell_morphology":
                curated.setdefault(row["Species"], set()).add(row["value"])
    notes = []
    for b in batches:
        for sp, meta in sorted(b["species"].items()):
            have = curated.get(sp)
            if not have:
                notes.append(f"no curated cell_morphology: {sp} (batch says {meta['growth_form']})")
            elif meta["growth_form"] not in have and not (
                    # a curated "dimorphic" agrees with a yeast-dominant dimorph here
                    "dimorphic" in have and meta.get("dominant_phase")):
                notes.append(f"CONFLICT {sp}: batch says {meta['growth_form']}, traits.csv says {sorted(have)}")
    return sorted(set(notes))


def make_short(species: str, used: set[str]) -> str:
    """Same rule as NovInvenio lib/master_pool.py make_short (genus[:4] + epithet[:4])."""
    parts = species.split()
    genus = re.sub(r"[^A-Za-z]", "", parts[0])[:4] if parts else "Sp"
    epithet = re.sub(r"[^A-Za-z]", "", parts[1])[:4] if len(parts) > 1 else ""
    base = (genus[:1].upper() + genus[1:].lower() + epithet.lower())[:8] or "Sp"
    short, n = base, 2
    while short in used:
        short = f"{base[:8 - len(str(n))]}{n}"
        n += 1
    used.add(short)
    return short


def known_sources(studies_root: Path, exclude: Path) -> dict[str, dict]:
    """Species -> one prefilled row, from tracked species.csv files that use a
    single UniProt proteome for it. Species with conflicting accessions across
    studies are dropped (left for bin/ni resolve)."""
    found: dict[str, dict[str, dict]] = {}
    for p in sorted(studies_root.glob("*/*/species.csv")):
        if exclude in p.parents:
            continue
        with open(p, newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("Protein_Source") != "uniprot" or not row.get("Protein_Accession"):
                    continue
                found.setdefault(row["Species"], {})[row["Protein_Accession"]] = row
    return {sp: next(iter(rows.values())) for sp, rows in found.items() if len(rows) == 1}


def render(batches, folder: Path, known: dict[str, dict], force: bool) -> list[dict]:
    index = []
    for b, s, ins, outs in iter_studies(batches):
        sid = s["id"]
        species_map = b["species"]
        used: set[str] = set()
        rows, prefilled = [], 0
        for group, members in (("IN", ins), ("OUT", outs)):
            for sp in members:
                k = known.get(sp)
                short = k["Short"] if k and k["Short"] not in used else make_short(sp, used)
                used.add(short)
                row = dict.fromkeys(SPECIES_FIELDS, "")
                row.update(Short=short, Species=sp, Group=group,
                           TaxonGroup=species_map[sp]["taxon_group"])
                if k:
                    prefilled += 1
                    for f in ("Strain", "Protein_Source", "Protein_Accession", "Taxon_ID",
                              "Genome_Source", "Genome_Accession"):
                        row[f] = k.get(f, "")
                rows.append(row)

        study_dir = folder / sid
        study_dir.mkdir(parents=True, exist_ok=True)
        species_csv = study_dir / "species.csv"
        if species_csv.exists() and not force:
            print(f"skip {species_csv} (exists; --force to overwrite)", file=sys.stderr)
        else:
            with open(species_csv, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=SPECIES_FIELDS, lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
        run_params = study_dir / "run_params.txt"
        if not run_params.exists() or force:
            run_params.write_text(RUN_PARAMS)
        index.append({
            "id": sid, "batch": b["batch"], "role": role_of(s["contrast"]),
            "origin": s["origin"], "contrast": s["contrast"], "compare_with": s.get("compare_with", ""),
            "n_in": len(ins), "n_out": len(outs),
            "n_prefilled": prefilled, "n_to_resolve": len(rows) - prefilled,
        })

    with open(folder / "study_index.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=INDEX_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(index)
    return index


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["check", "render"])
    ap.add_argument("--nii-root", type=Path, default=NII_ROOT)
    ap.add_argument("--master-pool", type=Path, help="NovInvenio config_support/master_pool.csv (check only)")
    ap.add_argument("--traits", type=Path, help="config_support/traits/traits.csv (check only; default: NII copy)")
    ap.add_argument("--force", action="store_true", help="render: overwrite existing species.csv/run_params.txt")
    args = ap.parse_args(argv)

    folder = args.nii_root / FOLDER_REL
    trait_defs = args.nii_root / "config_support/traits/trait_definitions.yaml"
    try:
        batches = load_batches(folder / "batches")
        problems = check(batches, load_growth_forms(trait_defs))
    except SpecError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    if problems:
        for p in problems:
            print(f"ERROR: {p}", file=sys.stderr)
        return 1

    if args.command == "check":
        n = sum(len(b["studies"]) for b in batches)
        print(f"ok: {len(batches)} batches, {n} studies")
        if args.master_pool:
            missing = check_master_pool(batches, args.master_pool)
            print(f"master pool: {len(missing)} species missing" + (": " + "; ".join(missing) if missing else ""))
        traits = args.traits or (args.nii_root / "config_support/traits/traits.csv")
        if traits.exists():
            for note in check_traits(batches, traits):
                print(f"traits: {note}")
        return 0

    known = known_sources(args.nii_root / "studies", exclude=folder)
    index = render(batches, folder, known, args.force)
    for r in index:
        print(f"{r['id']:<28} {r['role']:<9} IN={r['n_in']} OUT={r['n_out']} "
              f"prefilled={r['n_prefilled']} to_resolve={r['n_to_resolve']}")
    print(f"wrote {folder / 'study_index.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
