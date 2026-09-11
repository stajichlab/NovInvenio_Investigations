# Unified Study-Onboarding Source Dispatch Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `bin/build_study_config.py`'s UniProt/NCBI-only pipeline with a
unified per-species source dispatch (`uniprot | ncbi | local_faa | local_genome |
local_gff3`), migrate every existing study onto it (except `UHM_Akkermansia`, kept
study-specific), add an onboarding skill, and use it to build the new
`UHM_lachnoNovelclade` bacteria study.

**Architecture:** `species.csv` gains independent `Protein_Source`/`Protein_Accession`,
`Genome_Source`/`Genome_Accession`, `GFF3_Source`/`GFF3_Accession` column triplets.
`bin/build_study_config.py` dispatches each triplet to one of: a UniProt fetch, an
NCBI Datasets fetch, or a plain local-file copy with a provenance record — instead of
assuming every species is UniProt+NCBI. `Short` is the stem for anything not sourced
from UniProt (UniProt rows keep the existing `{Proteome_ID}_{Taxon_ID}` stem).

**Tech Stack:** Python 3.12, `pytest` (via `pixi run test` / `pixi run pytest tests/`),
existing `lib/provenance.py` helpers, no new dependencies.

**Spec:** `notes/superpowers/specs/2026-09-11-study-onboarding-design.md`

## Global Constraints

- No backwards-compat shim for the old `UniProt_Proteome_ID`/`GCA_Accession`
  `species.csv` schema — every study migrates in the same change (spec's migration
  scope section).
- `UHM_Akkermansia` is explicitly **not** migrated — its custom
  `build_akkermansia_config.py` stays as-is (spec's Akkermansia section).
- Annotation extraction (`extract_dat_annotations.py`) only runs for rows with
  `Protein_Source == "uniprot"`.
- `GFF3_Source` semantics (this plan's refinement of the spec, since the spec left
  the exact default behavior open): blank/unset means **auto** — if
  `Genome_Source == "ncbi"`, use whatever GFF3 that NCBI Datasets package included
  (today's behavior); otherwise leave the `config.csv` GFF3 cell empty.
  `GFF3_Source == "local_gff3"` always uses `GFF3_Accession`'s path, regardless of
  `Genome_Source`. `GFF3_Source == "none"` forces the cell empty even when
  `Genome_Source == "ncbi"` provided one.
- Every `local_*` copy gets a `build_record(...)` provenance entry
  (`lib/provenance.py`) folded into the study's `DATA_MANIFEST.yaml`, same as every
  other tracked file in this repo (`CLAUDE.md`'s provenance rule).
- Never commit a changed/new tracked data file without a provenance record — this
  applies to every `species.csv` edit in this plan (they're all mechanical
  column transforms of already-provenanced data, not new sourcing, so no new
  provenance entries are needed for the rename tasks — only Tasks 3 and 5, which
  introduce genuinely new local-file sourcing, need new `DATA_MANIFEST.yaml`
  entries, which `build_study_config.py` itself writes).
- Run tests via `pixi run pytest tests/test_build_study_config.py -v` (this repo's
  existing convention — see `pixi.toml`'s `[tasks] test`).

---

### Task 1: Rewrite `bin/build_study_config.py`'s dispatch logic

**Files:**
- Modify: `bin/build_study_config.py`
- Test: `tests/test_build_study_config.py` (new)

**Interfaces:**
- Produces: `main()` unchanged CLI surface (`--study-dir`, `--uniprot-cache`,
  `--ncbi-cache`, `--skip-fetch`) plus new `--local-license` (default
  `"Internal / unpublished (not yet released outside this project)"`).
- Produces (for later tasks): the new required `species.csv` header —
  `Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID,Genome_Source,Genome_Accession,GFF3_Source,GFF3_Accession`
  (`Taxon_ID` stays required-but-only-meaningful-for-uniprot-rows, exactly as today;
  blank for every other `Protein_Source`).
- Produces: `config.csv` schema unchanged
  (`GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup`).

This task rewrites the module end-to-end. Read the current file first
(`bin/build_study_config.py`) so you preserve its `run()`/`gunzip_to()`/
`load_provenance()` helpers unchanged — only the per-row body and `main()`'s argument
list change.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_build_study_config.py`:

```python
import csv
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bin"))
import build_study_config as bsc  # noqa: E402


def _write_species_csv(path, rows, header=None):
    header = header or [
        "Short", "Species", "Strain", "Group", "TaxonGroup",
        "Protein_Source", "Protein_Accession", "Taxon_ID",
        "Genome_Source", "Genome_Accession",
        "GFF3_Source", "GFF3_Accession",
    ]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        w.writerows(rows)


def test_local_faa_and_local_genome_row(tmp_path, monkeypatch):
    # Arrange: a study dir whose only row sources everything from local files.
    study_dir = tmp_path / "studies" / "bacteria" / "toy_study"
    study_dir.mkdir(parents=True)
    local_faa = tmp_path / "src" / "Sp1.faa"
    local_faa.parent.mkdir(parents=True)
    local_faa.write_text(">seq1\nMAAA\n")
    local_genome = tmp_path / "src" / "Sp1.fna"
    local_genome.write_text(">contig1\nACGT\n")

    _write_species_csv(study_dir / "species.csv", [{
        "Short": "Sp1", "Species": "Test species", "Strain": "T1",
        "Group": "IN", "TaxonGroup": "TestGroup",
        "Protein_Source": "local_faa", "Protein_Accession": str(local_faa),
        "Taxon_ID": "",
        "Genome_Source": "local_genome", "Genome_Accession": str(local_genome),
        "GFF3_Source": "", "GFF3_Accession": "",
    }])

    monkeypatch.setattr(sys, "argv", [
        "build_study_config.py",
        "--study-dir", str(study_dir),
        "--uniprot-cache", str(tmp_path / "data/uniprot"),
        "--ncbi-cache", str(tmp_path / "data/ncbi"),
    ])

    # Act
    rc = bsc.main()

    # Assert
    assert rc == 0
    config_csv = study_dir / "config.csv"
    with open(config_csv, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 1
    row = rows[0]
    assert row["Short"] == "Sp1"
    assert row["Protein"] == "Sp1.pep.fa"
    assert row["DNA"] == "Sp1.dna.fa"
    assert row["GFF3"] == ""
    assert (study_dir / "data_dir" / "pep" / "Sp1.pep.fa").read_text() == ">seq1\nMAAA\n"
    assert (study_dir / "data_dir" / "dna" / "Sp1.dna.fa").read_text() == ">contig1\nACGT\n"
    manifest = (study_dir / "DATA_MANIFEST.yaml").read_text()
    assert "local_path:" in manifest and str(local_faa) not in manifest  # path is the source, recorded as source_url not local_path
    assert "source_url:" in manifest


def test_local_faa_with_ncbi_genome_row(tmp_path, monkeypatch):
    # Arrange: protein already local, genome must be "fetched" -- mock the
    # subprocess call and pre-seed the ncbi cache the way fetch_genome_assembly.py
    # would, so we test build_study_config.py's own cache-reading logic, not the
    # network fetch itself.
    study_dir = tmp_path / "studies" / "bacteria" / "toy_study2"
    study_dir.mkdir(parents=True)
    local_faa = tmp_path / "src" / "Sp2.faa"
    local_faa.parent.mkdir(parents=True)
    local_faa.write_text(">seq2\nMBBB\n")

    ncbi_cache = tmp_path / "data/ncbi"
    genome_dir = ncbi_cache / "GCF_000000001.1" / "extracted" / "ncbi_dataset" / "data" / "GCF_000000001.1"
    genome_dir.mkdir(parents=True)
    (genome_dir / "genomic.fna").write_text(">contig2\nTTTT\n")

    _write_species_csv(study_dir / "species.csv", [{
        "Short": "Sp2", "Species": "Test species 2", "Strain": "T2",
        "Group": "OUT", "TaxonGroup": "TestGroup",
        "Protein_Source": "local_faa", "Protein_Accession": str(local_faa),
        "Taxon_ID": "",
        "Genome_Source": "ncbi", "Genome_Accession": "GCF_000000001.1",
        "GFF3_Source": "", "GFF3_Accession": "",
    }])

    calls = []

    def fake_run(cmd):
        calls.append(cmd)

    monkeypatch.setattr(bsc, "run", fake_run)
    monkeypatch.setattr(sys, "argv", [
        "build_study_config.py",
        "--study-dir", str(study_dir),
        "--uniprot-cache", str(tmp_path / "data/uniprot"),
        "--ncbi-cache", str(ncbi_cache),
    ])

    rc = bsc.main()

    assert rc == 0
    assert any("fetch_genome_assembly.py" in str(c) for call in calls for c in call)
    config_csv = study_dir / "config.csv"
    with open(config_csv, newline="") as fh:
        row = next(csv.DictReader(fh))
    assert row["DNA"] == "Sp2.dna.fa"
    assert (study_dir / "data_dir" / "dna" / "Sp2.dna.fa").read_text() == ">contig2\nTTTT\n"


def test_unknown_source_errors(tmp_path, monkeypatch):
    study_dir = tmp_path / "studies" / "bacteria" / "bad_study"
    study_dir.mkdir(parents=True)
    _write_species_csv(study_dir / "species.csv", [{
        "Short": "Sp3", "Species": "x", "Strain": "", "Group": "IN",
        "TaxonGroup": "x",
        "Protein_Source": "bogus", "Protein_Accession": "whatever",
        "Taxon_ID": "", "Genome_Source": "", "Genome_Accession": "",
        "GFF3_Source": "", "GFF3_Accession": "",
    }])
    monkeypatch.setattr(sys, "argv", [
        "build_study_config.py", "--study-dir", str(study_dir),
    ])
    try:
        bsc.main()
        assert False, "expected SystemExit"
    except SystemExit as e:
        assert "bogus" in str(e)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run pytest tests/test_build_study_config.py -v`
Expected: collection error or `AttributeError`/`ImportError` — `build_study_config.py`
doesn't have this row schema yet.

- [ ] **Step 3: Rewrite `bin/build_study_config.py`**

Replace the file's body (keep the module's `run()`, `gunzip_to()`,
`load_provenance()` helpers verbatim) with:

```python
#!/usr/bin/env python3
"""Drive a study's species.csv through per-species source dispatch and emit a
nf_NovInvenio-ready run config.

Input: <study-dir>/species.csv, columns:
    Short,Species,Strain,Group,TaxonGroup,
    Protein_Source,Protein_Accession,Taxon_ID,
    Genome_Source,Genome_Accession,
    GFF3_Source,GFF3_Accession

Protein_Source/Genome_Source/GFF3_Source are independent per row -- a species'
protein and genome can come from unrelated places (e.g. an already-local protein
FASTA paired with a genome that still needs an NCBI fetch). Recognized values:

  Protein_Source:
    uniprot    -- Protein_Accession = UniProt proteome ID (Taxon_ID also required).
                  Fetched via fetch_uniprot_proteome.py; also drives GO/Pfam/
                  InterPro/gene-name annotation extraction (extract_dat_annotations.py)
                  -- the only Protein_Source that does.
    local_faa  -- Protein_Accession = path to an existing protein FASTA. Copied in
                  directly, with a provenance record (no fetch).

  Genome_Source:
    ncbi          -- Genome_Accession = GCA/GCF assembly accession. Fetched via
                     fetch_genome_assembly.py (genome + GFF3 together).
    local_genome  -- Genome_Accession = path to an existing genome FASTA. Copied in
                     directly, with a provenance record (no fetch).
    (blank)       -- no genome for this row; DNA config.csv cell left empty.

  GFF3_Source (optional; blank means "auto" -- see below):
    local_gff3 -- GFF3_Accession = path to an existing GFF3. Always used, regardless
                  of Genome_Source.
    none       -- force the GFF3 config.csv cell empty even if Genome_Source=ncbi
                  would otherwise have supplied one.
    (blank)    -- auto: if Genome_Source=ncbi, use whatever GFF3 that NCBI package
                  included (today's behavior); otherwise empty.

Stem naming: a uniprot-sourced protein keeps the existing
"{Proteome_ID}_{Taxon_ID}" stem (self-documenting, collision-proof across
strains sharing a Short by mistake). Every other row uses Short as its stem --
Short is already required to be unique per study.

For each species this:
  1. Resolves Protein (fetch or local copy).
  2. Resolves Genome+GFF3 (fetch or local copy/copies).
  3. If Protein_Source == "uniprot", runs bin/extract_dat_annotations.py against
     the cached .dat.gz, writing <study-dir>/annotations/<Protein_Accession>.tsv
     (gitignored, regenerable, produced unconditionally -- even with --skip-fetch
     -- so bin/run_study.sh's post-run report sync always has current annotation
     to merge). Rows sourced any other way get no annotation file (unchanged from
     today's local-file studies).
  4. Writes <study-dir>/config.csv (GROUP,Species,Strain,Protein,DNA,GFF3,Short,
     TaxonGroup) with basenames only, and <study-dir>/DATA_MANIFEST.yaml
     aggregating every species' provenance record plus a record for this
     generation step itself.

Re-running with --skip-fetch rebuilds config.csv/data_dir from whatever's already
cached (uniprot/ncbi rows) without re-downloading; local_* rows are always
re-copied (there's nothing to "skip" -- they were never fetched).
"""
import argparse
import csv
import gzip
import shutil
import subprocess
import sys
import yaml
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from provenance import append_manifest, build_record, sha256_of  # noqa: E402

BIN = Path(__file__).parent
LOCAL_LICENSE_DEFAULT = "Internal / unpublished (not yet released outside this project)"


def run(cmd: list[str]) -> None:
    print(f"$ {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)


def gunzip_to(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def load_provenance(sidecar: Path) -> dict:
    with open(sidecar) as fh:
        return yaml.safe_load(fh)


def _local_copy_record(src: Path, dest: Path, license_: str) -> dict:
    return build_record(
        source_url=f"(internal -- {src})",
        source_release="local file, not independently versioned",
        license=license_,
        local_path=dest,
        checksum=sha256_of(dest),
        derived_by=f"copied from {src} via bin/build_study_config.py",
    )


def resolve_protein(row, args, pep_dir, manifest_records):
    """Returns (stem, protein_basename, dat_gz_path_or_None)."""
    short = row["Short"]
    source = row["Protein_Source"]
    accession = row["Protein_Accession"]

    if source == "uniprot":
        taxid = row["Taxon_ID"]
        if not args.skip_fetch:
            run([
                sys.executable, str(BIN / "fetch_uniprot_proteome.py"),
                "--proteome-id", accession, "--taxid", taxid,
                "--outdir", args.uniprot_cache, "--short", short,
            ])
        stem = f"{accession}_{taxid}"
        fasta_gz = Path(args.uniprot_cache) / accession / f"{stem}.fasta.gz"
        dat_gz = Path(args.uniprot_cache) / accession / f"{stem}.dat.gz"
        if not fasta_gz.exists():
            sys.exit(f"[{short}] ERROR: expected {fasta_gz} not found -- run without --skip-fetch first")
        pep_dir.mkdir(parents=True, exist_ok=True)
        pep_out = pep_dir / f"{stem}.pep.fa"
        gunzip_to(fasta_gz, pep_out)
        for sidecar in (fasta_gz, dat_gz):
            sc = sidecar.with_suffix(sidecar.suffix + ".provenance.yaml")
            if sc.exists():
                manifest_records.append(load_provenance(sc))
        return stem, pep_out.name, (dat_gz if dat_gz.exists() else None)

    if source == "local_faa":
        src = Path(accession)
        if not src.exists():
            sys.exit(f"[{short}] ERROR: expected local protein FASTA {src} not found")
        stem = short
        pep_dir.mkdir(parents=True, exist_ok=True)
        pep_out = pep_dir / f"{stem}.pep.fa"
        shutil.copyfile(src, pep_out)
        manifest_records.append(_local_copy_record(src, pep_out, args.local_license))
        return stem, pep_out.name, None

    sys.exit(f"[{short}] ERROR: unknown Protein_Source {source!r}")


def resolve_genome_and_gff3(row, args, stem, dna_dir, gff3_dir, manifest_records):
    """Returns (dna_basename, gff3_basename) -- either may be ''."""
    short = row["Short"]
    gsource = row["Genome_Source"]
    gaccession = row["Genome_Accession"]
    gffsource = row["GFF3_Source"]
    gffaccession = row["GFF3_Accession"]

    dna_out_name = ""
    gff3_out_name = ""

    if gsource == "ncbi":
        if not args.skip_fetch:
            run([
                sys.executable, str(BIN / "fetch_genome_assembly.py"),
                "--accession", gaccession, "--outdir", args.ncbi_cache, "--short", short,
            ])
        genome_dir = Path(args.ncbi_cache) / gaccession / "extracted" / "ncbi_dataset" / "data" / gaccession
        fna_candidates = sorted(genome_dir.glob("*.fna")) if genome_dir.exists() else []
        gff_candidates = sorted(genome_dir.glob("*.gff")) if genome_dir.exists() else []
        if not fna_candidates:
            sys.exit(f"[{short}] ERROR: expected genome under {genome_dir} not found -- run without --skip-fetch first")
        dna_dir.mkdir(parents=True, exist_ok=True)
        dna_out = dna_dir / f"{stem}.dna.fa"
        shutil.copyfile(fna_candidates[0], dna_out)
        dna_out_name = dna_out.name
        for f in (fna_candidates[:1] + gff_candidates[:1]):
            sc = f.with_suffix(f.suffix + ".provenance.yaml")
            if sc.exists():
                manifest_records.append(load_provenance(sc))
        if gffsource == "none":
            gff_candidates = []
        if gffsource == "" and gff_candidates:
            gff3_dir.mkdir(parents=True, exist_ok=True)
            gff3_out_name = f"{stem}.gff3"
            shutil.copyfile(gff_candidates[0], gff3_dir / gff3_out_name)
    elif gsource == "local_genome":
        src = Path(gaccession)
        if not src.exists():
            sys.exit(f"[{short}] ERROR: expected local genome FASTA {src} not found")
        dna_dir.mkdir(parents=True, exist_ok=True)
        dna_out = dna_dir / f"{stem}.dna.fa"
        shutil.copyfile(src, dna_out)
        dna_out_name = dna_out.name
        manifest_records.append(_local_copy_record(src, dna_out, args.local_license))
    elif gsource == "":
        pass
    else:
        sys.exit(f"[{short}] ERROR: unknown Genome_Source {gsource!r}")

    if gffsource == "local_gff3":
        src = Path(gffaccession)
        if not src.exists():
            sys.exit(f"[{short}] ERROR: expected local GFF3 {src} not found")
        gff3_dir.mkdir(parents=True, exist_ok=True)
        gff3_out = gff3_dir / f"{stem}.gff3"
        shutil.copyfile(src, gff3_out)
        gff3_out_name = gff3_out.name
        manifest_records.append(_local_copy_record(src, gff3_out, args.local_license))
    elif gffsource not in ("", "none"):
        sys.exit(f"[{short}] ERROR: unknown GFF3_Source {gffsource!r}")

    return dna_out_name, gff3_out_name


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study-dir", required=True, help="e.g. studies/fungi/pezizo_set1")
    ap.add_argument("--uniprot-cache", default="data/uniprot")
    ap.add_argument("--ncbi-cache", default="data/ncbi")
    ap.add_argument("--skip-fetch", action="store_true", help="Rebuild config/data_dir from an already-populated cache, no downloads")
    ap.add_argument("--local-license", default=LOCAL_LICENSE_DEFAULT)
    args = ap.parse_args()

    study_dir = Path(args.study_dir)
    species_csv = study_dir / "species.csv"
    if not species_csv.exists():
        sys.exit(f"ERROR: {species_csv} not found")

    link_dir = study_dir / "data_dir"
    pep_dir, dna_dir, gff3_dir = (link_dir / d for d in ("pep", "dna", "gff3"))

    config_rows = []
    manifest_records = []
    seen_stems: dict[str, str] = {}

    with open(species_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            short = row["Short"]
            stem, protein_name, dat_gz = resolve_protein(row, args, pep_dir, manifest_records)

            if stem in seen_stems:
                sys.exit(
                    f"ERROR: {species_csv} resolves two rows to the same stem {stem!r} "
                    f"(Short={seen_stems[stem]!r} and Short={short!r})"
                )
            seen_stems[stem] = short

            dna_name, gff3_name = resolve_genome_and_gff3(row, args, stem, dna_dir, gff3_dir, manifest_records)

            if row["Protein_Source"] == "uniprot" and dat_gz is not None:
                annotations_dir = study_dir / "annotations"
                annotations_dir.mkdir(parents=True, exist_ok=True)
                run([
                    sys.executable, str(BIN / "extract_dat_annotations.py"),
                    "--dat-gz", str(dat_gz),
                    "--output", str(annotations_dir / f"{row['Protein_Accession']}.tsv"),
                ])

            config_rows.append({
                "GROUP": row["Group"],
                "Species": row["Species"],
                "Strain": row["Strain"],
                "Protein": protein_name,
                "DNA": dna_name,
                "GFF3": gff3_name,
                "Short": short,
                "TaxonGroup": row["TaxonGroup"],
            })

    config_csv = study_dir / "config.csv"
    with open(config_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["GROUP", "Species", "Strain", "Protein", "DNA", "GFF3", "Short", "TaxonGroup"])
        w.writeheader()
        w.writerows(config_rows)

    manifest_records.append(build_record(
        source_url="(derived, no external source)",
        source_release="n/a",
        license="n/a",
        local_path=config_csv,
        derived_by=f"bin/build_study_config.py --study-dir {study_dir}",
    ))
    append_manifest(manifest_records, study_dir / "DATA_MANIFEST.yaml")

    print(f"\nWrote {config_csv} ({len(config_rows)} species)", file=sys.stderr)
    print(f"Wrote {link_dir} (pep/dna/gff3 -- pass as --data_dir to nf_NovInvenio)", file=sys.stderr)
    print(f"Wrote {study_dir / 'DATA_MANIFEST.yaml'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run pytest tests/test_build_study_config.py -v`
Expected: all 3 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add bin/build_study_config.py tests/test_build_study_config.py
git commit -m "Generalize build_study_config.py to a per-species source dispatch

Replaces the UniProt/NCBI-only pipeline with independent Protein_Source/
Genome_Source/GFF3_Source per species.csv row (uniprot/ncbi/local_faa/
local_genome/local_gff3), per notes/superpowers/specs/2026-09-11-study-onboarding-design.md.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 2: Migrate the 8 existing UniProt-only studies' `species.csv`

**Files:**
- Modify: `studies/fungi/pezizo_set1/species.csv`
- Modify: `studies/fungi/pezizo_set1_cluster/species.csv`
- Modify: `studies/fungi/mushrooms_tremella/species.csv`
- Modify: `studies/fungi/yeast_filamentous/species.csv`
- Modify: `studies/fungi/zoosporic_dikarya/species.csv`
- Modify: `studies/fungi/agaricomycetes_mmseqs/species.csv`
- Modify: `studies/fungi/agaricomycetes_novelty_discovery/species.csv`
- Modify: `studies/fungi/agaricomycetes_pairwise/species.csv`

**Interfaces:**
- Consumes: Task 1's new required header
  (`Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID,Genome_Source,Genome_Accession,GFF3_Source,GFF3_Accession`).
- Produces: nothing new for later tasks (this is pure data migration; `config.csv`/
  `data_dir` for these 8 studies are untouched, not rebuilt).

This is a mechanical column transform: `UniProt_Proteome_ID`→`Protein_Source=uniprot`
+ `Protein_Accession`, `GCA_Accession`→`Genome_Source=ncbi` + `Genome_Accession`,
`GFF3_Source`/`GFF3_Accession` both blank (today's behavior — auto from the NCBI
package). `Taxon_ID` column is kept verbatim. No fetch, no `config.csv`/`data_dir`
rebuild — these studies are already built and stay that way; only their input
recipe (`species.csv`) changes shape.

- [ ] **Step 1: Write the migration script**

Create a throwaway script (not committed — delete after use, per this repo's
YAGNI convention for one-time migrations) at
`/tmp/migrate_species_csv.py`:

```python
import csv
import sys
from pathlib import Path

STUDIES = [
    "studies/fungi/pezizo_set1",
    "studies/fungi/pezizo_set1_cluster",
    "studies/fungi/mushrooms_tremella",
    "studies/fungi/yeast_filamentous",
    "studies/fungi/zoosporic_dikarya",
    "studies/fungi/agaricomycetes_mmseqs",
    "studies/fungi/agaricomycetes_novelty_discovery",
    "studies/fungi/agaricomycetes_pairwise",
]

NEW_HEADER = [
    "Short", "Species", "Strain", "Group", "TaxonGroup",
    "Protein_Source", "Protein_Accession", "Taxon_ID",
    "Genome_Source", "Genome_Accession",
    "GFF3_Source", "GFF3_Accession",
]

REPO_ROOT = Path(sys.argv[1])

for rel in STUDIES:
    path = REPO_ROOT / rel / "species.csv"
    with open(path, newline="") as fh:
        old_rows = list(csv.DictReader(fh))
    new_rows = []
    for row in old_rows:
        upid = row["UniProt_Proteome_ID"]
        new_rows.append({
            "Short": row["Short"],
            "Species": row["Species"],
            "Strain": row["Strain"],
            "Group": row["Group"],
            "TaxonGroup": row["TaxonGroup"],
            "Protein_Source": "uniprot",
            "Protein_Accession": upid,
            "Taxon_ID": row["Taxon_ID"],
            "Genome_Source": "ncbi",
            "Genome_Accession": row["GCA_Accession"],
            "GFF3_Source": "",
            "GFF3_Accession": "",
        })
        if not upid:
            print(f"NOTE: {rel}: Short={row['Short']!r} has empty UniProt_Proteome_ID "
                  f"-- carried through as Protein_Source=uniprot with empty "
                  f"Protein_Accession; verify against this study's existing "
                  f"config.csv before relying on this row (see Task 2 Step 3).",
                  file=sys.stderr)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=NEW_HEADER)
        w.writeheader()
        w.writerows(new_rows)
    print(f"Migrated {path} ({len(new_rows)} rows)", file=sys.stderr)
```

- [ ] **Step 2: Run it**

Run: `python3 /tmp/migrate_species_csv.py /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations`
Expected: 8 "Migrated ..." lines, plus a NOTE line for any row whose
`UniProt_Proteome_ID` was already empty in the original file (e.g.
`agaricomycetes_pairwise`'s `Scom` row — see Step 3, this is a pre-existing data
question, not something this migration should paper over).

- [ ] **Step 3: Investigate and resolve any empty-`Protein_Accession` NOTE**

For each row the script flagged: open that study's existing `config.csv` and find
the same `Short`'s `Protein` column (e.g. `grep Scom studies/fungi/agaricomycetes_pairwise/config.csv`).
Compare its basename against the current `data_dir/pep/` contents
(`ls studies/fungi/agaricomycetes_pairwise/data_dir/pep/`). If the protein file
already on disk did **not** come from a UniProt fetch (i.e. no matching
`data/uniprot/<empty>/...` — which is impossible, confirming it wasn't UniProt),
determine its real source before finalizing this row: check
`DATA_MANIFEST.yaml` for that file's `source_url` (every tracked file must have a
provenance record per `CLAUDE.md`). Set `Protein_Source` to whatever that record
actually says (most likely `local_faa`, since `DATA_MANIFEST.yaml`'s `derived_by`
will name the real origin), and `Protein_Accession` to the path
`DATA_MANIFEST.yaml` names. Do not guess — if the manifest doesn't make the origin
clear, stop and ask rather than leaving a fabricated `local_faa` path that doesn't
resolve.

- [ ] **Step 4: Verify every migrated file parses and matches the old row count**

Run (per study):
```bash
python3 -c "
import csv
p = 'studies/fungi/pezizo_set1/species.csv'
with open(p) as fh:
    rows = list(csv.DictReader(fh))
assert all(r['Protein_Source'] in ('uniprot', 'local_faa') for r in rows)
assert all(r['Genome_Source'] in ('ncbi', 'local_genome', '') for r in rows)
print(f'{p}: {len(rows)} rows OK')
"
```
Expected: `OK` for all 8 files, same row count as `git show HEAD:<path> | wc -l`
minus 1 (header) for each.

- [ ] **Step 5: Clean up and commit**

```bash
rm /tmp/migrate_species_csv.py
git add studies/fungi/pezizo_set1/species.csv \
        studies/fungi/pezizo_set1_cluster/species.csv \
        studies/fungi/mushrooms_tremella/species.csv \
        studies/fungi/yeast_filamentous/species.csv \
        studies/fungi/zoosporic_dikarya/species.csv \
        studies/fungi/agaricomycetes_mmseqs/species.csv \
        studies/fungi/agaricomycetes_novelty_discovery/species.csv \
        studies/fungi/agaricomycetes_pairwise/species.csv
git commit -m "Migrate 8 UniProt-sourced studies' species.csv to the unified source schema

Mechanical column rename (UniProt_Proteome_ID/GCA_Accession -> Protein_Source=
uniprot/Genome_Source=ncbi + their *_Accession columns), per Task 2 of
notes/superpowers/plans/2026-09-11-study-onboarding-implementation.md. No
config.csv/data_dir rebuild -- these studies are already built.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 3: Migrate `UHM_Koxytoca` onto the unified dispatcher

**Files:**
- Modify: `studies/bacteria/UHM_Koxytoca/species.csv`
- Delete: `studies/bacteria/UHM_Koxytoca/bin/build_koxytoca_config.py`
- Delete: `studies/bacteria/UHM_Koxytoca/bin/fetch_koxytoca_outgroup_dna.sh`
- Modify (regenerated by Task 1's script): `studies/bacteria/UHM_Koxytoca/config.csv`,
  `studies/bacteria/UHM_Koxytoca/data_dir/`, `studies/bacteria/UHM_Koxytoca/DATA_MANIFEST.yaml`

**Interfaces:**
- Consumes: Task 1's `bin/build_study_config.py` (unmodified from Task 1).

Unlike Task 2, this study's `config.csv`/`data_dir` **are** regenerated — the old
custom script's file-naming convention (`<Accession>.faa`/`.fna` for OUT rows) will
change to the new unified stem convention (`Short`-based), so the fetch/copy has to
actually run once to produce filenames `config.csv` correctly points at.

- [ ] **Step 1: Read the current species.csv and confirm the mapping**

```bash
cat studies/bacteria/UHM_Koxytoca/species.csv
```
Every row has `Source` = `ncbi_refseq` (OUT rows) or `uhm_mag_metashot` (IN rows).
Per `studies/bacteria/UHM_Koxytoca/bin/build_koxytoca_config.py`'s own docstring,
**every** row's protein already comes from a local `.faa` — OUT from
`/bigdata/stajichlab/jpere468/klebsiella_story/koxytoca_outgroup_faa/<Accession>.faa`,
IN from `/bigdata/stajichlab/jpere468/klebsiella_story/koxytoca_ingroup_faa/<Accession>.faa`.
Genome differs: OUT genomes come from NCBI (`<Accession>` is a `GCF_*` accession),
IN genomes are local
(`/bigdata/stajichlab/jpere468/klebsiella_story/koxytoca_mags_full/<Accession>.fa`).

- [ ] **Step 2: Write the migration script**

Create `/tmp/migrate_koxytoca.py` (throwaway, not committed):

```python
import csv
from pathlib import Path

STUDY_DIR = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/bacteria/UHM_Koxytoca")
OUTGROUP_FAA_DIR = "/bigdata/stajichlab/jpere468/klebsiella_story/koxytoca_outgroup_faa"
INGROUP_FAA_DIR = "/bigdata/stajichlab/jpere468/klebsiella_story/koxytoca_ingroup_faa"
INGROUP_DNA_DIR = "/bigdata/stajichlab/jpere468/klebsiella_story/koxytoca_mags_full"

NEW_HEADER = [
    "Short", "Species", "Strain", "Group", "TaxonGroup",
    "Protein_Source", "Protein_Accession", "Taxon_ID",
    "Genome_Source", "Genome_Accession",
    "GFF3_Source", "GFF3_Accession",
]

with open(STUDY_DIR / "species.csv", newline="") as fh:
    old_rows = list(csv.DictReader(fh))

new_rows = []
for row in old_rows:
    accession = row["Accession"]
    source = row["Source"]
    if source == "ncbi_refseq":
        protein_path = f"{OUTGROUP_FAA_DIR}/{accession}.faa"
        genome_source, genome_accession = "ncbi", accession
    elif source == "uhm_mag_metashot":
        protein_path = f"{INGROUP_FAA_DIR}/{accession}.faa"
        genome_source, genome_accession = "local_genome", f"{INGROUP_DNA_DIR}/{accession}.fa"
    else:
        raise SystemExit(f"unknown Source {source!r} for Short={row['Short']!r}")
    new_rows.append({
        "Short": row["Short"], "Species": row["Species"], "Strain": row["Strain"],
        "Group": row["Group"], "TaxonGroup": row["TaxonGroup"],
        "Protein_Source": "local_faa", "Protein_Accession": protein_path, "Taxon_ID": "",
        "Genome_Source": genome_source, "Genome_Accession": genome_accession,
        "GFF3_Source": "", "GFF3_Accession": "",
    })

with open(STUDY_DIR / "species.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=NEW_HEADER)
    w.writeheader()
    w.writerows(new_rows)

print(f"Migrated {STUDY_DIR / 'species.csv'} ({len(new_rows)} rows)")
```

- [ ] **Step 3: Run it, then snapshot the old data_dir for comparison**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
python3 /tmp/migrate_koxytoca.py
sha256sum studies/bacteria/UHM_Koxytoca/data_dir/pep/*.faa \
          studies/bacteria/UHM_Koxytoca/data_dir/dna/*.fna \
          > /tmp/koxytoca_old_checksums.txt 2>/dev/null || true
```

- [ ] **Step 4: Regenerate config.csv/data_dir via the unified script**

```bash
rm -rf studies/bacteria/UHM_Koxytoca/data_dir
pixi run python bin/build_study_config.py --study-dir studies/bacteria/UHM_Koxytoca
```
Expected: no errors; `Wrote studies/bacteria/UHM_Koxytoca/config.csv (34 species)`
(or however many rows the study actually has — check `wc -l species.csv`).

- [ ] **Step 5: Verify content-identical, name-different files**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
python3 -c "
import csv
with open('studies/bacteria/UHM_Koxytoca/config.csv') as fh:
    rows = list(csv.DictReader(fh))
for r in rows:
    assert (open(f\"studies/bacteria/UHM_Koxytoca/data_dir/pep/{r['Protein']}\", 'rb').read()), r
    if r['DNA']:
        assert (open(f\"studies/bacteria/UHM_Koxytoca/data_dir/dna/{r['DNA']}\", 'rb').read()), r
print(f'{len(rows)} rows, all Protein/DNA files present and non-empty')
"
```
Expected: `<N> rows, all Protein/DNA files present and non-empty` with no
`AssertionError`.

- [ ] **Step 6: Delete the now-unused custom scripts**

```bash
rm studies/bacteria/UHM_Koxytoca/bin/build_koxytoca_config.py
rm studies/bacteria/UHM_Koxytoca/bin/fetch_koxytoca_outgroup_dna.sh
rmdir studies/bacteria/UHM_Koxytoca/bin 2>/dev/null || true
rm /tmp/migrate_koxytoca.py /tmp/koxytoca_old_checksums.txt
```

- [ ] **Step 7: Commit**

```bash
git add -A studies/bacteria/UHM_Koxytoca/
git status --short studies/bacteria/UHM_Koxytoca/  # confirm nothing unexpected staged
git commit -m "Migrate UHM_Koxytoca onto the unified build_study_config.py dispatcher

Both groups become Protein_Source=local_faa; Genome_Source=ncbi for the
outgroup (replaces the separate fetch_koxytoca_outgroup_dna.sh shim --
build_study_config.py now fetches directly), Genome_Source=local_genome
for the ingroup. Retires build_koxytoca_config.py/fetch_koxytoca_outgroup_dna.sh.
data_dir/ filenames change (Short-based instead of Accession-based for the
outgroup) but content is unchanged, per Task 3 of
notes/superpowers/plans/2026-09-11-study-onboarding-implementation.md.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 4: Add the `.claude/skills/new-study` onboarding skill

**Files:**
- Create: `.claude/skills/new-study/SKILL.md`

**Interfaces:**
- Consumes: Task 1's `species.csv` schema and `bin/build_study_config.py` CLI;
  `bin/run_study.sh`'s existing CLI (`bin/run_study.sh <domain>/<set_name> [args]`).

This is documentation/process, not code — no tests. It's the "AI-assisted judgment"
layer from the spec: classifying, per candidate species, which `Source` value
applies.

- [ ] **Step 1: Write the skill**

Create `.claude/skills/new-study/SKILL.md`:

```markdown
---
name: new-study
description: Use when onboarding a new study into studies/<domain>/<set_name>/ -- inventories candidate protein/genome directories, classifies each species' Protein_Source/Genome_Source/GFF3_Source, writes species.csv, and hands off to bin/build_study_config.py + bin/run_study.sh. Use when the user wants to start a new NII study, add a new species set, or bring in a new dataset (local FASTA directories, NCBI accessions, or UniProt proteome IDs).
---

# Onboarding a new study

Reference: `notes/superpowers/specs/2026-09-11-study-onboarding-design.md` for the
full schema rationale.

## Steps

1. **Pick `<domain>/<set_name>`.** Domain is one of `conf/domains.yaml`'s slugs
   (`fungi`, `animal`, `plant`, `bacteria`, `other`). `set_name` is a short,
   descriptive slug (matches existing studies like `pezizo_set1`,
   `UHM_lachnoNovelclade`). Create `studies/<domain>/<set_name>/`.

2. **Inventory every candidate data source the user gives you.** For each
   directory or accession list:
   - List files (`ls <dir>`). A directory of protein FASTA is a candidate
     `Protein_Source=local_faa` source; a directory of genome FASTA is a
     candidate `Genome_Source=local_genome` source.
   - For every species stem you find in a protein directory, check whether the
     *same stem* exists in a genome directory (exact filename match, not
     substring — a bin numbered `22` is not a substring match for `229`). If it
     does, that species is `local_genome`. If it doesn't, and the stem looks
     like an NCBI accession (`GCF_*`/`GCA_*`), that species is `Genome_Source=ncbi`
     with `Genome_Accession` = that accession. If neither, stop and ask the user
     where that species' genome comes from — do not guess.
   - A bare UniProt proteome ID (no local file at all) is `Protein_Source=uniprot`
     (needs a `Taxon_ID` too — look it up via UniProt if the user hasn't given
     one).

3. **Get real `Species`/`Strain`/`TaxonGroup` values from the user.** Filenames
   (MAG bin IDs, accessions) are not species names. Do not fabricate a
   Latin binomial or taxon group — ask, or point to whatever taxonomy/metadata
   file the user's source project already has (e.g. a GTDB-tk summary, a
   `groups.tsv`) and confirm your reading of it with the user before writing
   `species.csv`.

4. **Write `species.csv`** with the header:
   `Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID,Genome_Source,Genome_Accession,GFF3_Source,GFF3_Accession`
   `Group` is `IN` or `OUT`. Leave `GFF3_Source`/`GFF3_Accession` blank unless the
   user has an explicit local GFF3 to attach.

5. **Run** `pixi run python bin/build_study_config.py --study-dir studies/<domain>/<set_name>`.
   Fix any `ERROR:` it reports (missing file, unknown Source value) before moving on.

6. **Hand off**: `bin/run_study.sh <domain>/<set_name> [nextflow args]` runs the
   actual pipeline. See that script's own header comment for `NII_PIPELINE`/
   `NOVINVENIO_ROOT` local-checkout requirements.

## What this skill does NOT automate

Classifying `Source` values is judgment, not a fixed algorithm — directory
naming conventions vary per source project. This skill's job is to do that
classification carefully and show its reasoning, not to guess silently.
```

- [ ] **Step 2: Commit**

```bash
git add .claude/skills/new-study/SKILL.md
git commit -m "Add new-study onboarding skill

Wraps the source-classification judgment step (Protein_Source/Genome_Source/
GFF3_Source per species) that Task 1's unified build_study_config.py dispatch
can't do on its own, per notes/superpowers/specs/2026-09-11-study-onboarding-design.md's
onboarding-skill section.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 5: Build `UHM_lachnoNovelclade`

**Files:**
- Create: `studies/bacteria/UHM_lachnoNovelclade/species.csv`
- Create (via `bin/build_study_config.py`): `studies/bacteria/UHM_lachnoNovelclade/config.csv`,
  `studies/bacteria/UHM_lachnoNovelclade/data_dir/`,
  `studies/bacteria/UHM_lachnoNovelclade/DATA_MANIFEST.yaml`

**Interfaces:**
- Consumes: Task 1's `bin/build_study_config.py`, Task 4's onboarding skill process
  (this task is that skill's first real invocation).

**Step 1 is a hard checkpoint — do not fabricate species/strain/taxon names.** The
9 species below were identified purely by filename inventory (this plan's spec
document); none of `Species`, `Strain`, or `TaxonGroup` are known yet.

- [ ] **Step 1: Get real species identity from the user**

Ask the user (or check whether `/bigdata/stajichlab/jpere468/unknown_tree/` or
`/bigdata/stajichlab/jpere468/drep_herptile_95/` has a taxonomy/metadata file
covering these 9 stems — note `winners_with_metadata_plus_taxonomy 2.numbers`
exists at the `drep_herptile_95` root but is an Apple Numbers file; ask the user
for a plain-text export or the values directly) for `Species`, `Strain`, and
`TaxonGroup` for each of:

| Short (proposed) | stem |
|---|---|
| (assign) | `UHM1088.41098__UHM1088.41098_R.bin.56` |
| (assign) | `UHM1210.23070__UHM1210.23070_R.bin.85` |
| (assign) | `UHM896.23052__UHM896.23052_R.bin.22` |
| (assign) | `UHM904.23055__UHM904.23055_R.bin.119` |
| (assign) | `UHM207.23041__UHM207.23041_R.bin.71` |
| (assign) | `GCF_000687555.1` |
| (assign) | `GCF_002797975.1` |
| (assign) | `GCF_025149125.1` |
| (assign) | `GCF_900112885.1` |

Do not proceed to Step 2 until you have these values. A short, unique `Short`
code per row is also needed (existing studies use 4-8 character codes derived
from the genus/species, e.g. `Ccin` for *Coprinopsis cinerea* — ask the user to
confirm or supply one per row if the real species names don't make an obvious
short code).

- [ ] **Step 2: Write `species.csv`**

Once Step 1's answers are in hand, write
`studies/bacteria/UHM_lachnoNovelclade/species.csv` with header
`Short,Species,Strain,Group,TaxonGroup,Protein_Source,Protein_Accession,Taxon_ID,Genome_Source,Genome_Accession,GFF3_Source,GFF3_Accession`
and one row per stem above:

- 4 outgroup rows (from `clade1_faa`) and 1 ingroup row (`UHM207...`, from
  `lachno_ingroup_faa`): `Protein_Source=local_faa`,
  `Protein_Accession=/bigdata/stajichlab/jpere468/unknown_tree/<clade1_faa or lachno_ingroup_faa>/<stem>.faa`,
  `Genome_Source=local_genome`,
  `Genome_Accession=/bigdata/stajichlab/jpere468/drep_herptile_95/high_quality_genomes/<stem>.fa`
  (verified to exist in this plan's spec document's worked example).
- 4 ingroup rows (the `GCF_*` stems, from `lachno_ingroup_faa`):
  `Protein_Source=local_faa`,
  `Protein_Accession=/bigdata/stajichlab/jpere468/unknown_tree/lachno_ingroup_faa/<stem>.faa`,
  `Genome_Source=ncbi`, `Genome_Accession=<stem>` (the `GCF_*` accession itself).
- All 9 rows: `Taxon_ID` blank (not uniprot-sourced), `GFF3_Source`/
  `GFF3_Accession` blank.
- `Group`: `OUT` for the 4 `clade1_faa` rows, `IN` for the 5 `lachno_ingroup_faa`
  rows (matches the directory names given in the original request).

- [ ] **Step 3: Run the build**

```bash
cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
pixi run python bin/build_study_config.py --study-dir studies/bacteria/UHM_lachnoNovelclade
```
Expected: `Wrote studies/bacteria/UHM_lachnoNovelclade/config.csv (9 species)`, no
`ERROR:` lines. If any `local_faa`/`local_genome` path errors as not-found, re-check
the exact path against `ls` — do not adjust the path to "make it work" without
confirming the file it now points at is actually the right one.

- [ ] **Step 4: Verify**

```bash
cat studies/bacteria/UHM_lachnoNovelclade/config.csv
ls studies/bacteria/UHM_lachnoNovelclade/data_dir/pep studies/bacteria/UHM_lachnoNovelclade/data_dir/dna
```
Expected: 9 rows in `config.csv`, 9 files in `data_dir/pep/`, 9 files in
`data_dir/dna/` (all `Genome_Source` values resolved, none blank per this study's
design).

- [ ] **Step 5: Commit**

```bash
git add studies/bacteria/UHM_lachnoNovelclade/species.csv \
        studies/bacteria/UHM_lachnoNovelclade/config.csv \
        studies/bacteria/UHM_lachnoNovelclade/DATA_MANIFEST.yaml
git status --short studies/bacteria/UHM_lachnoNovelclade/  # data_dir/ is gitignored/regenerable -- confirm it's NOT staged
git commit -m "Add UHM_lachnoNovelclade study (9 species: 4 outgroup, 5 ingroup)

First study built via the new-study onboarding skill (Task 4) and the
unified build_study_config.py dispatch (Task 1): protein always local_faa
(already-pulled proteomes), genome mixed local_genome/ncbi per species.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

Do **not** run `bin/run_study.sh` as part of this task — that launches the actual
`nf_NovInvenio` pipeline (a real compute job), which is a separate decision the
user makes once they've reviewed `config.csv`.

---

## Self-Review Notes (for whoever executes this plan)

- Task 2's Step 3 exists because `agaricomycetes_pairwise`'s `Scom` row has an
  empty `UniProt_Proteome_ID` in the currently-committed `species.csv`, yet its
  `config.csv` shows a `Scom.pep.fa` that was clearly produced by *some* process
  — investigate before assuming it's safely `Protein_Source=uniprot` with a blank
  accession (it is very unlikely to be — `fetch_uniprot_proteome.py --proteome-id ""`
  would fail outright, so this row was evidently built by a different, unknown
  path). This is a pre-existing data-provenance gap this plan surfaces but does
  not resolve blindly.
- Tasks 2 and 3 deliberately do **not** touch `agaricomycetes_pairwise`,
  `agaricomycetes_mmseqs`, and `agaricomycetes_novelty_discovery`'s shared
  `config.csv`/`data_dir` (per `bin/run_study.sh`'s own note that these three
  share one `config.csv`/`data_dir`, symlinked) beyond the `species.csv` schema
  migration — do not regenerate their `config.csv` in this plan; only Task 3's
  `UHM_Koxytoca` gets a real regeneration.
