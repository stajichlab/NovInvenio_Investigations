"""File name of a species' annotation table: <study>/annotations/<name>.tsv.gz.

Why: the annotation file was named <Protein_Accession>.tsv.gz. That works for UniProt
proteome ids (UP000193719) and NCBI assembly accessions (GCA_002104975.1). For a
local_faa species, Protein_Accession is the path of the protein FASTA
(/bigdata/.../GCF_050916275.1.faa). `annotations_dir / f"{acc}.tsv.gz"` with an
absolute acc gives a path BESIDE THE SOURCE FILE (pathlib drops the left operand), in
a directory of another person or project. Every reader and writer of these files uses
the functions below, so the name is the same everywhere and is always a plain file
name inside the annotations directory.

Rules (annotation_names):
  1. Protein_Accession without "/": used unchanged (UniProt and NCBI names do not
     change), if it is a safe file name (SAFE_NAME_RE).
  2. Protein_Accession with "/" (a path): its basename, without a compression suffix
     (.gz, .zst, .bz2, .xz) and then without one proteome suffix (.pep.fa, .faa,
     .fasta, .fa, .pep, .fas).
     /x/GCA_032206225.1_FOAM13_A1_protein.faa.gz -> GCA_032206225.1_FOAM13_A1_protein
     /x/UHM1088.41098__UHM1088.41098_R.bin.56.faa -> UHM1088.41098__UHM1088.41098_R.bin.56
  3. Fall back to Short when the name from rule 1 or 2 is empty or not safe, or when
     two DIFFERENT Protein_Accession values give the same name (both rows fall back).
     Rows with the same Protein_Accession share one file (same proteome). A
     Protein_Source == uniprot row never falls back: bin/build_study_config.py writes
     its file as <Protein_Accession>.tsv.gz.
  4. The result must be unique per distinct Protein_Accession; otherwise ValueError.

annotation_path() and ensure_within() are the hard guards: they raise ValueError for a
name that is not a plain safe file name or a path outside the allowed directories.
"""
from __future__ import annotations

import re
from pathlib import Path

SAFE_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]*$")
COMPRESSION_SUFFIXES = (".gz", ".zst", ".bz2", ".xz")
PROTEOME_SUFFIXES = (".pep.fa", ".faa", ".fasta", ".fa", ".pep", ".fas")
ANNOTATION_SUFFIX = ".tsv.gz"


def _strip_suffix(name: str, suffixes: tuple[str, ...]) -> str:
    low = name.lower()
    for s in suffixes:  # longest first within each tuple (".pep.fa" before ".fa")
        if low.endswith(s) and len(name) > len(s):
            return name[: -len(s)]
    return name


def is_safe_name(name: str) -> bool:
    """A plain file name: no '/', not '.' or '..', no spaces or shell characters."""
    return bool(name) and bool(SAFE_NAME_RE.match(name)) and name not in (".", "..")


def raw_annotation_name(protein_accession: str) -> str:
    """Rules 1 and 2 only (no fallback). Return '' when nothing usable is left."""
    acc = (protein_accession or "").strip()
    if not acc:
        return ""
    if "/" not in acc:
        return acc
    base = acc.rstrip("/").rsplit("/", 1)[-1]
    base = _strip_suffix(base, COMPRESSION_SUFFIXES)
    base = _strip_suffix(base, PROTEOME_SUFFIXES)
    return base


def annotation_names(rows: list[dict]) -> dict[str, str]:
    """species.csv rows -> {Short: annotation file stem}. Use ALL rows of species.csv,
    not a subset: the collision rule depends on the other rows."""
    raw: dict[str, str] = {}
    acc_of: dict[str, str] = {}
    uniprot: set[str] = set()
    for r in rows:
        short = (r.get("Short") or "").strip()
        if not short:
            raise ValueError(f"species.csv row without Short: {r}")
        acc = (r.get("Protein_Accession") or "").strip()
        acc_of[short] = acc
        name = raw_annotation_name(acc)
        raw[short] = name if is_safe_name(name) else ""
        if (r.get("Protein_Source") or "").strip().lower() == "uniprot":
            uniprot.add(short)
    accs_by_name: dict[str, set[str]] = {}
    for short, name in raw.items():
        if name:
            accs_by_name.setdefault(name, set()).add(acc_of[short])
    out: dict[str, str] = {}
    for short, name in raw.items():
        collides = bool(name) and len(accs_by_name[name]) > 1
        if name and (not collides or short in uniprot):
            out[short] = name
            continue
        if short in uniprot:
            raise ValueError(f"{short}: UniProt Protein_Accession {acc_of[short]!r} is not a safe file name")
        if not is_safe_name(short):
            raise ValueError(f"{short!r}: fallback name (Short) is not a safe file name")
        out[short] = short
    # Rule 4: one name per distinct accession.
    owner: dict[str, str] = {}
    for short, name in out.items():
        acc = acc_of[short] or f"<Short {short}>"
        if name in owner and owner[name] != acc:
            raise ValueError(f"annotation name {name!r} is used by two different proteomes "
                             f"({owner[name]!r} and {acc!r})")
        owner[name] = acc
    return out


def ensure_within(path: Path, roots: list[Path]) -> Path:
    """Return path.resolve() if it is inside one of roots (resolved); else ValueError."""
    p = Path(path).resolve()
    for root in roots:
        r = Path(root).resolve()
        if p == r or r in p.parents:
            return p
    raise ValueError(f"refusing to write {p}: outside the allowed directories "
                     f"{[str(Path(r).resolve()) for r in roots]}")


def annotation_path(annotations_dir: Path, name: str, suffix: str = ANNOTATION_SUFFIX) -> Path:
    """<annotations_dir>/<name><suffix>, with the guard: name must be a safe plain file name
    and the result must be a direct child of annotations_dir."""
    if not is_safe_name(name):
        raise ValueError(f"unsafe annotation file name {name!r}")
    d = Path(annotations_dir).resolve()
    p = (d / f"{name}{suffix}").resolve()
    if p.parent != d:
        raise ValueError(f"refusing annotation path {p}: not directly in {d}")
    return p
