"""Which source made a per-species annotation file in studies/<d>/<s>/annotations/.

Two producers write files with the same column layout:
  - bin/extract_dat_annotations.py (via bin/build_study_config.py): UniProt .dat
    records. GO, Pfam and InterPro are filled.               -> SOURCE_UNIPROT
  - bin/hmm_domtbl_to_annotations.py: a Pfam-A HMMER search (hmmsearch/hmmscan).
    Pfam only; go_ids and interpro_ids are empty.            -> SOURCE_HMMER_PFAM

bin/merge_uniprot_annotations.py writes the source of each protein into the
`annotation_source` column, and bin/annotation_source_note.py writes a one-line
note into the regenerated report pages. Both use file_source() below.

Detection, in order:
  1. The file header has the column `pfam_best_ievalue` (only
     hmm_domtbl_to_annotations.py writes it)               -> hmmer-pfam
  2. A `<file>.provenance.yaml` sidecar exists and its `derived_by` names
     hmm_domtbl_to_annotations                              -> hmmer-pfam
  3. Otherwise                                              -> uniprot
"""
from __future__ import annotations

import csv
import re
from html import escape
from pathlib import Path

from compressed_io import open_text

SOURCE_UNIPROT = "uniprot"
SOURCE_HMMER_PFAM = "hmmer-pfam"
PFAM_MARKER_COLUMN = "pfam_best_ievalue"
HMMER_SCRIPT = "hmm_domtbl_to_annotations"
VOCABULARIES = {
    SOURCE_UNIPROT: ("GO", "Pfam", "InterPro"),
    SOURCE_HMMER_PFAM: ("Pfam",),
}


def sidecar_path(path: Path) -> Path:
    path = Path(path)
    return path.with_name(path.name + ".provenance.yaml")


def _sidecar_says_hmmer(path: Path) -> bool:
    side = sidecar_path(path)
    if not side.is_file():
        return False
    import yaml  # only needed when a sidecar exists

    with open(side) as fh:
        rec = yaml.safe_load(fh) or {}
    return HMMER_SCRIPT in str(rec.get("derived_by", ""))


def file_source(path: Path, header: list[str] | None = None) -> str:
    """SOURCE_UNIPROT or SOURCE_HMMER_PFAM for one annotation file. `header` is
    the file's column list when the caller has already read it."""
    if header is None:
        with open_text(path, "rt", newline="") as fh:
            header = next(csv.reader(fh, delimiter="\t"), [])
    if PFAM_MARKER_COLUMN in header:
        return SOURCE_HMMER_PFAM
    if _sidecar_says_hmmer(path):
        return SOURCE_HMMER_PFAM
    return SOURCE_UNIPROT


def file_accession(path: Path) -> str:
    """'annotations/GCA_002104975.1.tsv.gz' -> 'GCA_002104975.1' (the
    Protein_Accession the file is named after)."""
    name = Path(path).name
    for suffix in (".tsv.gz", ".tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def species_labels(species_csv: Path | None) -> dict[str, str]:
    """Protein_Accession and Short -> 'Species (Short)' from a study's species.csv.
    Missing file or missing columns -> {} (callers fall back to the key)."""
    if species_csv is None or not Path(species_csv).is_file():
        return {}
    out: dict[str, str] = {}
    with open(species_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            acc = (row.get("Protein_Accession") or "").strip()
            sp = (row.get("Species") or "").strip()
            short = (row.get("Short") or "").strip()
            label = f"{sp} ({short})" if sp and short else (sp or short or acc)
            for key in (acc, short):
                if key:
                    out[key] = label
    return out


def hmmer_proteomes_in_matrix(matrix: Path) -> list[str]:
    """Sorted source_proteome values of the rows whose annotation_source is
    hmmer-pfam in a merged matrix (bin/merge_uniprot_annotations.py output).
    [] when either column is missing."""
    found: set[str] = set()
    with open_text(matrix, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        if "annotation_source" not in cols or "source_proteome" not in cols:
            return []
        for row in reader:
            if row["annotation_source"] == SOURCE_HMMER_PFAM:
                found.add(row["source_proteome"])
    return sorted(found)


_NOTE_START = "<!-- nii-annotation-source:start -->"
_NOTE_END = "<!-- nii-annotation-source:end -->"
_NOTE_BLOCK_RE = re.compile(r"[ \t]*" + re.escape(_NOTE_START) + r".*?" + re.escape(_NOTE_END) + r"\n?", re.S)
_HEADER_END_RE = re.compile(r"</header>\n?")


def note_text(hmmer_labels: list[str]) -> str:
    """The one-line note, or '' when no species has HMMER-only annotation."""
    if not hmmer_labels:
        return ""
    return (
        "Annotation source: Pfam from a HMMER search (no GO/InterPro) for: "
        + ", ".join(hmmer_labels)
        + ". Other annotated species: UniProt (GO, Pfam, InterPro). "
        "For the species above, an empty GO or InterPro field means 'not searched', not 'no annotation'."
    )


def inject_note(html: str, note: str) -> tuple[str, bool]:
    """Put `note` right after the page's first </header>, between marker
    comments. Idempotent: an earlier block is replaced; an empty note removes it.
    Returns (new_html, anchor_found). anchor_found is False only when a note
    was given and the page has no </header>."""
    text = _NOTE_BLOCK_RE.sub("", html)
    if not note:
        return text, True
    m = _HEADER_END_RE.search(text)
    if m is None:
        return text, False
    block = (
        f"  {_NOTE_START}\n"
        f'  <p class="sub nii-annotation-source" style="margin:6px 0 10px">{escape(note)}</p>\n'
        f"  {_NOTE_END}\n"
    )
    return text[: m.end()] + block + text[m.end():], True
