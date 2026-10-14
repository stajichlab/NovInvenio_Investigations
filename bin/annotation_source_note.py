#!/usr/bin/env python3
"""Write a one-line annotation-source note into regenerated report pages.

bin/sync_reports.sh merges every annotations/*.tsv(.gz) into the uniprot_*
columns (bin/merge_uniprot_annotations.py), and nf_NovInvenio's make_report.py /
make_core_report.py / make_losses_report.py draw novelties/core/losses.html from
those columns. Files from bin/hmm_domtbl_to_annotations.py have Pfam only. This
script puts a note under each page's header that names those species:

  Annotation source: Pfam from a HMMER search (no GO/InterPro) for: <species>. ...

Which species: with --matrix (preferred), the source_proteome values of the
rows that bin/merge_uniprot_annotations.py labelled `hmmer-pfam` in the matrix
the page was drawn from. Without --matrix, every species with a Pfam-only file
in --annotations.

The page generators live in the nf_NovInvenio repo; this NII-side step edits
the generated HTML, between marker comments, after they run (same pattern as
lib/study_runs.py inject_run_report_links). Idempotent. When no species
qualifies, an earlier note is removed and nothing else changes. A missing page
is skipped.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from annotation_sources import (  # noqa: E402
    SOURCE_HMMER_PFAM, file_accession, file_source, hmmer_proteomes_in_matrix, inject_note,
    note_text, species_labels,
)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--matrix", type=Path,
                     help="Merged matrix the pages were drawn from (has annotation_source + source_proteome)")
    src.add_argument("--annotations", nargs="+", type=Path,
                     help="The annotation files given to merge_uniprot_annotations.py")
    ap.add_argument("--species", type=Path, default=None,
                    help="Study species.csv (Short / Protein_Accession -> 'Species (Short)'). Optional; "
                         "without it the note names the Short or file accession.")
    ap.add_argument("--html", nargs="+", required=True, type=Path,
                    help="Report pages to edit in place (missing files are skipped)")
    args = ap.parse_args(argv)

    labels = species_labels(args.species)
    if args.matrix is not None:
        keys = hmmer_proteomes_in_matrix(args.matrix)
    else:
        keys = [file_accession(p) for p in sorted(args.annotations) if file_source(p) == SOURCE_HMMER_PFAM]
    hmmer = [labels.get(k, k) for k in keys]
    note = note_text(hmmer)
    for page in args.html:
        if not page.is_file():
            print(f"{page}: not found, skipped", file=sys.stderr)
            continue
        old = page.read_text(encoding="utf-8")
        new, ok = inject_note(old, note)
        if not ok:
            print(f"WARNING: {page} has no </header>; annotation-source note not added", file=sys.stderr)
            continue
        if new == old:
            print(f"{page}: annotation-source note already current", file=sys.stderr)
            continue
        page.write_text(new, encoding="utf-8")
        print(f"{page}: annotation-source note {'set' if note else 'removed'} "
              f"({len(hmmer)} HMMER Pfam-only species)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
