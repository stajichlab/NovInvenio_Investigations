#!/usr/bin/env python3
"""Re-run TBLASTN for one candidate protein against one outgroup genome, printing a
human-readable pairwise alignment -- for manually reviewing a novelty candidate that
already has a TBLASTN hit (i.e. a candidate the protein-level search missed but the
genome-level search found -- worth eyeballing before trusting it as a real novelty).

The pipeline's own TBLASTN (nf_NovInvenio modules/tblastn.nf) only ever writes
-outfmt 6 (tabular: qseqid,sseqid,evalue,bitscore,pident,length,qstart,qend,sstart,
send) -- no alignment text, by design (that's what feeds tblastn_summary.tsv's
protein x genome presence/absence calls, alignment text would be wasted for that).
This script does NOT re-run makeblastdb -- it reuses the per-genome nucleotide BLAST
DBs the pipeline already built and storeDir-cached under
results/UHM_Koxytoca/tblastn/db/, so it's fast (single tblastn call, no DB rebuild).

Usage:
  # See every outgroup genome this protein hit, sorted by e-value
  studies/bacteria/UHM_Koxytoca/bin/show_tblastn_alignment.py --protein-id k141_119261_30 --list

  # Show the alignment against its best (lowest e-value) hit
  studies/bacteria/UHM_Koxytoca/bin/show_tblastn_alignment.py --protein-id k141_119261_30

  # Show the alignment against one specific outgroup genome
  studies/bacteria/UHM_Koxytoca/bin/show_tblastn_alignment.py --protein-id k141_119261_30 --genome KoxO_017310465

  # Same, but for a loss candidate (queries loss_candidates.fa + tblastn/loss_*
  # once that side of the pipeline has been (re-)run with ingroup DNA -- see
  # build_koxytoca_config.py's ingroup DNA fix)
  studies/bacteria/UHM_Koxytoca/bin/show_tblastn_alignment.py --protein-id <id> --loss
"""
from __future__ import annotations

import argparse
import glob
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

NII_ROOT = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
RESULTS_DIR = NII_ROOT / "results" / "UHM_Koxytoca"
BLAST_MODULE = "ncbi-blast/2.17.0+"


def ensure_tblastn_on_path() -> str:
    exe = shutil.which("tblastn")
    if exe:
        return exe
    sys.exit(
        f"ERROR: \`tblastn\` not found on PATH.\n"
        f"Run: source /etc/profile.d/modules.sh && module load {BLAST_MODULE}"
    )


def load_member_to_rep(cluster_tsv: Path) -> dict[str, str]:
    """{member_id: rep_id} from an mmseqs easy-cluster TSV (rep\\tmember per line,
    including a rep\\trep self-line) -- see summarize_tblastn.py's own parser."""
    member_to_rep: dict[str, str] = {}
    if not cluster_tsv.exists():
        return member_to_rep
    with open(cluster_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                member_to_rep[parts[1]] = parts[0]
    return member_to_rep


def find_hits(protein_id: str, tblastn_dir: Path, member_to_rep: dict[str, str] | None = None) -> list[dict]:
    """Every outgroup-genome hit for this protein, across all *.tblastn.tsv files,
    sorted by e-value ascending (best first).

    TBLASTN only ever runs on cluster REPRESENTATIVE proteins (modules/tblastn.nf's
    reps_fa input), not every candidate -- most candidates are non-representative
    cluster members and never appear as a qseqid in these files directly (that's
    how tblastn_summary.tsv's per-protein presence/absence expansion works: a hit
    on the rep gets propagated to every member of its cluster). So if protein_id
    itself isn't found, fall back to its cluster representative's hits, when
    member_to_rep is given -- otherwise a non-rep candidate would wrongly look like
    it has zero TBLASTN hits even though tblastn_summary.tsv says it does.
    """
    def _search(qid: str) -> list[dict]:
        found = []
        for tsv in sorted(tblastn_dir.glob("*.tblastn.tsv")):
            genome_short = tsv.stem.replace(".tblastn", "")
            with open(tsv) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 10 or parts[0] != qid:
                        continue
                    found.append({
                        "genome": genome_short,
                        "sseqid": parts[1], "evalue": float(parts[2]), "bitscore": float(parts[3]),
                        "pident": float(parts[4]), "length": int(parts[5]),
                        "qstart": parts[6], "qend": parts[7], "sstart": parts[8], "send": parts[9],
                    })
        return found

    hits = _search(protein_id)
    if not hits and member_to_rep:
        rep_id = member_to_rep.get(protein_id)
        if rep_id and rep_id != protein_id:
            hits = _search(rep_id)
    hits.sort(key=lambda h: h["evalue"])
    return hits


def extract_query_fasta(protein_id: str, candidates_fa: Path, dest: Path) -> None:
    with open(candidates_fa) as fh, open(dest, "w") as out:
        writing = False
        for line in fh:
            if line.startswith(">"):
                header_id = line[1:].split()[0]
                writing = header_id == protein_id
            if writing:
                out.write(line)
    if dest.stat().st_size == 0:
        sys.exit(f"ERROR: protein_id {protein_id!r} not found in {candidates_fa}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--protein-id", required=True, help="Query protein id, e.g. k141_119261_30")
    ap.add_argument("--genome", help="Outgroup Short (e.g. KoxO_017310465); default: best (lowest e-value) hit")
    ap.add_argument("--loss", action="store_true", help="Look in loss_candidates.fa / loss-side TBLASTN output instead of the novelty side")
    ap.add_argument("--list", action="store_true", help="Just list every genome this protein hit (no alignment run)")
    ap.add_argument("--outfmt", default="0", help="BLAST -outfmt value (default 0 = pairwise text; try 4 for a flat query/subject/match block)")
    ap.add_argument("--evalue", default="10", help="BLAST -evalue ceiling for this ad-hoc re-run (default generous: 10, so a real hit always shows even if just below the pipeline's own cutoff)")
    args = ap.parse_args()

    candidates_fa = RESULTS_DIR / ("loss_candidates.fa" if args.loss else "candidates.fa")
    tblastn_dir = RESULTS_DIR / ("loss_tblastn" if args.loss else "tblastn")
    cluster_tsv = RESULTS_DIR / "clusters" / ("loss_clusters_cluster.tsv" if args.loss else "clusters_cluster.tsv")
    db_dir = tblastn_dir / "db"

    if not tblastn_dir.exists():
        sys.exit(f"ERROR: {tblastn_dir} not found -- has the {'loss' if args.loss else 'novelty'} TBLASTN step run yet?")

    member_to_rep = load_member_to_rep(cluster_tsv)
    hits = find_hits(args.protein_id, tblastn_dir, member_to_rep)
    if not hits:
        sys.exit(f"No TBLASTN hits found for {args.protein_id} in {tblastn_dir}")

    print(f"# {len(hits)} outgroup hit(s) for {args.protein_id}:", file=sys.stderr)
    for h in hits:
        print(f"#   {h['genome']:20s} {h['sseqid']:25s} evalue={h['evalue']:.2e}  bitscore={h['bitscore']:.0f}  "
              f"pident={h['pident']:.1f}%  len={h['length']}  q:{h['qstart']}-{h['qend']}  s:{h['sstart']}-{h['send']}",
              file=sys.stderr)
    if args.list:
        return 0

    chosen = next((h for h in hits if h["genome"] == args.genome), hits[0]) if args.genome else hits[0]
    print(f"\n# Aligning against {chosen['genome']} (evalue={chosen['evalue']:.2e}) ...\n", file=sys.stderr)

    tblastn = ensure_tblastn_on_path()
    genome_db = db_dir / f"{chosen['genome']}.genome_db"
    if not any(glob.glob(f"{genome_db}.*")):
        sys.exit(f"ERROR: BLAST db {genome_db}.* not found under {db_dir}")

    with tempfile.TemporaryDirectory() as td:
        query_fa = Path(td) / "query.fa"
        extract_query_fasta(args.protein_id, candidates_fa, query_fa)
        subprocess.run(
            ["tblastn", "-query", str(query_fa), "-db", str(genome_db),
             "-outfmt", args.outfmt, "-evalue", args.evalue],
            check=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
