#!/usr/bin/env python3
"""Batch-generate TBLASTN alignment text files for every novelty (or loss) candidate
that has at least one outgroup TBLASTN hit -- so the whole set can be reviewed
offline instead of one protein at a time via show_tblastn_alignment.py.

For UHM_Koxytoca, nearly every novelty candidate (696/697) has *some* TBLASTN hit
against the 16 outgroup genomes (see results/UHM_Koxytoca/tblastn_summary.tsv) --
worth eyeballing each one's actual alignment before trusting the "novelty" call,
since a clean full-length hit usually means the protein-level search just missed an
already-present gene, not a real lineage-specific gain.

Reuses show_tblastn_alignment.py's helpers (same script, same directory) -- no
separate BLAST DB build, one `tblastn` call per candidate against its single best
(lowest e-value) outgroup hit.

Output: one text file per candidate under --outdir (default:
results/UHM_Koxytoca/tblastn_alignments/<protein_id>.txt), each containing the full
hit list (every outgroup genome this protein matched, sorted by e-value) followed by
the actual pairwise alignment against the best hit -- plus an index.tsv summarizing
every candidate's best hit (genome, e-value, %identity) for quick sorting/scanning
without opening 696 files.

This is a derived, regenerable archive (like everything else under results/ --
gitignored, DESIGN.md Sec 4) -- re-run any time after a pipeline re-run changes the
candidate set or TBLASTN results.

Usage:
  studies/bacteria/UHM_Koxytoca/bin/archive_tblastn_alignments.py
  studies/bacteria/UHM_Koxytoca/bin/archive_tblastn_alignments.py --loss
  studies/bacteria/UHM_Koxytoca/bin/archive_tblastn_alignments.py --outdir /path/to/archive --evalue 1e-5
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from show_tblastn_alignment import (  # noqa: E402
    RESULTS_DIR, ensure_tblastn_on_path, extract_query_fasta, find_hits, load_member_to_rep,
)


def load_candidate_ids(candidates_txt: Path) -> list[str]:
    """candidates.txt lines are '<Short>::<protein_id>' -- protein_id alone is
    what candidates.fa/tblastn tsvs key on (verified collision-free for this study:
    697 candidates, 697 unique bare ids)."""
    ids = []
    with open(candidates_txt) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            ids.append(line.split("::", 1)[-1])
    return ids


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loss", action="store_true", help="Archive loss candidates instead of novelty candidates")
    ap.add_argument("--outdir", help="Default: results/UHM_Koxytoca/{tblastn,loss_tblastn}_alignments/")
    ap.add_argument("--evalue", default="10", help="BLAST -evalue ceiling for the re-run (default generous: 10)")
    ap.add_argument("--outfmt", default="0", help="BLAST -outfmt value (default 0 = pairwise text)")
    args = ap.parse_args()

    candidates_txt = RESULTS_DIR / ("loss_candidates.txt" if args.loss else "candidates.txt")
    candidates_fa = RESULTS_DIR / ("loss_candidates.fa" if args.loss else "candidates.fa")
    tblastn_dir = RESULTS_DIR / ("loss_tblastn" if args.loss else "tblastn")
    cluster_tsv = RESULTS_DIR / "clusters" / ("loss_clusters_cluster.tsv" if args.loss else "clusters_cluster.tsv")
    db_dir = tblastn_dir / "db"
    outdir = Path(args.outdir) if args.outdir else RESULTS_DIR / ("loss_tblastn_alignments" if args.loss else "tblastn_alignments")

    if not candidates_txt.exists():
        sys.exit(f"ERROR: {candidates_txt} not found")
    if not tblastn_dir.exists():
        sys.exit(f"ERROR: {tblastn_dir} not found -- has the {'loss' if args.loss else 'novelty'} TBLASTN step run yet?")

    ensure_tblastn_on_path()
    outdir.mkdir(parents=True, exist_ok=True)
    member_to_rep = load_member_to_rep(cluster_tsv)

    protein_ids = load_candidate_ids(candidates_txt)
    n_hit, n_nohit, n_err = 0, 0, 0
    index_rows = []

    with tempfile.TemporaryDirectory() as td:
        query_fa = Path(td) / "query.fa"
        for i, protein_id in enumerate(protein_ids, 1):
            hits = find_hits(protein_id, tblastn_dir, member_to_rep)
            if not hits:
                n_nohit += 1
                continue

            best = hits[0]
            out_txt = outdir / f"{protein_id}.txt"
            try:
                extract_query_fasta(protein_id, candidates_fa, query_fa)
                with open(out_txt, "w") as out_fh:
                    out_fh.write(f"# {len(hits)} outgroup hit(s) for {protein_id}:\n")
                    for h in hits:
                        out_fh.write(
                            f"#   {h['genome']:20s} {h['sseqid']:25s} evalue={h['evalue']:.2e}  "
                            f"bitscore={h['bitscore']:.0f}  pident={h['pident']:.1f}%  len={h['length']}  "
                            f"q:{h['qstart']}-{h['qend']}  s:{h['sstart']}-{h['send']}\n"
                        )
                    out_fh.write(f"\n# Aligning against {best['genome']} (evalue={best['evalue']:.2e}) ...\n\n")
                    out_fh.flush()
                    subprocess.run(
                        ["tblastn", "-query", str(query_fa), "-db", str(db_dir / f"{best['genome']}.genome_db"),
                         "-outfmt", args.outfmt, "-evalue", args.evalue],
                        check=True, stdout=out_fh,
                    )
                n_hit += 1
                index_rows.append((protein_id, best["genome"], best["sseqid"], best["evalue"], best["pident"], best["length"]))
            except subprocess.CalledProcessError as e:
                print(f"[{protein_id}] ERROR: tblastn failed ({e}) -- skipping", file=sys.stderr)
                n_err += 1

            if i % 100 == 0:
                print(f"...{i}/{len(protein_ids)} candidates processed", file=sys.stderr)

    index_rows.sort(key=lambda r: r[3])  # best evalue first, worst (least novel) last
    index_tsv = outdir / "index.tsv"
    with open(index_tsv, "w") as fh:
        fh.write("protein_id\tbest_genome\tbest_sseqid\tbest_evalue\tbest_pident\tbest_length\n")
        for r in index_rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.2e}\t{r[4]:.1f}\t{r[5]}\n")

    print(f"\nWrote {n_hit} alignment file(s) + {index_tsv} to {outdir}", file=sys.stderr)
    print(f"({n_nohit} candidate(s) had no TBLASTN hit at all -- nothing to align; {n_err} tblastn error(s))", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
