#!/usr/bin/env python3
"""Reference-sequence-based HAC/hacA screen across all 293 strains --
notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md,
component 7, run directly against per-strain proteomes via diamond/
hmmsearch rather than through this study's own tier-1 clustering (which
has not been run on the full 293-strain set yet; see
PANGENOME_CLUSTER_PROFILE_NOTES.md's open items).

Three signals, deliberately kept separate:

1. hacA (Afu3g04070) presence/identity/coverage per strain -- single-copy
   reference locus, diamond blastp best hit per strain.
2. hrmA (Afu5g14900) presence/identity/coverage per strain -- same method.
3. PF11001 family COPY NUMBER per strain (hmmsearch --tblout hit count) --
   context only, NOT a presence/absence signal for "the HAC locus". Real
   data (2026-09-13) showed PF11001 hits ~3-9 strongly-scoring (>170
   bitscore) paralogs in literally ALL 293 strains: it is a broader,
   ancient, non-mobile paralogous domain family ("YDR124W-like, helical
   bundle domain"), of which the specific hrmA copy is the one Starship-
   mobile, subtelomeric member the paper calls "HAC". Screening by raw
   PF11001 presence would call every strain "positive" and hide the real
   variability -- use signal 2 (the hrmA locus itself) for that question,
   and report PF11001 copy number alongside it as context (e.g. to spot a
   strain with an unusual expansion/contraction of the whole family).

Usage:
  hac_reference_screen.py \\
      --diamond_tsv hacA_hrmA_vs_study.tsv \\
      --hmmsearch_tblout PF11001_vs_study.tblout \\
      --config config.csv \\
      --haca_query Afu3g04070-T-p1 --haca_qlen 433 \\
      --hrma_query Afu5g14900-T-p1 --hrma_qlen 362 \\
      --output hac_reference_screen.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from id_crosswalk import (  # noqa: E402
    parse_diamond_blastp_hits_per_strain,
    parse_hmmsearch_tblout_hits,
)

DEFAULT_MIN_PIDENT = 50.0
DEFAULT_MIN_QCOV = 0.5


def is_present(hit: dict, qlen: int, min_pident: float, min_qcov: float) -> bool:
    qcov = hit["align_len"] / qlen
    return hit["pident"] >= min_pident and qcov >= min_qcov


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diamond_tsv", required=True)
    ap.add_argument("--hmmsearch_tblout", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--haca_query", default="Afu3g04070-T-p1")
    ap.add_argument("--haca_qlen", type=int, default=433)
    ap.add_argument("--hrma_query", default="Afu5g14900-T-p1")
    ap.add_argument("--hrma_qlen", type=int, default=362)
    ap.add_argument("--min_pident", type=float, default=DEFAULT_MIN_PIDENT)
    ap.add_argument("--min_qcov", type=float, default=DEFAULT_MIN_QCOV)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
    from novinvenio_path import add_novinvenio_lib_to_path  # noqa: E402

    add_novinvenio_lib_to_path()
    from config_parser import parse_config  # noqa: E402

    samples = [s for s in parse_config(args.config) if s.group == "IN"]
    all_shorts = [s.short for s in samples]

    with open(args.diamond_tsv) as fh:
        diamond_lines = fh.readlines()
    diamond_hits_raw = {}
    for line in diamond_lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12:
            continue
        diamond_hits_raw.setdefault(parts[0], []).append(line)

    haca_hits = parse_diamond_blastp_hits_per_strain(diamond_hits_raw.get(args.haca_query, []))
    hrma_hits = parse_diamond_blastp_hits_per_strain(diamond_hits_raw.get(args.hrma_query, []))

    with open(args.hmmsearch_tblout) as fh:
        pf11001_hits = parse_hmmsearch_tblout_hits(fh.readlines())

    with open(args.output, "w") as out:
        out.write(
            "Short\thacA_present\thacA_pident\thacA_qcov\t"
            "hrmA_present\thrmA_pident\thrmA_qcov\tPF11001_copy_number\n"
        )
        for short in sorted(all_shorts):
            haca = haca_hits.get((args.haca_query, short))
            hrma = hrma_hits.get((args.hrma_query, short))
            pf_n = len(pf11001_hits.get(short, []))

            haca_present, haca_pident, haca_qcov = "N", "", ""
            if haca:
                qcov = haca["align_len"] / args.haca_qlen
                haca_present = "Y" if is_present(haca, args.haca_qlen, args.min_pident, args.min_qcov) else "N"
                haca_pident, haca_qcov = f"{haca['pident']:.1f}", f"{qcov:.2f}"

            hrma_present, hrma_pident, hrma_qcov = "N", "", ""
            if hrma:
                qcov = hrma["align_len"] / args.hrma_qlen
                hrma_present = "Y" if is_present(hrma, args.hrma_qlen, args.min_pident, args.min_qcov) else "N"
                hrma_pident, hrma_qcov = f"{hrma['pident']:.1f}", f"{qcov:.2f}"

            out.write(
                f"{short}\t{haca_present}\t{haca_pident}\t{haca_qcov}\t"
                f"{hrma_present}\t{hrma_pident}\t{hrma_qcov}\t{pf_n}\n"
            )

    n_haca = sum(1 for s in all_shorts if haca_hits.get((args.haca_query, s)))
    n_hrma = sum(1 for s in all_shorts if hrma_hits.get((args.hrma_query, s)))
    print(
        f"{len(all_shorts)} strains screened: hacA hit in {n_haca}, "
        f"hrmA hit in {n_hrma} (pre-threshold hit counts; see --output for "
        "identity/coverage-filtered present/absent calls)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
