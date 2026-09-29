#!/usr/bin/env python3
"""Build one genome-subset input set from a finished pangenome run.

Sampling unit: representative ingroup strains (strain_inventory.tsv
is_representative == 1). The pipeline's co-occurrence test uses only
representatives, so N here is the number of genomes the test sees.
- Af293 is always kept, so the Af293 Starship coordinate check works in
  every subset. The other N-1 are drawn at random (seeded).
- All outgroup (GROUP=OUT) strains are kept, matching the full run.
- Non-representative strains are dropped.

Writes into --out_dir: presence_matrix.rescued.tsv, samplesheet.with_clades.csv,
strain_inventory.tsv, family_positions.tsv.zst, subset_strains.txt.
Rows of the matrix are kept even when a family is absent from every kept
strain; the frequency-bin step labels those, and co-occurrence needs
>= --min_strain_count carriers, so they do not enter the test.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

import numpy as np

FORCED_STRAIN = "Asfu_Af293"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pangenome_dir", required=True, type=Path)
    ap.add_argument("--n", required=True, type=int, help="representative ingroup strains to keep")
    ap.add_argument("--seed", required=True, type=int)
    ap.add_argument("--out_dir", required=True, type=Path)
    args = ap.parse_args()
    pdir, out = args.pangenome_dir, args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    with (pdir / "samplesheet.with_clades.csv").open() as fh:
        sheet = list(csv.DictReader(fh))
    group_of = {r["Short"]: r["GROUP"] for r in sheet}
    inventory = read_tsv(pdir / "strain_inventory.tsv")

    reps = sorted(r["Short"] for r in inventory
                  if r["is_representative"] == "1" and group_of.get(r["Short"]) == "IN")
    outgroups = sorted(s for s, g in group_of.items() if g == "OUT")
    if FORCED_STRAIN not in reps:
        raise SystemExit(f"{FORCED_STRAIN} is not a representative ingroup strain")
    if not 1 <= args.n <= len(reps):
        raise SystemExit(f"--n must be in 1..{len(reps)}")

    rng = np.random.default_rng(args.seed)
    pool = [s for s in reps if s != FORCED_STRAIN]
    drawn = list(rng.choice(pool, size=args.n - 1, replace=False)) if args.n > 1 else []
    keep_in = sorted([FORCED_STRAIN, *drawn])
    keep = set(keep_in) | set(outgroups)
    (out / "subset_strains.txt").write_text("\n".join(keep_in + outgroups) + "\n")

    with (out / "samplesheet.with_clades.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sheet[0]))
        w.writeheader()
        w.writerows(r for r in sheet if r["Short"] in keep)

    with (out / "strain_inventory.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(inventory[0]), delimiter="\t")
        w.writeheader()
        w.writerows(r for r in inventory if r["Short"] in keep)

    with (pdir / "presence_matrix.rescued.tsv").open() as src, \
            (out / "presence_matrix.rescued.tsv").open("w") as dst:
        header = src.readline().rstrip("\n").split("\t")
        missing = keep - set(header[1:])
        if missing:
            raise SystemExit(f"strains missing from presence matrix: {sorted(missing)}")
        idx = [0] + [i for i, s in enumerate(header) if s in keep]
        dst.write("\t".join(header[i] for i in idx) + "\n")
        for line in src:
            parts = line.rstrip("\n").split("\t")
            dst.write("\t".join(parts[i] for i in idx) + "\n")

    # family_positions: Short, family, contig, rank -- keep rows of kept strains.
    fp_in = pdir / "family_positions.tsv.zst"
    fp_out = out / "family_positions.tsv.zst"
    awk = ("NR==FNR{k[$1]=1;next} FNR==1||($1 in k)")
    cmd = (f"zstd -dc {fp_in} | awk -F'\\t' '{awk}' {out / 'subset_strains.txt'} - "
           f"| zstd -q -T0 -o {fp_out} -f")
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)
    print(f"subset n={args.n} seed={args.seed}: {len(keep_in)} IN + {len(outgroups)} OUT -> {out}")


if __name__ == "__main__":
    main()
