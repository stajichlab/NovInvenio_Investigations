#!/usr/bin/env python3
"""Add (or replace) a "Species tree" card in a study run's report.html.

The card holds a PNG of the tree, embedded as a base64 data URI so report.html stays
one self-contained file. The PNG is rendered from species_tree.svg (written by
nf_NovInvenio bin/species_tree.py) with ImageMagick `convert`.

The caption states the method, the number of markers kept and how the tree was rooted
(read from species_tree.rooting.json). Running the script again replaces the card.
A pipeline rerun that regenerates report.html removes the card; run this script again.

Usage:
  bin/inject_species_tree.py --report docs/fungi/<study>/<run>/report.html \\
      --tree-dir studies/fungi/<study>/analysis/species_tree
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

CARD_RE = re.compile(r'\n?[ \t]*<section class="card" id="species-tree">.*?</section>\n?', re.S)
ANCHOR = '<section class="card">\n    <h2 class="card-title">Run summary</h2>'


def render_png(svg: Path, density: int = 150) -> bytes:
    with tempfile.TemporaryDirectory() as tmp:
        png = Path(tmp) / "tree.png"
        subprocess.run(
            ["convert", "-density", str(density), "-background", "white", str(svg), str(png)],
            check=True,
        )
        return png.read_bytes()


def rooting_text(info: dict) -> str:
    method = info.get("method")
    n_out, outside = info.get("n_out", 0), info.get("n_out_outside", 0)
    if method == "outgroup":
        return "Rooted on the OUT group (monophyletic)."
    if method == "outgroup_best_split":
        return (f"Rooted on the best OUT/non-OUT split; {outside} of {n_out} OUT tips "
                "are outside the OUT side.")
    return "No OUT tips; midpoint rooted."


def build_card(png: bytes, n_markers: int, info: dict) -> str:
    b64 = base64.b64encode(png).decode("ascii")
    return (
        '<section class="card" id="species-tree">\n'
        '    <h2 class="card-title">Species tree</h2>\n'
        f'    <img alt="Species tree" style="max-width:100%;height:auto" '
        f'src="data:image/png;base64,{b64}">\n'
        f'    <p class="card-note">phyling (FastTree, consensus of {n_markers} single-marker trees; '
        f'top {n_markers} markers by treeness/RCV), fungi_odb12. {rooting_text(info)}</p>\n'
        '  </section>\n'
    )


def inject(html: str, card: str) -> str:
    html = CARD_RE.sub("\n", html)
    if ANCHOR not in html:
        raise SystemExit("anchor 'Run summary' card not found in report.html; not modified")
    return html.replace(ANCHOR, card + "  " + ANCHOR, 1)


def update_run_json(meta_path: Path, old_bytes: bytes, new_bytes: bytes) -> str:
    """Keep run.json's report_sha256 in step with an edited report.html, as
    lib/study_runs.py does for the links and notes edits: update it only when it
    matched the report before the edit. Returns a status suffix."""
    if not meta_path.is_file():
        return ""
    meta = json.loads(meta_path.read_text())
    if meta.get("report_sha256") == hashlib.sha256(old_bytes).hexdigest():
        meta["report_sha256"] = hashlib.sha256(new_bytes).hexdigest()
        meta_path.write_text(json.dumps(meta, indent=2) + "\n")
        return "; run.json report_sha256 updated"
    if "report_sha256" in meta:
        return "; WARNING: run.json report_sha256 did not match report.html before the edit, not changed"
    return ""


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--tree-dir", required=True, type=Path)
    args = ap.parse_args(argv)

    svg = args.tree_dir / "species_tree.svg"
    rooting = args.tree_dir / "species_tree.rooting.json"
    n_markers = len(list((args.tree_dir / "filter").glob("*.mfa")))
    for p in (svg, rooting, args.report):
        if not p.exists():
            raise SystemExit(f"missing: {p}")
    if n_markers == 0:
        raise SystemExit(f"no marker alignments (*.mfa) in {args.tree_dir / 'filter'}")

    card = build_card(render_png(svg), n_markers, json.loads(rooting.read_text()))
    old_bytes = args.report.read_bytes()
    new_bytes = inject(old_bytes.decode("utf-8"), card).encode("utf-8")
    args.report.write_bytes(new_bytes)
    note = update_run_json(args.report.parent / "run.json", old_bytes, new_bytes)
    print(f"{args.report}: species tree card written ({n_markers} markers){note}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
