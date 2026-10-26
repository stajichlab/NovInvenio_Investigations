#!/usr/bin/env python3
"""Write the cross-contrast comparison page for yeast_convergence.

Reads analysis/convergence/ (convergence_compare.py output + PROVENANCE.yaml) and writes
<docs>/fungi/yeast_convergence/comparison.html plus data/ with the small tables (class 2).
Numbers come from the tables, not from typed text.

Usage: python3.12 bin/yc_comparison_page.py --repo-root <worktree of main with docs/>
Absolute NII path for the study folder (study-specific script, NII CLAUDE.md).
"""
import argparse, csv, gzip, shutil
from html import escape
from pathlib import Path

SRC = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/yeast_convergence")
CONV = SRC / "analysis/convergence"
DATA_FILES = ["summary.md", "PROVENANCE.yaml", "term_recurrence.tsv.gz", "pairwise_jaccard.tsv",
              "reversal_overlap.tsv", "focus_terms.tsv", "annotation_coverage.tsv"]


def tsv(p):
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def lichen_screen_item() -> str:
    """Limits-box line from the lichen proteome and candidate screens (analysis/lichen_check,
    yc_pseudo_lichen/analysis/candidate_taxcheck), read from the summary tables. Empty if absent."""
    p1 = CONV.parent / "lichen_check/summary.tsv"
    p2 = SRC / "yc_pseudo_lichen/analysis/candidate_taxcheck/summary.tsv"
    if not (p1.is_file() and p2.is_file()):
        return ""
    prot = {(r["proteome"], r["group_of_best_hit"]): r for r in tsv(p1)}
    lup = prot.get(("Lethlupi", "Fungi"))
    cand_tot = cand_nonfungal = cand_nohit = 0
    for r in tsv(p2):
        n = int(r["n_proteins"])
        if r["group_of_best_hit"] == "no_hit":
            cand_nohit += n
        elif r["group_of_best_hit"] != "Fungi":
            cand_nonfungal += n
        if r["group_of_best_hit"] != "no_hit":
            cand_tot += n
    if not lup:
        return ""
    plant = int(prot.get(("Lethlupi", "Viridiplantae"), {"n_proteins": 0})["n_proteins"])
    return (f"<li>Lichen proteomes: <em>Letharia lupina</em> WasteWater1 was screened against NCBI nr "
            f"(diamond, best hit by bit score): {100 * float(lup['fraction_of_all_proteins']):.1f}% of its "
            f"{int(lup['n_proteins_total']):,} proteins have a fungal best hit, and {plant} have a green-plant or algal one. "
            f"Of {cand_tot + cand_nohit:,} novelty and loss candidate proteins of the lichen contrast, {cand_nonfungal} have a "
            f"non-fungal best hit and {cand_nohit} have no hit. The screen is a check, not an assignment. "
            f"<em>C. rubrotincta</em> was not screened as a whole proteome.</li>\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", type=Path, required=True)
    a = ap.parse_args()
    out = a.repo_root / "docs/fungi/yeast_convergence"
    (out / "data").mkdir(parents=True, exist_ok=True)
    for f in DATA_FILES:
        shutil.copyfile(CONV / f, out / "data" / f)
    idx = tsv(SRC / "study_index.tsv")
    cov = {r["id"]: r for r in tsv(CONV / "annotation_coverage.tsv")}
    rec = tsv(CONV / "term_recurrence.tsv.gz")
    hits = [r for r in rec if r["tier"] == "primary" and r["class"] in ("yeast_convergent", "streamlining_generic")]
    n_secondary = sum(1 for r in rec if r["tier"] == "secondary" and r["class"] in ("yeast_convergent", "streamlining_generic"))
    # Numbers quoted in the limits box come from the files, not from typed text.
    busco = {r["name"]: float(r["complete_pct"]) for r in tsv(SRC / "analysis/busco/summary.tsv")}
    # BUSCO rows are keyed by manifest name; match each species.csv proteome by its accession
    # (a replaced proteome keeps the old Short, so the Short is not a safe key).
    by_acc = {r["accession"]: r["name"] for r in tsv(SRC / "analysis/busco/manifest_all.tsv")}
    used = {}
    for sp in SRC.glob("yc_*/species.csv"):
        for r in csv.DictReader(open(sp)):
            acc = r["Protein_Accession"]
            used[acc] = by_acc.get(acc) or Path(acc).name.replace(".proteins.fa.gz", "")
    n_low = sum(1 for k in used.values() if k in busco and busco[k] < 85.0)
    amf = SRC.parent.parent.parent / "results/yc_pseudo_arbuscular_mycorrhizal/candidates.txt"
    n_amf = sum(1 for _ in open(amf)) if amf.is_file() else None
    roles = {}
    for r in idx:
        roles.setdefault(r["role"], []).append(r["id"])
    have_tables = [r for r in idx if (a.repo_root / f"docs/fungi/{r['id']}/pairwise-diamond/enrichment.html").is_file()]
    n_role = {k: sum(1 for r in have_tables if r["role"] == k) for k in roles}
    rows = "".join(
        f"<tr><td>{escape(r['direction'])}</td><td>{escape(r['vocab'])}</td><td>{escape(r['term'])}</td>"
        f"<td>{escape(r['name'])}</td><td>{r['n_origin_sig']}/{r['n_origin_tested']}</td>"
        f"<td>{escape(r['origins_sig'].replace(';', ', '))}</td><td>{float(r['q_recur']):.3g}</td>"
        f"<td>{r['n_pseudo_sig']}/{r['n_pseudo_tested']}</td><td>{escape(r['class'])}</td></tr>"
        for r in sorted(hits, key=lambda r: (r["direction"], float(r["q_recur"]))))
    n_conv = sum(1 for r in hits if r["class"] == "yeast_convergent")
    contrast_rows = "".join(
        f"<tr><td>{escape(r['id'])}</td><td>{escape(r['role'])}</td><td>{escape(r['origin'])}</td>"
        f"<td>{r['n_in']} / {r['n_out']}</td>"
        + (f"<td><a href=\"../{r['id']}/pairwise-diamond/enrichment.html\">novelty</a> &middot; "
           f"<a href=\"../{r['id']}/pairwise-diamond/enrichment_loss.html\">loss</a> &middot; "
           f"<a href=\"../{r['id']}/pairwise-diamond/report.html\">report</a></td>"
           if r in have_tables else "<td>not run</td>") + "</tr>"
        for r in idx)
    not_run = [r["id"] for r in idx if r not in have_tables]
    not_run_item = ("<li>Not run, so not in the comparison: " + ", ".join(f"<code>{escape(i)}</code>" for i in not_run)
                    + f". {n_role.get('pseudo', 0)} controls are included.</li>") if not_run else ""
    lichen_item = lichen_screen_item()
    files = "".join(f"<li><a href=\"data/{f}\">{f}</a></li>" for f in DATA_FILES)
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>yeast_convergence comparison</title>
<style>
:root {{ color-scheme: light dark; --bg:#f7f7f5; --fg:#1a1a1a; --muted:#666; --line:#d8d8d4; --card:#fff; --link:#0b5cad; --warn-bg:#fff6e0; --warn-line:#e0b84c; --head:#efefeb; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#16181c; --fg:#e6e6e6; --muted:#a0a4ab; --line:#33373f; --card:#22252b; --link:#7db4ff; --warn-bg:#2e2716; --warn-line:#8a6d1f; --head:#2a2e35; }} }}
body {{ margin:0; padding:1.5rem 1rem 3rem; background:var(--bg); color:var(--fg); font-family:system-ui,-apple-system,sans-serif; line-height:1.45; }}
main {{ max-width:1100px; margin:0 auto; }} a {{ color:var(--link); }} h1 {{ font-size:1.5rem; margin:.2rem 0; }} h2 {{ font-size:1.15rem; margin:2rem 0 .5rem; }}
.note {{ font-size:.88rem; color:var(--muted); }} .scroll {{ overflow-x:auto; }}
table {{ border-collapse:collapse; font-size:.82rem; background:var(--card); }} th,td {{ border:1px solid var(--line); padding:.25rem .5rem; text-align:left; vertical-align:top; }} th {{ background:var(--head); }}
.caveats {{ background:var(--warn-bg); border:1px solid var(--warn-line); border-radius:8px; padding:.6rem 1rem; margin:1rem 0; font-size:.9rem; }} .caveats ul {{ margin:.3rem 0; padding-left:1.2rem; }}
</style></head><body><main>
<p class="note"><a href="../index.html">fungi</a> / <a href="index.html">yeast convergence</a> / comparison. Design and run log: <a href="https://github.com/stajichlab/NovInvenio_Investigations/tree/main/studies/fungi/yeast_convergence">studies/fungi/yeast_convergence</a> (PR #43).</p>
<h1>yeast_convergence: cross-contrast comparison</h1>
<p>Sixteen contrasts test whether yeast growth form evolved independently with the same functional gene families.
{n_role.get('origin', 0)} origin contrasts compare a yeast lineage with its filamentous sister. {n_role.get('pseudo', 0)} pseudo-origin controls
(filamentous lineages with large gene loss and no change in growth form) and {n_role.get('reversal', 0)} reversal contrasts have enrichment tables.
A term counts as recurrent when at least 3 origins call it significant (Pfam and Pfam clan; BH q &le; 0.05 against a Poisson-binomial null).
Terms that a pseudo-origin control also calls are labelled <em>streamlining_generic</em>.</p>
<h2>Result</h2>
<p><strong>{n_conv} terms are recurrent and absent from all controls</strong> (class <em>yeast_convergent</em>). {len(hits) - n_conv} terms are recurrent and also called in at least one control.
This does not show that convergence is absent. The analysis has few origins, small ingroups for some contrasts, and a null model that ignores family size.</p>
<div class="scroll"><table><thead><tr><th>direction</th><th>vocab</th><th>term</th><th>name</th><th>origins sig/tested</th><th>origins</th><th>q</th><th>controls sig/tested</th><th>class</th></tr></thead><tbody>{rows}</tbody></table></div>
<div class="caveats"><strong>Limits</strong><ul>
{not_run_item}{lichen_item}
<li>Pfam is the primary vocabulary for every contrast: from UniProt records, or from a HMMER search for proteomes without UniProt annotation. GO and InterPro are used only for contrasts whose species are all UniProt proteomes. {n_secondary} terms are recurrent in those secondary vocabularies.</li>
<li>{n_low} proteomes score under 85% BUSCO complete (fungi_odb10) and are kept. Loss calls for those species may include missing genes, not only losses.</li>
<li>The arbuscular mycorrhizal control has {f"{n_amf:,}" if n_amf else "many"} novelty candidates. About half of its proteins have no diamond hit in the Mucoromycota outgroups. Whether this is real expansion is not tested. Several of its top Pfam bars are transposon families.</li>
<li><code>yc_rev_neolecta</code> has one ingroup species. The pipeline builds the presence matrix from proteins with a hit to another proteome, so it cannot call novelty candidates there. This contrast adds loss tables only.</li>
<li><code>yc_saccharomycotina</code> and <code>yc_taphrinomycotina</code> share two outgroup species. Neither makes a recurrent call, so the overlap does not change these results.</li>
<li>Growth-form codings for some species are predicted or unsourced (<a href="https://github.com/stajichlab/NovInvenio_Investigations/issues/44">issue 44</a>).</li>
</ul></div>
<h2>Definitions</h2>
<p><strong>Novelty</strong>: a gene family present in at least 75% of the ingroup (IN) species and absent from all outgroup (OUT) species. <strong>Loss</strong>: the reverse, present in at least 75% of the OUT species and absent from all IN species. Calls come from pairwise diamond searches (E-value &lt; 1e-5). Enrichment is a one-sided hypergeometric test per Pfam family, BH-corrected within each vocabulary.</p>
<h2>Contrasts</h2>
<div class="scroll"><table><thead><tr><th>id</th><th>role</th><th>lineage</th><th>ingroup / outgroup</th><th>pages</th></tr></thead><tbody>{contrast_rows}</tbody></table></div>
<h2>Files and provenance</h2>
<p class="note">Inputs, script commit, parameters and SHA-256 of every input and output are in PROVENANCE.yaml.</p>
<ul>{files}</ul>
</main></body></html>
"""
    (out / "comparison.html").write_text(html, encoding="utf-8")
    print(f"wrote {out/'comparison.html'}: {len(hits)} recurrent terms, {n_conv} yeast_convergent, contrasts with tables {len(have_tables)}")


if __name__ == "__main__":
    main()
