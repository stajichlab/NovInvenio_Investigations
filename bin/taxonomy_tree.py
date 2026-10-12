"""Cladogram from NCBI taxonomy lineages for the species in a study's species.csv.
Usage: taxtree.py <species.csv> <outdir>
Branching nodes are labelled with their NCBI taxon name. Unary nodes are collapsed.
Branch lengths are not divergence times (all 1 per collapsed step)."""
import csv, json, sys, urllib.request
from pathlib import Path
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

API = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{}/dataset_report"
sp = list(csv.DictReader(open(sys.argv[1])))
out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)

def fetch(ids):
    r = urllib.request.urlopen(API.format(",".join(map(str, ids))), timeout=60)
    return json.load(r)["reports"]

lin = {}
for r in fetch([s["Taxon_ID"] for s in sp]):
    t = r["taxonomy"]
    lin[str(t["tax_id"])] = [str(p) for p in t["parents"]] + [str(t["tax_id"])]

# build nested dict tree
root = {}
for s in sp:
    node = root
    for tid in lin[s["Taxon_ID"]]:
        node = node.setdefault(tid, {})
    node["__leaf__"] = s

internal = set()
def walk(n):
    for k, v in n.items():
        if k != "__leaf__":
            internal.add(k); walk(v)
walk(root)
names = {}
ids = sorted(internal)
for i in range(0, len(ids), 50):
    for r in fetch(ids[i:i + 50]):
        t = r["taxonomy"]; names[str(t["tax_id"])] = t["current_scientific_name"]["name"]

def build(n, label):
    kids = [(k, v) for k, v in n.items() if k != "__leaf__"]
    if "__leaf__" in n and not kids:
        s = n["__leaf__"]
        return Clade(branch_length=1, name=f'{s["Species"]} [{s["Group"]}]')
    subs = [build(v, names.get(k, k)) for k, v in kids]
    # collapse unary nodes
    if len(subs) == 1:
        c = subs[0]; c.branch_length = (c.branch_length or 0) + 1
        return c
    return Clade(branch_length=1, name=label, clades=subs)

kids = [(k, v) for k, v in root.items()]
tree = Tree(root=Clade(clades=[build(v, names.get(k, k)) for k, v in kids]), rooted=True)
# drop the unary chain from the top: use the first branching node as root
c = tree.root
while len(c.clades) == 1:
    c = c.clades[0]
tree = Tree(root=c, rooted=True)
tree.root.branch_length = 0
Phylo.write(tree, out / "taxonomy_tree.nwk", "newick")

def lab(c):
    return c.name if c.name else None
for cl in tree.find_clades():
    cl.branch_length = None
fig, ax = plt.subplots(figsize=(9, 5))
Phylo.draw(tree, axes=ax, do_show=False, label_func=lab)
for t in ax.texts:
    if "[IN]" in t.get_text():
        t.set_color("#b34700"); t.set_fontweight("bold")
    elif "[OUT]" in t.get_text():
        t.set_color("#1f5fa8")
ax.set_title("NCBI taxonomy cladogram (topology only; branch lengths not scaled)\norange = IN, blue = OUT", fontsize=10)
ax.set_xlabel(""); ax.set_ylabel("")
fig.tight_layout()
fig.savefig(out / "taxonomy_tree.png", dpi=160)
fig.savefig(out / "taxonomy_tree.svg")
print(open(out / "taxonomy_tree.nwk").read())
