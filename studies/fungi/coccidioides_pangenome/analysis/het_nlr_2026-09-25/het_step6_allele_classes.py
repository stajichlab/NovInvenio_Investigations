#!/usr/bin/env python3
"""Step 6: allele classes per site. Families at one site whose reps share >= 95%
identity (het_rep_pairwise_annotated.tsv / het_locus_allele_pairwise.tsv) are
one class (length or gene-model variants). Classes below 95% are treated as
divergent allele classes. Class membership is written out explicitly below.
Frequencies use the own-family site occupancy from step 5."""
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
A = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/analysis/het_nlr_2026-09-25")
CLASSES = {
    "L005": {"L005-A (PGAP1-type)": ["H050", "H059", "H069", "H072", "H077"],
             "L005-B (Abhydrolase_6/DUF676-type)": ["H052", "H057", "H041"],
             "L005-C (ANAPC4_WD40-type)": ["H054"]},
    "L006": {"L006-A": ["H047", "H065"], "L006-B (Abhydrolase_6/DUF676-type)": ["H051", "H075"]},
    "L008": {"L008-fused Patatin-NB-ARC-TPR": ["H025"], "L008-NB-ARC-TPR only": ["H028"]},
    "L004": {"L004-fused (+Ank)": ["H026", "H062"], "L004-no Ank": ["H027", "H037"]},
    "L002": {"L002-G1 Goodbye": ["H009", "H015"], "L002-G2 Goodbye": ["H010"], "L002-G3 Goodbye": ["H011", "H013"],
             "L002-G4 Goodbye": ["H012"], "L002-H014 (unassigned fragment)": ["H014"]},
    "L001a": {"L001a NB-ARC fragments (one class)": ["H029", "H030", "H031", "H032", "H033", "H034", "H035", "H039", "H042", "H043", "H044"]},
    "L007": {"L007 HET-Het6_barrel (one class)": ["H001", "H003", "H004", "H005", "H006"]},
}
occ = pd.read_csv(A / "het_sites_ownfamily_occupancy.tsv.gz", sep="\t")
rows = []
for site, cl in CLASSES.items():
    o = occ[(occ.site == site) & occ.state.isin(["occupied", "empty_site"])]
    n = o.groupby("species").size()
    fams = o.families.fillna("").str.split(",")
    for name, members in cl.items():
        has = fams.apply(lambda l: bool(set(l) & set(members)))
        ci = int(has[o.species == "Ci"].sum()); cp = int(has[o.species == "Cp"].sum())
        nci, ncp = int(n.get("Ci", 0)), int(n.get("Cp", 0))
        orr, p = fisher_exact([[ci, nci - ci], [cp, ncp - cp]])
        rows.append(dict(site=site, allele_class=name, families=",".join(members), Ci_carriers=ci, Ci_resolved=nci,
                         Ci_frac=round(ci / nci, 3), Cp_carriers=cp, Cp_resolved=ncp, Cp_frac=round(cp / ncp, 3), fisher_p=p))
d = pd.DataFrame(rows)
d["bh_q"] = multipletests(d.fisher_p, method="fdr_bh")[1]
inr = lambda x: (x >= 0.05) & (x <= 0.95)
d["poly_within"] = inr(d.Ci_frac) | inr(d.Cp_frac)
d["poly_between"] = d.bh_q < 0.05
d["in_both_species_ge5pct"] = (d.Ci_frac >= 0.05) & (d.Cp_frac >= 0.05)
d.to_csv(A / "het_allele_classes.tsv", sep="\t", index=False)
print(d.to_string())
