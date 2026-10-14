#!/usr/bin/bash
# Convert the four Pfam domain tables in this folder into annotation files
# (studies/fungi/neocallimastigomycota_vs_early_fungi/annotations/<Protein_Accession>.tsv.gz
# + .provenance.yaml) with bin/hmm_domtbl_to_annotations.py. Runs in about 1 min on a login node.
# The two UniProt files (UP000193719, UP000193944) are not touched.
# Absolute paths only (no BASH_SOURCE).
set -euo pipefail
NII=/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations
S=$NII/studies/fungi/neocallimastigomycota_vs_early_fungi
H=$S/analysis/hmmscan
PFAM=/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm
export PATH=$NII/.pixi/envs/default/bin:$PATH
# Short|domain table|Protein_Accession|Taxon_ID (from species.csv)
for entry in \
    "Neocalif|Neocalif.pfam.hmmsearch.domtbl.zst|GCA_002104975.1|1754190" \
    "Neolanat|Neolanat.pfam.hmmsearch.domtbl.zst|GCA_016946835.1|2767002" \
    "Neoconst|Neoconst.pfam.hmmsearch.domtbl.zst|GCA_050613775.1|3108450" \
    "PiroE2sp|PiroE2sp.pfam.domtbl.zst|GCA_002157105.1|73868"; do
  IFS='|' read -r SHORT TBL ACC TAX <<< "$entry"
  python3 -I "$NII/bin/hmm_domtbl_to_annotations.py" --domtbl "$H/$TBL" \
      --fasta "$S/data_dir/pep/$SHORT.pep.fa" --taxon-id "$TAX" --pfam-hmm "$PFAM" \
      --out "$S/annotations/$ACC.tsv.gz"
done
