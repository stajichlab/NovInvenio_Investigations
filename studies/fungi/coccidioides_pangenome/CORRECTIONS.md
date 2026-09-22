# Data corrections log

## 2026-09-21: two mislabeled strains (species swap)

**What**: `485B-1_L_OLD_CPA0023` was labeled *Coccidioides immitis*; corrected to
*Coccidioides posadasii*. `B3476` was labeled *Coccidioides posadasii*; corrected
to *Coccidioides immitis*.

**Why**: `species.csv`'s `Species`/`TaxonGroup` columns were populated directly
from the species word embedded in the shared annotation-freeze directory's
filenames (`Coccidioides_<species>_<strain>.*`, see
`bin/build_coccidioides_species_csv.py`), which is itself the wrong label for
these two strains. Mash-distance clustering analysis
(`notes/pangenome-method-investigations/2026-09-20-phylogrouping-and-relatedness.md`)
found:
- `485B-1_L_OLD_CPA0023` sits at mash distance 0.000072 from a *C. posadasii*
  sibling in the same isolate series -- essentially identical, not a distinct
  species.
- `B3476` sits centrally in the *C. immitis* mash-distance cloud on a
  clean assembly (not a fragmented/ambiguous case).
- Across all 529 strains, these are the ONLY two that fall outside the
  0.9-1.5 mash-distance "hybrid window" that would suggest genuine
  uncertainty (they sit at 0.16 and 0.22) -- consistent with a data-entry
  swap upstream, not biological ambiguity.

**What changed**: `Species`/`TaxonGroup`/`Taxon_ID` in `species.csv` (the
source of truth), and the same two columns in the derived `config.csv`,
`config_immitis.csv`/`config_posadasii.csv` (rows moved between files to
match), and the study-specific `config_genus_vs_ureesii.csv` /
`config_immitis_in_posadasii_out.csv` / `config_posadasii_in_immitis_out.csv`
(GROUP column also flipped where that file's IN/OUT split is species-keyed).
`Protein_Accession`/`Genome_Accession`/`GFF3_Accession` paths are UNCHANGED --
they still point at the real sequence files on the shared annotation-freeze
mount; only the species/taxon LABEL is corrected, not which file is used.

**Not fixed by this commit**: the upstream shared annotation-freeze
directory's filenames themselves
(`/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/annotation_freeze/20260112/`)
still carry the old (wrong) species word. That's outside this repo and
outside this study's provenance boundary -- flagging for whoever owns that
directory, not renaming it here.

**Blocked, separately**: `bin/build_study_config.py --study-dir
studies/fungi/coccidioides_pangenome --skip-fetch` (the normal way to
regenerate `config.csv` from `species.csv`) currently fails for ALL local
strains, not just these two -- the shared freeze directory's DNA symlinks
resolve to `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/...`,
which does not exist (the real path is
`.../Population_Genomics/Coccidioides/...` vs the shared project having been
reorganized under a different root than the symlinks assume). This is an
unrelated shared-storage issue -- flag to whoever manages that mount before
relying on `build_study_config.py` again for this study; in the meantime the
derived config CSVs above were hand-corrected to match `species.csv` instead
of being regenerated.
