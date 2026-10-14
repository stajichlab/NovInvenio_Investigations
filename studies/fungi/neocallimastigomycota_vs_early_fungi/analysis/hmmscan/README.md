# Pfam annotation of the 4 NCBI/JGI IN proteomes

Why: 4 of the 6 IN species (Neocalif, Neolanat, Neoconst, PiroE2sp) come from NCBI
GenBank assemblies (JGI annotation). They have no UniProt record, so
`bin/build_study_config.py` wrote no `annotations/<Protein_Accession>.tsv.gz` for them.
This folder holds the Pfam-A search of these proteomes and the scripts that made it.

Database: Pfam-A 38.2, `/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm`
(release from the directory name). HMMER 3.4 (Aug 2023), `--cut_ga`.

## Files

| File | What | Made by |
|---|---|---|
| `run_hmmscan.sh` | pilot: hmmscan of one proteome | job 29887866 (PiroE2sp, 14,612 proteins, 16 cores, 9,371 s) |
| `PiroE2sp.pfam.domtbl.zst` | hmmscan domtblout (col 1 Pfam name, col 2 Pfam accession, col 4 protein) | `run_hmmscan.sh` |
| `bench_hmmsearch_vs_hmmscan.sh` | 500 PiroE2sp proteins, 4 cores: hmmscan 362 s, hmmsearch 81 s, same 525 (protein, Pfam) pairs | job 29943954 |
| `run_hmmsearch.sh` | hmmsearch of one proteome, 8 cores | jobs 29944241-29944243 |
| `Neocalif.pfam.hmmsearch.domtbl.zst` | hmmsearch domtblout (col 1 protein, col 4 Pfam name, col 5 Pfam accession); 20,214 proteins, 1,911 s | job 29944241 |
| `Neolanat.pfam.hmmsearch.domtbl.zst` | same; 27,675 proteins, 2,765 s | job 29944242 |
| `Neoconst.pfam.hmmsearch.domtbl.zst` | same; 23,663 proteins, 2,448 s | job 29944243 |
| `convert_to_annotations.sh` | converts the 4 tables into `../../annotations/<Protein_Accession>.tsv.gz` + `.provenance.yaml` with `bin/hmm_domtbl_to_annotations.py` (about 6 s, login node) | 2026-10-14 |
| `validation_PiroE2sp_hmmsearch/` | test of `bin/annotate_proteomes_pfam.py --submit` on PiroE2sp, split in 2 chunks of 7,306 proteins (jobs 29956183/29956184: 635 s and 717 s; finalize 29956185) | 2026-10-14 |

## Commands

```bash
# convert the existing tables (does not touch the UniProt files UP000193719, UP000193944)
bash studies/fungi/neocallimastigomycota_vs_early_fungi/analysis/hmmscan/convert_to_annotations.sh

# the general tool for a new study (writes jobs under <study>/analysis/pfam_hmmsearch/)
bin/annotate_proteomes_pfam.py --study fungi/neocallimastigomycota_vs_early_fungi --dry-run
```

## Results (2026-10-14)

| Short | Protein_Accession | Proteins | With >= 1 Pfam | % |
|---|---|---|---|---|
| Neocalif | GCA_002104975.1 | 20,214 | 13,069 | 64.7 |
| Neolanat | GCA_016946835.1 | 27,675 | 14,712 | 53.2 |
| Neoconst | GCA_050613775.1 | 23,663 | 14,423 | 61.0 |
| PiroE2sp | GCA_002157105.1 | 14,612 | 8,158 | 55.8 |

For comparison, the UniProt Pfam cross-references of the 2 UniProt IN species give
67.0% (Pirofinn, 7,219 of 10,768) and 65.1% (Anaerobu, 8,234 of 12,656). These come from
a different pipeline (UniProt/InterPro) and Pfam release, so the numbers are not strictly
comparable.

hmmscan vs hmmsearch check (PiroE2sp, full proteome): the hmmscan table and the
validation hmmsearch table have the same 22,484 domain rows (protein, Pfam accession,
envelope start, envelope end). The two annotation files are identical in every column
except `pfam_best_ievalue`: the i-Evalue depends on the search-space size, which differs
between hmmscan, hmmsearch and a chunked hmmsearch. Bit scores are identical.

## Limitation

These files have Pfam only. GO and InterPro columns are empty. In the enrichment
analysis Pfam and Pfam-clan tests cover 6 of 6 IN species; GO and InterPro tests cover
2 of 6 (Pirofinn, Anaerobu). The enrichment page states this.

## Tracking

`annotations/` is gitignored (CLAUDE.md: regenerable). The `.domtbl.zst` files in this
folder are not tracked either; they are the inputs to regenerate the annotation files
and cost about 40 min of compute per proteome to remake.
