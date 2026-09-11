# UHM_lachnoNovelclade

Lachnospiraceae novelty/loss study: 5 UHM metagenome-assembled genomes
(an outgroup clade plus one ingroup MAG of uncertain taxonomic placement,
hence "novel clade") vs. 4 NCBI reference ingroup genomes for named
*Lachnospiraceae*/*Lacrimispora*/*Enterocloster* species.

## Datasets

- **Outgroup (`Group=OUT`, `TaxonGroup=Lachnospiraceae`)** — 4 UHM
  metagenome-assembled genomes (`Ub56`, `Ub85`, `Ub22`, `Ub119`), unclassified
  below family level (`Species=Lachnospiraceae sp.`, `Strain=`<original MAG
  bin id>). Both protein (`Protein_Source=local_faa`, from
  `/bigdata/stajichlab/jpere468/unknown_tree/clade1_faa/`) and genome
  (`Genome_Source=local_genome`, from
  `/bigdata/stajichlab/jpere468/drep_herptile_95/high_quality_genomes/`) are
  local files, copied in by `bin/build_study_config.py` (no fetch).
- **Ingroup (`Group=IN`, `TaxonGroup=Lachnospiraceae`)** — 5 species:
  - `Ub71` — one more unclassified UHM MAG bin, same sourcing as the outgroup
    rows above (`local_faa`/`local_genome`).
  - `Lsp1` (*Lacrimispora* sp.), `Lsph` (*Lacrimispora sphenoides*), `Easp`
    (*Enterocloster asparagiformis*), `Cami` (*[Clostridium] aminophilum*) —
    4 NCBI reference genomes. Protein is a local file
    (`Protein_Source=local_faa`, already pulled into
    `/bigdata/stajichlab/jpere468/unknown_tree/lachno_ingroup_faa/<GCF_*>.faa`
    ahead of this study); genome is fetched live from NCBI
    (`Genome_Source=ncbi`, keyed by the same `GCF_*` accession as the protein
    file's name) since no local genome copy exists for these 4.

Species/strain identity for the 5 MAG-bin rows and the exact NCBI accession
→ species mapping for the 4 reference genomes were supplied directly by the
study author (2026-09-11) — no automated classification (e.g. GTDB-tk) was
run as part of onboarding this study.

Full per-file provenance (source URL/path, release, checksum, license):
`DATA_MANIFEST.yaml`.

## Building

```bash
pixi run python bin/build_study_config.py --study-dir studies/bacteria/UHM_lachnoNovelclade
```

Already run once (2026-09-11) to produce `config.csv`/`DATA_MANIFEST.yaml` and
this study's local `data_dir/` (gitignored, regenerable from `species.csv` via
the command above). Safe to re-run any time `species.csv` changes.

## Running

No `run_params.txt`/`run.sh` yet — this study hasn't been run through
`nf_NovInvenio` yet. See `studies/bacteria/UHM_Koxytoca/run_Koxytoca.sh` for
the pattern once run parameters (`--cluster_tool`, `--run_tool`, etc.) are
decided for this study.
