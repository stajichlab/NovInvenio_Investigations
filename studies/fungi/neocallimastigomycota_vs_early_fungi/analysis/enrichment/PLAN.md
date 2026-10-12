# Enrichment of novelty candidates: strategy and prototype (draft, 2026-10-12)

Status: prototype only. Nothing here is committed or published. Numbers come from
`outputs/neocallimastigomycota_vs_early_fungi/` (run of `enrich_novelties.py` on
`results/neocallimastigomycota_vs_early_fungi`, 6 IN species, pairwise-diamond).

## Inputs (what exists, checked)

- `presence_matrix.uniprot.tsv`: one row per IN protein that the pipeline tested
  (103,108 rows). Columns `uniprot_go_ids` (`GO:nnnnnnn:EVIDENCE|...`),
  `uniprot_pfam_ids` (`PFnnnnn|...`), `uniprot_pfam_names` (parallel list),
  `uniprot_interpro_ids` (`IPRnnnnnn|...`). Presence/absence per protein only. No
  domain coordinates, no copy number.
- `candidates.txt`: `<Short>::<protein_id>`, 32,717 lines.
- `novelties.<Short>.tsv`: `family_id` (mmseqs representative). Empty `family_id`
  means singleton. 11,296 families; 6,470 singletons; largest family 468 proteins.
- `studies/.../annotations/<UPID>.tsv.gz`: same UniProt fields, keyed by bare
  accession. Only the 15 UniProt proteomes have one.
- `data/go/go-basic.obo`: release 2026-07-26, provenance sidecar present.
  `data/` is gitignored; re-fetch with `bin/fetch_go_obo.py`.
- Pfam clans: `/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat`
  (`#=GF CL`) and `Pfam-C` (clan names). Not in the repo. UniProt's own Pfam release
  for these entries is not recorded in the annotation TSV (not checked).
- InterPro names/types: `/srv/projects/db/interproscan/6.0.0/110.0/interpro/110.0/entries.json`.
  Labels only.

Existing tools: `bin/domain_enrichment.py` (Pfam/InterPro, one-sided hypergeometric,
BH) and `bin/go_enrichment.py` (goatools). The prototype reproduces
`domain_enrichment.py` exactly for Pirofinn protein-level Pfam (same 404 terms, same
p and q values). For GO the p-values match goatools on the top terms, but see
finding F1.

## Plan

1. **Background, scope A (whole set).** Universe = all IN proteins that the pipeline
   tested (rows of `presence_matrix.uniprot.tsv`), restricted per vocabulary to
   proteins with >=1 term. Do not use OUT proteins. Reason: candidates are drawn
   from IN proteomes only. OUT proteins can never be candidates. Adding them mixes
   lineage composition into the test. Note: the matrix has fewer proteins than the
   input FASTA (Pirofinn 10,284 of 10,768; Anaerobu 11,375 of 12,656; PiroE2sp
   12,863 of 14,612). The cause of the filter is not verified. Use the matrix rows,
   because those are the proteins that could become candidates.
2. **Background, scope B (per species).** Universe = that species' rows in the
   matrix, restricted to annotated proteins. Candidates = that species' candidates.
   This removes proteome-size and annotation-depth differences between species.
3. **Redundancy.** Run two levels and report both:
   - `protein`: one unit per protein.
   - `family`: candidates collapsed to one unit per mmseqs family (scope A: across
     species; scope B: within the species). Unit terms = union of annotated
     members' terms. The background is not clustered by the pipeline, so
     background proteins stay single units. This is conservative only on the
     candidate side.
   Effect here: whole-set Pfam terms at FDR<0.05 drop from 66 (protein) to 42
   (family); GO from 97 to 26. Treat a term as robust only if it passes at the
   family level. Later step: cluster the background too (mmseqs on the IN
   proteomes, or the existing `<sp>_vs_<sp>.diamond.tsv.gz` self hits) so both
   sides use the same unit.
4. **Annotation coverage.** Measured (coverage.tsv):
   - Pirofinn and Anaerobu (UniProt): 76-78% of all proteins have >=1 GO/Pfam/IPR
     term; candidates only 43% and 50%.
   - Neocalif, Neoconst, Neolanat, PiroE2sp (NCBI): 0% annotated (0 of 81,449
     proteins). `uniprot_match=seq_other` rows (137) also carry no terms.
   - Whole set: 3,015 of 32,717 candidates (9.2%) have any term.
   - 13,501 candidates (41%) sit in 6,934 families with no Pirofinn/Anaerobu
     member, so family transfer cannot reach them.
   Handling: (a) restrict each universe to annotated units (done); (b) print
   coverage next to every result; (c) skip species with no annotated candidates
   and say so (done: status `skipped_no_annotated_candidates`); (d) family-level
   transfer gives NCBI proteins the terms of their UniProt family members (done at
   scope A family level); (e) later, optional: run Pfam hmmscan on the 4 NCBI
   proteomes (81,449 proteins) plus unannotated UniProt candidates. Cost not
   measured. The only local timing is `ANNOTATE_MATRIX` in agaricomycetes_pairwise
   (9,961 candidates, 11-16 min realtime); whether that task ran hmmscan is not
   verified. InterProScan 6 is installed (`/srv/projects/db/interproscan/6.0.0`) but
   no local timing exists. The 4 NCBI proteomes lack UniProt-style IEA GO; GO would
   need InterPro2GO on the InterProScan output. Restricting to annotated proteins
   does not remove bias: annotated candidates are a non-random subset (known
   domains only). Say this in every report.
5. **Statistics.** One-sided hypergeometric P(X>=k). Test only terms with >=5
   background units (`--min-bg`, filter on background count only). BH within
   each test (scope x level x vocabulary) = `fdr_bh`. For scope B also BH across
   all species tests of the same level x vocabulary = `fdr_bh_all_species`. Report
   the number of tested terms per test (`tests.tsv`). Report a term only with
   k>=3 in the summary page (display filter, not a test filter).
6. **GO.** True-path propagation over `is_a` (goatools `get_all_parents`), drop the
   3 namespace roots. All 45,748 GO assignments in this run are IEA. So GO is
   derived from InterPro/Pfam here and is not independent evidence. Do not count a
   GO hit and its source domain as two findings. `--exclude-iea` would remove all
   GO for this study.
7. **Pfam clans and InterPro levels.** Test Pfam families and Pfam clans as
   separate vocabularies. Clans merge split families (Ank_2/Ank/Ank_4 -> CL0465).
   Many Pfam families have no clan (CBM_10 has none), so the clan test does not
   replace the family test. InterPro entries nest (Family/Domain/Homologous
   superfamily); report the entry type and expect redundant hits.
8. **Domain-level vs whole-protein.** The candidate call is whole-protein. UniProt
   annotation is per protein. So the test asks "are proteins that carry domain X
   over-represented among whole-protein novelties". It does not find novel domain
   architectures. Evidence from this run: 487 IN proteins carry CBM_10
   (dockerin); 183 are candidates; the other 304 have an OUT hit (most to
   Umbevina 227, Rhizdele 164, Mucocirc 159, Neurcras 141), probably through the
   fused catalytic domain (not verified). A domain-level test needs domain
   coordinates (hmmscan) and per-domain presence calls. That is a pipeline change.
9. **Outputs.** Per run: `coverage.tsv`, `tests.tsv`, `enrichment.tsv.gz` (all
   tested terms; candidate IDs only for terms with fdr_bh<0.05). Publishing (design
   only): a script like `bin/sync_pangenome_report.py` stages
   `docs/<domain>/<study>/<run>/enrichment.html` (sortable table, one tab per
   scope x level x vocabulary, coverage banner at top) and
   `archive/enrichment.tsv.gz`. Because the candidate ID lists and the number of
   significant terms grow with candidate count, both are class 3 (release asset).
   The coverage numbers (constant size) go into `report.html` as one table with a
   link to `enrichment.html`.
10. **Placement.** The prototype takes a results dir and is not study-specific.
    If adopted, move it to `bin/enrich_novelties.py` and wire it after
    `bin/sync_reports.sh`. Add tests (fixture with a known hypergeometric answer;
    a family-collapse case; a species with 0 annotation).
11. **Next studies** (coverage measured with `--coverage-only`, outputs in
    `outputs/coverage_other_studies/`):
    - `chytrid_vs_mucoromycota`: 4,868 candidates, 86% annotated, 11 IN all
      UniProt. Best next run; per-species tests possible (337-548 candidates each).
    - `zoosporic_dikarya`: 2,927 candidates, 88% annotated, 4 IN. Good.
    - `serinales_vs_pichiales`: 2,995 candidates, 48% annotated, 18 IN. Scope A
      fine; per-species tests have 23-205 candidates each, low power.
    - `malassezia_vs_ustilaginomycotina`: 530 candidates, 10% annotated (about 54).
      Scope A only; expect few or no terms.
    - `metschnikowiaceae_cah`: 206 candidates, 18% annotated (about 38). Scope A
      only; probably underpowered.

## Findings about existing code (not fixed; tracked files)

- F1. `bin/go_enrichment.py` says "one-sided" but goatools' default
  `fisher_scipy_stats` calls `scipy.stats.fisher_exact` with the default
  two-sided alternative (`goatools/pvalcalc.py`). Then the script keeps only
  `enrichment == 'e'` terms and runs BH on that subset only. Two-sided p is
  conservative; BH on a direction-filtered subset is anti-conservative. Example
  (Pirofinn, protein GO): its q for GO:0005576 is 4.1e-9; the prototype (BH over
  all 2,358 tested terms) gives 2.0e-8.
- F2. `bin/go_enrichment.py` / `domain_enrichment.py` take bare accessions and the
  per-UPID annotation TSV. They cannot read NCBI proteins (no annotation file) and
  cannot do family collapse. The prototype reads the merged matrix instead.

## Not verified

- Why the matrix has fewer proteins than the FASTA inputs.
- Whether OUT hits for GH48 (all 17 hit Umbevina only) and PFL (all 38 hit Cateangu
  only) are true homologs or spurious; this decides whether these known anaerobic
  fungal genes should be candidates at all.
- Biological interpretation of any term. No literature check was done here.
- hmmscan / InterProScan runtime.
