# Afumigatus_pangenome

Pangenome cluster-profile analysis for 293 *Aspergillus fumigatus* strains (+2
outgroup: *A. lentulus*, *A. fischeri*) — detecting Starship-driven gene
gain/loss, including a targeted screen for the HAC (hrmA-Associated Cluster)
gene family.

- **Method/design**: `notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md`
  (in the NII repo root) — read this for the *why*.
- **Findings/decisions/open items for this study specifically**:
  `PANGENOME_CLUSTER_PROFILE_NOTES.md` — read this for what's actually been run
  and what it found.
- **This file**: the *how* — copy-paste commands to reproduce or continue the
  analysis, and where things currently stand.

Every command below sets `STUDY=/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome`
(a symlink to this study's real location in the NII repo) so paths work
whether your cwd is the NII repo root or `$RUN`. `pixi run` (add
`--manifest-path /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/pixi.toml`
if you're not in that directory) resolves mmseqs2/diamond/blast/mash/hmmer
without an activated shell. `$RUN` is wherever you're writing output
(`$STUDY/results/<name>/` — `results/` is gitignored, never commit run
output). A few steps (mash, mmseqs) write scratch files to cwd and need
`cd "$RUN"` first — noted per step below.

## Current state (as of 2026-09-13)

Done, on the real 293+2-strain data:

- [x] `config.csv` GROUP/TaxonGroup columns fixed and filled (clade labels —
  see notes item 6)
- [x] Mash+PCoA clade assignment (`bin/assign_clades.py`)
- [x] HAC/hacA reference-sequence screen (`bin/hac_reference_screen.py`)
- [x] Tier-1 mmseqs clustering (`results/full_293run/tier1_cluster.tsv`)
- [x] Presence matrix (`results/full_293run/presence_matrix.tsv`)
- [x] Strain dereplication (`results/full_293run/strain_inventory.tsv`)
- [x] Frequency binning (`results/full_293run/frequency_table.tsv`)

Not yet done:

- [ ] Isoform collapse (skipped so far as a checked, documented assumption —
  see notes item 6; not a proper fix)
- [ ] Genome-level tblastn rescue pass
- [ ] Co-occurrence — **blocked on implementing the scaling mitigations the
  design spec's component 10/Opus review call for** (bitset presence vectors,
  streaming BH, analytic prefilter) before running at 47,983-family scale; do
  not run `cooccurrence.py` as-is against the real data yet
- [ ] Synteny / accessory-island detection
- [ ] Pair classification (component 8)
- [ ] `pangenome.html` report (component 9)
- [ ] ParSNP full run — deliberately deferred (Mash+PCoA judged sufficient for
  now); `run_parsnp_ingroup293.sh` is ready if revisited

## 1. Clade assignment (Mash + PCoA + k-means)

```bash
STUDY=/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome
RUN=$STUDY/results/clade_assignment
mkdir -p "$RUN"
cd "$RUN"   # mash writes sketch files to cwd
export NOVINVENIO_ROOT=/bigdata/stajichlab/jstajich/projects/NovInvenio
pixi run --manifest-path /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/pixi.toml \
  python3 "$STUDY/bin/assign_clades.py" \
  --config "$STUDY/config.csv" \
  --data_dir "$STUDY/data_dir" \
  --groups IN \
  --k_range 2,15 \
  --output "$RUN/clade_assignments.tsv"
```

`--groups IN` (default) is required — mixing in the outgroups swamps the
within-species clade signal (a 294-vs-1 split, confirmed on this data). Use
`--k_fixed N` to skip the silhouette sweep and force a specific k.

## 2. Build the Short-prefixed, all-strain concatenated proteome

No dedicated script for this yet — a one-liner per the ID convention in
`bin/build_presence_matrix.py`'s docstring:

```bash
STUDY=/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome
RUN=$STUDY/results/<run_name>
mkdir -p "$RUN"
cd "$STUDY"
: > "$RUN/all_strains.fa"
awk -F, 'NR>1 && ($1=="IN" || $1=="OUT") {print $7","$4}' config.csv | \
while IFS=, read -r short protein; do
  awk -v s="$short" '/^>/{sub(/^>/,">"s"|")}1' "data_dir/pep/$protein" >> "$RUN/all_strains.fa"
done
```

No isoform collapsing is applied here yet (see "Not yet done" above) — a
checked, not resolved, simplification for *A. fumigatus*'s low alt-splicing
rate.

## 3. Tier-1 clustering

```bash
STUDY=/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome
RUN=$STUDY/results/<run_name>
cd "$RUN"
export NOVINVENIO_ROOT=/bigdata/stajichlab/jstajich/projects/NovInvenio
pixi run --manifest-path /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/pixi.toml \
  python3 "$STUDY/bin/cluster_backend.py" mmseqs-tier1 \
  --fasta "$RUN/all_strains.fa" --out_prefix "$RUN/tier1"
```

Real timing on all 295 strains (2,788,402 proteins): ~14 minutes wall-clock,
47,983 families. Delete `tmp_mmseqs/` and `tier1_all_seqs.fasta` afterward —
both are large, regenerable mmseqs intermediates not needed downstream (only
`tier1_cluster.tsv` and `tier1_rep_seq.fasta` are).

## 4. Presence matrix

Run from the NII repo root (`cd /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations`),
not from `$RUN` — this script and the ones below don't write to cwd:

```bash
pixi run python3 "$STUDY/bin/build_presence_matrix.py" \
  --cluster_tsv "$RUN/tier1_cluster.tsv" \
  --config "$STUDY/config.csv" \
  --output "$RUN/presence_matrix.tsv"
```

Default `--groups IN,OUT` — needed so `cooccurrence.py`'s outgroup gain/loss
polarization has outgroup columns to read (see the design spec's component 8
must-fix M7 for a real caveat about this at 90%-identity tier-1 clustering).

## 5. Strain dereplication (Mash/ANI)

Writes mash sketch files to cwd, so `cd "$RUN"` first, using the same absolute
`$STUDY` variable set in step 1/3 (don't reuse a relative path here):

```bash
cd "$RUN"
pixi run --manifest-path /bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations/pixi.toml \
  python3 "$STUDY/bin/dereplicate_strains.py" \
  --config "$STUDY/config.csv" \
  --data_dir "$STUDY/data_dir" \
  --output "$RUN/strain_inventory.tsv"
```

Real result on this data: 123 representative strains from 295 at the default
Mash threshold (0.001) — a substantial reduction, worth sanity-checking against
a stricter/looser threshold before trusting it blindly.

## 6. Frequency binning

Back in the NII repo root (not `$RUN`):

```bash
pixi run python3 "$STUDY/bin/frequency_bins.py" \
  --matrix "$RUN/presence_matrix.tsv" \
  --config "$STUDY/config.csv" \
  --inventory "$RUN/strain_inventory.tsv" \
  --output "$RUN/frequency_table.tsv"
```

Real result: 6,600 core / 582 soft-core / 3,943 shell / 9,238 cloud / 27,620
singleton (of 47,983 families, 121 dereplicated ingroup strains). The 58%
singleton fraction needs the completeness/annotation-source check (design
spec must-fix M6) before trusting it as pure biology rather than partly an
annotation-heterogeneity artifact.

## 7. Co-occurrence — DO NOT RUN YET at full scale

`bin/cooccurrence.py` exists and is tested, but its all-pairs enumeration over
shell/cloud/singleton-adjacent families is O(n²) in family count (design
spec's Opus-review must-fix M8). At this study's real scale (~40k eligible
families) that is untested and likely to hang or exhaust memory. Implement the
mitigations in the design spec's component 10 step 7 first (bitset presence
vectors, an analytic prefilter, streaming BH, permutation only on
FDR-survivors), or run on a deliberately restricted family subset first to
characterize actual runtime before the real run.

## 8. Synteny / accessory islands, and pair classification

Not yet run against real data. `bin/synteny_windows.py`'s CLI is intentionally
a stub (its functions are meant to be called directly for a real run — see its
own module docstring); component 8 (`pair_classification.tsv`) needs a new
captain-gene (DUF3435) hmmsearch step that doesn't exist yet either — see the
design spec's component 8 for the full, revised procedure.

## HAC/hacA reference screen (already run, independent of the chain above)

This doesn't need tier-1 clustering — it's a direct diamond/hmmsearch screen
against the raw proteomes:

```bash
STUDY=/rhome/jstajich/projects/NII/studies/fungi/Afumigatus_pangenome
RUN=$STUDY/results/hac_crosswalk
mkdir -p "$RUN"

# 1. Build all_ingroup293.fa (step 2 above, IN only, output named all_ingroup293.fa)
# 2. Extract the two reference query sequences (hacA=Afu3g04070, hrmA=Afu5g14900)
#    from the local FungiDB AF293 proteome into $RUN/hac_reference_queries.fa
#    (headers Afu3g04070-T-p1 / Afu5g14900-T-p1 in
#    NovInvenio/db/modelorgs/FungiDB-68_AfumigatusAf293_AnnotatedProteins.fasta)
# 3. Build the diamond database and run both searches:
pixi run diamond makedb --in "$RUN/all_ingroup293.fa" -d "$RUN/all_ingroup293" -p 8
pixi run diamond blastp -q "$RUN/hac_reference_queries.fa" -d "$RUN/all_ingroup293" \
  -o "$RUN/hacA_hrmA_vs_study.tsv" --outfmt 6 --evalue 1e-5 --max-target-seqs 400 -p 8
pixi run hmmfetch /bigdata/stajichlab/jstajich/projects/NovInvenio/db/pfam/Pfam-A.hmm \
  AFUB_07903_YDR124W_hel > "$RUN/PF11001.hmm"   # PF11001's NAME field, not its accession
pixi run hmmsearch --tblout "$RUN/PF11001_vs_study.tblout" -E 1e-3 --cpu 8 \
  "$RUN/PF11001.hmm" "$RUN/all_ingroup293.fa"
pixi run python3 "$STUDY/bin/hac_reference_screen.py" \
  --diamond_tsv "$RUN/hacA_hrmA_vs_study.tsv" \
  --hmmsearch_tblout "$RUN/PF11001_vs_study.tblout" \
  --config "$STUDY/config.csv" \
  --output "$RUN/hac_reference_screen.tsv"

# delete all_ingroup293.fa / .dmnd afterward -- large, regenerable, not needed once
# the tblout/tsv outputs above exist
```

Real result: hacA present in all 293 strains; the true hrmA ortholog (not just
any PF11001-family member — see the notes for why that distinction matters)
in only 8/293. Full writeup with Starship cross-referencing in the notes file.

## Tests

From the NII repo root:

```bash
pixi run pytest studies/fungi/Afumigatus_pangenome/tests/ -v
```

78 tests as of this writing, all passing.
