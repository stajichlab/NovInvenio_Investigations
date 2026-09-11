# mmseqs cluster.tsv ID-restoration backfill status

Date: 2026-09-09

See `NovInvenio/docs/superpowers/specs/2026-09-09-mmseqs-cluster-id-restoration-design.md`
and `NovInvenio/docs/superpowers/plans/2026-09-09-mmseqs-cluster-id-restoration.md`
for the full bug and fix. Short version: `mmseqs easy-cluster` silently
collapses UniProt-style (`sp|`/`tr|`) and several NCBI-style FASTA headers to
a shorter field in its own `*_cluster.tsv` output, breaking gene-family
grouping, TBLASTN rep→member hit expansion, and alignment shards for any
UniProt-sourced study, until `NovInvenio/bin/restore_mmseqs_cluster_ids.py`
is run against the affected `results/<study>/clusters/*_cluster.tsv` file.

Any future Nextflow run of the pipeline (post `NovInvenio` commit `c862554`
and later) produces already-correct `*_cluster.tsv` files automatically —
this note only tracks **already-computed, pre-fix `results/` output** that
needed a one-time manual correction.

## Backfilled (both `clusters_cluster.tsv` and `loss_clusters_cluster.tsv` corrected)

- `pezizo_set1` — 3354 / 124 headers restored
- `zoosporic_dikarya` — 2734 / 1221 headers restored
- `mushrooms_tremella` — 5036 / 2395 headers restored
- `yeast_filamentous` — 335 / 5637 headers restored

Command used per study (run from a `NovInvenio` checkout, using `pixi run
python` — the bare system `python3` is too old for `lib/fasta.py`'s union-type
syntax):

```bash
pixi run python bin/restore_mmseqs_cluster_ids.py \
    --input-fasta <NII_ROOT>/results/<study>/candidates.fa \
    --cluster-tsv <NII_ROOT>/results/<study>/clusters/clusters_cluster.tsv
pixi run python bin/restore_mmseqs_cluster_ids.py \
    --input-fasta <NII_ROOT>/results/<study>/loss_candidates.fa \
    --cluster-tsv <NII_ROOT>/results/<study>/clusters/loss_clusters_cluster.tsv
```

## Not affected (confirmed MAG/prodigal headers, no `sp|`/`tr|`/NCBI-style prefix — `mmseqs_id()` is a no-op)

- `UHM_Akkermansia`
- `UHM_Koxytoca`

## Deferred (out of scope for this fix, blocks re-verifying TBLASTN rep→member expansion)

`pezizo_set1`'s `results/pezizo_set1/tblastn/*.tblastn.tsv` files are stale
(10 columns; the pipeline's TBLASTN module now emits 13) and have already
been removed from disk — `bin/summarize_tblastn.py` needs fresh TBLASTN
output regenerated before its rep→member hit-expansion fix (a consequence of
this same cluster.tsv correction) can be re-verified end to end for that
study.
