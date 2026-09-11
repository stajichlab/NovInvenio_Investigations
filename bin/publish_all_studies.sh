#!/usr/bin/bash
# Publish every study's alignment + report release assets in one pass, then
# trigger exactly one Pages deploy at the end.
#
# Why not just call bin/publish_alignment_release.sh / bin/publish_report_release.sh
# once per study: each of those, by default, triggers a full static.yml deploy run
# (`gh workflow run static.yml`) after uploading. static.yml's single deploy job
# re-downloads and re-merges EVERY alignments-*/reports-* release on every run,
# regardless of which one triggered it, and then rebuilds+redeploys the whole
# Pages site from scratch -- there is no partial/incremental merge and nothing is
# cached. So triggering it once per study when publishing several in a row is
# pure waste: `concurrency: {group: "pages", cancel-in-progress: false}` means
# the triggered runs queue and execute sequentially in full, not overlap, and
# only the LAST run's output is ever actually visible. This script uploads every
# study's release(s) with NII_SKIP_DEPLOY_TRIGGER=1 (see both scripts' own
# headers) and triggers the workflow itself, once, after the loop.
#
# Discovers studies by the presence of docs/<domain>/<set>/report.html (the same
# "this is a real published gallery entry" signal bin/sync_reports.sh's own
# domain/set index generation uses) -- not studies/*/*/, which also contains
# comparison/dev studies with no docs/ entry at all, and not raw docs/*/*/
# existence, which would also catch non-study directories (docs/ no longer
# has any such directories -- Mycelium artifacts like superpowers notes live
# under notes/superpowers/ instead -- but the same reasoning applies to any
# future non-study docs/ entry with no report.html).
#
# Usage: bin/publish_all_studies.sh [--dry-run]
#   --dry-run  list which studies would be published, upload nothing, trigger nothing.

set -euo pipefail

DRY_RUN=0
if [ "${1:-}" = "--dry-run" ]; then
    DRY_RUN=1
fi

REPO_ROOT="$(cd "$(dirname "${0}")/.." && pwd)"
cd "$REPO_ROOT"

published_any=0

for report in docs/*/*/report.html; do
    [ -e "$report" ] || continue
    study_dir="$(dirname "$report")"                  # docs/<domain>/<set>
    study="${study_dir#docs/}"                          # <domain>/<set>

    echo "== $study =="
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  (dry-run: would publish alignments-${study//\//-} and reports-${study//\//-})"
        continue
    fi

    NII_SKIP_DEPLOY_TRIGGER=1 bin/publish_alignment_release.sh "$study"
    NII_SKIP_DEPLOY_TRIGGER=1 bin/publish_report_release.sh "$study"
    published_any=1
done

if [ "$DRY_RUN" -eq 1 ]; then
    exit 0
fi

if [ "$published_any" -eq 0 ]; then
    echo "== no studies with docs/<domain>/<set>/report.html found -- nothing published, no deploy triggered =="
    exit 0
fi

echo "== triggering one Pages deploy for all studies above =="
gh workflow run static.yml
echo "== done =="
