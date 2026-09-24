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
# Discovers runs by docs/<domain>/<study>/<run>/run.json (study/run layout,
# notes/superpowers/specs/2026-09-24-study-run-site-layout-design.md). run.json
# is written only by the sync scripts for a folder that has a publish.yaml, so
# unpublished comparison/dev folders and old flat docs/<domain>/<set>/ folders
# are never picked up.
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

# Study/run layout: one release pair per run (docs/<domain>/<study>/<run>/run.json).
# Old flat docs/<domain>/<set>/ folders without run.json are not published here.
for meta in docs/*/*/*/run.json; do
    [ -e "$meta" ] || continue
    run_dir="$(dirname "$meta")"                       # docs/<domain>/<study>/<run>
    study="${run_dir#docs/}"                           # <domain>/<study>/<run>
    IFS=/ read -r d s r <<< "$study"

    echo "== $study =="
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  (dry-run: would publish alignments-$d-$s--$r and reports-$d-$s--$r)"
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
    echo "== no runs with docs/<domain>/<study>/<run>/run.json found -- nothing published, no deploy triggered =="
    exit 0
fi

echo "== triggering one Pages deploy for all studies above =="
gh workflow run static.yml
echo "== done =="
