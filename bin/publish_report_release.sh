#!/usr/bin/bash
# Package a study's report pages (docs/<domain>/<set>/{novelties,core,losses}.html,
# summary.pdf) as a GitHub Release asset, so a rebuilt report reaches the
# published site without ever being git-committed -- the same mechanism
# bin/publish_alignment_release.sh already uses for TBLASTN alignment shards
# (see DESIGN.md Sec 8 and issue #3's .gitignore rule), extended here to cover
# the report HTML/PDF files themselves.
#
# Usage: bin/publish_report_release.sh <domain>/<set_name>
#        bin/publish_report_release.sh <domain>/<set_name>/<run>
#
# The 3-part form publishes one pangenome.nf run's release-only files
# (docs/<domain>/<set>/<run>/{figures,figures_pdf,archive}/, island_synteny.html,
# assembly_quality.html -- staged by bin/sync_pangenome_report.py). Tag:
# reports-<domain>-<set>--<run>; manifest.json then also carries "run", and
# static.yml merges the files into docs/<domain>/<set>/<run>/.
#
# Tag: reports-<domain>-<set> (one release per study, overwritten in place on
# every rerun via --clobber -- matches how alignments-<domain>-<set> and
# docs/<domain>/<set>/report.html itself are already always-overwritten with
# no run history).
#
# The tarball carries a manifest.json ({"domain": ..., "set": ...}) at its root
# alongside whichever of novelties.html/core.html/losses.html/summary.pdf exist,
# so the deploy workflow (.github/workflows/static.yml) can recover the
# destination path by reading the manifest instead of parsing the tag name --
# robust regardless of how a domain/set name is spelled.
#
# --clobber does NOT re-fire the release:published webhook (that only fires once,
# at creation) -- this script explicitly triggers the Pages deploy workflow
# afterward via `gh workflow run` (workflow_dispatch), since a rerun would
# otherwise silently never redeploy.
#
# This is a deliberate, standalone publish step -- NOT wired into
# bin/sync_reports.sh's automatic post-run chain, since it creates a real GitHub
# Release and triggers a real Pages deploy (external, visible side effects).
# Run it by hand once a study's regenerated reports are ready to actually go live.
#
# Set NII_SKIP_DEPLOY_TRIGGER=1 to upload the release without triggering
# static.yml -- see bin/publish_alignment_release.sh's matching note /
# bin/publish_all_studies.sh, which uses this to publish every study first and
# trigger exactly one deploy at the end.

set -euo pipefail

STUDY="${1:?Usage: bin/publish_report_release.sh <domain>/<set_name>[/<run>]}"
REPO_ROOT="$(cd "$(dirname "${0}")/.." && pwd)"
cd "$REPO_ROOT"  # so `gh`'s repo auto-detection (git remote) works regardless of caller's cwd

IFS=/ read -r DOMAIN SET_NAME RUN EXTRA <<< "$STUDY"
if [ -z "$DOMAIN" ] || [ -z "$SET_NAME" ] || [ -n "${EXTRA:-}" ]; then
    echo "ERROR: expected <domain>/<set_name>[/<run>], got '$STUDY'" >&2
    exit 1
fi
DOCS_DIR="$REPO_ROOT/docs/$STUDY"

REPORT_FILES=()
if [ -n "$RUN" ]; then
    for f in figures figures_pdf archive island_synteny.html assembly_quality.html; do
        [ -e "$DOCS_DIR/$f" ] && REPORT_FILES+=("$f")
    done
else
    for f in novelties.html core.html losses.html summary.pdf; do
        [ -f "$DOCS_DIR/$f" ] && REPORT_FILES+=("$f")
    done
fi

if [ ${#REPORT_FILES[@]} -eq 0 ]; then
    echo "== no report files found at $DOCS_DIR -- nothing to publish ==" >&2
    exit 0
fi

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

if [ -n "$RUN" ]; then
    TAG="reports-${DOMAIN}-${SET_NAME}--${RUN}"
    printf '{"domain": "%s", "set": "%s", "run": "%s"}\n' "$DOMAIN" "$SET_NAME" "$RUN" > "$WORKDIR/manifest.json"
    TARBALL="$WORKDIR/${SET_NAME}--${RUN}-reports.tar.gz"
else
    TAG="reports-${DOMAIN}-${SET_NAME}"
    printf '{"domain": "%s", "set": "%s"}\n' "$DOMAIN" "$SET_NAME" > "$WORKDIR/manifest.json"
    TARBALL="$WORKDIR/${SET_NAME}-reports.tar.gz"
fi
TAR_ARGS=(-czf "$TARBALL" -C "$WORKDIR" manifest.json -C "$DOCS_DIR" "${REPORT_FILES[@]}")
tar "${TAR_ARGS[@]}"

echo "== publishing $TAG (asset: $(basename "$TARBALL"), $(du -h "$TARBALL" | cut -f1), files: ${REPORT_FILES[*]}) ==" >&2
if gh release view "$TAG" >/dev/null 2>&1; then
    gh release upload "$TAG" "$TARBALL" --clobber
else
    gh release create "$TAG" "$TARBALL" \
        --title "$TAG" \
        --notes "Report pages for $STUDY -- data-only release asset (not a software release), see DESIGN.md Sec 8. Downloaded and merged into docs/ at Pages-deploy time; never git-committed."
fi

if [ -n "${NII_SKIP_DEPLOY_TRIGGER:-}" ]; then
    echo "== done: $TAG published (deploy trigger skipped, NII_SKIP_DEPLOY_TRIGGER set) ==" >&2
else
    echo "== triggering Pages deploy ==" >&2
    gh workflow run static.yml
    echo "== done: $TAG published, deploy triggered ==" >&2
fi
