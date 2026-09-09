#!/usr/bin/bash
# Package a study's TBLASTN alignment shards (docs/<domain>/<set>/alignments/,
# .../loss_alignments/ -- produced by nf_NovInvenio's BUILD_ALIGNMENT_SHARDS,
# issue #73, and relocated into the two-tier gallery path by bin/sync_reports.sh,
# issue #5) as a GitHub Release asset, so they reach the published site without
# ever being git-committed (see DESIGN.md Sec 8 and issue #3's .gitignore rule).
#
# Usage: bin/publish_alignment_release.sh <domain>/<set_name>
#
# Tag: alignments-<domain>-<set> (one release per study, overwritten in place on
# every rerun via --clobber -- matches how docs/<domain>/<set>/report.html itself
# is already always-overwritten with no run history, per nf_NovInvenio's own
# convention).
#
# The tarball carries a manifest.json ({"domain": ..., "set": ...}) at its root
# alongside alignments/ and loss_alignments/, so the deploy workflow
# (.github/workflows/static.yml, issue #2) can recover the destination path by
# reading the manifest instead of parsing the tag name -- robust regardless of
# how a domain/set name is spelled.
#
# --clobber does NOT re-fire the release:published webhook (that only fires once,
# at creation) -- this script explicitly triggers the Pages deploy workflow
# afterward via `gh workflow run` (workflow_dispatch), since a rerun would
# otherwise silently never redeploy.
#
# This is a deliberate, standalone publish step -- NOT wired into
# bin/run_study.sh's automatic post-run chain, since it creates a real GitHub
# Release and triggers a real Pages deploy (external, visible side effects).
# Run it by hand once a study's reports are ready to actually go live.

set -euo pipefail

STUDY="${1:?Usage: bin/publish_alignment_release.sh <domain>/<set_name>}"
REPO_ROOT="$(cd "$(dirname "${0}")/.." && pwd)"
cd "$REPO_ROOT"  # so `gh`'s repo auto-detection (git remote) works regardless of caller's cwd

SET_NAME="$(basename "$STUDY")"
DOMAIN="$(dirname "$STUDY")"
DOCS_DIR="$REPO_ROOT/docs/$STUDY"

if [ ! -d "$DOCS_DIR/alignments" ] && [ ! -d "$DOCS_DIR/loss_alignments" ]; then
    echo "== no alignment shards found at $DOCS_DIR -- nothing to publish ==" >&2
    exit 0
fi

TAG="alignments-${DOMAIN}-${SET_NAME}"
WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

printf '{"domain": "%s", "set": "%s"}\n' "$DOMAIN" "$SET_NAME" > "$WORKDIR/manifest.json"

TARBALL="$WORKDIR/${SET_NAME}-alignments.tar.gz"
TAR_ARGS=(-czf "$TARBALL" -C "$WORKDIR" manifest.json)
[ -d "$DOCS_DIR/alignments" ] && TAR_ARGS+=(-C "$DOCS_DIR" alignments)
[ -d "$DOCS_DIR/loss_alignments" ] && TAR_ARGS+=(-C "$DOCS_DIR" loss_alignments)
tar "${TAR_ARGS[@]}"

echo "== publishing $TAG (asset: $(basename "$TARBALL"), $(du -h "$TARBALL" | cut -f1)) ==" >&2
if gh release view "$TAG" >/dev/null 2>&1; then
    gh release upload "$TAG" "$TARBALL" --clobber
else
    gh release create "$TAG" "$TARBALL" \
        --title "$TAG" \
        --notes "TBLASTN alignment shards for $STUDY -- data-only release asset (not a software release), see DESIGN.md Sec 8. Downloaded and merged into docs/ at Pages-deploy time; never git-committed."
fi

echo "== triggering Pages deploy ==" >&2
gh workflow run static.yml

echo "== done: $TAG published, deploy triggered ==" >&2
