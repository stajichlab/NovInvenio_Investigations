# Study-group cards on the domain gallery: design

Status: DRAFT for review. No code is written. Date: 2026-10-26.
Request: the 16 `yc_*` contrasts of `yeast_convergence` fill the Fungal gallery as 16 separate cards. They are one specialized study. Show them under one card.

This document names every file the change touches and every decision that needs a reviewer. Facts below come from `bin/generate_docs.py`, `lib/site_pages.py`, `lib/study_runs.py` and `conf/domains.yaml` on `origin/main` (550866f, 2026-10-24). I did not run the generator for this draft.

## 1. Current behavior

- `bin/generate_docs.py` `find_studies(domain)` returns one dict per study name. Sources: `docs/<domain>/<study>/` (published, has `*/run.json`) and `studies/<domain>/<folder>/publish.yaml` (not yet published).
- `study_folders()` already treats a folder without `species.csv` as a grouping folder and reads its children (`studies/fungi/yeast_convergence/<contrast>/`). The grouping is not kept after this step. Each child becomes its own study.
- `lib/site_pages.py` `render_domain_index()` draws one `.card` per study dict, then one card per hand-made `pages:` entry from `conf/domains.yaml`.
- `render_top_level()` shows `n_studies` per domain: `len(domain_data[slug])`.
- The comparison page is at `docs/fungi/yeast_convergence/index.html`. It is written by `studies/fungi/yeast_convergence/bin/yc_comparison_page.py`. It is also a hand-made `pages:` card in `conf/domains.yaml`.
- Each `yc_*` `publish.yaml` has a `related:` link to `yeast_convergence/index.html` (written by `yc_batches.py`, line 68).
- Member URLs today: `docs/fungi/yc_<name>/<run>/...`.

## 2. Goal

1. The domain page shows one card, "yeast_convergence", in place of the 16 contrast cards.
2. The card opens a group page. The group page shows the comparison page link and one card per contrast.
3. The mechanism is generic. Another set of related studies can use it without code changes.

## 3. Decisions for the reviewer

### D1. Where the membership comes from

Option A (recommended): derive it. A study folder under a grouping folder belongs to the group named by the parent folder (`yeast_convergence`). `study_folders()` already finds these. Add an optional `group:` key to `publish.yaml` that overrides the derived value (for a study that lives elsewhere, or to opt out with `group: null`).

Option B: an explicit `group:` key in each `publish.yaml`. 16 files carry the key. `yc_batches.py` writes it.

Option C: list members in `conf/domains.yaml`. One place to edit, but the list drifts from the study folders.

Recommendation: A. It needs no change to the 16 files and cannot drift. Cost: the grouping depends on folder layout. The override key covers exceptions.

### D2. Where the group's title and description come from

`conf/domains.yaml` gets an optional `groups:` list per domain: `slug`, `title`, `desc`, optional `page` (see D4). A group with no entry gets a default (title = slug, desc = empty) and a warning.

### D3. Do member URLs move?

Option 1 (recommended for the first PR): no move. Only the gallery presentation changes. Member pages stay at `docs/fungi/yc_*/`. No redirects. No change to the sync scripts, `publish_report_release.sh` asset paths, release-asset names, or relative links.

Option 2: move members to `docs/fungi/yeast_convergence/yc_*/`. The URL shows the grouping. You said URL changes are acceptable. The cost: every script that builds a `docs/<domain>/<study>/` path (`rebuild_study_pages`, `sync_reports.sh`, the publish scripts, `yc_*` generators, enrichment page generator, relative `../` links, `related:` hrefs, release-asset naming) must learn a two-level path, and 16 old URLs need redirect stubs. I have not counted the call sites. A reviewer should decide whether the URL benefit is worth it.

Recommendation: Option 1 now. Option 2 as a separate, later change if wanted.

### D4. The comparison page and the group page

Today `docs/fungi/yeast_convergence/index.html` is the comparison page. The group page needs the same address, so one of them must move.

Option I (recommended): the generic group page is written to `docs/<domain>/<group>/index.html`. The comparison page moves to `docs/<domain>/<group>/comparison.html`. `yc_comparison_page.py` writes the new name. The group page links to it as a card (title and href from the group's `page:` entry in `conf/domains.yaml`). The old address now shows the group page, so no stub is needed there. The 8+ `related:` links and the `pages:` card change to `yeast_convergence/comparison.html`.

Option II: keep the comparison page at `index.html` and have `generate_docs.py` add a member list into it. This edits a hand-built page by string insertion. It is fragile. Not recommended.

### D5. What the group card counts

- Domain page card meta: `N contrasts` (the group's member count) and, when members carry a role, a short breakdown.
- `render_top_level()` `n_studies`: count groups as one study or count members? Recommendation: count members, so the number matches the research output (it stays 16 for this group). The label stays "studies". A reviewer may prefer groups.
- Role breakdown (9 origin, 5 control, 2 reversal): the member `publish.yaml` `notes:` text is free text. I would not parse it. Options: omit the breakdown, or add `role:` to `publish.yaml`. Recommendation: omit in the first PR; show `N contrasts`.

## 4. Proposed change

### 4.1 `conf/domains.yaml`

```yaml
- slug: fungi
  groups:
    - slug: yeast_convergence
      title: yeast_convergence
      desc: >-
        Independent yeast origins in fungi: 9 origin contrasts, 5 filamentous
        controls and 2 reversals.
      page:
        title: comparison across contrasts
        href: yeast_convergence/comparison.html
```

The existing `pages:` card for the comparison page is removed (the group page carries it).

### 4.2 `bin/generate_docs.py`

- `find_studies()` adds `"group": <slug or None>` to each study dict. Source: `publish.yaml` `group:` if present, else the parent folder name when the parent has no `species.csv`, else `None`. Published studies (from `docs/`) are matched to their group through the `studies/` folder with the same study name. A published study with no `studies/` folder stays ungrouped.
- New `load_domain_groups()` reads `groups:`.
- `main()`: split the study list into ungrouped studies and groups. Call `render_domain_index()` with both. For each group with at least one member, write `docs/<domain>/<group>/index.html` with `render_group_index()`.
- A group slug that equals an existing study name is an error (the address is taken). The script stops with a message.
- A group with exactly one member renders as that member's own card (no group page), with a warning.

### 4.3 `lib/site_pages.py`

- `render_domain_index(..., groups=[...])`: draw one group card per group, between study cards and `pages:` cards. The card links to `<group>/index.html`. Meta: `N contrasts`.
- New `render_group_index(domain_name, group, members, site_name)`: title, desc, breadcrumb (`← <domain>`), the `page:` card, then one card per member. Member card markup is the existing study-card markup, factored out and reused so the two pages cannot differ. Member hrefs are written relative to the group page (`../yc_x/<run>/report.html`).
- `render_top_level()`: no change except the `n_studies` value (D5).

### 4.4 `studies/fungi/yeast_convergence/`

- `bin/yc_comparison_page.py`: write `comparison.html` instead of `index.html`. Fix any self-links.
- `bin/yc_batches.py` line 68 and the 16 `publish.yaml` `related:` links: `yeast_convergence/index.html` → `yeast_convergence/comparison.html`. Regenerate the files with `yc_batches.py render`.
- Add a test or check that no `related:` href points at a missing file.

### 4.5 Tests

- `tests/` for `site_pages`: group card appears; member cards do not appear on the domain page; group page lists all members; a one-member group falls back; unknown `group:` key warns; slug collision fails.
- `generate_docs.py` on a small fixture tree: derived group, `group: null` opt-out, published-only study.

### 4.6 Rollout

1. PR 1: code, tests, `conf/domains.yaml`, comparison rename, regenerated `docs/` index pages. The `docs/` index pages are class 2 (small, do not grow with data).
2. Docs-only follow-up commit to `main` for the regenerated gallery, as in the NII commit workflow.
3. Check the live site: domain page, group page, comparison page, one member page, the old comparison address.

## 5. Risks

- External links to `yeast_convergence/index.html` (for example in issue comments and PR #43) will show the group page, not the comparison. The comparison link is one click away. I have not searched outside the repo for other links.
- Without a built `docs/` tree on the reviewer's machine, the layout cannot be checked by eye. The PR should include rendered screenshots or the generated HTML.
- The derived-group rule (D1 option A) would group any future nested study folder automatically. A reviewer should confirm that is wanted for `studies/fungi/Afumigatus_pangenome/`-style folders (they have `species.csv` at the top, so they are not grouping folders today).

## 6. Not decided here

- Role labels per contrast (origin, control, reversal) on member cards.
- Moving member URLs (D3, Option 2).
- A second group, for example the pangenome runs of one genus.

## 7. Questions for the reviewing model

1. Is deriving the group from the parent folder (D1-A) safe, or is an explicit key better?
2. Is the comparison-page rename (D4-I) the least risky way to free `index.html`?
3. Does counting members in `n_studies` (D5) fit the top-level page wording?
4. Which call sites does D3 Option 2 touch that I did not list?
5. Is any step above missing a test?
