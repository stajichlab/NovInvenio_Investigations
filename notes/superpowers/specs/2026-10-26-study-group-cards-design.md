# Study-group cards on the domain gallery: design

Status: revision 2 (2026-10-26). Revised after an Opus review of revision 1. Implemented in the code PR linked from PR #52.

Request: the 16 `yc_*` contrasts of `yeast_convergence` filled the Fungal gallery as 16 cards. They are one specialized study. Show them under one card.

Facts come from `bin/generate_docs.py`, `lib/site_pages.py`, `lib/study_runs.py` and `conf/domains.yaml` on `origin/main` (550866f). The reviewer checked section 1 against the code and found no false statement, with the corrections noted there.

## 1. Behavior before the change

- `find_studies(domain)` returns one dict per study. Entries are keyed by the `study:` value of `publish.yaml`, not by folder name. One study can come from several folders (`agaricomycetes` comes from three). A `studies/` folder needs `species.csv`.
- `study_folders()` reads one level of grouping folder (a folder with no `species.csv`). Today 4 folders in fungi qualify: `yeast_convergence` (16 children) and three with none.
- `render_domain_index()` draws one card per study, then one card per `pages:` entry of `conf/domains.yaml`.
- `render_top_level()` shows `len(domain_data[slug])` as the study count (33 for fungi).
- The comparison page was `docs/fungi/yeast_convergence/index.html`, written by `yc_comparison_page.py`, and also a `pages:` card.
- Each `yc_*` `publish.yaml` had a `related:` link to `yeast_convergence/index.html`. So do 32 published enrichment pages (`enrichment.html` and `enrichment_loss.html` of each contrast, line 46), written by `bin/generate_enrichment_page.py` from `related:`.
- Member pages live at `docs/fungi/yc_<name>/<run>/`.

## 2. Goal

1. The domain page shows one card, "yeast convergence", in place of the 16 contrast cards.
2. The card opens a group page with the comparison link and one card per contrast.
3. The mechanism is generic: another set of studies can use it with configuration only.

## 3. Decisions

### D1. Membership: derived, and enabled by configuration

A study's group is, in order:
1. the `group:` key of its `publish.yaml`, if the key is present (`group: null` opts out; the code tests `"group" in data`, because a missing key and `null` differ);
2. else the name of the grouping folder that holds its folder (`studies/<domain>/<group>/<folder>/`);
3. else none.

The group is resolved per study, not per folder. If a study's folders give different groups, `generate_docs.py` stops with an error that names them.

A group takes effect only when the domain lists it under `groups:` in `conf/domains.yaml`. An unlisted group prints a warning and the studies stay ungrouped. This keeps grouping from switching on by accident when someone adds a nested folder.

### D2. Group metadata: `conf/domains.yaml`

```yaml
groups:
  - slug: yeast_convergence
    title: yeast convergence      # default: slug with spaces
    member_label: contrasts       # default: studies
    desc: ...
    page: {title: ..., href: yeast_convergence/comparison.html, desc: ...}
```

The slug must pass `validate_slug()`.

### D3. Member URLs do not move

Member pages stay at `docs/fungi/yc_*/`. Only the gallery presentation changes. No redirects are needed.

Moving members to `docs/fungi/yeast_convergence/yc_*/` is not part of this change. The reviewer counted about 17 places that build `docs/<domain>/<study>/` paths: `static.yml` (`dest`), `generate_docs.py`, `site_pages.py`, `study_runs.py`, `generate_enrichment_page.py`, `study_pages.py`, `sync_reports.sh`, `sync_pangenome_report.py`, `lib/ni_pangenome.py`, `lib/pangenome_site.py`, both release publish scripts, `publish_all_studies.sh`, `yc_comparison_page.py`, `.gitignore` and the guard test.

One hard-rule risk if it is ever done: the `.gitignore` patterns for class-3 files (`novelties.html`, `alignments/`, `figures/`, `clinker/` and others) have a fixed depth (`docs/*/*/` and `docs/*/*/*/`). At a third level they would not match and release-only files could be committed. The guard test also checks only 4-part paths. Both must become depth-independent first.

### D4. Comparison page and group page

The group page takes `docs/<domain>/<group>/index.html`. The comparison page moves to `comparison.html` in the same folder. `yc_comparison_page.py` writes the new name. The group page links to it through the group's `page:` entry. The old address now shows the group page, so no stub is needed.

Link updates: `yc_batches.py` (the `PUBLISH_YAML` template), the 16 `publish.yaml` files (edited directly), `conf/domains.yaml`, and the 32 published enrichment pages.

`yc_batches.py render` is not used for the 16 files. It rewrites `publish.yaml` only when the file is missing or with `--force`, and `--force` also overwrites every `species.csv`, which holds hand edits.

### D5. Counts and labels

- The top-level count stays the number of member studies (33 for fungi). A group does not hide its members.
- The group card and group page say "N <member_label>".
- Role labels per contrast are not shown. `study_index.tsv` has a structured `role` column (`origin`, `pseudo`, `reversal`), so a later change can add an optional `role:` key to `publish.yaml`.

## 4. Change

### 4.1 `lib/study_runs.py`
`PublishTarget` gets `group` and `group_set`. `group` goes through `validate_slug()`.

### 4.2 `bin/generate_docs.py`
- `load_domain_groups()` reads and validates `groups:`.
- `find_studies()` resolves `group` per study (D1) and fails on a conflict.
- `split_groups()` decides what the domain page shows:
  - listed group with two or more members: one group card and one group page;
  - listed group with one member: the member's own card, plus the group's `page:` as a normal page card, with a warning;
  - listed group with no members: the group's `page:` as a page card, with a warning;
  - unlisted group: ungrouped, with a warning.
- A group slug equal to a study name stops with an error.
- `check_group_target()` refuses to overwrite `docs/<domain>/<group>/index.html` unless it is a generated group page (it carries `<!-- ni-group-page -->`). The one exception is the old comparison page, and only after the group's `page:` file exists at its new address.
- `warn_stale_group_pages()` warns about a marked group page whose group is no longer configured. The generator never deletes files.
- Card order: ungrouped studies by name, then group cards, then page cards.

### 4.3 `lib/site_pages.py`
- The study card is a function `_study_card(s, prefix)`. The prefix (`../` on the group page) goes in front of the main link, the run-list link and every extra link.
- `render_domain_index(..., groups=[...])` draws group cards (`.card.group`, a thicker left border; no new colors, so dark mode needs no new rule).
- `render_group_index()` writes the group page: breadcrumb to the domain, the `page:` card, then member cards.

### 4.4 `studies/fungi/yeast_convergence/` and `conf/domains.yaml`
As in D4.

### 4.5 Tests
`tests/test_study_runs.py`: group card and page; relative links on the group page (all three kinds); unlisted group stays ungrouped; `group: null` and override; one-member group; conflicting folders; slug equal to a study name; refusal to overwrite a foreign index, and the allowed overwrite after the move; stale group page warning; `group` key parsing. `test_gallery_lists_nested_folders` still passes, because its config lists no group.
`studies/fungi/yeast_convergence/tests`: each committed `publish.yaml` equals the `PUBLISH_YAML` template; no committed `publish.yaml` links to `yeast_convergence/index.html`.

### 4.6 Rollout
1. PR 2 (code): items 4.1 to 4.5 and the 16 `publish.yaml` edits. No `docs/` files.
2. After merge, one docs-only commit to `main` with: the regenerated domain page and `docs/index.html`; `yeast_convergence/comparison.html` (from `yc_comparison_page.py`); the group page `yeast_convergence/index.html` (same commit, so the comparison never disappears); the 32 enrichment pages (rerun `generate_enrichment_page.py` for each contrast and direction; needs the enrichment outputs on disk).
3. Check the live site: domain page, group page, comparison page, one member page, one enrichment page link.

## 5. Risks and limits

- Links outside the repo to `yeast_convergence/index.html` now show the group page. The comparison is one click away. Not searched for outside the repo.
- Member run lists and enrichment pages still link back to the domain page, not the group page. Accepted for this change.
- The group title is shown with spaces ("yeast convergence"), matching how study cards show names.
- The rewrite of the 32 enrichment pages depends on enrichment outputs that I have not confirmed are on disk.

## 6. Open items

- A `role:` key for per-contrast labels.
- Moving member URLs (D3).
- A second group, for example the pangenome runs of one genus.
