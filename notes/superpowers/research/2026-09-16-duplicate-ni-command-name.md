# `bin/ni` exists in both NovInvenio and NovInvenio_Investigations (2026-09-16)

Investigation only, no changes made. Triggered by: "I see a ~/projects/NII/bin/ni
and ~/projects/NI/bin/ni - I think there should only be one implementation."

## Finding: these are NOT duplicate/out-of-sync copies of the same tool

Diffed both files directly (`~/projects/NI/bin/ni` vs `~/projects/NII/bin/ni`,
i.e. `NovInvenio/bin/ni` vs `NovInvenio_Investigations/bin/ni` — `~/projects/NI`
and `~/projects/NII` are just symlinks to the full repo names). They are
completely different scripts that happen to share the short command name `ni`:

| | `NovInvenio/bin/ni` | `NovInvenio_Investigations/bin/ni` |
|---|---|---|
| Subcommands | `init` only | `resolve`, `fetch`, `discover`, `run` |
| Purpose | Scaffold a **brand-new** NII-style repo from a reference checkout | Operate **within an existing** NII checkout: fill in `species.csv` accessions, fetch data, discover NCBI genome populations, launch the pipeline |
| History | 1 commit, `711a36f` (2026-09-09), closes NovInvenio issue #76 | 7 commits, 2026-09-11 → 2026-09-15, design doc `notes/superpowers/specs/2026-09-11-ni-dataset-resolution-design.md` |
| Why it lives where it does | Must exist somewhere OTHER than the not-yet-created target repo — can't live inside a repo it's responsible for creating | Needs `lib/ni_resolve.py`, `lib/ni_discover.py`, `bin/build_study_config.py`, `bin/run_study.sh` — all NII-local code |

Each placement is individually well-justified and was a deliberate design choice
(NovInvenio's commit message and `.living/decisions.md:917` both explain the
`init`-must-live-outside-the-target-repo reasoning explicitly). **This is not
implementation drift or an accidental fork** — nothing needs reconciling at the
code level.

## But it IS a real naming/documentation problem worth fixing

1. **Same short name, disjoint subcommands, no cross-reference.** `ni init`
   (NovInvenio) and `ni resolve|fetch|discover|run` (NII) are two CLIs with zero
   subcommand overlap. Neither repo's docs mention the other tool at all —
   confirmed by grepping both READMEs/DESIGN.md for `bin/ni` and `scaffold`.
   NII's own README (line 17) claims `bin/ni` is "the one entrypoint" for the
   study-building workflow, which isn't quite true once the earlier
   repo-bootstrapping stage is counted — a real, if minor, documentation gap.
2. **Not currently colliding in practice.** Checked `~/.bashrc`, `~/.bash_profile`,
   `~/.profile` — neither `NI/bin` nor `NII/bin` is on `$PATH`, so today each is
   only invoked via an explicit path inside its own repo (`bin/ni resolve` run
   from inside a NII checkout, `bin/ni init` run from inside a NovInvenio
   checkout). No live ambiguity right now.
3. **Latent risk.** If either `bin/` directory is ever added to `$PATH` (plausible
   — both are literally named to invite a global `ni` command), whichever comes
   first on `$PATH` silently shadows the other, and `ni resolve` typed while the
   wrong one is first would fail with a confusing "invalid choice" rather than
   "command not found." Also a real onboarding-confusion risk independent of
   `$PATH`: someone reading both repos' docs would reasonably expect a single
   `ni` tool with `init` alongside `resolve`/`fetch`/`discover`/`run`, not two
   separate binaries.

## Recommendation (not implemented — investigation only)

Keep the functional split (real architectural constraint: the scaffolding tool
genuinely cannot live inside a repo it hasn't created yet), but resolve the name
collision and the missing cross-reference:

1. **Rename NovInvenio's tool** to something that doesn't collide, e.g.
   `bin/ni-init` or `bin/scaffold_ni_repo.py` — it's the less-frequently-invoked
   of the two (run once per new deploy repo, vs. NII's `ni` used continuously
   for every study), so it should yield the short name. Keep NII's `bin/ni` as
   the canonical `ni` — it already has the deeper investment (7 commits, its own
   design doc, real subcommands people use every day) and NII's README already
   calls it "the one entrypoint."
2. **Add a one-line cross-reference in both repos' docs**: NovInvenio's
   scaffolding tool's docstring/README should say "for building/running studies
   inside an already-scaffolded NII repo, see that repo's `bin/ni`"; NII's
   README should have a "before this: scaffolding a brand-new repo" note
   pointing at NovInvenio's tool, so "the one entrypoint" framing becomes
   accurate (or gets qualified) instead of silently omitting the earlier stage.
3. **Do not attempt to merge them into one script/one repo.** That would require
   either vendoring NII's study-operation code into NovInvenio (defeats the
   "always copy from a live, working reference" design principle the
   scaffolding tool is built around) or vendoring the scaffolding tool's
   copy-manifest logic into NII (meaningless there, since NII already exists by
   definition once you're running its own `bin/ni`). The repo split itself is
   correct; only the naming and documentation need to change.

No code changes made in this investigation, per instruction to recommend only.

## Update: implemented (2026-09-16, same day, on user approval)

User confirmed the renaming recommendation ("yes that renaming seems
justified") and approved implementing it now despite the coccidioides
session having parallel uncommitted work on the same `NovInvenio` branch
(`pangenome-profiling-module`). Verified before touching anything that
`bin/ni` itself had no uncommitted edits from that session, then made only
the two renames (`git mv bin/ni bin/ni-init`, `git mv tests/test_ni.py
tests/test_ni_init.py`) plus their docstring cross-reference updates,
staged and committed only those two files (`git status` confirmed the other
session's 9 modified files were untouched before and after), and left the
commit unpushed per that branch's established local-commit convention.
`NovInvenio` commit: `adc1f10`. Ran `pixi run pytest tests/test_ni_init.py`
first (9/9 passed) to confirm the rename didn't break anything.

Also updated this repo's own `README.md` with a one-line cross-reference to
`NovInvenio`'s `bin/ni-init` for the earlier bootstrapping stage, so the
"one entrypoint" claim about this repo's `bin/ni` is now accurately scoped.
