# Growth-form literature search: plan for an agent set

Purpose: code the growth form of the 106 yeast_convergence species from published text, with a source
for every coding (issue #44). This plan records what a pilot did on 10 species (2026-10-25), what
worked, what failed, and how to run it better. A new agent set can follow it without the chat history.

Status words used here: **tested** = I ran it and saw it work. **used** = the pilot agents used it.
**not tested** = proposed only.

## 1. The question and the output

For each species, answer two questions from the literature:

1. Growth-form class: `yeast`, `filamentous` or `dimorphic`.
2. If dimorphic: which phase is the normal form in standard rich medium at the usual culture
   temperature: `yeast`, `hyphal` or `unclear`.

Never infer from the genus. If no source states it for this species, the answer is `unsourced`.
Dimorphism often depends on conditions (temperature, CO2, medium, host). Record the condition. Do not
reduce a condition-dependent result to one word without saying so.

Output is an append-only evidence table, one row per quote (section 6), plus one adjudicated row per
species. Nothing goes to `config_support/traits/` until a person approves a sample. Repo rules apply
(`CLAUDE.md`): tracked data needs a provenance record; no guessing; Simplified Technical English.

## 2. What the pilot did

Ten species: *Thamnidium elegans*, *Yarrowia lipolytica*, *Eremothecium cymbalariae* (the three open
conflicts in #44), *Rhizopus delemar* and *Mucor circinelloides* (known answers, tested by Tahiri et al.
2026), and *Aureobasidium pullulans*, *Exophiala dermatitidis*, *Taphrina deformans*, *Protomyces
inouyei*, *Neolecta irregularis*.

Five reader agents (Sonnet), two species each, in parallel. They did not see our codings. Each ran 3 to 4
searches per species, read at most 8 papers, and returned records with a verbatim quote of at most
40 words. A second step (done by hand) re-checked 17 quotes or key phrases against the sources (8 abstract quotes,
9 full-text key phrases; see `analysis/growth_form_pilot/verification.tsv`). I checked no other record.

Results: 17 of 17 checked items matched the sources (full-text items are key phrases, not whole quotes). 2 of 2 known answers were right. Coverage: 7 of 10
species had a direct species-level statement on the class; a dominant phase was supported for 2 of 6
dimorphic species. Cost: about 516,000 tokens in total (about 52,000 per species), 43 to 70 seconds per
agent. The five reports, the prompt and the checks are in `analysis/growth_form_pilot/`; issue #44 holds the species list.

## 3. Tools and resources

### 3.1 Used or tested in the pilot

| Resource | How | Status | Notes |
|---|---|---|---|
| PubMed search (MCP `plugin_bio-research_pubmed`: `search_articles`) | fielded query, returns PMIDs | used | Load tools first with `ToolSearch` (`select:<tool name>`). Rate limits hit in the pilot. |
| PubMed metadata (`get_article_metadata`) | PMIDs in, abstract and IDs out | used, tested | Reliable for abstract quotes. Cite PubMed and give DOI links in any user-facing text (tool notice requires it). |
| PMC full text (`get_full_text_article`) | PMCIDs in | used | Returned an empty body for several papers. Strips italics, so species names vanish from sentences. Do not rely on it. |
| NCBI E-utilities `efetch` for PMC XML | `curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=<numeric PMCID>&rettype=xml"` | tested (6 papers) | Full text with italic tags kept. Strip tags, collapse whitespace, then match strings. This is the best full-text route. |
| NCBI E-utilities `esearch`, db `pmc` | `...esearch.fcgi?db=pmc&retmode=json&retmax=N&term="<species>"[All Fields] AND (dimorphic OR yeast OR hyphae)` | tested | Searches full text, not only abstracts. *T. elegans* gave 105 hits. |
| NCBI E-utilities `esearch`, db `pubmed` | `term="<species>"[Title/Abstract] AND dimorphism` | tested | *Y. lipolytica* gave 144 hits. |
| Local NCBI taxonomy | `/srv/projects/db/ncbi/preformatted/20260128/names.dmp` (tax id, name, unique name, class) | tested | Gives synonyms. *Y. lipolytica* (4952): *Candida lipolytica*, *Endomycopsis lipolytica*, *Mycotorula lipolytica*. *E. dermatitidis* (5970): *Wangiella dermatitidis*, *Hormiscium dermatitidis*. Not every taxon has synonyms (*T. elegans* has none). |
| `species.csv` `Taxon_ID` column | key into names.dmp | tested | The 16 `yc_*/species.csv` files hold the IDs. |
| WebSearch / WebFetch | web pages | used | WebFetch returns a summary, not text. Never quote from it. |
| Species list and our codings | `studies/fungi/yeast_convergence/batches/*.yaml` (`species:` map), `config_support/traits/traits.csv` (NovInvenio repo) | used | Keep these hidden from readers; compare after. |

### 3.2 Seen, not used or failed

- Europe PMC REST: HTTP 500 in the pilot. Retry later; not tested since.
- bioRxiv and Consensus MCP tools exist in this environment: not tested for this task.
- Paywalled books and monographs (Kurtzman et al., *The Yeasts: A Taxonomic Study*; Benny 1992 on
  *Thamnidium*): not readable by agents. A person with library access can read the species entries and
  type the codings in.

### 3.3 Proposed, not tested

- PMC Open Access bulk files (FTP) as a local corpus. With the full-text XML on disk, search by `grep`
  and no network limit.
- NCBI API key (raises the E-utilities rate limit from 3 to 10 requests per second; stated in NCBI
  documentation, not verified here).
- PubTator or Europe PMC annotations for species mentions; OpenAlex or Crossref for citation metadata.
- SLURM `scrontab` for scheduling (see section 8). Not checked whether it is enabled here.

## 4. Pipeline (what to build)

Stages A to C are plain scripts. Only stage D needs a model.

**A. Name table (script).** For each species, read its `Taxon_ID`, list all names from `names.dmp`
(scientific, synonym, authority forms), and add names from the batch YAML and from the literature
(for example *Mucor lusitanicus* for *M. circinelloides*; also *Rhizopus oryzae* for *R. delemar*,
*Ashbya gossypii* for the *Eremothecium* clade). Output: `names.tsv` (species, name, source).

**B. Candidate retrieval (script).** For each name, run `esearch` on pubmed (Title/Abstract) and pmc
(Text Word) in three modes: `dimorph` (the name AND dimorphism terms such as `dimorphic`,
`yeast-to-hypha`, `pseudohyphal growth`), `growth` (the name AND the generic growth terms `yeast`,
`hyphae`, `mycelium`, `budding`, ...), and `name` (the name alone). Priority order: dimorph, growth, name.
The `dimorph` mode is needed because for well-studied species the generic growth query matches
thousands of papers, and the relevance-sorted top 200 can miss the dimorphism papers (the pilot's
*Yarrowia* source, Liang 2017, was missed this way). Save the PMIDs and PMCIDs, the query, the date and
the hit count. Keep every query (section 6).

**B2. Article selection (script, in `fetch`).** Do not take the first N rows by search rank. Rank ties
among the many name queries of one species fill a small cap with arbitrary papers (the pilot's
*Mucor* source, Lee 2013, was rank 2 and still cut). Merge the PubMed and PMC rows of one article by
PMCID, then order by (priority, number of queries that returned it, best rank). Use caps large enough
that selection is not the limiting step: the sentence stage costs about 1,200 tokens per species, so
60 articles (40 full texts) per species is cheap. `gf_lit.py fetch --selection score` does this;
`--selection rank` is the earlier order.

For species with very large literatures, add a second step: `--selection two-step`. It takes the top
`--screen-articles` (300) candidates in the score order, fetches only their PubMed abstracts, and scores
each title and abstract (dimorphism terms +3 each, growth terms +1 each up to 5, species name in the title
+5, in the abstract +2). It then fills half of the slots (`--base-share`, 0.5) from the score order and the
rest from the abstract-screen order. Neither order alone is enough: in the 19-row gold test the score order
missed a dimorphism paper that the abstract screen found (Liang 2017), and the abstract screen missed a
paper whose abstract has no dimorphism term (Tahiri 2026) that the score order found. The hybrid found both.

**C. Text fetch and sentence filter (script).**
1. For each PMCID, `efetch` the XML once. Cache it on `/bigdata`, compressed (`zstd`; CLAUDE.md prefers
   compressed text). Record the SHA-256 of the XML.
2. For abstract-only PMIDs, store the abstract from `get_article_metadata` or `efetch db=pubmed`.
3. Strip tags, split into sentences (keep the section title and table captions).
4. Keep sentences that contain a species name from stage A *and* a growth term, plus the one
   sentence before and after. Names are matched in three ways: the full name, the abbreviation
   (`G. species`) when the full name is in the same article, and a name written with both genera
   (`Wangiella (Exophiala) dermatitidis`, built from the species' synonyms).
   **Article-level rule** (`--article-level abstract`, the default): when the full name is in the
   article's title or abstract, the article is about the species, and its growth sentences that do not
   name the species ("the fungus ...", "these cells ...") are kept too, labelled `match_level article`,
   at most `--article-cap` (12) per article, sentences with a dimorphism term first. Without this rule a
   gold sentence such as "the switch to sustained filamentous growth" in a paper about *R. delemar* is
   lost. A genus-level sentence ("*Aureobasidium* is a genus of dimorphic yeast-like
   fungi") is kept and labelled genus-level.

**D. Classification (model).** Give the model only the kept sentences of one species, not the papers.
Ask for one JSON row per supporting sentence (section 6). The model decides: what is said, which
organism it is about (species, genus, other), whether it is tested, described or predicted, and the
condition. It copies the quote. It does not summarize.

**E. Scripted quote check (script, no model).** For each row, confirm that `quote` appears verbatim in
the cached text of that `pmcid`/`pmid` after whitespace normalization. Reject rows that fail. This
replaces the manual check done in the pilot and should run on 100% of rows.

**F. Adjudication (model, once per species).** Input: the verified rows for one species. Output: the
class, the dominant phase and a confidence, with the conflicts listed (no majority vote). Rules in
section 5.

**G. Review table (script).** Join with our current codings and `traits.csv`. Mark: agree, differs,
new, unsourced. Write `growth_form_evidence.tsv` and `growth_form_adjudicated.tsv`. A person reviews
every `differs` row and a random 10% of `agree` rows before anything is promoted.

## 5. Decision rules for the adjudicator

- Evidence ranks: tested at species level > described at species level (review, table, database) >
  genus-level statement > predicted. A genus-level statement alone gives `unsourced` for the species
  plus a note, not a coding.
- Abstract-only evidence is allowed for the class, but flag it. Prefer one full-text source.
- Dominant phase = the form in the species' standard laboratory culture (rich medium, usual
  temperature). If only a condition-dependent result exists, write `unclear` and give the conditions.
  The pilot showed that sources often name the conditions for a switch and not the default.
- If an organism name in a source differs (synonym, strain of another species), the row stays with the
  name the source uses, mapped through `names.tsv`. A different species (for example *A. melanigenum*
  for *A. pullulans*) is not evidence.
- Absence of a dimorphism report is not evidence of filamentous growth. Say `unsourced` unless a source
  states the growth form.
- Two sources that disagree: report both quotes. A person decides.

## 6. Record formats

Evidence row (one per quote, TSV or JSONL):
`species, taxon_id, name_used_in_source, pmid, pmcid, doi, read_level (full_text|abstract),
evidence_type (tested|described|predicted), organism_level (species|genus|other), condition, location,
quote (max 40 words), supports (yeast|filamentous|dimorphic|other), phase_stated (yeast|hyphal|none),
source_text_sha256, retrieved_utc, query_id, model, prompt_version`

Query log (one row per search): `species, name, database, query string, date, hit_count`.
Adjudicated row: `species, class, dominant_phase, confidence, n_records, n_full_text, conflicts,
note`. Provenance record for any file that is later tracked: use `lib/provenance.py`
(`build_record`), with the script commit, the NCBI taxonomy file version (`20260128`), the E-utilities
query date and the SHA-256 of every cached source text.

Store short quotes only. Do not copy full text or whole abstracts into the repo.

## 7. Failures in the pilot and the fix

| Failure | Fix |
|---|---|
| PMC text tool returned an empty body; italics (species names) stripped | stage C uses `efetch` XML and keeps the raw tags until the filter has run |
| Agents ran fewer searches than planned; synonyms skipped (*Wangiella*, *Candida lipolytica*, *M. lusitanicus*, *Pullularia*) | stage A builds names; stage B runs every name; the agent does not decide |
| Rate limits and HTTP 500 | stage B is a script with retries, back-off and a cache; one request every 0.35 s without an API key |
| Known answer (Tahiri 2026) not found by search | stage B also searches the species name with each of the three Mucorales papers' typical terms; keep a gold list (section 9) and test recall on it |
| Dominant phase rarely stated | accept `unclear`; add the culture-default rule; ask for a person to supply book tables |
| Genus-level statements read as species evidence (*Protomyces*, *Aureobasidium*) | `organism_level` is a column; the adjudicator may not use genus rows for a species coding |
| Manual quote check | stage E is a script |

## 8. Simplifications

1. Move retrieval, text fetch, sentence filtering and quote verification out of the model (stages A, B,
   C, E). The pilot spent its tokens on search and reading whole papers. A sentence filter should cut
   the text that the model reads by a large factor. I have not measured the factor; measure it on the
   10 pilot species first.
2. Use one cheap-model call per species (stage D) and one stronger-model call per species (stage F)
   instead of an open-ended agent per species pair.
3. Cache everything (XML, abstracts, query results). A later run changes only the new items.
4. Drop the species whose coding is already sourced and unchanged.
5. Keep an `unsourced` list as a first-class result. Do not spend tokens trying to force an answer.

## 9. Calibration and quality control

- **Gold set.** Tahiri et al. 2026 (Nat Commun 17:2148, PMC12957363): *M. lusitanicus* and
  *M. circinelloides* dimorphic; *R. microsporus* and *R. delemar* not dimorphic (tested). Add species
  from Nagy et al. 2014 (Nat Commun 5:4471) and its supplement, once read. Report recall and precision
  on the gold set at every run. A change to a prompt or a model needs a gold-set rerun.
- **Metrics** per run: quote-verification rate (target 100% after stage E), fraction of species with a
  species-level source, fraction with a stated phase, number of conflicts, number of `differs` rows.
- **Human review**: every `differs` row, every conflict, and a random 10% of `agree` rows.
- **Blind comparison**: readers never see our codings.

## 10. If tokens were not limited

- Run stage D with two different models independently on every species and let a third adjudicate.
  Disagreement between the two readers is a flag for review.
- Read the full text of every hit, not only the filtered sentences, for species that stay `unsourced`
  after the first pass. Include figure legends and supplementary tables.
- Build the local Open Access corpus (PMC OA bulk files) and grep it for every name. Then recall no
  longer depends on how a query is written.
- Follow citations: for each species, take the references of the best review papers and search those.
- Add bioRxiv and Europe PMC (when its API works) as further sources.
- Ask for tables explicitly: parse the XML table elements for rows that name a species and a morphology
  column (yeast, hypha, mycelium, dimorphic). Tables in reviews are the likeliest source of a clean answer.
- Spend human time where agents cannot: books and monographs. Give the person a worksheet with the
  species, the book and the page to look up.

## 11. Small, slow, periodic background process

A design that needs few tokens and no one watching:

1. **State.** One SQLite or TSV ledger: per species and name, `last_searched_utc`, `n_hits`,
   `status` (new, pending_model, verified, adjudicated, needs_review, reviewed). Append only; never
   overwrite a reviewed row.
2. **Schedule.** One SLURM job, for example weekly, that does: (a) stage B for names whose last search
   is older than N days, limited to papers added since `last_searched_utc` (the E-utilities `edat` date
   filter; not tested), (b) stages C and E for new papers, (c) stage D only for species with new
   sentences, (d) stage F only when a species has new verified rows, (e) a short report of changes.
   Use `scrontab` if SLURM enables it (not checked), or a self-resubmitting `sbatch --begin=now+7days`
   job. Use a wall time near the 1 to 1.5 hours of real work that CLAUDE.md recommends. A first full
   pass over 106 species is a few waves, for example 20 species per run.
3. **Rate and cost control.** One request every 0.35 s (or faster with an API key). Cap model calls per
   run. A cheap model for stage D, a stronger one only for stage F. Stop the run and write the ledger if
   the cap is hit.
4. **Safety.** The job never writes outside its own output folder, never edits `config_support/`, and
   never opens a PR. It writes `needs_review` rows. A person (or an agent that a person starts) promotes
   rows after review.
5. **Alerts.** One line per run in a log (counts of new papers, new verified rows, conflicts, failures).
   Failures (HTTP errors, empty text) do not mark a species as searched.
6. **Rebuilds.** Keep the cached XML and the query log under `/bigdata/.../growth_form/cache/` with
   SHA-256 values, so a later run can re-verify old quotes without network access.

Expected first-pass cost with the simplified design: unmeasured. The pilot's unsimplified figure was
about 52,000 tokens per species; the simplified design should be lower, but I have not shown it.

## 12. Prompts (starting point)

Reader prompt used in the pilot (adapt for stage D by replacing the search instructions with "Here are
the sentences"): see the pilot task text in the session of 2026-10-25. Keep these rules in any prompt:

- Quote verbatim, maximum 40 words. If you cannot quote it, do not report it.
- Say whether the sentence is about this species, the genus, or another organism.
- Evidence type: tested, described or predicted.
- Report conflicts. Do not choose between sources.
- No answer for a species means `unsourced`. Do not use the genus, and do not guess.
- Never invent a citation. Give a PMID, DOI or PMCID for every row.

## 13. First steps for the next agent set

1. Build `names.tsv` for the 106 species from `names.dmp` (stage A).
2. Run stages B, C and E for the 10 pilot species. Compare to the pilot results above. This measures the
   token saving and shows whether the scripts find the same sources.
3. Run the gold set (Tahiri 2026 species) and report recall.
4. If recall and the quote-verification rate are acceptable, run the 6 low-BUSCO species and the three
   #44 conflicts, then the rest in waves.
5. Hand the unsourced list and the book worksheet to a person.

## 14. Working with more than one agent

Several agents, possibly from different vendors and with different costs, can share this work. They do
not need to talk to each other. They need one file contract, a script that checks every output, and a
rule for who touches which file. A message to a peer session is not reliable for this: a session list
shows only sessions started from the same tool, and an agent from another tool is not on it.

### 14.1 Roles

| Role | Does | Model |
|---|---|---|
| Retrieval | stages A to C (`gf_lit.py species, names, search, fetch, sentences`): scripts, no model | none |
| Extractor | stage D: reads `sentences/<species>.jsonl`, writes evidence rows | cheap model; two independent passes with different models |
| Verifier | `gf_lit.py verify`: checks each quote against the cached source text | none (script) |
| Adjudicator | stage F: one decision per species, from the verified rows of all passes | strong model; only for species where passes disagree, or where evidence is conflicting or genus-level |
| Reviewer | a person: reads every `differs` and `conflict` row | human |

An agent may hold several roles. A cheap extractor is safe because the verifier rejects any quote that is
not verbatim in the source. A model cannot make an unverifiable claim pass.

### 14.2 File contract

All shared files live under `studies/fungi/yeast_convergence/growth_form/`. Agents write only there. No
agent writes to `config_support/traits/`, a batch YAML, or any file outside `growth_form/`.

| Path | Writer | Content |
|---|---|---|
| `sentences/<species_slug>.jsonl` | retrieval | input for extractors (gitignored; holds source text) |
| `prompts/<prompt_version>.txt` | whoever changes a prompt | the exact prompt; never edited after use, a new version is a new file |
| `evidence/<pass_tag>/<species_slug>.jsonl` | one extractor | evidence rows (section 6 of this plan). `<pass_tag>` names the model and prompt, for example `haiku-v1`. One directory per pass, so passes stay independent. |
| `evidence/<pass_tag>/<species_slug>.verified.tsv`, `.rejected.tsv` | the verifier | output of `gf_lit.py verify` |
| `adjudicated/<species_slug>.json` | adjudicator | class, dominant phase, confidence, conflicts, the rows used |
| `work/<stage>/<species_slug>.claim/` | any agent | claim marker (14.3) |
| `cost_log.tsv` | any agent | one row per call batch: pass_tag, species, model, input tokens, output tokens, wall seconds, date |
| `calibration.tsv` | whoever runs a calibration | per model: quote-verified rate, gold sources found, class agreement (14.5) |
| `STATUS.md` (study folder) | each agent, own section only | what is done, blocked or unclear |

Every evidence row carries `model` and `prompt_version`. An agent never edits another pass's files.

### 14.3 Claiming work

1. List the species with no `evidence/<pass_tag>/<species_slug>.jsonl`.
2. Claim one with `mkdir work/D/<species_slug>.claim`. If `mkdir` fails, another agent has it: take the
   next species. (`mkdir` is atomic on `/bigdata`.)
3. Write `owner`, `model` and the UTC time into `work/D/<species_slug>.claim/info.txt`.
4. Write the output file, then remove the claim directory.
5. A claim older than 6 hours with no output file is stale. Any agent may remove it and must add a line
   to its STATUS section.

### 14.4 The check that gates everything

`gf_lit.py verify --evidence <file>` is the only way a row becomes `verified`. It accepts a row when the
quote is at most 40 words and appears verbatim in the cached source (after whitespace and punctuation
normalization). It works for any model, so the cheap and the strong model are held to the same rule.
Rows that fail go to `.rejected.tsv`. Report the rejected fraction per pass in `calibration.tsv`: a
high fraction means the model paraphrases, and that model should not be used for extraction.

### 14.4a Comparing passes

`gf_stage_d.py compare` is written (section 16.1). It reads the verified rows of all passes for a species and
report: the class and phase each pass supports, whether they agree, and the quotes behind any
disagreement. The adjudicator reads only species that are flagged. Build it before running two passes.

### 14.5 Choosing the cheap model by measurement

1. Pick the gold-set species (`gold_sources.tsv`, 19 rows from the pilot's verified items plus Tahiri 2026).
2. Run stage D on these species with each candidate model and the same prompt version.
3. For each model record in `calibration.tsv`: share of rows that verify; share of gold sources that give
   a row; agreement of class and phase with the verified pilot result; tokens used.
4. Use the cheapest model whose numbers match the strong model on these four. If none does, use the strong
   model for extraction and keep the cheap one as the second pass only.

This sample is small (about 10 species). State that in any conclusion drawn from it.

### 14.6 What an outside agent needs to join

- This plan and `STATUS.md`.
- `bin/gf_lit.py` and its `--help` (run with `pixi run python bin/gf_lit.py ...`).
- Read access to `growth_form/sentences/` and the cache, write access to `growth_form/evidence/<pass_tag>/`
  and `growth_form/work/`.
- An API key for its own model, if it uses one. Never write a key to a file in the repo.
- No network access to NCBI is needed for stage D: the sentences are local.

### 14.7 Rules that protect the data

- Source text stays out of git: `cache/` and `sentences/` are gitignored. Evidence rows hold quotes of at
  most 40 words.
- Nothing is promoted to `config_support/traits/` until a person has reviewed a sample.
- A change to a prompt or a model is a new `pass_tag`. Results of different tags are never mixed in one file.

## 15. Gold sets and what they show (2026-10-25)

Two gold sets test retrieval and selection (`gf_lit.py recall --gold <file> --recall-out <file>`):

- `growth_form/gold_sources.tsv` (19 rows, 12 sources): from the pilot's verified items. A method that
  copies the pilot's search scores well here, so it is a regression test, not a measure of coverage.
- `growth_form/gold_sources_independent.tsv` (82 rows, 29 species): built by `bin/gf_gold_tables.py` from
  the tables and body text of 336 review and comparative papers found by seven searches that are not the
  pipeline's. Rows have no phrase: they test whether the paper is retrieved and selected. Each row has a
  `strength`: `table`, `about` (species in the title or abstract), or `passing` (a sentence only).

Result of the hybrid selection (caps 60 articles, 40 full texts) on the independent set:

| strength | rows | retrieved | selected |
|---|---|---|---|
| table | 13 | 12 | 4 (31%) |
| about | 40 | 39 | 28 (70%) |
| passing | 29 | 23 | 6 (21%) |

Excluding *C. albicans*, *C. neoformans* and *S. cerevisiae*: about 28 rows, 26 selected (93%); table 8 rows,
4 selected, all four misses being one review table that lists many species once. For the three
species: table 0 of 5, about 2 of 12 selected. With caps of 300 articles (200 full texts) and a screen of
800 for these three: table 3 of 5, about 7 of 12. So the cap limits heavily studied species, and a larger
cap does not remove the problem. A species-specific `review[pt]` query mode is a candidate fix; not tried.
Do not read `passing` recall as a failure: a review that names 30 species once is not a source to select for each.

### 15.1 Review query and the per-species measure

A `review` query mode (the name AND dimorphism terms AND `review[pt]`, PubMed, same priority as `dimorph`)
returns a small, complete set even for heavily studied species (162 reviews for *C. albicans* against 1,677
dimorphism papers). Run on the whole pipeline (hybrid selection, article-level sentences), it changed little
per paper (independent gold, the three species: about 2 of 12 to 3 of 12 selected; table rows still 0 of 5;
pilot gold unchanged at 17, 17, 17 and 13 of 19). It did rescue *S. cerevisiae* by one review.

The per-paper measure understates what matters. The missed papers of the three species are mostly
tangential research papers (an antifungal compound that blocks hyphal growth, an enzyme paper), at rank 20
to 170 with one agreeing query; the authoritative reviews of dimorphism were selected. The useful measure is
per species: is at least one gold source selected? On the independent set that is 17 of 20 species
(16 before the review query). The three left are *A. pullulans* and *H. werneckii*, whose only gold row is a
multi-species table in a review about smut fungi, and *T. rubrum*, whose only gold row is a paper on the
mating-type locus. None of the three is a clear selection failure. Report this per-species number next to
the per-paper one.

## 16. Stage D calibration (2026-10-26): Haiku and Sonnet on ten species

Ten species with known answers (the pilot's), candidate sentences from the hybrid selection with article-level
sentences, prompt `prompts/extract-v1.txt`, one call per chunk of 200 sentences (21 chunks), the evidence checked by
`gf_lit.py verify` and scored by `gf_stage_d.py score`. Opus was not run.

| Pass | Rows | Verified | Class match | Gold sources found | Agent tokens (total) |
|---|---|---|---|---|---|
| haiku-v1 (unchunked input; flawed baseline) | 368 | 94% | 10 of 10 | 10 of 13 | about 1.15 million |
| haiku-v2 (chunked) | 789 | 96% | 9 of 10 | 12 of 13 | about 1.29 million |
| sonnet-v2 (chunked) | 543 | 94.5% | 9 of 10 | 10 of 13 | about 0.98 million |

Token totals are the agents' own reports (input, output and tool reads together); I did not convert them to cost.

What it shows:

- **Quotes are not the problem.** 367 of 368 first-pass quotes were verbatim in the source. Nearly all rejections were
  two gate rules, now fixed: the quote had to hold the full species name even when the sentence used "M. circinelloides",
  and article-level and genus-level sentences were rejected outright.
- **Reading protocol matters.** One unchunked run read 350 of 775 lines. Input is now cut into chunks.
- **Haiku found more and as accurately as Sonnet, on this sample.** It gave more rows and found more gold sources
  (12 of 13 against 10 of 13); both had a verified rate near 95%. Its first pass also gave different counts from its
  second on the same sentences (for example 5 against 10 rows for *Protomyces*), so one pass is not a stable answer.
- **The crude class call was the weak part, not extraction.** My first rule let one `predicted` row (Tahiri 2026: *T.
  elegans* behaves like the dimorphic *Mucor*) flip *Thamnidium* to "dimorphic" against ten descriptive "filamentous"
  rows, and one loose "the yeast E. cymbalariae" row flip *Eremothecium*. The rule now ignores `predicted` rows and needs
  a minority form to reach 20%. **I changed the rule after seeing these results**, so the 9 of 10 match is not an
  independent test; both models gain, haiku-v1 is unaffected. The adjudicator (stage F) must do this job properly.
- ***Protomyces* is the one miss.** Both models call it "yeast" from a sentence that reports *P. inouyei* cell length
  in YPD culture; my reference (the pilot: genus-level only) may be the thing that is wrong.
- **Caps still bind.** Chunks that return 55 to 60 objects (the prompt allows 60 per chunk) were common for *Taphrina*,
  *Exophiala* and *Mucor*, so the number of rows per species is bounded by the prompt, not by the evidence.
- **Row precision is not measured.** I read 12 verified rows by hand for haiku-v1: about 10 clearly stated the growth form
  of the species, 2 were marginal. A verified row means the quote is real and sits in a sentence that names the species,
  not that the model's label is right.

Open: two independent passes with different models would show disagreements for the adjudicator; `gf_lit.py compare`
is not written; a larger test than ten species is needed before any model is chosen for the 106.

### 16.1 `compare`: where two passes disagree

`gf_stage_d.py compare --out <dir> --pass-tags haiku-v2 sonnet-v2` reads the verified rows of each pass and reports, per
species, the class and phase of each pass, whether they agree, the sentences each found, and `needs_adjudication`:
class_disagree, phase_disagree, thin_evidence (a pass with one or two counted rows), no_species_rows_in. A species
whose rows show both yeast and hyphal growth is not flagged: that is what a dimorphic species looks like (an earlier
version flagged all five dimorphic species; it was changed). It also writes the rows an adjudicator should read.

On the ten species: Haiku and Sonnet agree on the class for 9; three species need a look: *Neolecta* and *Protomyces*
(thin evidence: one row or none) and *Taphrina* (phase). Haiku against its own first pass: class agrees for 8; the two
that flip are the same thin species. So thin evidence is where a single model run is unstable.

**The models differ on phase.** For *Taphrina*, Sonnet filled the `phase_stated` field on 12 of 32 counted rows (mostly
"yeast": the yeast state grows on culture media and filamentous growth appears in the host), Haiku on 1 of 56. Phase is
the weakest field of the existing codings (the batch YAML codes *Taphrina* `dominant_phase: hyphal`), and a class pass
on Haiku would leave it mostly blank. Options: ask for phase in a separate prompt, or use Sonnet for it. Not tested.
The quotes are verified as real; my reading of them is not reviewed.

### 16.2 Phase pass (a second step on verified rows)

Extraction on Haiku left `phase_stated` almost empty, so a second prompt (`prompts/phase-v1.txt`) works on the
rows that already passed `verify` (quotes only, about 1,600 to 5,600 tokens per species). Per quote it asks the
context (culture, induced, host, natural, unspecified) and the phase; per species it asks for
`dominant_phase_in_culture` (yeast, hyphal, unclear, not_stated), up to five basis sentence ids, and a confidence.
`gf_stage_d.py phase-input` writes the input and `phase-compare` compares passes and prints the cited quotes.
Tested on the five dimorphic species with Haiku and Sonnet:

| Species | Haiku | Sonnet | Batch YAML `dominant_phase` |
|---|---|---|---|
| *Mucor circinelloides* | hyphal (medium) | hyphal (medium) | hyphal |
| *Taphrina deformans* | yeast (medium) | yeast (high) | hyphal |
| *Aureobasidium pullulans* | unclear (low) | unclear (low) | yeast |
| *Yarrowia lipolytica* | unclear (low) | unclear (medium) | yeast |
| *Exophiala dermatitidis* | unclear (low) | yeast (medium) | yeast |

- The two models agree on 4 of 5 and all cited sentence ids were in the input. The disagreement (*Exophiala*) is about
  the same quotes: "in most rich media ... budding yeast", "exists predominantly as a yeast form", and "exclusively as a
  filamentous morphotype in liquid media". Sonnet reads the rich-medium statements as the culture default; Haiku reads the
  conflict as unclear. The prompt defines the question as a rich medium, which favours Sonnet's reading; I did not review
  the source papers.
- Only *Mucor* has a firm reference (hyphal in aerobic culture, tested by Lee et al. 2013), and both models match it.
  The others show agreement between models and quotes a person can read, not accuracy.
- If the quotes are right, three of the five existing `dominant_phase` codings are not supported: *Taphrina* is yeast in
  culture (hyphae in the host), and *Aureobasidium* and *Yarrowia* depend on the medium. This is from five species and has
  not been applied to any file.
- Agent tokens: about 0.37 million (Haiku) and 0.34 million (Sonnet) for five species. Most of this is the fixed cost of
  starting an agent (each used 60,000 to 80,000 tokens on a 2,000 to 6,000 token input), so one agent per species is
  wasteful; batching species into one call would cut it. Not tried.

