# Pangenome visualization tooling: reuse existing tools vs. build a custom report

Date: 2026-09-14

## Context

`studies/fungi/Afumigatus_pangenome` has a real, already-computed 295-strain (293
ingroup + 2 outgroup), 47,983-gene-family presence/absence matrix from mmseqs2
tier-1 clustering at ~90% identity (see
`notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md`) — not a
Roary/Panaroo/Prokka-derived pangenome, so there is no GFF3-derived
Roary-native input already sitting on disk. Component 8 of that design classified
1,259,432 FDR-significant co-occurring gene-family pairs into 522,933 "trans"
(physically unlinked, candidate interaction/co-evolution pairs), 2,416
"starship_explained", 5,242 "unexplained_physical", plus smaller categories.
Component 9 speced a custom `pangenome.html` report (canvas heatmap + a
`trans`-pair network panel) built as a structural twin of NovInvenio's existing
`novelties.html`/`core.html`/`losses.html` reports. This note investigates whether
an existing pangenome-visualization tool could be reused/adapted instead of, or
alongside, that custom build, and how to visualize the 522,933-edge co-occurrence
network specifically. Every claim below is sourced against the tool's own
docs/repo/paper — see inline citations; anywhere a primary source did not confirm
a claim, that is stated explicitly rather than asserted.

Maintenance was checked directly against each tool's GitHub API `pushed_at`
timestamp on 2026-09-14 (not a changelog page, which can lag actual commit
activity):

| Tool | Last repo push (2026-09-14 check) | Verdict |
|---|---|---|
| Roary (`sanger-pathogens/Roary`) | 2025-11-03 | repo receives pushes, but see caveat below — this is `pushed_at`, not confirmed to be substantive commits from the original maintainers |
| Panaroo (`gtonkinhill/panaroo`) | 2026-07-02 | actively maintained |
| PPanGGOLiN (`labgem/PPanGGOLiN`) | 2026-09-14 (same day) | actively maintained |
| Phandango (`jameshadfield/phandango`) | 2024-02-26 | >2.5 years stale — fails the "real commit activity in the last 2-3 years" bar |
| PIRATE (`SionBayliss/PIRATE`) | 2022-07-08 | >4 years stale — fails the bar |
| Coinfinder (`fwhelan/coinfinder`) | 2023-10-20 | ~3 years stale, borderline; included below because it is the one directly pangenome-specific co-occurrence-network precedent found |

## 1. Existing pangenome visualization tools

### Roary / `roary_plots.py`

- **Input format**: `gene_presence_absence.csv`. Metadata columns (14, matching
  the script's `--skipped-columns` default) per Roary's own docs page
  (https://sanger-pathogens.github.io/Roary/): Gene name, Non-unique gene name,
  Functional annotation, Number of isolates in cluster, Number of sequences in
  cluster, Average sequences per isolate, Genome fragment, Order within fragment,
  Accessory fragment, Accessory order with fragment, QC comments, Min/Max/Avg
  sequence length (nt) — followed by one column per genome. Roary also emits a
  simpler `gene_presence_absence.Rtab`, a plain tab-delimited binary matrix where,
  per the same docs page, "a 1 indicates the gene is present in the sample, a 0
  indicates it is absent." The plotting script itself
  (https://github.com/sanger-pathogens/Roary/blob/master/contrib/roary_plots/roary_plots.py)
  reads the CSV directly: `Gene` as index, skips the first 14 metadata columns,
  and treats a non-NaN cell in a genome column as "present" via regex matching (it
  does not require the Rtab file).
  - **Conversion effort for our data**: our TSV is already family × strain with
    present/absent/genome_only values (three-state, not Roary's implicit
    present/absent). Producing a Roary-shaped CSV is a small script: synthesize
    the 14 metadata columns (most can be placeholder/derived — e.g. "Number of
    isolates," "Min/Max/Avg length" are computable from our own data;
    "annotation" can carry our `gene_name`/Pfam/SwissProt columns), map
    `genome_only` to either "present" or a separate flag (the format has no native
    third state), and write one column per strain with either a gene ID or blank.
    Low-to-moderate effort (a single conversion script, a few hours), with the
    caveat that the three-state distinction would need to be dropped or bolted on
    as a side annotation since Roary's format is binary.
- **Figures produced out of the box**
  (https://github.com/sanger-pathogens/Roary/blob/master/contrib/roary_plots/roary_plots.py):
  (1) pangenome frequency histogram, (2) pangenome presence/absence matrix
  plotted against a phylogenetic tree (`pangenome_matrix.<fmt>`), (3) pangenome
  composition pie chart (core/soft-core/shell/cloud). No rarefaction/accumulation
  curve and no interactive component — everything is a static image
  (png/tiff/pdf/svg).
- **Scale**: the script has no hardcoded genome/gene ceiling, but per a direct
  read of the script's own logic, matplotlib font-size scaling is done by
  decrementing per genome, which degrades to unreadable labels well before 300
  genomes; the matrix figure's fixed 17×10 inch canvas does not rescale for
  ~48,000 rows either. It would very likely run without crashing on 295 strains ×
  47,983 families but produce an illegible `pangenome_matrix` figure at this row
  count — this is an inference from reading the script's plotting logic, not a
  documented limit stated by the authors.

### Panaroo

- **Input/native operation**: Panaroo is a full clustering pipeline (takes GFF3 +
  annotated genomes as input), not a presence/absence-matrix visualizer taking
  arbitrary input. Its own docs
  (https://github.com/gtonkinhill/panaroo/blob/master/docs/gettingstarted/output.md)
  state its `gene_presence_absence.csv` is "the same format as that given by
  Roary" — so the conversion path above would work equally for feeding
  Panaroo-format consumers, but Panaroo itself is not something you point at an
  existing matrix; you'd be reusing its *output format*, not running the tool.
- **Figures**: per its own docs and the companion tool `panstripe`
  (https://github.com/gtonkinhill/panstripe, described in Panaroo's own repo
  listing as "post processing of bacterial pangenome gene presence/absence
  matrices"), Panaroo's own visualization story is (a) viewing its pangenome graph
  in Cytoscape (external desktop tool, not embeddable in a self-contained HTML
  report) and (b) presence/absence heatmaps. Panstripe adds regression-style
  plots for gene gain/loss dynamics against a phylogeny, not a frequency
  histogram or pie chart.
- **Reusable for our data?** Only via Roary-format-compatible ingestion (same
  caveats/effort as Roary above); Panaroo itself is a clustering pipeline we are
  not using (we already have mmseqs2 tier-1 clusters), so there is nothing to
  "run" here beyond adopting its CSV schema.

### PPanGGOLiN

- **Input**: builds its own pangenome graph from raw genome/annotation input; I
  did not find, in the fetched docs pages
  (https://ppanggolin.readthedocs.io/en/latest/user/PangenomeAnalyses/pangenomeFigures.html),
  any statement that it accepts an externally-computed presence/absence matrix as
  a starting point — its own docs page did not address this question at all
  (explicitly unverified rather than assumed absent).
- **Figures**: (1) an interactive HTML "U-shaped plot" — family count (y) vs.
  number of genomes each family appears in (x), with a dashed line marking the
  soft-core threshold (default 95%) and a Heaps'-law fit; (2) an interactive HTML
  "tile plot" — a genuine heatmap (families × genomes, colored present/absent,
  genomes clustered by Jaccard distance, families ordered by partition/
  prevalence). Both are described in
  https://ppanggolin.readthedocs.io/en/latest/user/PangenomeAnalyses/pangenomeFigures.html.
- **Scale — the single most directly relevant finding in this section**: the same
  PPanGGOLiN docs page states explicit numeric ceilings: **the tile/U-shaped plot
  is not drawn at all when using the standard workflow subcommands once a
  pangenome exceeds 60,000 gene families**; **above 32,767 families, only the
  'shell' and 'persistent' partitions are drawn, cloud is dropped**; and **above
  65,535 families or genomes, the docs state no browser can render the tile plot
  at all** (a stated browser limitation, not a PPanGGOLiN algorithmic one). Our
  47,983 families sits **under** the 60,000 draw-cutoff but **over** the 32,767
  full-partition cutoff — i.e., PPanGGOLiN's own authors, at a very comparable
  scale to ours, found it necessary to start dropping data from their own
  browser-rendered heatmap. This is strong, directly-comparable, primary-sourced
  evidence that a plain per-cell heatmap render (whether ours or a borrowed one)
  needs real row/column reduction at this scale — reinforcing, from an independent
  tool's real engineering experience, the design spec's already-planned bitset
  packing + column virtualization for `pangenome.html` (component 9,
  "should-consider S1/S2").
- **Reusable for our data?** Not directly — no confirmed external-matrix ingestion
  path found, and it is designed around its own graph-building pipeline. Its
  documented scale ceilings are more useful to us as **engineering evidence for
  our own report's design** than as a tool to adopt.

### Phandango

- Interactive viewer (phylogeny + metadata + Roary-format pan-genome data,
  BRATNextGen, PLINK, SEER) per its own repo
  (https://github.com/jameshadfield/phandango) and its paper (Hadfield et al.,
  bioRxiv, https://www.biorxiv.org/content/10.1101/119545v1.full.pdf). **Excluded
  from further evaluation**: last repository push 2024-02-26, over 2.5 years
  before this check (2026-09-14) — fails the "real commit activity in the last
  2-3 years" bar the task set, even though its GitHub Issues page still shows
  recent user-filed issues (activity from users, not from maintainers).

### PIRATE

- Per its own README (https://github.com/SionBayliss/PIRATE), ships optional
  R-based plotting (`-r|--rplots` flag, ggplot2/dplyr/ggtree/phangorn-based
  `PIRATE_plots.pdf`) and documents Roary-format-compatible output for use with
  Phandango, Scoary, and PanX. **Excluded from further evaluation**: last
  repository push 2022-07-08, over 4 years before this check — clearly fails the
  "real commit activity in the last 2-3 years" bar.

### Coinfinder (Whelan, Rusilowicz, McInerney 2020)

Not a presence/absence-matrix visualizer in the Roary/ppanggolin sense — a
lineage-aware co-occurrence *tester*, already methodologically relevant to
component 4/8 of the design spec (see section 2 below for its network-output
details). Included here only to note it is **not** an alternative to
Roary/ppanggolin/Panaroo for the core/shell/cloud matrix-heatmap question; it
answers a different part of our pipeline (co-occurrence testing and, separately,
network export), not the presence/absence heatmap.

### Summary for section 1

No tool surveyed accepts our exact three-state family×strain TSV natively; all
would need a conversion script into Roary's CSV/Rtab schema (moderate, one-time
effort, with the caveat that the three-state present/genome_only/absent
distinction has no native slot in that schema). Of the actively-maintained
options, **PPanGGOLiN is the most informative primary source for scale
engineering** (concrete family-count ceilings directly comparable to our 47,983),
but it does not appear to ingest external matrices per its own docs. **Roary's own
`roary_plots.py`** is the most mechanically reusable — a real, if unpolished,
static-figure generator that would run on a converted matrix — but its plotting
code was written assuming far smaller genome/gene counts than ours (inferred from
reading the script, not a documented limit) and produces only static images, not
an interactive report. None of the surveyed tools is a drop-in replacement for the
speced `pangenome.html`; at best, Roary's static-figure path is a legitimate
low-effort supplement (see section 3).

## 2. Visualizing large-scale gene co-occurrence networks

(522,933 "trans" pairwise associations — a measured number from the real 293-
strain run — far too many to render as a single graph.)

### 2.1 Community/module detection methods

- **Louvain** (Blondel, Guillaume, Lambiotte, Lefebvre 2008, "Fast unfolding of
  communities in large networks," J. Stat. Mech. — cited here via its role as
  Leiden's predecessor in the Leiden paper, not independently re-verified against
  the 2008 paper directly in this pass). Widely implemented
  (`networkx.algorithms.community.louvain_communities`, `python-louvain`,
  igraph's `community_multilevel`), but has a documented structural defect.
- **Leiden** (Traag, Waltman, van Eck 2019, "From Louvain to Leiden: guaranteeing
  well-connected communities," *Scientific Reports* 9:5233,
  https://www.nature.com/articles/s41598-019-41695-z, open preprint
  https://arxiv.org/abs/1810.08473). The paper's own abstract states Louvain "may
  yield arbitrarily badly connected communities. In the worst case, communities
  may even be disconnected, especially when running the algorithm iteratively,"
  with an empirical finding of "up to 25% of the communities are badly connected
  and up to 16% are disconnected." The authors "prove that the Leiden algorithm
  yields communities that are guaranteed to be connected," and state "the Leiden
  algorithm runs faster than the Louvain algorithm" while producing better
  partitions. This is a direct, citable reason to prefer Leiden over Louvain here,
  not a fashion choice.
  - Reference implementation: `leidenalg`
    (https://github.com/vtraag/leidenalg, docs
    https://leidenalg.readthedocs.io/en/stable/intro.html), C++ core exposed to
    Python via `python-igraph`. Its own docs state: "The implementation scales
    well, and can be run on graphs of millions of nodes (as long as they can fit
    in memory)." Our graph (≤47,983 nodes, 522,933 edges) is well inside this
    stated envelope. Supports modularity, CPM, RBConfiguration, Significance,
    Surprise as optimizable objectives.
- **Other methods** (Infomap, label propagation) are standard alternatives in the
  general community-detection literature, but I did not independently verify
  Infomap's own current-maintenance or scale claims against its primary docs in
  this pass — flagged as unverified rather than asserted. Leiden's explicit,
  primary-sourced superiority claim over Louvain plus its stated million-node
  scalability makes it the better-justified default here without needing to reach
  for these alternatives.

### 2.2 Python graph libraries at ~500K-edge scale

Our graph is ≤47,983 nodes / 522,933 edges.

- **networkx**: pure Python, dict-of-dicts internal structure. NetworkX's own docs
  (https://networkx.org/documentation/stable/reference/introduction.html) make no
  quantitative scale claim at all — only that MultiGraph's flexibility causes
  "some degradation in performance, though usually not significant," with no
  numeric threshold given anywhere on that page. NetworkX 3.x's own backend-dispatch
  docs (https://networkx.org/documentation/stable/reference/backends.html)
  describe optional GPU (`nx-cugraph`, https://github.com/rapidsai/nx-cugraph) and
  joblib-parallel (`nx-parallel`, https://github.com/networkx/nx-parallel)
  backends that can accelerate specific algorithms — but these are opt-in extra
  dependencies, not networkx's baseline. graph-tool's own performance page (next
  bullet) states NetworkX's slowdown "can scale up quickly" once graphs reach
  "hundreds of thousands, or millions of vertices/edges" — this is graph-tool's
  characterization of networkx, not networkx's own claim about itself, so treat it
  as a competitor's framing, not a primary NetworkX statement.
- **python-igraph**: C-core library (~1,100 user-facing functions per the
  library's own paper, Csárdi et al., "igraph enables fast and robust network
  analysis across programming languages," https://arxiv.org/pdf/2311.10260), with
  64-bit integer support throughout. That paper reports constructing a
  3.2-billion-edge Bethe lattice in 4 minutes on a single core, and ~30 billion
  edges / 100M nodes fitting in 1TB RAM. Our 522,933-edge graph is trivially
  within igraph's own demonstrated envelope, and it composes directly with
  `leidenalg` for Leiden clustering.
- **graph-tool**: C++/Boost-backed, OpenMP-parallel for supported algorithms. Its
  own performance page (https://graph-tool.skewed.de/performance.html) benchmarks
  a real graph of **N=39,796 vertices, E=301,498 edges** (the PGP web-of-trust
  network) — the single most scale-comparable primary-source data point found in
  this research pass (our graph is ~1.6-1.7x that benchmark's edge count). The
  page states graph-tool is "the most steady performer" across the tasks it
  benchmarks and beats NetworkX by "40 to 250 times." It does not publish
  benchmarks at larger scale than the ~40K/300K graph, so behavior at millions of
  edges is not confirmed from graph-tool's own docs — only that it clearly handles
  a scale close to ours well.
- **scikit-network**: SciPy-CSR-sparse-matrix-backed, Cython-compiled, per its own
  JMLR paper (Bonald et al. 2020, https://jmlr.org/papers/v21/20-412.html /
  https://arxiv.org/abs/2009.07660) and docs
  (https://scikit-network.readthedocs.io/). States it targets large sparse
  graphs via SciPy CSR + Cython + parallelism, but I did not find an explicit
  numeric node/edge benchmark figure in the docs pages surfaced during this
  search — weaker sourcing than igraph's or graph-tool's numeric claims, flagged
  explicitly.
- **leidenalg**: see 2.1 — million-node scale stated directly in its own docs.

**Library recommendation basis**: igraph + leidenalg is the best-evidenced pairing
— both have primary-source scale claims that comfortably exceed our 522,933-edge
graph, and they compose directly. Plain networkx is fine for loading/inspecting a
*post-clustering, collapsed* graph (a few hundred community "super-nodes"), but
its own docs give no assurance for operating on the full 522,933-edge graph
directly, and graph-tool's competitor-framed comparison (not a primary networkx
claim) puts it far behind at this scale.

### 2.3 Published precedent: pangenome-specific vs. adjacent microbiome/WGCNA literature

- **Coinfinder** (Whelan, Rusilowicz, McInerney 2020, "Coinfinder: detecting
  significant associations and dissociations in pangenomes," *Microbial Genomics*
  6(3),
  https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000338,
  preprint https://www.biorxiv.org/content/10.1101/859371v1.full, code
  https://github.com/fwhelan/coinfinder) is genuine, directly on-topic
  **pangenome-specific** precedent, not microbiome/OTU literature. It is
  lineage-aware (uses a user-supplied phylogenetic tree to avoid calling
  co-occurrence that is really shared clade membership — directly analogous to
  our design's stratified-permutation-null approach), tests pangenome gene-family
  presence/absence for statistically significant association *and* dissociation,
  and was validated on 534 real *Streptococcus pneumoniae* genomes.
  - **On visualization specifically**, per a direct read of the tool's own GitHub
    README: Coinfinder does **not** render its own network plot. It outputs
    `.gephi` (GEXF) network files plus text/pdf summaries and presence/absence
    heatmaps against the phylogeny; the documented workflow states "the
    association network displayed in part A was made by inputting the coinfinder
    output .gephi file into the Gephi software" — visualization is explicitly
    delegated to Gephi, an external desktop application, not rendered by
    Coinfinder itself. No scale/node-count claim was found in Coinfinder's own
    docs (neither a stated maximum nor a benchmark number).
  - This is close conceptual precedent for our co-occurrence *testing*
    methodology (already reflected in the design spec's stratified permutation
    null) but is **not** a solved answer for large-scale, self-contained,
    `file://`-openable network visualization — it hands that problem to an
    external, non-web, non-self-contained tool, which is directly incompatible
    with `pangenome.html`'s existing no-network-fetch/self-contained-file
    constraint (per `CLAUDE.md`'s "Constraints to preserve when editing the page").
- No published, pangenome-specific paper was found describing a methodology for
  *visualizing* (as opposed to statistically testing) hundreds of thousands of
  gene-family co-occurrence edges. Searches for "pangenome accessory gene
  co-occurrence network," "gene presence absence network bacteria," and
  "Roary/Panaroo co-occurrence network" surfaced Coinfinder as the only
  pangenome-specific co-occurrence-network software found.
- **Explicitly distinguished, not blurred**: the well-established
  collapse-to-modules / cap-displayed-edges / fixed-layout methodology for very
  large networks comes from the **general network-science and adjacent
  microbiome co-occurrence-network literature** (SparCC, SPIEC-EASI, WGCNA-style
  module analysis of taxon/OTU or expression-correlation networks) — this
  research pass did not verify those specific tools' own docs (out of scope for
  the pangenome-focused question asked), and they are named here only to mark
  that this is where the general methodology originates, not to claim they were
  primary-sourced in this pass. Our planned approach (Leiden-collapse the
  522,933-edge network + cap displayed edges to a fixed top-N + mandatory
  edge-list table twin, as already speced in component 9's "should-consider S5")
  is methodologically sound and precedented **in that adjacent literature**, but
  is not itself validated against a published pangenome-accessory-gene case study
  at comparable scale — treat it as adapted general practice, not as reproducing
  an established pangenome-specific protocol.

## 3. Recommendation

### The two real numbers that drive this decision

- **47,983 families × 295 strains** for the presence/absence matrix. This is
  *directly* comparable to a documented scale problem another actively-maintained
  tool (PPanGGOLiN) has already hit and solved for: at 47,983 families we are
  under their 60,000-family "don't even draw a plain HTML tile plot" cutoff but
  over their 32,767-family "drop the cloud partition" cutoff
  (https://ppanggolin.readthedocs.io/en/latest/user/PangenomeAnalyses/pangenomeFigures.html).
  That is independent, primary-sourced confirmation — from a different team
  solving the same rendering problem — that the design spec's already-planned
  bitset-packed payload + column virtualization for `pangenome.html` (component 9,
  S1/S2) is not over-engineering; it is the same order of mitigation another real
  tool needed at this scale.
- **522,933 "trans" edges.** No tool surveyed in section 1 renders a network at
  all (Roary/Panaroo/PPanGGOLiN don't have a co-occurrence-network feature;
  Coinfinder does, but delegates rendering entirely to external Gephi, which is
  incompatible with our self-contained-report constraint). Whatever we build here
  is *de novo* work regardless of whether we reuse an existing matrix-heatmap
  tool for section 1's problem.

### Can any existing tool's plotting code be adapted directly?

- **Not a full substitute, but a legitimate cheap supplement**: converting our
  matrix to Roary's `gene_presence_absence.csv`/`.Rtab` schema and running the
  unmodified `roary_plots.py`
  (https://github.com/sanger-pathogens/Roary/blob/master/contrib/roary_plots/roary_plots.py)
  is real, low-effort reuse — it is a single, dependency-light Python script (not
  a full pipeline invocation) that would produce the frequency histogram and
  core/soft-core/shell/cloud pie chart essentially for free once the conversion
  script exists. Its presence/absence-vs-tree matrix figure is likely to render
  poorly at 295 strains × 47,983 rows (inferred from the script's own fixed
  figure-size and per-genome font-scaling logic, not a documented limit) and
  should not be relied on as the primary heatmap.
- **PPanGGOLiN and Panaroo are not adaptable this way** — no confirmed path to
  feed either an externally-computed matrix (PPanGGOLiN: unconfirmed either way
  from its own docs; Panaroo: is itself a from-genomes clustering pipeline we've
  already bypassed by using mmseqs2 tier-1 directly).
- **No existing tool renders a co-occurrence network at our scale in a
  self-contained way.** Coinfinder's own approach (hand off to Gephi) is the only
  pangenome-specific precedent and is structurally incompatible with our report
  format.

### Minimum-effort path to a useful, honest visualization

Given the two numbers above, and that the design spec's component 9 has already
been revised once by independent review (Opus, 2026-09-13) to fix real technical
gaps (bitset payload encoding, column virtualization, capped/tabular network
panel) — the honest minimum-effort path is a **staged commitment**, not an
all-or-nothing choice between "just run Roary's script" and "build the full
canvas report right now":

1. **Immediately, cheaply, and in parallel with anything else**: write the
   small Roary-format conversion script and run the unmodified
   `roary_plots.py` against it. This gets a real frequency histogram and
   core/shell/cloud pie chart for near-zero net-new code, using a tool whose
   plotting logic is already written and tested by someone else. Treat its
   presence/absence-matrix figure as unreliable at this scale (per the
   font/figure-size caveat above) and do not use it for anything beyond the
   histogram/pie chart. This is squarely a "do this first, it's nearly free"
   step, not a replacement for anything else.
2. **Build the interactive `pangenome.html` matrix view (component 9's first
   bullet) as already speced**, because no existing tool provides a
   self-contained, `file://`-openable, three-state (present/genome_only/absent)
   heatmap at this scale — PPanGGOLiN's own team hit the same scale wall we
   would hit with a naive per-cell encoding and mitigated it the same way
   (partition-dropping/row limits) that the design spec's bitset-packing
   solution avoids more gracefully (we keep all three states and all families,
   just packed). This is not optional/deferrable work substitutable by an
   existing tool — build it as speced.
3. **For the `trans`-pair network panel specifically**: do not attempt to render
   522,933 edges directly in any tool, existing or custom — no primary source
   surveyed (including graph-tool's own real ~300K-edge benchmark) claims
   pleasant *interactive browser* rendering at this edge count; graph-tool's
   own benchmark is a computational, not visual-rendering, benchmark. Run Leiden
   clustering (igraph + leidenalg, both with primary-sourced scale headroom well
   past 522,933 edges) offline as a batch step to collapse the graph into
   communities, then have `pangenome.html` render only the collapsed
   community-graph (fixed, non-simulated layout, hard-capped top-N edges by
   effect size) plus the mandatory full edge-list table — exactly as component
   9's S5 already specs, now with a concrete, primary-sourced algorithm/library
   choice (Leiden via igraph/leidenalg) behind the previously-abstract "collapse
   to communities" step.
4. **Sequencing**: steps 1 and 2 can proceed in parallel; step 3 depends on
   component 8's pair-classification output already existing (it does, per the
   design spec) and on picking the Leiden resolution parameter empirically
   against the real 522,933-edge graph — this is new, un-precedented tuning work
   (no pangenome-specific published precedent exists for what resolution/edge-cap
   values are "right" for this data, per section 2.3) and should be validated
   against a few known `starship_explained` pairs (do they land in sensible
   communities?) before trusting the collapsed view for `trans` pairs, mirroring
   the design spec's own "do not build this against synthetic data" and
   benchmark-suite discipline (components 6 and 9's last bullet).

**Net recommendation**: commit to the full interactive `pangenome.html`
(component 9) as speced — it is not replaceable by an existing tool at this
scale/shape — but do not treat that as blocking step 1's near-free Roary-plot
supplement, and treat the network panel's community-collapse step as the one
genuinely open, unprecedented design choice in this whole pathway (no published
pangenome-specific network-visualization precedent exists to validate against;
Leiden/igraph is the best-evidenced algorithm/library choice, but the specific
resolution parameter and edge-cap value need empirical tuning against the real
data, not literature lookup).
