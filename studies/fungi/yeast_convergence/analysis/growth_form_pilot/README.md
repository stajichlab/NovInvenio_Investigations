# Growth-form literature pilot (2026-10-25)

Record of a pilot that tested an agent-based literature search for the growth form (yeast, filamentous,
dimorphic) of 10 yeast_convergence species. Plan for a full run: `../../GROWTH_FORM_LITERATURE_PLAN.md`.
Issue: #44.

## Files

| File | Content |
|---|---|
| `reports/*.txt` | The five agent reports, as returned. Two species per report. Copied from the chat session, not edited. |
| `prompt.txt` | The reader prompt. Only the species line differed between agents. |
| `verification.tsv` | What I checked afterwards: 17 quotes or key phrases from 11 sources. All matched. |
| `PROVENANCE.yaml` | How the files were made, tools used, and the SHA-256 of each file. |

## How far to trust the reports

- The reports are model output. A quote is only as good as its check. Only the 17 items in
  `verification.tsv` were checked. Every other record is unchecked.
- The full-text checks searched a key phrase of the quote in the PMC XML, not the whole quote.
- The PMC text tool strips italic species names. Where a report says "(species name stripped)", the
  agent put the name in from context. The PMC XML keeps the names.
- A record that supports `other` or is labelled genus-level is not evidence for the species.

## Result in one table (agents blind to our codings; compared afterwards)

| Species | Supported by the literature | Batch coding | Note |
|---|---|---|---|
| *Eremothecium cymbalariae* | filamentous | filamentous | contradicts `traits.csv` "yeast" |
| *Rhizopus delemar* | filamentous | filamentous | known answer (Tahiri 2026) |
| *Mucor circinelloides* | dimorphic, hyphal | dimorphic, hyphal | known answer |
| *Yarrowia lipolytica* | dimorphic; phase depends on medium | dimorphic, yeast | phase not supported |
| *Aureobasidium pullulans* | dimorphic; phase unclear | dimorphic, yeast | abstracts only |
| *Exophiala dermatitidis* | dimorphic, yeast-like cells | dimorphic, yeast | abstracts only |
| *Taphrina deformans* | dimorphic; phase unclear | dimorphic, hyphal | abstracts only |
| *Protomyces inouyei* | genus-level only | dimorphic, hyphal | no species-level source |
| *Thamnidium elegans* | unsourced | dimorphic, hyphal | only incidental hyphal culture |
| *Neolecta irregularis* | unsourced for culture | filamentous | never cultured |

Cost: about 516,000 tokens in total for 10 species, 43 to 70 seconds per agent.
