"""Normalize a UniProt-derived protein_id to its bare accession.

nf_NovInvenio's presence_matrix.tsv/novelties.<Short>.tsv protein_id column carries
the full FASTA header token UniProt writes -- "sp|P12345|NAME_ORG" (reviewed) or
"tr|P12345|NAME_ORG" (unreviewed) -- since that's the first whitespace-delimited
token of each sequence's FASTA header, and that's what the pipeline's own FASTA
parsing uses as protein_id. bin/extract_dat_annotations.py, by contrast, keys its
output by the bare accession alone (from the .dat.gz AC line). Every script that
joins the two (bin/merge_uniprot_annotations.py, bin/go_enrichment.py,
bin/domain_enrichment.py) needs both sides normalized to the same key -- this is
that one normalizer, used everywhere so the two representations don't quietly
drift back out of sync.
"""
import re

_SP_TR_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|")


def bare_accession(protein_id: str) -> str:
    """'sp|P12345|NAME_ORG' or 'tr|P12345|NAME_ORG' -> 'P12345'. Already-bare
    accessions (e.g. from a non-UniProt source, or already-normalized input)
    pass through unchanged."""
    m = _SP_TR_RE.match(protein_id)
    return m.group(1) if m else protein_id
