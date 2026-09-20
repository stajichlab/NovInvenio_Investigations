# Design decisions agreed during grilling (2026-09-20)
Q1 genome_only != present analytically. AGREED (user: "no they aren't equivalent").
Q2 assembly-quality-vs-content: emit regression as QC every run. AGREED -> issue #130.
    No statistical correction, no auto-exclusion.
Q3 split min_clades into two jobs: keep k=2-style stratification for the null
    (FPR 0.050 vs 0.327 unstratified); stop using min_clades as a `trans` FILTER
    when the structure diagnostics say the partition is noise; emit diagnostics
    either way. AGREED ("I agree with all of this").
Q4 gain/loss: KEEP direction, ADD asymmetry_a, GATE like #130 with FIXED thresholds
    to start + docs to revisit per dataset; support >2 groups. AGREED.
    Decisive evidence: A. fumigatus gain:loss 7.1:1 with 11.6% ambiguous and
    5,871,769 rescued cells; Coccidioides 143:1 with 0% ambiguous and 0 rescued.
    The method is sound; Coccidioides is the broken case.
Q5 assumptions register: cover NON-PARAMETER assumptions too. AGREED.
    Standard = publication-defensible: each cutoff justified by statistics,
    empirical result, or theory+practical constraint. Anything else = open liability.
    User will start a separate session developing pangenome TEST SETS.
Q5b sensitivity analysis on clustering/cutoffs: APPROVED as needed for publication.
Q6 rescue: NOT a threshold problem. 77.6% of rescuable cells overlap a predicted
    gene of a DIFFERENT family; 71.9% of query reps are <80% the length of the gene
    they hit, with a cliff at exactly mmseqs -c 0.8. Redundant-locus 74.2%,
    true paralog 3.3%, genuine dropout ~1.3%. Tightening pident/qcov does NOT help
    (class mix flat 19-23% from (90,80) to (99.9,99)) because the artifact hits are
    100%-identical real DNA. AGREED SEQUENCING: fix clustering first (#132), add the
    structural overlap criterion second. Leave rescue thresholds alone meanwhile.
    Independent confirmation: fragmented assemblies carry 2.7x more strain-private
    families (median 54 vs 20; rho +0.29 contigs / -0.27 N50) - same root cause as
    the accessory-inflation finding.
Q7 enforcement: ADVISORY by default + --pangenome_strict to promote warnings to
    failures. AGREED. Warnings MUST land at the top of report.md AND the report
    HTML pages (not buried in stderr - the rescue failure emitted 530 warnings
    per run for a whole study and nobody saw them). Each warning must also offer
    CONCRETE PRUNING OPTIONS with measured trade-offs, not just state a problem.
