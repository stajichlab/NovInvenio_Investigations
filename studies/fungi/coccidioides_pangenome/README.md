## Pangenome analysis of Coccidioides with possible sub-split to C. immitis and C. posadasii runs


See the studies/fungi/Afumigatus_pangenome/PANGENOME_CLUSTER_PROFILE_NOTES.md and NEXTFLOW_MIGRATION_NOTES.md in that folder or any notes on how the workflow was constructed, assumptions that might need to be addressed.

## Data
* gives inventory of strains and their species designation /bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/samples.csv
* /bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/annotation_freeze/20260112 are the genome annotation sets, proteins, genome, GFF are all in there
* Previous Orthofinder results are in /bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome/ for some comparisons later on numbers gotten from mmseqs and the NII pipeline
* A phylogeny of strains based on ASTRAL consensus of gene trees /bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Phylogeny/results/tree_cds_ascomycota/astral/consensus_aster2.nw
* a second ML phylogeny is in here /bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Phylogeny/results/msa_filter_cds_ascomycota-buildtree/Cocci_cds.488taxa_ascomycota.fa.raxml.rba.raxml.bestTree noting that not all strains are present
* note the GFF3 for this dataset is from funannotate and will be coded different in how gene locus IDs are found from the NCBI GFF3 in the Afumigatus_pangenome data
* BUSCO results named by strainid.AAFTF are in /bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/BUSCO

The goals 
1. run the new ni (nf_NovInvenio) tool with --pipeline pangenome on Coccidioides perhaps on the whole set and then also on species groups alone (C. immitis and C. posadasii)
2. Develop figures to examine pangenome properties, open/closed 
3. Look for co-gain and co-lost families
4. Though no starship / starfish analysis has been run yet, examine candidates for clusters of lost/gained genes to see if there are any patterbns
5. Also look for significant correlations of gain/loss events of 2 or more genes, even in trans, to suggest co-necessary genes, potential secondary metabolite processes

Develop markdown reports with embedded PNG figures and links to high quality PDFs. Table which represent the summary statistics about genomes, gene count. Consider conditioning data on quality of assemblies and removing low quality (by BUSCO scores or possible other metrics) that could be outliers in quality

Propose other investigations based on experts in population genomics and bioinformatics to identify patterns of gene, gene family, and gene organization evolution and change that might distinguish geographically different isolates, or isolates from environmental vs clinical origins.

Assess any data problems or quality first.  Review sensibility of results with Fable where appropriate and propose any workflow code improvements to generalize.

## Additional consideration
* We are using mmseqs clustering as the primary clustering. Name the result folders in a way that it would be possible to replace it with diamond clustering to see how it makes a difference
