It should be noted that some HPV genotypes may not be well amplified with this primer set.

For example HPV52 has 5 mutations from the forward primer and 3 mutations from the reverse primer.

A small faction of reads can be aligned to the HPV52 genome reference.  The BWA-MEM alignment is shown below:
![BWA-MEM L1 Alignment](HPV_genotype_by_Age.png "BWA-MEM L1 Alignment")

If Bowtie1 is used instead of BWA-MEM, then the HPV52 false positives art lower:
![Bowtie1 L1 Alignment](HPV_genotype_by_Age.png "Bowtie1 L1 Alignment")

However, we believe that we underestimate the percent of human off-target reads if we use Bowtie1 instead of BWA-MEM.

I should also acknowledge that I learned about this from an anonymous reviewer, citing the [Matsukura and Sugase 2004](https://pubmed.ncbi.nlm.nih.gov/15207629/) publication.