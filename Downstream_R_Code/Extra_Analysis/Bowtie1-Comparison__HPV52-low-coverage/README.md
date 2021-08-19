It should be noted that some HPV genotypes may not be well amplified with this primer set.

For example, HPV52 has 5 mutations from the forward primer and 3 mutations from the reverse primer.

A small faction of reads can be aligned to the HPV52 genome reference.  The BWA-MEM alignment is shown below:
![BWA-MEM L1 Alignment](IGV_BWA-MEM_HPV52-zoom.png "BWA-MEM L1 Alignment")

If Bowtie1 is used instead of BWA-MEM, then the HPV52 false positives art lower:
![Bowtie1 L1 Alignment](IGV_Bowtie1-HPV52-zoom.png "Bowtie1 L1 Alignment")

However, we believe that we **underestimate** the percent of *human off-target reads* if we use Bowtie1 instead of BWA-MEM.

I should also acknowledge that I learned about this from an anonymous reviewer, citing the [Matsukura and Sugase 2004](https://pubmed.ncbi.nlm.nih.gov/15207629/) publication.

The Bowtie1 alignment were run outside of the Docker image for this project.  However, Bowtie version 1.2.2 was used with the following command:

```
bowtie --best -p $THREADS -S $REF -1 $R1cutadapt -2 $R2cutadapt $SAMout
```
