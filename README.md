# HPV_typing
Code for Assigning HPV Types (tested with L1 amplicons)

### Reference Preparation ###

Please download full hg38 reference.  Code uses reference with UCSC chromosome names (starting with "chr"), so I would recommend downloading sequence from ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

Will provide script for combining and indexing reference, along with HPV reference set with paper submission.

### Order to Run Scripts ###

Scripts use a parameter file (set as `parameters.txt` in code).  Please see below for description of parameters.

Dependencies are installed in this [Docker image](https://hub.docker.com/r/cwarden45/hpv-project/).

1) `python downsample_reads.py` You may need to modify regular expression to extract sample ID from your reads.

2) `python cutadapt_filter.py`  You may need to modify regular expression to extract sample ID from your reads.

3) `python BWA_MEM_alignment.py`

4) `python call_genotype_mixed_ref.py`

5) `PE_HPVtype_counts.R`

Within an R session, you can run `source("PE_HPVtype_counts.R")`. Running `Rscript PE_HPVtype_counts.R` via command line should also do the trick.

### Dependencies (some optional) ###

R: https://cran.r-project.org/

Biopython: http://biopython.org/wiki/Biopython

Cutadapt: http://cutadapt.readthedocs.io/en/stable/index.html

PEAR: http://sco.h-its.org/exelixis/web/software/pear/

BWA: http://bio-bwa.sourceforge.net/

samtools: http://samtools.sourceforge.net/

### Recommendations for Downstream Analysis ###

* Remember that the final line in the count table produced by `PE_HPVtype_counts.R` is the total joint aligned fragments from idxstats (they can be used to define library sizes or define abundances, but should not be treated like a HPV type count)

* Likewise, the 2nd to last line in the count table produced by `PE_HPVtype_counts.R` is the number of human-aligned reads, which you may want to take into consideration.

* limma can be used for analysis of abundances (percent HPV type, values between 0 and 1, for example) or counts (with limma-voom)

### Parameter Values ###
| Parameter | Value|
|---|---|
|Reads_Folder|Path to Raw Reads|
|Max_Reads|Maximum Number of Down-Sampled Reads|
|Forward_Primer|Sequence for Forward Primer for Amplicon|
|Reverse_Primer|Sequence for Reverse Primer for Amplicon|
|Alignment_Folder|Output folder for BWA-MEM Alignment and Associated Files|
|BWA_Ref|Path to hg38 + HPV Indexed Reference|
|Threads|Number of Threads for BWA|
|HPV_Freq_Cutoff|Minimum Frequency to Call Sample HPV+ (Percent, 0-100)|
|Type_Freq_Cutoff|Minimum Frequency to Call Specific HPV Type (Percent, 0-100)|
|Summary_File|HPV Type Table, for all samples in *Alignment_Folder*|
|Count_File|Table of Counts for HPV types (with total idxstats fragment counts at bottom of table), for all samples in *Alignment_Folder*|
