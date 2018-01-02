# HPV_typing
Code for Assigning HPV Types (tested with L1 amplicons)

### Reference Preparation ###

Please download full hg38 reference.  Code uses reference with UCSC chromosome names (starting with "chr"), so I would recommend downloading sequence from ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

Will provide script for combining and indexing reference, along with HPV reference set with paper submission.

### Order to Run Scripts ###

Scripts use a parameter file (set as `parameters.txt` in code).  Please see below for description of parameters.

Dependencies are installed in this [Docker image](https://hub.docker.com/r/cwarden45/hpv-project/).

1) `python downsample_reads.py` - assumes you have uncompressed reads with .fastq extension

2) `python cutadapt_filter.py` - uses sample description to extract / re-name samples (in format used for next step)

3) `python BWA_MEM_alignment.py`

4) `python call_genotype_mixed_ref.py`

### Dependencies (some optional) ###

R: https://cran.r-project.org/

Biopython: http://biopython.org/wiki/Biopython

Cutadapt: http://cutadapt.readthedocs.io/en/stable/index.html

BWA: http://bio-bwa.sourceforge.net/

samtools: http://samtools.sourceforge.net/


### Parameter Values ###
| Parameter | Value|
|---|---|
|Sample_Description_File|Name of Tab-Delimited Sample Description File; **SampleID** column for sample ID, **Foward_Read** column for forward (R1) read, and **Reverse_Read** column for reverse (R2) read.  Used by `python cutadapt_filter.py`.  You only need to provide the basenames (e.g. Sample1_R1.fastq and Sample1_R2.fastq), with the full path for each step of read processing determined using *Reads_Folder*.|
|Reads_Folder|Path to Raw Reads|
|Max_Reads|Maximum Number of Down-Sampled Reads|
|Forward_Primer|Sequence for Forward Primer for Amplicon|
|Reverse_Primer|Sequence for Reverse Primer for Amplicon|
|Alignment_Folder|Output folder for BWA-MEM Alignment and Associated Files|
|BWA_Ref|Path to hg38 + HPV Indexed Reference|
|Threads|Number of Threads for BWA|
|HPV_Human_Ref_FC|Minimum Fold-Change above Off-Target Human Counts, for overall HPV status (if value is 1.5, human counts must be below 40%, since 1.5x of 40% is 60%).  Most likely over-rides *HPV_Freq_Cutoff*.|
|HPV_Freq_Cutoff|Minimum Frequency to Call Sample HPV+ (Percent, 0-100)|
|Type_Human_Ref_FC|Minimum Fold-Change above Off-Target Human Counts, for specific HPV type.  Should be less than or equal to *HPV_Human_Ref_FC*|
|Type_Freq_Cutoff|Minimum Frequency to Call Specific HPV Type (Percent, 0-100)|
|Summary_File|HPV Type Table, for all samples in *Alignment_Folder*|
