# HPV_typing
Code for Assigning HPV Types (tested with L1 amplicons)

### Reference Preparation ###

Please download full hg38 reference.  Code uses reference with UCSC chromosome names (starting with "chr"), so I would recommend downloading sequence from ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

Will provide script for combining and indexing reference, along with HPV reference set with paper submission.

### Order to Run Scripts ###

Scripts use a parameter file (set as `parameters.txt` in code).  Please see below for description of parameters.

Most dependencies are available in this [Docker image](https://hub.docker.com/r/cwarden45/hpv-project/).

1) `downsample_reads.py` You may need to modify regular expression to extract sample ID from your reads.

2) `cutadapt_filter.py`  You may need to modify regular expression to extract sample ID from your reads.

3) `BWA_MEM_alignment.py`

4) `call_genotype_mixed_ref.py`

### Dependencies (some optional) ###

Biopython: http://biopython.org/wiki/Biopython

Cutadapt: http://cutadapt.readthedocs.io/en/stable/index.html

PEAR: http://sco.h-its.org/exelixis/web/software/pear/

BWA: http://bio-bwa.sourceforge.net/

samtools: http://samtools.sourceforge.net/

limma: http://bioconductor.org/packages/release/bioc/html/limma.html

heatmap.3: https://github.com/obigriffith/biostar-tutorials/blob/master/Heatmaps/heatmap.3.R

heatmap.3 example: https://www.biostars.org/p/18211/

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
