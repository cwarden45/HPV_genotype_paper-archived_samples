# HPV_typing
Code for Assigning HPV Types (tested with L1 amplicons)

### Reference Preparation ###

Please download full hg38 reference.  Code uses reference with UCSC chromosome names (starting with "chr"), so I would recommend downloading sequence from ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

Will provide script for combining and indexing reference, along with HPV reference set with paper submission.

### Order to Run Scripts ###

Most dependencies are available in this [Docker image](https://hub.docker.com/r/cwarden45/hpv-project/).

1) `downsample_reads.py`


### Dependencies (some optional) ###

Cutadapt: http://cutadapt.readthedocs.io/en/stable/index.html

Biopython: http://biopython.org/wiki/Biopython

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
