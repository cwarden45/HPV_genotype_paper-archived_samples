P-values for ancestry (**EUR vs AMR**) analysis were calculated using the following strategies:

**1)** ***Fisher’s exact test*** (using tentative **genotype** assignments, for 15% read fractions), using the `fisher.test()` R-base function.

**2a)** ***limma-voom*** ([Law et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)) on **read counts** with 1-variable (the ancestry assignment) 

**2b)** ***limma-voom*** ([Law et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)) on **read counts** with 2-variables (*adjusting for sample type*, along with predicted ancestry).

Ancestry predictions required the QC Array, which was only available for archived DNA and frozen tissue samples (so, only those two types could be combined).
