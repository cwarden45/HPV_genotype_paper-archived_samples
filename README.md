```diff
+ This code was developed for a project within the lab of Javier Gordon Ogembo.
```

To improve readabilty, we focus on one particular genotyping strategy in the paper.

However, we want to be more careful about providing the right emphasis for the HPV genotyping code, since **we are not sure those exact parameters are actually the best to use across sample types and/or batches**.  In fact, for at least some of our samples, hypothesize that the more likely scenario is that a more conservative strategy may need to be used (*especially if users did not have access to something like the qPCR concentrations used to flag samples with low amounts of DNA amplified*).

Nevertheless, this is the subfolder with code that may be more easy for others to use:

```
Paried-End_FASTQ_Code
```

and this is the subfolder with R code to help reproduce figures for our study (but which we don't expect others to directly use for new samples):

```
Downstream_R_Code
```
