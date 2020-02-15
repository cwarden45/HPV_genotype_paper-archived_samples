Statistical significance of the effects of sample type on HPV genotype assignments in HPV+ tumor samples.
-----------------

Limma-voom uses a multivariate model to detect differences between sample types, adjusting for the percentage of human reads. Regular linear regression using linear percentages or log2(percent + 1) values was conducted for 1 variable (sample type) to provide a p-value without knowing the effect of human reads, then for 2 variables (sample type + percentage of human reads) to match the other comparisons. Log-transformed linear regression is therefore different from the binomial generalized linear model (GLM), which represents logistic regression for the status of the discrete HPV type.

Fisher’s exact test (FE) and logistic regression (binomial GLM) require discrete genotype assignments. All other p-values were calculated without making genotype assignments. 

Unless the analysis directly uses the 15% read fraction genotypes, it should be assumed that genotypes and HPV+ samples were defined with 5% read fractions (to be more similar to not using a genotyping step).
