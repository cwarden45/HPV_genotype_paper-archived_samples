Statistical significance of the effects of sample type on HPV genotype assignments in HPV+ tumor samples.
-----------------

### Overview for Table

Unless the analysis directly uses the 15% read fraction genotypes, it should be assumed that genotypes and HPV+ samples were defined with 5% read fractions (to be more similar to not using a genotyping step).

Fisher’s exact test (FE) and logistic regression (binomial GLM) require discrete genotype assignments. All other p-values were calculated without making genotype assignments. False Discovery Rates (FDRs) were calculated using the method of [Benjamini and Hochberg](https://www.jstor.org/stable/2346101).

### P-value Methods Summary

For analysis of sample type (FFPE versus DNA, FFPE versus frozen, and frozen versus DNA) samples were compared using the following strategies:

**1)** *Fisher’s exact test* (using tentative ***genotype*** assignments, for both 5% and 15% read fractions)

**2)** *limma-voom* ([Law et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)) on ***read counts*** (adjusting for percent human reads), 

**3)** *logistic regression* (binomial generalized linear model, or binomial GLM) was performed on tentative ***genotype*** assignments using the `glm()` R-base function with the parameter “family=’binomial’”

**4)** *linear regression* on ***linear read fractions*** from ***read counts*** (values between 0 and 100, with and without a second variable for percent human reads) using the R-base `lm()` function

**5)** *linear regression* on ***log2(read fraction+ 1) values*** from ***read counts*** (with and without a second variable for percent human reads) using the R-base `lm()` function.

Log-transformed linear regression is therefore different from the binomial generalized linear model (GLM), which represents logistic regression for the status of the discrete HPV type.
