If you follow the similar analysis to [Human-Read_Adjusted_Divergence](https://github.com/cwarden45/HPV_type_paper-archived_samples/tree/master/Downstream_R_Code/Extra_Analysis/Human-Read_Adjusted_Divergence), you can see that applying the qPCR filter improves read fraction concordance in the sense that it removes the HPV58+ samples that had a lower correlation coefficient.

### Read frequencies for common HPV genotypes among tumor samples from the same patient.

![paired sample correlations](qPCR_flag_removes_HPV58-pos_pairs.png "paired sample correlations")

**(TOP)** Three types of tumor-tumor pairs were compared to evaluate consistency in HPV16, HPV18, and HPV58 genotyping: pairs of FFPE tissue samples from the same patient, as reported in patient records (“FFPE:Both”); pairs of frozen and FFPE tissue samples from the same patient, as reported in patient records (“Mixed:Reported”); and pairs of archived DNA and frozen tissue samples, matched via QC Array data (not reported in sample records, “Mixed:QCarray”). Correlations between the read frequencies for each pair were lower for HPV58 than for HPV16 and HPV18. 

**(BOTTOM)** Same as (TOP), but with qPCR flagged samples removed.  Notice that HPV58+ sample pairs can no longer be defined, so we it is no longer reasonable to calculate a correlation coefficient for those samples.  However, correlation coefficient for HPV16 is slightly higher (now matching HPV18).