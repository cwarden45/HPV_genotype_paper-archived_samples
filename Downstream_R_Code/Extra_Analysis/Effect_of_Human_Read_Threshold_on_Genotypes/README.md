Without the qPCR threshiold, the human read threshold helps for some conditions but not others (for 5% read fraction versus 20% read fraction results).  However, the combination of the 20% read fraction and qPCR filter produces the most consistent results for HPV18 versus HPV58.

These tables are also similar to what is presented in the [Human-Read_Adjusted_Divergence](https://github.com/cwarden45/HPV_genotype_paper-archived_samples/tree/master/Downstream_R_Code/Extra_Analysis/Human-Read_Adjusted_Divergence) subfolder.

The values are generated from [Table1_Supplemental_TableS1_S4-LIMITED-INPUT.R](https://github.com/cwarden45/HPV_genotype_paper-archived_samples/blob/master/Downstream_R_Code/Table1_Supplemental_TableS1_S4-LIMITED-INPUT.R).

**A) Effect of Read Fraction Thresholds *without* qPCR Filter with *5% Read Fractions***

**B) Effect of Read Fraction Thresholds *without* qPCR Filter with *20% Read Fractions***

<table>
  <tbody>
    <tr>
	<th align="center" colspan="2"></th>
	<th align="center" colspan="4">Parameters for Assigning HPV Genotypes<sup>a</sup></th>
    </tr>
    <tr>
	<td align="center" colspan="2">Overall HPV Threshold</td>
	<td align="center"><b>&gt1.5x Human</b></td>
	<td align="center"><b>&gt1.2x Human</b></td>
	<td align="center"><b>&gt1.0x Human</b></td>
	<td align="center"><b>&gt0.8x Human</b></td>
    </tr>
    <tr>
	<td align="center" colspan="2">HPV Genotype Threshold</td>
  	<td align="center"><b>&gt1.2x Human</b></td>
	<td align="center"><b>&gt1.0x Human</b></td>
	<td align="center"><b>&gt0.8x Human</b></td>
	<td align="center"><b>&gt0.6x Human</b></td>
    </tr>
    <tr>
	<td align="center" rowspan="3">Samples Positive for HPV Genotypes </td>
  	<td align="center">HPV16</td>
	<td align="center">74</td>
	<td align="center">74</td>
	<td align="center">74</td>
	<td align="center">76</td>
    </tr>
    </tr>
    <tr>
  	<td align="center">HPV18</td>
	<td align="center">24</td>
	<td align="center">24</td>
	<td align="center">24</td>
	<td align="center">24</td>
    </tr>
    </tr>
    <tr>
  	<td align="center">HPV58</td>
	<td align="center">17</td>
	<td align="center">20</td>
	<td align="center">20</td>
	<td align="center">26</td>
    </tr>
    <tr>
	<td align="center" rowspan="3">HPV58/HPV18 Ratio<sup>b</sup></td>
  	<td align="center">Archived DNA</td>
	<td align="center">0.25 (1 HPV58)</td>
	<td align="center">0.25 (1 HPV58)</td>
	<td align="center">0.25 (1 HPV58)</td>
	<td align="center">0.25 (1 HPV58)</td>
    </tr>
    </tr>
    <tr>
  	<td align="center">Frozen Tissue</td>
	<td align="center">0.22 (2 HPV58)</td>
	<td align="center">0.56 (5 HPV58)</td>
	<td align="center">0.56 (5 HPV58)</td>
	<td align="center">1.11 (10 HPV58)</td>
    </tr>
    </tr>
    <tr>
  	<td align="center">FFPE Tissue</td>
	<td align="center">1.27 (14 HPV58)</td>
	<td align="center">1.27 (14 HPV58)</td>
	<td align="center">1.27 (14 HPV58)</td>
	<td align="center">1.36 (15 HPV58)</td>
    </tr>
    <tr>
	<th align="center" colspan="2">Samples with “Unclear” Genotype Assignments<sup>c</sup></th>
	<td align="center"><b>Overall:</br></b>2 DNA</br>7 frozen</br>1 FFPE</br></br></br><b>Genotype-specific:</b></br>1 frozen</td>
	<td align="center"><b>Overall:</br></b>2 DNA</br>3 frozen</br>1 FFPE</br></br></br><b>Genotype-specific:</b></br>2 frozen</td>
	<td align="center"><b>Overall:</br></b>2 DNA</br>2 frozen</br></br><b>Genotype-specific:</b></br>3 frozen</br>1 FFPE></br></td>
	<td align="center"><b>Overall:</br></b>2 DNA</br>1 frozen</br></br><b>Genotype-specific:</b></br>[none]</br></br></td>
    </tr>
</tbody>
</table>

**C) Effect of Read Fraction Thresholds *with* qPCR Filter with *20% Read Fractions***

<sup>a<sup>All strategies require at least 20% HPV reads overall and for specific genotypes. See Table S3 for the effects of varying the read threshold. Samples with low amplified DNA concentrations are included in the calculations above. Exclusion of samples with amplified DNA concentrations <2-nM threshold reduces the number of HPV58+ samples, as well as the HPV58/HPV18 ration.
<sup>b<sup>Given that 1) the numbers of HPV18+ and HPV58+ samples (including co-infections) are similar, 2) there is more variation among HPV18 sequences than among HPV58 sequences, and 3) off-target human reads positively correlate with HPV58 reads, we believe that a high HPV58/HPV18 ratio likely indicates a relatively high false positive rate for HPV58 assignments. 
<sup>c<sup>Samples with an “unclear” overall HPV genotype have >20% HPV reads but cannot be assigned a specific HPV genotype due to a high percentage of human reads. Samples with an “unclear” specific genotype meet the read requirement for assignment as overall HPV+ but do not meet the read requirement for any specific HPV genotype. 
