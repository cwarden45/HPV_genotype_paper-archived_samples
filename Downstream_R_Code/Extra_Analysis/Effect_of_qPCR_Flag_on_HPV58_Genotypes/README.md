The more stringent read fraction filters reduce the variation in HPV58 genotype assignments in FFPE samples.  However, the qPCR filter essentially removes the increased frequency of HPV58 assignments in FFPE samples (at any read fraction threshold).

**A) Effect of Read Fraction Thresholds *without* qPCR Filter**

<table>
  <tbody>
    <tr>
	<th align="center"></th>
	<th align="center">Archived</br>DNA</th>
	<th align="center">Frozen</br>Tissue</th>
	<th align="center">FFPE</br>Tissue</th>
	<th align="center" rowspan="2">FFPE vs DNA</br>FE P-value</th>
    </tr>
    <tr>
	<td align="left">Number of Samples</td>
	<td align="center">28*</td>
	<td align="center">40*</td>
	<td align="center">57</td>
    </tr>
    <tr>
	<td align="left">HPV58 Detected</br>(>5% Reads)</td>
	<td align="center">1 (3.6%)</td>
	<td align="center">6 (15.0%)</td>
	<td align="center">22 (38.6%)</td>
	<td align="center"><b>5.0 x 10<sup>-4</sup></b></td>
    </tr>
    <tr>
	<td align="left">HPV58 Alternative Detection Status</br>(>20% Reads)</td>
	<td align="center">1 (3.6%)</td>
	<td align="center">4 (10.0%)</td>
	<td align="center">14 (24.6%)</td>
	<td align="center"><b>0.017</b></td>
    </tr>
    <tr>
	<td align="left">HPV58 Plurality</br>(Most Abundant HPV-Type)</td>
	<td align="center">1 (3.6%)</td>
	<td align="center">4 (10.0%)</td>
	<td align="center">10 (17.5%)</td>
	<td align="center"><b>0.042</b></td>
    </tr>
    <tr>
	<td align="left">HPV58 Dominant / Majority</br>(>50% Reads)</td>
	<td align="center">1 (3.6%)</td>
	<td align="center">2 (5.0%)</td>
	<td align="center">8 (14.0%)</td>
	<td align="center">0.078</td>
    </tr>
</tbody>
</table>

**B) Effect of Read Fraction Thresholds *WITH* qPCR Filter**


\*There are 36 Archived DNA samples, but 28 is the count of samples excluding 6 prostate cancer negative controls (1 of which would have been excluded if 2 nM is used as the qPCR threshold for a QC flag).  Frozen adjacent normal samples were also excluded

P-values are calculating using a Fisher’s exact test (FE) for a given set of genotype assignments.  These p-values are slightly different than reported in the folder [Sample_Type_Full_Statistical_Analysis](https://github.com/cwarden45/HPV_genotype_paper-archived_samples/tree/master/Downstream_R_Code/Extra_Analysis/Sample_Type_Full_Statistical_Analysis) folder because that other analysis including more sample filtering for HPV+ samples (looking at 81 rather than 85 samples, where all 4 filtered samples were HPV58-).  More specifically, these results are slightly more significant (since you have fewer HPV58- samples).

All criteria for HPV58 detection require that total HPV reads > 1.2x human reads, HPV58 read count was greater than 1.0x human reads, and >5% HPV58 reads.  Notice that a plurality requirement is stricter than increasing detection threshold from 5% to 20%.
