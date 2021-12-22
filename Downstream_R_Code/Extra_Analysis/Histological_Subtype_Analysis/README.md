P-values for histological subtype (**SCC vs Adeno Types**) analysis were calculated using the following strategies:

For analysis, **"Adeno Types"** represents both "*Adenocarcinoma*" as well as "*Adenosquamous*".

The following comparisons were run on the samples with **>20% read fractions and qPCR filter** as well as **>5% read fractions (and no qPCR filter)**:

**1)** ***Fisher’s exact test*** (using tentative **genotype** assignments, for 20% read fractions) for **separate archive types**, using the `fisher.test()` R-base function.

**2a)** ***limma-voom*** ([Law et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)) for **separate archive types** on **read counts** with 1-variable (the ancestry assignment) 

**2b)** ***limma-voom*** ([Law et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)) for **overall samples** on **read counts** with 2-variables (*adjusting for archive type*, along with histological subtype).

**3a)** binomial GLM for logistic regression for **separate archive types** on **read counts** with 1-variable (the ancestry assignment), using the `glm()` R-base function with the parameter `family="binomial"`.

**3b)** binomial GLM for logistic regression for **overall samples** on **read counts** with 2-variables (*adjusting for archive type*, along with histological subtype), using the `glm()` R-base function with the parameter `family="binomial"`.

**5% Minimum Read Fraction Assignments *without* qPCR Filter:**

<table>
  <tbody>
    <tr>
	<th align="center">Histological Subtype</th>
	<th align="center">Archived</br>DNA</th>
	<th align="center">Frozen</br>Tissue</th>
	<th align="center">FFPE</br>Tissue</th>
    </tr>
    <tr>
	<td align="left">Squamous Cell Carcinoma (SSC)</td>
	<td align="center"><u><b>HPV16</b></u>: 18 / 23 (<b><i>78%</b></i>)</br><u>HPV18</u>: 2 / 23 (<i>9%</i>)</br><u>HPV16 + HPV18</u>: 0 / 23 (<i>0%</i>)</td>
	<td align="center"><u><b>HPV16</b></u>: 15 / 25 (<b><i>60%</b></i>)</br><u>HPV18</u>: 4 / 25 (<i>16%</i>)</br><u>HPV16 + HPV18</u>: 0 / 25 (<i>0%</i>)</td>
	<td align="center"><u><b>HPV16</b></u>: 15 / 37 (<b><i>41%</b></i>)</br><u>HPV18</u>: 5 / 37 (<i>14%</i>)</br><u>HPV16 + HPV18</u>: 8 / 37 (<i>22%</i>)</td>
    </tr>
    <tr>
	<td align="left">Adenocarcinoma (Adeno)</td>
	<td align="center"><u>HPV16</u>: 0 / 3 (<i>0%</i>)</br><u><b>HPV18</b></u>: 1/ 3 (<b><i>33%</b></i>)</br><u>HPV16 + HPV18</u>: 0 / 3 (<i>0%</i>)</td>
	<td align="center"><u>HPV16</u>: 1 / 7 (<i>14%</i>)</br><u><b>HPV18</b></u>: 4 / 7 (<b><i>57%</b></i>)</br><u>HPV16 + HPV18</u>: 0 / 7 (<i>0%</i>)</td>
	<td align="center"><u>HPV16</u>: 7 / 13 (<i>54%</i>)</br><u>HPV18</u>: 4 / 13 (<i>31%</i>)</br><u>HPV16 + HPV18</u>: 1 / 13 (<i>8%</i>)</td>
    </tr>
    <tr>
	<td align="left">Adenosquamous Carcinoma (ASC)</td>
	<td align="center"><u>HPV16</u>: 0 / 1 (<i>0%</i>)</br><u>HPV18</u>: 1/ 1 (</b><i>100%</b></i>)</br><u>HPV16 + HPV18</u>: 0 / 1 (<i>0%</i>)</td>
	<td align="center"><u>HPV16</u>: 0 samples </br><u>HPV18</u>: 0 samples </br><u>HPV16 + HPV18</u>: 0 samples </td>
	<td align="center"><u>HPV16</u>: 2 / 5 (<i>40%</i>)</br><u>HPV18</u>: 1 / 5 (<i>20%</i>)</br><u>HPV16 + HPV18</u>: 1 / 5 (<i>20%</i>)</td>
    </tr>
</tbody>
</table>

**20% Minimum Read Fraction Assignments *with* qPCR Filter:**

<table>
  <tbody>
    <tr>
	<th align="center">Histological Subtype</th>
	<th align="center">Archived</br>DNA</th>
	<th align="center">Frozen</br>Tissue</th>
	<th align="center">FFPE</br>Tissue</th>
    </tr>
    <tr>
	<td align="left">Squamous Cell Carcinoma (SSC)</td>
	<td align="center"><u><b>HPV16</b></u>: 18 / 23 (<b><i>78%</b></i>)</br><u>HPV18</u>: 2 / 23 (<i>9%</i>)</br><u>HPV16 + HPV18</u>: 0 / 23 (<i>0%</i>)</td>
	<td align="center"><u><b>HPV16</b></u>: 14 / 25 (<b><i>56%</b></i>)</br><u>HPV18</u>: 4 / 25 (<i>16%</i>)</br><u>HPV16 + HPV18</u>: 0 / 25 (<i>0%</i>)</td>
	<td align="center"><u><b>HPV16</b></u>: 12 / 27 (<b><i>44%</b></i>)</br><u>HPV18</u>: 5 / 27 (<i>19%</i>)</br><u>HPV16 + HPV18</u>: 0 / 37 (<i>0%</i>)</td>
    </tr>
    <tr>
	<td align="left">Adenocarcinoma (Adeno)</td>
	<td align="center"><u>HPV16</u>: 0 / 3 (<i>0%</i>)</br><u><b>HPV18</b></u>: 1/ 3 (<b><i>33%</b></i>)</br><u>HPV16 + HPV18</u>: 0 / 3 (<i>0%</i>)</td>
	<td align="center"><u>HPV16</u>: 1 / 7 (<i>14%</i>)</br><u><b>HPV18</b></u>: 4 / 7 (<b><i>57%</b></i>)</br><u>HPV16 + HPV18</u>: 0 / 7 (<i>0%</i>)</td>
	<td align="center"><u>HPV16</u>: 6 / 11 (<i>55%</i>)</br><u>HPV18</u>: 4 / 11 (<i>36%</i>)</br><u>HPV16 + HPV18</u>: 1 / 11 (<i>9%</i>)</td>
    </tr>
    <tr>
	<td align="left">Adenosquamous Carcinoma (ASC)</td>
	<td align="center"><u>HPV16</u>: 0 / 1 (<i>0%</i>)</br><u>HPV18</u>: 1/ 1 (</b><i>100%</b></i>)</br><u>HPV16 + HPV18</u>: 0 / 1 (<i>0%</i>)</td>
	<td align="center"><u>HPV16</u>: 0 samples </br><u>HPV18</u>: 0 samples </br><u>HPV16 + HPV18</u>: 0 samples </td>
	<td align="center"><u>HPV16</u>: 2 / 3 (<i>67%</i>)</br><u>HPV18</u>: 1 / 3 (<i>33%</i>)</br><u>HPV16 + HPV18</u>: 0 / 3 (<i>0%</i>)</td>
    </tr>
</tbody>
</table>

In the tables above, replicate samples are counted more than once.  However, samples must have a histological subtype, meaning that these are all cervical cancer samples.
**HPV16 has been previously been published to be more common in SSC, and HPV18 has previously been published to be more common in Adenocarcinoma.**  As expected, HPV16 + HPV18 co-infections decrease with the qPCR filter.


Here is a barplot of the samples with **20%** HPV genotype assignments (**with** the qPCR filter), for the 3 most common histological subtypes:

![20% qPCR Filtered Histological Subtype Distribution](HistSubtype-qPCR20.png "20% qPCR Filtered Histological Subtype Distribution")

Here is a barplot of the samples with **5%** HPV genotype assignments (***without*** the qPCR filter), for the 3 most common histological subtypes:

![5% Histological Subtype Distribution](HistSubtype-all5.png "%Histological Subtype Distribution")

Among the FFPE HPV58+ samples, you can see the percent off-target human reads per sample below (when separated by histological subtype):

![FFPE HPV58 Human Reads by Histological Subtype](FFPE_HPV58_HumanReads_by_HistologicalSubtype.png "FFPE HPV58 Human Reads by Histological Subtype")

Points are colored based upon whether they are also kep in the 20% assignments (with the qPCR filter).

The ANOVA p-value for percent off-target human reads varying as a function of histological subtype is 0.23 (among the above set of HPV58+ FFPE samples).
