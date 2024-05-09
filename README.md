# Goals

The goal of this project is to understand what is going on in a GWAS, understand compuational problems with large data, and understrand how to interpret and cross-reference GWAS results in practice.

# Introduction

Running a GWAS on HDL Cholesterol Levels where we wish to establish relationships between SNP (Single Nucelotide Polymoforisms) and gene expression in HDL leves on 2,504 individuals from the the 1000 Genomes Project. A GWAS - Genome Wide Association Test - is used to perform an association test on many different SNPs at once to determine if their effect on the phenotype - in this case simulated HDL Cholesterol levels that have been transformed - is significant. A GWAS typically occurs in 5 steps : Sample collection → genotyping → quality control → association analysis → supplemental analysis. The data received has already been quality controlled and thus the focus of this paper will be on the association analysis of this process. 

Our goal is to identify SNPs that correlate with HDL Cholesterol levels so that we can motivate possible medicines to regulate gene expression and reduce high HDL Cholesterol levels.

# Data

If you would like the simulated data, please email me at steventurnbull22@gmail.com

The data comes from four different files:

.bed : Primary Representation of genotype calls at biallelic variants. The format of the .bed file is in binary. This file must be accompanied by the .bim and .fam files

.bim : Extended variant information to accompany .bed binary genotype table. Represented in a text file with no header, each line is one variant with the following field:
  1. Chromosome Code (integer, X,Y, XY or 0(unknown))
  2. Variant key
  3. Morgan Position : A dummy value of 0 for our postions
  4. Base-pair coordinate
  5. Allele 1
  6. Allele 2

.fam : Sample information file. Again represented as a text file with no head, and one variant per line:
  1. Family ID
  2. Within-family ID (individual ID if no family info)
  3. Within-family ID of father (0 if father is missing)
  4. Within-family ID of mother (0 if mother is missing)
  5. Sex code (1 = male, 2 = female, 0 = unknown)
  6. Phenotype value (1 = control, 2 = case, -9/0 = missing)

phenotype.txt : Phenotype values for each individual: 
  1. Family ID
  2. Individual ID
  3. Phenotype Value
     
# Methods 

Our PLINK GWAS returned promising results in regards to association between many SNPs and our simulated HDL Cholesterol levels.

After re-coding the binary PLINK files to a transposed text file and reading the file into R, linear regression was used to find significance. The results are equivalent to those found via plink!

We perform association tests for each SNP. For each SNP, we test the following Hypotheses:

$H_0 : \beta_{variant} = 0$

$H_A : \beta_{variant} \neq 0$

![alt text](https://github.com/stevturn3/Simulated_GWAS/blob/main/Turnbull_manhattan.png?raw=true)

We create a Manhattan plot to see how chromosomes are correlated to HDL Cholesterol. The red line indicates the negative log-probability of the significance level. Consequently, any SNP above this line is labeled significant in terms of its respective association test - thereby rejecting the respective null hypotheses.  In total, 65,875 SNPs(7.9% of SNPs) are labeled significant. Chromosome two has the highest percentage (ratio of significant/total) of significant SNPs, about 9%, and chromosome 11 has the lowest percentage of significant SNPs, about 7%. However, all of the chromosomes are within two percent of chromosome two’s ratio, thus we can reason that our results show that each chromosome has a significant effect on the HDL cholesterol levels of the individuals. Chromosomes 15,16 and 18 have SNPs that extremely significant (very tiny p-values) and would be interesting SNPs to study in an experimental fashion. It is also important to note that valleys, such as around 9 and 1, found in the Manhattan plot are regions of the chromosomes that very likely have no association with HDL Cholesterol levels and likely would not be labeled significant for any reasonable significance level.

# Studies

We compare our results to hose of a ukBB GWAS on the same subject. We expect some differences in results do to our lack of data. Here we display the minimum values on each chromsome.

![alt text]([https://github.com/stevturn3/Simulated_GWAS/blob/main/Turnbull_MinP.png]?raw=true)

Additionally, we take a look at three grouping of these SNPs. We take the UK bio bank p-values as the true p-values of our experiment as it was performed in a similar way but including more data.

a.UK bio bank p-value is not significant and no studies found:

(We note that this is relatively unsurprising since the UK bio bank p-values
weren’t suggestive of an association, verifying it’s results for these SNPs - at least given the studies that have been conducted)
rs4846920, rs1112403, rs6532041,rs1036172, rs2819974, rs834793, rs4300315, rs10846772, rs9519977, rs35772501, rs9983496, rs732381

b.Both p-values are significant and no studies found:

(The p-values computed in both our simulated GWAS and the UK bio bank GWAS are very significant so it is surprising there are no studies that show an association. This mean either these results were false positives or the association has not been confirmed within a study)
rs35617716, rs6073958, rs2740488(Only mention in a study was about age-related muscular degeneration), rs429358 (Not mentioned in any studies related to HDL but appears to be associated with many other phe- notypes/diseases)

c.Both p-values and significant and mentioned in studies on HDL Cholesterol levels:

rs6754295: From Boes E et al., Experimental gerontology, 2009 showed an overall null effect of APOB(the gene that this SNP is located on) on HDL Cholesterol levels. This result was confirmed by Weissglas-Volkov D et al.,Journal of lipid research,2010, as they showed the association is insignificant with a p value of 4.4 × 10−8.

rs10750097: From Boes E et al.,Experimental gerontology,2009: ”rs10750097 (APOA5) never showed an as- sociation in any [of their] studies.” However, Shirts BH et al.,Atherosclerosis, 2009, found an association with ”25OHD-dependent changes in APOA5 promoter activity in HEP3B and HEK293 cells.” These researchers identified an interaction between rs10750097 and low HDL-C levels in individuals with low winter dietary 25OHD.

rs1077835: Consistently found to have an association to HDL Cholesterol levels.Brinkley TE et al.,Journal of applied physiology, 2011, found rs1077835 and other SNPs on the LIPC promoter were found to account for approximately ”20-30%” of HL variation and Hodoglugil U et al., Journal of lipid research,2010, found an association with plasma HDL-C levels.

rs11076175: Consistently found to show ”convincing” association with low HDL-C levels. Carlquist JF et al., Journal of clinical experimental cardiology, 2011, show it was the only SNP (out of a selection of SNPs from lipid-related genes such as ApoF, LIPC, LPL, SCARBI, LCAT, CETP - where rs11076175 is found) found to be significant under the Bonferroni adjusted significance level. GG homozygotes had a significant, 15% reduction in HDL levels.

rs4939883: Found to be significant in many different papers (including Shirts BH et al.,Atherosclerosis,2012) however, the association and effect was not discussed in depth.

