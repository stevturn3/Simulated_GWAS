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
  4. 
# Methods 

Our PLINK GWAS returned promising results in regards to association between many SNPs and our simulated HDL Cholesterol levels.

After re-coding the binary PLINK files to a transposed text file and reading the file into R, linear regression was used to find significance. The results are equivalent to those found via plink!

We perform association tests for each SNP. For each SNP, we test the following Hypotheses:

$H_0 : \Beta_{variant} = 0$

$H_A : \Beta_{variant} \neq 0$

![alt text](https://github.com/stevturn3/Simulated_GWAS/blob/main/Turnbull_manhattan.png?raw=true)

We create a Manhattan plot to see how chromosomes are correlated to HDL Cholesterol. The red line indicates the negative log-probability of the significance level. Consequently, any SNP above this line is labeled significant in terms of its respective association test - thereby rejecting the respective null hypotheses.  In total, 65,875 SNPs(7.9% of SNPs) are labeled significant. Chromosome two has the highest percentage (ratio of significant/total) of significant SNPs, about 9%, and chromosome 11 has the lowest percentage of significant SNPs, about 7%. However, all of the chromosomes are within two percent of chromosome two’s ratio, thus we can reason that our results show that each chromosome has a significant effect on the HDL cholesterol levels of the individuals. Chromosomes 15,16 and 18 have SNPs that extremely significant (very tiny p-values) and would be interesting SNPs to study in an experimental fashion. It is also important to note that valleys, such as around 9 and 1, found in the Manhattan plot are regions of the chromosomes that very likely have no association with HDL Cholesterol levels and likely would not be labeled significant for any reasonable significance level.

# Studies

We compare our results to hose of a ukBB GWAS on the same subject. 

