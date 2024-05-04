#!/usr/bin/bash

$1/plink --noweb --bfile $1/snpdata --linear --pheno $1/phenotype.txt --out $1/plink

$1/plink --noweb --bfile $1/snpdata --recode --transpose --out $1/snp_dat_trans

py File_split.py $1

Rscript Project_Script.R $1 $2 $3 $4 $5 $6
