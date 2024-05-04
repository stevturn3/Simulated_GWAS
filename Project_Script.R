# Access command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Load and download required packages
list_of_packages <- c("data.table", "gridExtra", 'qqman')
links <- c('https://Rdatatable.gitlab.io/data.table',"http://cran.univ-lyon1.fr",'https://github.com/stephenturner/qqman/')
packages_to_install <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]

if (length(packages_to_install)) {
  install.packages(packages_to_install, links[which(list_of_packages %in% packages_to_install)])
}

library('data.table')
library('gridExtra')
library('qqman')
source("functions.R")

path <- args[[1]]

# Perform GWAS in R to confirm Plink Results
p_vals <- compare_to_plink_results(path)

# Save P values if specified in arguments
if (args[[3]]) {
  save(file = paste(path, "/Turnbull_p_vals.R", sep=''), p_vals)
}

# Match corresponding SNPs for comparison
matches <- match(unname(unlist(plink_regression[,2])), p_vals[,1])
p_vals_reordered <- p_vals[matches,]
cat("Difference Between plink and R results: ", mean((p_vals_reordered[,2] - plink_regression[,1])^2))

build_manhattan(path, plink_regression)
compare_to_gwas(plink_regression, path)

# Build chromosome summary table if specified in arguments
if (args[[6]]) {
  chromosome_significance(plink_regression, path)
}