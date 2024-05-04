#' Reads in a certain subset of a list of files
#' 
#' @param files_list A list of file paths
#' @param i Index indicating which file to read from the list
#' @return A data frame containing the read-in data
read_in <- function(files_list, i) {
  ped1 <- as.data.frame(fread(files_list[i]))
  ped1
}

#' Merges the  adjacent columns to form allele form columns
#' 
#' @param ped_mat Input matrix with data to be fixed and augmented
#' @return Augmented and corrected data frame
fix_columns_and_augment <- function(ped_mat) {
  i <- which(snpdata_bim$SNP_id == unname(unlist(ped_mat[,2]))[1])
  ped_mat <- data.frame(ped_mat[1:4], snpdata_bim$Allele_1[c(seq(i,(i + 9744)))],snpdata_bim$Allele_2[c(seq(i,(i + 9744)))], mapply(paste0, ped_mat[-c(seq(1,4))][rep(c(T,F), 2504)],ped_mat[-c(seq(1,4))][rep(c(F,T), 2504)]),stringsAsFactors = FALSE)
  ped_mat
}

#' Encodes genetic data based on reference allele
#' 
#' @param ped_mat Input data frame containing genetic information
#' @param which_file Index indicating the reference file
#' @return Encoded matrix with genetic data
encode <- function(ped_mat, which_file) {
  output <- matrix(nrow = 9745, ncol = 2504)
  for(i in seq(1,9745)) {
    allele_1 = ped_mat[i,5]
    allele_2 = ped_mat[i,6]
    hold<-as.character(gsub(paste(allele_1,allele_1,sep =''),'0',as.character(unname(ped_mat[i,c(7:2510)]))))
    hold<-gsub(paste(allele_1, allele_2,sep = ''),'1',hold)
    hold<-gsub(paste(allele_2, allele_1,sep = ''),'1',hold)
    hold<-gsub(pattern = paste(allele_2, allele_2,sep = ''),'2',hold)
    output[i,] <- as.numeric(hold)
  }
  output <- cbind(ped_mat[c(1:4)],output)
  output
}

#' Loads PLINK regression results from a specified path
#' 
#' @param path The directory path where PLINK regression results are located
#' @return A data frame containing the PLINK regression results
load_plink <- function(path) {
  plink_regression <- as.data.frame(fread(input = paste(path,'/plink.assoc.linear',sep='')))
  plink_regression
}

#' Runs linear regression models on genetic data for each SNP
#' 
#' This function reads in genetic and expression data, fixes data format issues, encodes SNPs,
#' and then runs linear regression models using R's built-in lm function.
#' 
#' @param files List of file paths containing genetic data
#' @param i Index indicating which file to process from the list
#' @param ped Genetic data frame
#' @param gene_expr Expression data to be regressed against genetic data
#' @return List of linear regression models for each SNP
run_lms <- function(files, i, gene_expr) {
  # First, read in and preprocess the genetic data
  ped <- read_in(files, i)
  ped <- fix_columns_and_augment(ped)
  ped <- encode(ped, i)
  
  # Run linear regression models for each SNP
  lms <- apply(ped[, c(seq(5, 2508))], 1, function(X) lm(gene_expr ~ X))
  return(list(lms, ped))
}

compare_to_plink_results <- function(path)
{
  #The output of the linear regression ran using plink:
  #plink --bfile snpdata --linear --pheno phenotype.txt
  #We load this in form comparison:
  plink_regression = load_plink(path)
  if(args[[2]])
  {
    head(plink_regression)
  }
  
  #Read in functional files
  snpdata_bim <- as.data.frame(fread(paste(path,'/snpdata.bim',sep='')))
  snpdata_phenotype <- read.table(paste(path,'/phenotype.txt',sep=''))
  colnames(snpdata_bim) <- c("Chromosome", "SNP_id", "Gen_distance", "Base_Pair_pos(bp)", "Allele_1", "Allele_2")
  #We will load the .bed using the files achieved by converting it to a .ped.
  #To speed this up, we split the files into 85 different .ped files.
  files <-paste(paste(path, '/files_ts/', sep=''),list.files(pattern ="*.txt", path = paste(path, '/files_ts', sep='')), sep ='')
  #Functions to help with analysis/data modification
  p_vals <- matrix(c(rep("A",9745*85),rep(1,9745*85)),nrow = 9745 * 85, ncol = 2)
  gene_expr <- snpdata_phenotype[,3]
  for(i in seq(1,85))
  {
    res = run_lms(files, i, gene_expr)
    lms = res[[1]]
    ped = res[[2]]
    p_vals[seq((((i - 1) * 9745) + 1), i*9745),2] <- unname(unlist(lapply(lms, function(X) summary(X)$coefficients[2,4])))
    p_vals[seq((((i - 1) * 9745) + 1), i*9745),1] <- unname(unlist(ped[,2]))
    print(paste(i,"done"))
    rm(ped)
  }
  return(list(p_vals))
}

#' Analyzes significance of SNPs per chromosome based on PLINK regression results
#' 
#' This function calculates the number of significant SNPs (P-value <= 5e-8) per chromosome
#' and generates visualizations and saves the results.
#' 
#' @param plink_regression Data frame containing PLINK regression results
#' @param path Directory path for saving the results
chromosome_significance <- function(plink_regression, path) {
  sig <- plink_regression$P <= 5e-8
  temp <- split(sig, plink_regression$CHR)
  len <- unname(unlist(lapply(temp, length)))
  temp <- unname(unlist(lapply(temp, sum)))
  temp <- data.frame("CHR" = seq(1, 22), "Significant SNPs" = temp, "Total SNPs" = len, "Ratio" = round(temp / len, 3))
  save(temp, file = paste(path, "/Turnbull_sig_per_chrom.R", sep = ''))
  png(paste(path, '/sig_per_chrom.png', sep = ''), height = 1000, width = 500)
  p <- tableGrob(temp)
  grid.arrange(p)
  dev.off()
}

#' Builds a Manhattan plot for PLINK regression results
#' 
#' This function generates a Manhattan plot for PLINK regression results and saves it as an image.
#' 
#' @param path Directory path for saving the Manhattan plot image
#' @param plink_regression Data frame containing PLINK regression results
build_manhattan <- function(path, plink_regression) {
  # Build Manhattan plot for PLINK results
  if (args[[4]]) {
    png(file = paste(path, '/manhattan.png', sep=''), width = 800, height = 480)
    p <- manhattan(plink_regression, col = colors()[c(600:604, 608:644)], ylim = c(0, 350), cex = 0.5)
    dev.off()
  }
}

#' Compares PLINK regression results with GWAS data
#' 
#' This function compares the minimum P-values per chromosome from PLINK regression 
#' results with corresponding data from a GWAS file and saves the comparison results.
#' 
#' @param plink_regression Data frame containing PLINK regression results
#' @param path Directory path for accessing files and saving results
compare_to_gwas <- function(plink_regression, path) {
  splits <- split(plink_regression, plink_regression$CHR)
  temp <- lapply(splits, function(x) x[x$P == min(x$P), ])
  min_per_chrom <- do.call(rbind.data.frame, temp)
  
  variants <- paste(min_per_chrom$CHR, min_per_chrom$BP, sep=':')
  rm(splits)
  rm(temp)
  
  # Load and compare with GWAS data
  gwas_ukBB <- as.data.frame(fread(paste(path, "/30760_irnt.gwas.imputed_v3.both_sexes.tsv", sep = '')))
  matches <- data.frame()
  for (i in variants) {
    matches <- rbind(matches, gwas_ukBB[which(startsWith(gwas_ukBB$variant, i)), c(1, 11)])
  }
  
  min_per_chrom <- cbind(min_per_chrom[, c(1, 2)], matches$variant, min_per_chrom$P)
  min_per_chrom$ukBB_P <- matches[, 2]
  names(min_per_chrom)[c(3, 4)] <- c("variant", "Plink_P")
  
  if (args[[5]]) {
    save(min_per_chrom, file = paste(path, "/Turnbull_Minimum_P.R", sep=''))
    png(paste(path, '/Minimum_P.png', sep=''), height=1000, width=500)
    p <- tableGrob(min_per_chrom)
    grid.arrange(p)
    dev.off()
  }
}