## Given a multiVCF, count SNP differences between each pair of samples
## The difference counting I did for the Acineto paper only applied to biallelic sites
## Want to develop a better tool that counts SNPs at multiallelic sites
## Started 17 December 2015


###################################################
################# IMPORTANT NOTE ##################
###################################################

# The multiVCF that we feed this script
# must have "#CHROM" replaced with "CHROM".
# So run this: sed 's/#CHROM/CHROM/' in.vcf > out.vcf

###################################################
################# LOAD LIBARIES ###################
###################################################

library(stringr)

###################################################
################ DEFINE FUNCTIONS #################
###################################################

processer <- function(dataset) {

  ## REMOVE FIRST NINE COLUMNS FROM THE MULTIVCFs
  data <- dataset[-c(1:9)]
  genos <- as.data.frame(sapply(data, function(x) str_extract(x, "[0123456789]")))
  
  mat <- as.matrix(genos)
  obj <- NULL
  
  # Count pairwise differences
  # Alleles are coded as 0,1,2, or 3 depending on the nucleotide.
  # So find the differences first, but then convert to as.logical to just a T or F answer
  # Then convert back to numeric to turn all the trues into 1s
  # So essentially this just takes any number > 0 and turns it into a 1
  # Then sum that to identify the total number of pairwise differences.
  for (i in 1:ncol(mat)) {
    for (j in 1:ncol(mat)) {
      obj <- append(obj, sum(as.numeric(as.logical(as.numeric(mat[,i]) - as.numeric(mat[,j]))), na.rm = TRUE))
    } 
  }
  
  res <- matrix(obj, nrow = ncol(mat), ncol = ncol(mat)) # make a matrix out of it
  
  rownames(res) <- colnames(genos)
  colnames(res) <- colnames(genos)
  
  return(res)
}

###################################################
################## READ IN DATA ###################
###################################################

## READ IN THE MULTIVCF

st131 <- read.table("st131.biallelic.vcf", comment.char="#", header=TRUE)
nonst131 <- read.table("nonst131.biallelic.vcf", comment.char="#", header=TRUE)

results_st131 <- processer(st131)
results_nonst131 <- processer(nonst131)

write.table(results_st131, "results_biallelic_st131.txt", sep="\t")
write.table(results_nonst131, "results_biallelic_nonst131.txt", sep="\t")


