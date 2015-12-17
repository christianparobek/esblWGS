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

## Given a multiVCF, get the genotypes
genoGet <- function(dataset) {
  
  ## REMOVE FIRST NINE COLUMNS FROM THE MULTIVCFs
  data <- dataset[-c(1:9)]
  geno <- as.data.frame(sapply(data, function(x) str_extract(x, "[0123456789]")))
  
}


###################################################
################## READ IN DATA ###################
###################################################

## READ IN THE Pf AND Pv MULTIVCFs
pf <- read.table("esbl_UG.vcf", comment.char="#", header=TRUE)

genos <- genoGet(pf)
mat <- as.matrix(genos)
class(mat)
genos$ESBL01 - genos$ESBL02

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

obj[order(obj)] # find lowest numer of SNPs

