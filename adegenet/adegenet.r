# For genetic analysis of our ESBL surveillance data
# Started 16 Dec 2015
# Basics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf
# Genomics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf
# Extra Commands: http://www.inside-r.org/packages/cran/adegenet/docs/.rmspaces


######################
### LOAD LIBRARIES ###
######################

library(adegenet)
library(pegas)
library(RColorBrewer)
#display.brewer.all()

########################
### DEFINE FUNCTIONS ###
########################

## Function to create genlight from VCF.
genlight.maker <- function(infile) {
  loci <- read.vcf(infile, from = 1, to = 1000000) # reads first million sites
  genlight <- new("genlight", loci) # convert data frame into genlight object
  ploidy(genlight) <- as.integer(1) # add back population information
  return(genlight)
}


### READ IN THE DATA ###
gl <- genlight.maker("../variants/esbl_UG.pass.vcf")
indNames(gl) # check to see if indnames are in there
ploidy(gl) # check to make sure it's haploid

st131 <- as.numeric(indNames(gl) %in% c("ESBL06", "ESBL04", "ESBL27", "ESBL70", "ESBL02", "ESBL05", "ESBL10", "ESBL64", "ESBL16", "ESBL58", "ESBL33", "ESBL39", "ESBL68", "ESBL31", "ESBL61", "ESBL21", "ESBL63", "ESBL51", "ESBL09", "ESBL32", "ESBL29", "ESBL13", "ESBL26", "ESBL24", "ESBL03", "ESBL23", "ESBL36", "ESBL71", "ESBL38", "ESBL25", "ESBL60", "ESBL62", "ESBL37", "ESBL49", "ESBL42", "ESBL18", "ESBL59", "ESBL34", "ESBL07"))*131
st405 <- as.numeric(indNames(gl) %in% c("ESBL17", "ESBL44", "ESBL45", "ESBL57", "ESBL65", "ESBL67"))*405
st648 <- as.numeric(indNames(gl) %in% c("ESBL35", "ESBL41", "ESBL46"))*648
st1284 <- as.numeric(indNames(gl) %in% c("ESBL40", "ESBL66", "ESBL20"))*1284
st224 <- as.numeric(indNames(gl) %in% c("ESBL50", "ESBL53"))*224
st12 <- as.numeric(indNames(gl) %in% c("ESBL54", "ESBL55"))*12
st38 <- as.numeric(indNames(gl) %in% c("ESBL12", "ESBL48"))*38
st69 <- as.numeric(indNames(gl) %in% c("ESBL01", "ESBL30"))*69

pops <- st131 + st405 + st648 + st1284 + st224 + st12 + st38 + st69
pops <- pops + 1

### DO THE PCA CALCULATIONS ###
pca <- glPca(gl) # for genlight

### PLOT PCA EIGENVALUES ###
barplot(pca1$eig, xlab = "", ylab = "Variance", main = "P. vivax Eigenvalues") # barplot

### ADD JITTER ###
pca_jit <- jitter(pca$scores, factor = 300)

### PLOT PCA PICTURE ###
palette(brewer.pal(length(levels(as.factor(pops))), "Set1"))

svg("pca.svg", width=5, height=5.5)
tiff("pca.tiff", width=5, height=5.5, res=300, unit="in", compression="lzw")
plot(pca_jit[,1], pca_jit[,2], 
     pch=19,
     cex=1,
     col=as.factor(pops),
     axes=FALSE, 
     xlab=paste("PC1 - ", round(pca$eig[1]/sum(pca$eig)*100), "% of total variance", sep = ""),
     ylab=paste("PC2 - ", round(pca$eig[2]/sum(pca$eig)*100), "% of total variance", sep = ""),
     main=""
)
axis(1)
axis(2)
legend(-95, 152, 
       legend = c("Unique", "ST-12", "ST-38", "ST-69", "ST-131", "ST-224", "ST-405", "ST-648", "ST-1284"), 
       col = palette(), pch = 16, ncol = 3, cex = 0.85, pt.cex=1.5, bty = "n")
dev.off()