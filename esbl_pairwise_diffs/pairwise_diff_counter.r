## Christian Parobek
## December 9, 2016

## We want to calculate number of SNP diffs
## between ESBL isolates, excluding recombinant regions
## that were identified by Gubbins.
## Use the FASTA output from Gubbins and 
## analyze using Adegenet.


#####################################
####### LOAD USEFUL LIBRARIES #######
#####################################

library(adegenet)
library(ape)
library(stats)


#####################################
###### DEFINE USEFUL FUNCTION #######
#####################################

diff.calc <- function(in_fasta, out_table) {
  bin <- fasta2DNAbin(in_fasta)
  dist <- dist.dna(bin, model = "N", as.matrix = TRUE, pairwise.deletion = FALSE)
    ## pairwise.deletion = FALSE is the default, and how MEGA makes it's trees too
  write.table(dist, file = out_table, sep = "\t")
}


#####################################
######## CALCULATE SNP DIFFS ########
#####################################

diff.calc("gubbins_output_fastas/esbl_all_multifasta.filtered_polymorphic_sites.fasta", "diff_tables/esbl_all.csv")
diff.calc("gubbins_output_fastas/st131multi.filtered_polymorphic_sites.fasta", "diff_tables/st131.csv")
diff.calc("gubbins_output_fastas/nonst131multi.filtered_polymorphic_sites.fasta", "diff_tables/nonst131.csv")


diff.calc("gubbins_output_fastas/esbl_all_multifasta.filtered_polymorphic_sites.fasta", "diff_tables/pairwise_del_true.csv")
diff.calc("gubbins_output_fastas/esbl_all_multifasta.filtered_polymorphic_sites.fasta", "diff_tables/pairwise_del_false.csv")
