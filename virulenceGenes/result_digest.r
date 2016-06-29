## Script to digest reuslts of blasting
## Oxford virulence genes blasted against ESBL de novo assemblies
## This distills those results into a binary presence/absence table
## Christian Parobek
## Started 28 June 2016


##############################
####### LOAD LIBRARIES #######
##############################

library(stringr)


##############################
######## READ IN DATA ########
##############################

blastfiles = list.files(path="res", full.names = TRUE)
  # get paths

myfiles = lapply(blastfiles, read.table)
  # read in tables into a list


##############################
######## HOUSEKEEPING ########
##############################

sampNames <- str_extract(blastfiles, "P*ESBL[0-9]+")
  # extract sample names

for (i in 1:length(sampNames)){
  myfiles[[i]]$Sample <- sampNames[i]
} # add sample names to each table within the liat

megatable <- do.call("rbind", myfiles)
  # turn the list of tables into one big table

summarytable <- as.data.frame(with(megatable, table(Sample, V1)) > 0L) + 0L
  # make a binary summary table

summarytable["Total" ,] <- colSums(summarytable)
summarytable[, "Total"] <- rowSums(summarytable)
  # add column and row sums
  # gives an idea of how virulent samples are and how common given factors are

write.table(x = summarytable, file = "final_blast_results.txt", sep = "\t", quote = FALSE)

##############################
##### FIND NO-MATCH GENES ####
##############################

virNames <- as.character(read.table("virulence_names")$V1)

summarytable <- summarytable[row.names(summarytable) %in% virNames,]
  # get rid of that wierd "nothing" that came from the ESBL28 file

virNames[!virNames %in% row.names(summarytable)]
  # get the virulence genes which had no matches