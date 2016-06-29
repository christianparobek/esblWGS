## Quick script to blast Oxford virulence genes
## Against our ESBL denovo assemblies
## Started 28 June 2016

for sample in `ls ../denovo/aln/ | grep fa | sed 's/.fa//'`
do

blastn -subject ../denovo/aln/$sample.fa -query virulence_genes.fa -perc_identity 80 -qcov_hsp_perc 80 -outfmt "6 std qcovs" > res/$sample\_blast.txt

done
