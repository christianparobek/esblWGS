## A de novo assembler for Hajime's ESBL/CRE reads
## Started 23 December 2015
## Make sure and rename CRE15 -> CRE015, for example




### Shuffle the reads; required input for docker
#for i in `ls | grep -v shuf | grep R1 | cut -f1 -d "_"`
#do

#echo $i

#r1=`ls | grep $i | grep R1`
#r2=`ls | grep $i | grep R2`

#echo $r1
#echo $r2

#paste $r1 $r2 | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' > $i.shuf.fq

#done


## Assemble the reads

for i in `ls shortreads | cut -f1 -d "."`
do

echo $i


./assemble koadman/docker-a5-miseq default shortreads/$i.short.shuf.fq.gz aln/
	# runs the assemble executable file from nucleotid.es
	# which in turn runs the a5-miseq docker
	# which had by far the best performance for de novo assembly
	# in my hands, and according to the online benchmarks

mv aln/contigs.fa aln/$i.fa #rename output fasta
mv aln/log.txt aln/$i.log # rename log

done
