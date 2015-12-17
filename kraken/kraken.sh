## Script to run Kraken to identify species of our genomes
## 13 December 2015

reads=/proj/julianog/users/ChristianP/esblWGS/symlinks
krakenwd=/proj/julianog/users/ChristianP/esblWGS/kraken

#for i in {1..9}
#do
#echo ESBL0$i\_R1.fq.gz
#kraken \
#	--preload \
#	--db /proj/julianog/src/kraken-0.10.4-beta/minikraken_20141208/ \
#	--fastq-input \
#	--gzip-compressed $reads/ESBL0$i\_R1.fq.gz \
#	> $krakenwd/results/ESBL0$i.kraken
#done

#for i in {10..71}
#do
#echo ESBL$i\_R1.fq.gz
#kraken \
#	--preload \
#	--db /proj/julianog/src/kraken-0.10.4-beta/minikraken_20141208/ \
#	--fastq-input \
#	--gzip-compressed $reads/ESBL$i\_R1.fq.gz \
#	> $krakenwd/results/ESBL$i.kraken
#done

#./proj/julianog/src/kraken-0.10.4-beta/scripts/kraken-translate \
#	--db /proj/julianog/src/kraken-0.10.4-beta/minikraken_20141208/ \
#	/proj/julianog/users/ChristianP/esblWGS/kraken/results/ESBL01.kraken > ESBL01.labels


## what percent of reads are E coli?
for i in {1..9}
do
echo ESBL0$i
wc -l $krakenwd/results/ESBL0$i.kraken
cut -f3 $krakenwd/results/ESBL0$i.kraken | grep -c "585056\|562\|543" # ecoli ids
echo "\n"
done

for i in {10..71}
do
echo ESBL$i
wc -l $krakenwd/results/ESBL$i.kraken
cut -f3 $krakenwd/results/ESBL$i.kraken | grep -c "585056\|562\|543"
echo "\n"
done


