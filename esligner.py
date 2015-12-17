## Likely have to realign for LUMPY variant detection
## SVTyper is having major trouble with the headers of my
## PICARD & GATK pipeline BAMs
## So realign all using bwa mem and samtools

workdir: '/proj/julianog/users/ChristianP/esblWGS'
REF = '/proj/julianog/refs/Ec985_HG941718/Ec958_HG941718.fa'
readWD = '/proj/julianog/users/ChristianP/esblWGS/'
SAMPLES, = glob_wildcards('/proj/julianog/users/ChristianP/esblWGS/symlinks/{sample}_R1.fq.gz')
TMPDIR = '/netscr/prchrist/tmp_for_picard/'
PICARD = '/nas02/apps/picard-1.88/picard-tools-1.88'
GATK = '/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar'


####### Target #######
rule all:
	#input: expand('aln/{sample}.bam', sample = SAMPLES)
	#input: expand('aln/{sample}.sorted.bam', sample = SAMPLES)
	#input: expand('aln/{sample}.dedup.bam', sample = SAMPLES)
	#input: expand('aln/{sample}.dedup.bai', sample = SAMPLES)
	#input: expand('aln/{sample}.realn.bam', sample = SAMPLES)
	#input: 'coverage/cov_plot.pdf'
	#input: 'variants/esbl_UG.vcf'
	input: 'variants/esbl_UG.pass.vcf'

rule select_variants :
	input: 'variants/esbl_UG.qual.vcf'
	output: 'variants/esbl_UG.pass.vcf'
	shell: 'java -jar {GATK} \
		-T SelectVariants \
		-R {REF} -V {input} -o {output} \
		-select "vc.isNotFiltered()" \
		-restrictAllelesTo BIALLELIC'
			# keeps only unfiltered sites

rule filter_variants :
	input: vcf = 'variants/esbl_UG.vcf'
	output: 'variants/esbl_UG.qual.vcf'
	shell: 'java -jar {GATK} \
		-T VariantFiltration \
		-R {REF} \
		-V {input.vcf} \
		--filterExpression "QD < 25.0" \
		--filterName "QD" \
		--filterExpression "MQ < 59.0" \
		--filterName "MQ" \
		--filterExpression "FS > 4.0" \
		--filterName "FS" \
		--filterExpression "MQRankSum < -4.0" \
		--filterName "MQRankSum" \
		--filterExpression "ReadPosRankSum < -4.0" \
		--filterName "ReadPosRankSum" \
		--logging_level ERROR \
		-o {output}'

rule filter_min_depth :
	input: 'variants/{list}_UG.vcf'
	output: 'intervals/{list}_UG_05xAT100%.intervals'
	shell: 'java -Xmx2g -jar {GATK} \
		-T CoveredByNSamplesSites \
		-R {REF} -V {input} -out {output} \
		-minCov 05 -percentage 0.99999'
		 #Output interval file contains sites that passed
		 #Would be more elegant to use 1.0 instad of 0.99999, but that doesn't work

rule unified_genotyper :
	input: 'names/esbl_5x@75%.list'
	output: 'variants/esbl_UG.vcf'
	shell: 'java -jar {GATK} -T UnifiedGenotyper \
		-R {REF} -I {input} -nt 8 \
		-ploidy 1 -o {output}'

rule plot_coverage:
	input: 'coverage/covPlotter.r'
	output: 'coverage/cov_plot.pdf'
	shell: 'Rscript {input}'

rule make_R_cov_script:
	input: 'coverage/coverage.txt'
	output: 'coverage/covPlotter.r'
	shell: """echo '## Read in data
		coverage <- read.table("coverage/coverage.txt", header=FALSE)

		## Order data acendingly
		coverage <- coverage[with(coverage, order(V2)), ]

		## Add a third column with numbers in it
		coverage$V6 <- 1:length(coverage$V1)

		##Plot the coverage graph
		pdf(file="coverage/cov_plot.pdf", width = 9, height = 5)
		plot(coverage$V2 ~ coverage$V6, axes=FALSE, xlab="", ylab="Frac Genome Covered", ylim=c(0,1), col="black", pch=20, type = "n")
		points(coverage$V5 ~ coverage$V6, col="grey", pch=20)
		points(coverage$V2 ~ coverage$V6, col="blue", pch=20)
		points(coverage$V3 ~ coverage$V6, col="red", pch=20)
		points(coverage$V4 ~ coverage$V6, col="black", pch=20)
		axis(1, at=1:length(coverage$V1), labels=coverage$V1, cex.axis=.7, las=3, cex = 0.5)
		axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), las = 2)
		legend(30,0.4, legend=c("5x Coverage", "10x Coverage", "25x Coverage", "50x Coverage"), col=c("blue", "red", "black", "grey"), pch=20)
		dev.off()' > coverage/covPlotter.r
		"""

rule digest_coverage:
	input: expand('coverage/data/{sample}.{cov}', sample = SAMPLES, cov = 'cov05 cov10 cov25 cov50'.split())
	output: 'coverage/coverage.txt'
	shell: 'for name in `ls coverage/data/ | grep cov05 | sed "s/\.cov..//"`; \
		do \
		cov05=$(tail -1 coverage/data/$name.cov05 | cut -f 5); \
		cov10=$(tail -1 coverage/data/$name.cov10 | cut -f 5); \
		cov25=$(tail -1 coverage/data/$name.cov25 | cut -f 5); \
		cov50=$(tail -1 coverage/data/$name.cov50 | cut -f 5); \
		echo -e $name"\t"$cov05"\t"$cov10"\t"$cov25"\t"$cov50 >> coverage/coverage.txt; \
		done'

rule calculate_50x_cov:
	input: 'aln/{sample}.realn.bam'
	output: 'coverage/data/{sample}.cov50'
	shell: 'bedtools genomecov \
		-ibam {input} -max 50 | grep genome \
		> {output}'

rule calculate_25x_cov:
	input: 'aln/{sample}.realn.bam'
	output: 'coverage/data/{sample}.cov25'
	shell: 'bedtools genomecov \
		-ibam {input} -max 25 | grep genome \
		> {output}'

rule calculate_10x_cov:
	input: 'aln/{sample}.realn.bam'
	output: 'coverage/data/{sample}.cov10'
	shell: 'bedtools genomecov \
		-ibam {input} -max 10 | grep genome \
		> {output}'

rule calculate_05x_cov:
	input: 'aln/{sample}.realn.bam'
	output: 'coverage/data/{sample}.cov05'
	shell: 'bedtools genomecov \
		-ibam {input} -max 5 | grep genome \
		> {output}'

rule realn_indels:
	input: bam = 'aln/{sample}.dedup.bam', targets = 'aln/{sample}.realigner.intervals', 
	output: 'aln/{sample}.realn.bam'
	shell: 'java -jar {GATK} -T IndelRealigner \
		-R {REF} -I {input.bam} \
		-targetIntervals {input.targets} \
		-o {output}' 

rule find_indels:
	input: bam = 'aln/{sample}.dedup.bam', index = 'aln/{sample}.dedup.bai'
	output: 'aln/{sample}.realigner.intervals'
	shell: 'java -jar {GATK} -T RealignerTargetCreator \
		-R {REF} -I {input.bam} \
		-o {output}'

rule index_merged: 
	input: 'aln/{sample}.dedup.bam'
	output: 'aln/{sample}.dedup.bai'
	shell: 'java -jar {PICARD}/BuildBamIndex.jar INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'

rule mark_dups:
	input: 'aln/{sample}.sorted.bam'
	output:'aln/{sample}.dedup.bam','aln/{sample}.dedup.metrics'
	shell: 'java -jar {PICARD}/MarkDuplicates.jar \
		I={input} O={output[0]} \
		METRICS_FILE={output[1]} \
		TMP_DIR={TMPDIR} REMOVE_DUPLICATES=TRUE \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000'

rule sort_bam:
	input: 'aln/{sample}.bam'
	output: 'aln/{sample}.sorted.bam'
	shell: 'java -jar {PICARD}/SortSam.jar \
		INPUT={input} OUTPUT={output} SO=coordinate'

rule fastq_to_bam:
	input: 'symlinks/{sample}_R1.fq.gz', 'symlinks/{sample}_R2.fq.gz'
	output: 'aln/{sample}.bam'
	shell: 'bwa mem {REF} {readWD}{input[0]} {readWD}{input[1]} \
			-M -R "@RG\tID:bwa\tSM:{wildcards.sample}\tLB:lib1\tPL:illumina"|\
			samtools view -Sb - > {output}'
