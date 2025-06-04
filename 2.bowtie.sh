#!/bin/bash
##use bsub on HPC for this bash file
##bsub -oo output.log -q long -W 24:00 -n 8 -R rusage[mem=8G] -R span[hosts=1] "sh alignment.sh index"
##index is the name of individual sample for the analysis and equals to arg1 in this bash script
##conda environment install all required packages including bowtie2, trimmomatic, samtools, picard, HOMER, deeptools

arg1="$1"

cd /path/to/sample/fastq/folder/${arg1}/
## For CUT&RUN (TruSeq adapter)
java -Xmx32g -jar /path/to/trimmomatic.jar PE ./${arg1}_all_1.fq.gz ./${arg1}_all_2.fq.gz ${arg1}_trim_1_paired.fq.gz ${arg1}_trim_1_unpaired.fq.gz ${arg1}_trim_2_paired.fq.gz ${arg1}_trim_2_unpaired.fq.gz ILLUMINACLIP:/path/to/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:25:7:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 LEADING:3 TRAILING:3 MINLEN:10
## For CUT&Tag (Nextera adapter)
java -Xmx32g -jar /path/to/trimmomatic.jar PE ./${arg1}_all_1.fq.gz ./${arg1}_all_2.fq.gz ${arg1}_trim_1_paired.fq.gz ${arg1}_trim_1_unpaired.fq.gz ${arg1}_trim_2_paired.fq.gz ${arg1}_trim_2_unpaired.fq.gz ILLUMINACLIP:/path/to/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:25:7:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 LEADING:3 TRAILING:3 MINLEN:10
#the percentage mapping is in ${arg1}_map.txt.
bowtie2 -q --threads 8 -N 1 -X 1000 -x /path/to/bowtie_gemone/mm10 -1 ./${arg1}_trim_1_paired.fq.gz -2 ./${arg1}_trim_2_paired.fq.gz -S ./${arg1}.bowtie.sam &> ${arg1}_map.txt
java -Xmx4g -jar /path/to/picard.jar SortSam INPUT=${arg1}.bowtie.sam OUTPUT=${arg1}_sort.sam VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp SORT_ORDER=coordinate
java -Xmx4g -jar /path/to/picard_2.10.9/picard.jar MarkDuplicates INPUT=${arg1}_sort.sam OUTPUT=${arg1}_dup.sam VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp METRICS_FILE=dup.txt REMOVE_DUPLICATES=true
# flag 99 and 147 is the forward unique mapping read pairs and flag 83 and 163 is the reverse unique mapping read pairs
# This filter is also optional since HOMER would automatically keep unique mapped reads for the analysis
samtools view -h -f 99 ${arg1}_dup.sam > ${arg1}_C99.sam
samtools view -h -f 147 ${arg1}_dup.sam > ${arg1}_C147.sam
samtools view -h -f 83 ${arg1}_dup.sam > ${arg1}_C83.sam
samtools view -h -f 163 ${arg1}_dup.sam > ${arg1}_C163.sam
java -Xmx4g -jar /pi/thomas.fazzio-umw/Tong/Tools/picard_2.10.9/picard.jar MergeSamFiles I=${arg1}_C83.sam I=${arg1}_C99.sam I=${arg1}_C147.sam I=${arg1}_C163.sam OUTPUT=${arg1}_split.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
#filter MAPQ quality larger than 10
## All CUT&Tag, CUT&RUN and MapR would used ${arg1}_filter.bam for downstream analysis
samtools view -bq 10 ./${arg1}_split.bam > ./${arg1}_filter.bam
samtools index ./${i}_filter.bam
## For ATAC-seq analysis, excluded mitochondria reads and subnucleosome fragments
samtools view --no-PG -h ./${arg1}_filter.bam | grep -vwE "(chrM)" | samtools view -S -T /path/to/GRCm38.primary_assembly.genome.fa -b -o ./${arg1}_noMT.bam
samtools view --no-PG -h ./${arg1}_noMT.bam | awk ' /^@/ { next } ($9 <= 120 && $9 >= 1)||($9 >= -120 && $9 <= -1)' | samtools view -S -T /path/to/GRCm38.primary_assembly.genome.fa -b -o ./${arg1}_sub.bam
samtools index ./${arg1}_sub.bam
## For MNase-seq analysis, spliting read into subnucleosome fragment and mononucleosome fragment.
samtools view --no-PG -h ./${arg1}_filter.bam | awk ' /^@/ { next } ($9 <= 120 && $9 >= 1)||($9 >= -120 && $9 <= -1)' | samtools view -S -T /path/to/GRCm38.primary_assembly.genome.fa -b -o ./${arg1}_sub.bam
samtools view --no-PG -h ./${arg1}_filter.bam | awk ' /^@/ { next } ($9 <= 200 && $9 >= 150)||($9 >= -200 && $9 <= -150)' | samtools view -S -T /path/to/GRCm38.primary_assembly.genome.fa -b -o ./${arg1}_mono.bam
samtools index ./${arg1}_sub.bam
samtools index ./${arg1}_mono.bam

## All CUT&Tag, CUT&RUN and MapR
mkdir -p ./HOMER_${arg1}/
## HOMER used RPM (read per 10 million for the normalization)
makeTagDirectory ./HOMER_${arg1} ./${arg1}_filter.bam
## generate bedgraph file
makeUCSCfile ./HOMER_${arg1}/ -o auto

## To generate bigwig file, either converting from bedgraph file or using bamCoverage from deeptool to generate.
#option1:
gzip -d ./HOMER_${arg1}/HOMER_${arg1}.ucsc.bedGraph.gz
##mm10_size.txt contains chromosome name and size information
bedGraphToBigWig ./HOMER_${arg1}/HOMER_${arg1}.ucsc.bedGraph /path/to/mm10_size.txt ./${arg1}_RPM.bw
gzip ./HOMER_${arg1}/HOMER_${arg1}.ucsc.bedGraph
#option2: RPGC for normalization
bamCoverage --bam ./${i}_filter.bam -o ./${i}_filter_RG.bw --binSize 10 -p max --normalizeUsing RPGC --effectiveGenomeSize 2494787188 --ignoreForNormalization chrX chrM --extendReads

## ATAC-seq
mkdir -p ./HOMER_${arg1}_sub/
## HOMER used RPM (read per 10 million for the normalization)
makeTagDirectory ./HOMER_${arg1}_sub ./${arg1}_sub.bam
## generate bedgraph file
makeUCSCfile ./HOMER_${arg1}_sub/ -o auto
## To generate bigwig file, either converting from bedgraph file or using bamCoverage from deeptool to generate.
#option1:
gzip -d ./HOMER_${arg1}_sub/HOMER_${arg1}_sub.ucsc.bedGraph.gz
##mm10_size.txt contains chromosome name and size information
bedGraphToBigWig ./HOMER_${arg1}_sub/HOMER_${arg1}_sub.ucsc.bedGraph /path/to/mm10_size.txt ./${arg1}_sub_RPM.bw
gzip ./HOMER_${arg1}_sub/HOMER_${arg1}_sub.ucsc.bedGraph
#option2: RPGC for normalization
bamCoverage --bam ./${arg1}_sub.bam -o ./${i}_sub_RG.bw --binSize 10 -p max --normalizeUsing RPGC --effectiveGenomeSize 2494787188 --ignoreForNormalization chrX chrM --extendReads
