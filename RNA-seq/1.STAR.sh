#!/bin/bash

ml rsem/1.3.3
ml star/2.7.10a
ml samtools/1.16.1

cd /path/to/reference/genome/build/up/
#to build up the reference genome for STAR
rsem-prepare-reference --gtf gencode.vM25.annotation.gtf --star --star-path /share/pkg/star/2.7.10a/ -p 8 GRCm38.primary_assembly.genome.fa ./mm10/mm10_STAR/

#ID.txt contains all fastq file names
for i in `cat ID.txt`; do
    # each folder contains read1 and read2 fastq file as *_1.fq.gz and *_2.fq.gz
	cd /path/to/fastq/file/${i}/
    # trimming fastq file is optional
	java -Xmx4g -jar /path/to/trimmomatic.jar PE ./${i}*_1.fq.gz ./${i}*_2.fq.gz ${i}_trim_1_paired.fq.gz ${i}_trim_1_unpaired.fq.gz ${i}_trim_2_paired.fq.gz ${i}_trim_2_unpaired.fq.gz ILLUMINACLIP:/path/to//trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:25:7:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 LEADING:3 TRAILING:3 MINLEN:10
	#--quantMode TranscriptomeSAM would be critical for using rsem-calculate-expression to get raw counts, FPMK ...etc.
    STAR --runThreadN 8 --genomeDir /path/to/reference/genome/build/up/mm10 --sjdbGTFfile /path/to/reference/genome/build/up/gencode.vM25.annotation.gtf --readFilesIn ./${i}_trim_1_paired.fq.gz ./${i}_trim_2_paired.fq.gz --readFilesCommand zcat --outFileNamePrefix ./ -outFilterType BySJout --outMultimapperOrder Random --outFilterMultimapNmax 200 --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMstrandField None --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts --outReadsUnmapped Fastx
	rsem-calculate-expression --alignments --bam --strandedness reverse --paired-end -p 16 --no-bam-output ./Aligned.toTranscriptome.out.bam /path/to/reference/genome/build/up/mm10/mm10_STAR/ ./${i}.RSEM
    #two output file with ${i}.RSEM.genes.results and ${i}.RSEM.isoforms.results in the folder
    #use ${i}.RSEM.genes.results for the downstream DESeq2 analysis
done
