#!/bin/bash

#check peak folder for peak files
## -d list all tag directory for each marks
annotatePeaks.pl tss mm10 -size 4000 -hist 20 -d /path/to/tagdirectory/control_rep1/ /path/to/tagdirectory/control_rep2/ /path/to/tagdirectory/Dox_rep1/ /path/to/tagdirectory/Dox_rep2/ ... > ATAC_TSS.txt
annotatePeaks.pl /path/to/DHS.bed mm10 -size 4000 -hist 20 -d /path/to/tagdirectory/control_rep1/ /path/to/tagdirectory/control_rep2/ /path/to/tagdirectory/Dox_rep1/ /path/to/tagdirectory/Dox_rep2/ ... > ATAC_DHS.txt

##ATAC signal on different gene expression class
##class1, class2 and class3 were classified based on FPMK (see bulk-RNA-seq_downstream.R)
computeMatrix reference-point --numberOfProcessors max -S /list/all/bigwig.bw ... -R ./class1.bed ./class2.bed ./class3.bed --referencePoint TSS -b 2000 -a 2000 --sortRegions keep --missingDataAsZero --skipZeros -o ./ATAC_expressionlevel.gz
#ATAC_expressionlevel_enrich_matrix.tab was used to plot in R 
plotProfile -m ./ATAC_expressionlevel.gz -out ./ATAC_expressionlevel.pdf --numPlotsPerRow 3 --refPointLabel TSS --regionsLabel high medium low --outFileNameData ./ATAC_expressionlevel_enrich_matrix.tab

##identify peaks
findPeaks /path/to/Tag/directory/control_rep1_sub/ -style factor > ./control_rep1_peak.txt
findPeaks /path/to/Tag/directory/control_rep2_sub -style factor > ./control_rep2_peak.txt
findPeaks /path/to/Tag/directory/control_rep3_sub -style factor > ./control_rep3_peak.txt
findPeaks /path/to/Tag/directory/control_rep4_sub -style factor > ./control_rep4_peak.txt
findPeaks /path/to/Tag/directory/Dox_rep1_sub -style factor > ./Dox_rep1_peak.txt
findPeaks /path/to/Tag/directory/Dox_rep2_sub -style factor > ./Dox_rep2_peak.txt
findPeaks /path/to/Tag/directory/Dox_rep3_sub -style factor > ./Dox_rep3_peak.txt
findPeaks /path/to/Tag/directory/Dox_rep4_sub -style factor > ./Dox_rep4_peak.txt

pos2bed.pl control_rep1_peak.txt | sort -k1,1 -k2,2n > control_rep1_peak.bed
pos2bed.pl control_rep2_peak.txt | sort -k1,1 -k2,2n > control_rep1_peak.bed
pos2bed.pl control_rep3_peak.txt | sort -k1,1 -k2,2n > control_rep1_peak.bed
pos2bed.pl control_rep4_peak.txt | sort -k1,1 -k2,2n > control_rep1_peak.bed
pos2bed.pl Dox_rep1_peak.txt | sort -k1,1 -k2,2n > Dox_rep1_peak.txt.bed
pos2bed.pl Dox_rep2_peak.txt | sort -k1,1 -k2,2n > Dox_rep2_peak.txt.bed
pos2bed.pl Dox_rep3_peak.txt | sort -k1,1 -k2,2n > Dox_rep3_peak.txt.bed
pos2bed.pl Dox_rep4_peak.txt | sort -k1,1 -k2,2n > Dox_rep4_peak.txt.bed

mergePeaks -d given control_rep1_peak.bed control_rep2_peak.bed control_rep3_peak.bed control_rep4_peak.bed -prefix control -matrix control -venn control > N_ATAC_overlap.bed
mergePeaks -d given Dox_rep1_peak.txt.bed Dox_rep2_peak.txt.bed Dox_rep3_peak.txt.bed Dox_rep4_peak.txt.bed -prefix Dox -matrix Dox -venn Dox > Dox_overlap.bed
##only keep peak that overlap across all replicates in each condition

cat control_ATAC_overlap.bed Dox_ATAC_overlap.bed | sort -k1,1 -k2,2n | bedtools merge -d 200 | bed2pos.pl > ATAC_overlap_merge.txt

annotatePeaks.pl ATAC_overlap_merge.txt mm10 -d /list/all/tag/directory/for/ATAC/control/ /list/all/tag/directory/for/ATAC/Dox/ ... -raw > ATAC_rawcount_peak.txt
getDiffExpression.pl ATAC_rawcount_peak.txt control control control control Dox Dox Dox Dox -peak -rlog > ATAC_diff.txt