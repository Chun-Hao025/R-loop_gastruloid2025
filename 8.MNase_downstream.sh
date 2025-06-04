#!/bin/bash

##pool data from all titers together (10U, 25U, 100U)
##ID=[control_rep1, control_rep2, Dox_rep1, Dox_rep2]
for i in `cat ID.txt`; do
    java -Xmx4g -jar /path/to/picard_2.10.9/picard.jar MergeSamFiles I=${i}_10U_mono.bam I=${i}_25U_mono.bam I=${i}_100U_mono.bam OUTPUT=${i}_mono.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
    java -Xmx4g -jar /path/to/picard_2.10.9/picard.jar MergeSamFiles I=${i}_10U_sub.bam I=${i}_25U_sub.bam I=${i}_100U_sub.bam OUTPUT=${i}_sub.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
    samtools index ${i}_mono.bam
    samtools index ${i}_sub.bam
    mkdir -p ./HOMER_${i}_mono/
    mkdir -p ./HOMER_${i}_sub/
    makeTagDirectory ./HOMER_${i}_mono/ ./${i}_mono.bam
    makeTagDirectory ./HOMER_${i}_sub/ ./${i}_sub.bam
    makeUCSCfile ./HOMER_${i}_mono/ -o auto
    makeUCSCfile ./HOMER_${i}_sub/ -o auto
    gzip -d ./HOMER_${i}_mono/HOMER_${i}_mono.ucsc.bedGraph.gz
    gzip -d ./HOMER_${i}_sub/HOMER_${i}_sub.ucsc.bedGraph.gz
    bedGraphToBigWig ./HOMER_${i}_mono/HOMER_${i}_mono.ucsc.bedGraph /path/to/mm10_size.txt ./${i}_mono_RPM.bw
    bedGraphToBigWig ./HOMER_${i}_sub/HOMER_${i}_sub.ucsc.bedGraph /path/to/mm10_size.txt ./${i}_sub_RPM.bw
    gzip ./HOMER_${i}_mono/HOMER_${i}_mono.ucsc.bedGraph
    gzip ./HOMER_${i}_sub/HOMER_${i}_sub.ucsc.bedGraph
done

##MNase-seq signal on different gene expression class
##class1, class2 and class3 were classified based on FPMK (see bulk-RNA-seq_downstream.R)
computeMatrix reference-point --numberOfProcessors max -S /list/to/all/bigwig_mono.bw ... -R ./class1.bed ./class2.bed ./class3.bed --referencePoint TSS -b 2000 -a 2000 --sortRegions keep --missingDataAsZero --skipZeros -o ./MN_pool_mono_expressionlevel.gz
computeMatrix reference-point --numberOfProcessors max -S /list/to/all/bigwig_sub.bw ... -R ./class1.bed ./class2.bed ./class3.bed --referencePoint TSS -b 2000 -a 2000 --sortRegions keep --missingDataAsZero --skipZeros -o ./MN_pool_sub_expressionlevel.gz
#*_expressionlevel_enrich_matrix.tab was used to plot in R 
plotProfile -m ./MN_pool_mono_expressionlevel.gz -out ./MN_mono_expressionlevel.pdf --numPlotsPerRow 3 --refPointLabel TSS --regionsLabel high medium low --outFileNameData ./MN_mono_expressionlevel_enrich_matrix.tab
plotProfile -m ./MN_pool_sub_expressionlevel.gz -out ./MN_sub_expressionlevel.pdf --numPlotsPerRow 3 --refPointLabel TSS --regionsLabel high medium low --outFileNameData ./MN_sub_expressionlevel_enrich_matrix.tab

##chromHMM aggregation plot
computeMatrix scale-regions --numberOfProcessors max -S /list/to/all/bigwig_mono.bw ... -R /list/to/all/chromHMM/state.bed ... --startLabel "start" --endLabel "end" -b 500 -a 500 --missingDataAsZero --skipZeros --smartLabels -p max/2 -o ./MN_mono_chromHMM_matrix.gz
computeMatrix scale-regions --numberOfProcessors max -S /list/to/all/bigwig_sub.bw ... -R /list/to/all/chromHMM/state.bed ... --startLabel "start" --endLabel "end" -b 500 -a 500 --missingDataAsZero --skipZeros --smartLabels -p max/2 -o ./MN_sub_chromHMM_matrix.gz
#*_matrix_enrich.tab was used to plot aggregation plot in R
plotProfile -m ./MN_mono_chromHMM_matrix.gz -out ./MN_mono_chromHMM.pdf --numPlotsPerRow 4 --perGroup --outFileNameData ./MN_mono_chromHMM_matrix_enrich.tab
plotProfile -m ./MN_sub_chromHMM_matrix.gz -out ./MN_sub_chromHMM.pdf --numPlotsPerRow 4 --perGroup --outFileNameData ./MN_sub_chromHMM_matrix_enrich.tab
