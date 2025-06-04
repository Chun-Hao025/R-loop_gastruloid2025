#!/bin/bash

#check peak folder for peak files
## -d list all tag directory for each marks
##H4ac, CTCF, H3.3, MapR 
annotatePeaks.pl tss mm10 -size 4000 -hist 20 -d /path/to/tagdirectory/control_rep1/ /path/to/tagdirectory/control_rep2/ /path/to/tagdirectory/Dox_rep1/ /path/to/tagdirectory/Dox_rep2/ > hitone_mark_TSS.txt

#MapR peak call
getDifferentialPeaksReplicates.pl -t /path/to/tagdirectory/MapR/replicate/ ... -b /path/to/tagdirectory/MapR_RNaseA/replicate/ ... -i /path/to/tagdirectory/pA-MNase/replicate/ ... -genome mm10 -style factor -tagThreshold 40 -region > ./MapR_Peaks2.txt

##RNAP2S2P, H3K36me3
##To profile meta profile from TSS to TTS
makeMetaGeneProfile.pl rna mm10 -d /path/to/tagdirectory/control_rep1/ /path/to/tagdirectory/control_rep2/ /path/to/tagdirectory/Dox_rep1/ /path/to/tagdirectory/Dox_rep2/ ... > histone_marker_meta.txt

##H3K4me3, H3K9me3, H3K27me3, H3K27ac, ZX1_H2A.Z
annotatePeaks.pl /path/to/peakfile.bed mm10 -size 4000 -hist 20 -d /path/to/tagdirectory/control_rep1/ /path/to/tagdirectory/control_rep2/ /path/to/tagdirectory/Dox_rep1/ /path/to/tagdirectory/Dox_rep2/ ... > hitone_mark_peak.txt

#H2A.Z
annotatePeaks.pl tss mm10 -size 4000 -hist 20 -d /path/to/tagdirectory/control_rep1/ /path/to/tagdirectory/control_rep2/ /path/to/tagdirectory/Dox_rep1/ /path/to/tagdirectory/Dox_rep2/ ... > hitone_mark_TSS.txt
annotatePeaks.pl /path/to/peakfile.bed mm10 -size 4000 -hist 20 -d /path/to/tagdirectory/control_rep1/ /path/to/tagdirectory/control_rep2/ /path/to/tagdirectory/Dox_rep1/ /path/to/tagdirectory/Dox_rep2/ ... > hitone_mark_peak.txt
annotatePeaks.pl /path/to/DHS.bed mm10 -size 4000 -hist 20 -d /path/to/tagdirectory/control_rep1/ /path/to/tagdirectory/control_rep2/ /path/to/tagdirectory/Dox_rep1/ /path/to/tagdirectory/Dox_rep2/ ... > hitone_mark_DHS.txt
##peakcalling from control and Dox separately
##this is based on DESeq2 method to identify peaks from replicates
getDifferentialPeaksReplicates.pl -t /list/all/tag/directory/for/H2AZ/control/ ... -i /list/all/tag/directory/for/IgG/control/and/Dox/ ... -genome mm10 -style histone -tagThreshold 40 > N_H2AZ_peaks.txt
getDifferentialPeaksReplicates.pl -t /list/all/tag/directory/for/H2AZ/Dox/ ... -i /list/all/tag/directory/for/IgG/control/and/Dox/ ... -genome mm10 -style histone -tagThreshold 40 > D_H2AZ_peaks.txt

pos2bed.pl N_H2AZ_Peaks_input.txt > N_H2AZ_peak_input.bed
pos2bed.pl D_H2AZ_Peaks_input.txt > D_H2AZ_peak_input.bed

cat N_H2AZ_Peaks.bed D_H2AZ_Peaks.bed > H2AZ_peak.bed
sort -k1,1 -k2,2n H2AZ_peak.bed > H2AZ_sort.bed
bedrolls merge -i H2AZ_sort.bed -d 2000 > H2AZ_bedtool.bed

#differential peak analysis within the peak file identified above
getDifferentialPeaksReplicates.pl -t /list/all/tag/directory/for/H2AZ/control/ ... -b /list/all/tag/directory/for/H2AZ/Dox/ ... -p H2AZ_bedtool.bed > H2AZ_down2.txt
getDifferentialPeaksReplicates.pl -t /list/all/tag/directory/for/H2AZ/Dox/ ... -b /list/all/tag/directory/for/H2AZ/control/ ... -p H2AZ_bedtool.bed > H2AZ_UP2.txt

pos2bed.pl H2AZ_down2.txt > H2AZ_down2.bed
pos2bed,pl H2AZ_UP2.txt > H2AZ_UP2.bed
mergePeaks -d 2000 H2AZ_bedtool.bed H2AZ_down2.bed H2AZ_UP2.bed -prefix overlap

pos2bed.pl overlap_H2AZ_bedtool.bed > H2AZ_share.bed
pos2bed.pl overlap_H2AZ_bedtool.bed_H2AZ_down2.bed.bed > H2AZ_down3.bed
pos2bed.pl overlap_H2AZ_bedtool.bed_H2AZ_UP2.bed > H2AZ_UP3.bed

#heatmap
computeMatrix reference-point --numberOfProcessors max -S /list/representative/bigwig.bw ... -R ./peak/H2AZ_down3.bed ./peak/H2AZ_UP3.bed --referencePoint center -b 3000 -a 3000 --missingDataAsZero --skipZeros -o ./H2AZ_difpeak_center.gz 
plotHeatmap -m ./H2AZ_difpeak_center.gz -out ./H2AZ_heatmap.pdf --colorMap inferno --whatToShow 'heatmap and colorbar' --zMin auto --zMax auto --heatmapHeight 15 --refPointLabel "peak center"

#genome metadata analysis and the adjancent gene to the differential peaks
#these files would be used to do GO analysis and plot in R 
annotatePeaks.pl /path/to/differential/H2AZ_down3.bed mm10 -annStats > ./down_annotation.txt
annotatePeaks.pl /path/to/differential/H2AZ_UP3.bed mm10 -annStats > ./up_annotation.txt

#combine output from -annStats of UP peaks and down peaks to get meta_analysis.csv


##chromHMM aggregation plot
computeMatrix scale-regions --numberOfProcessors max -S /list/all/bigwig.bw ... -R /list/to/all/chromHMM/state.bed ... --startLabel "start" --endLabel "end" -b 500 -a 500 --missingDataAsZero --skipZeros --smartLabels -p max/2 -o ./H2AZ_chromHMM_matrix.gz
#H2AZ_chromHMM_matrix_enrich.tab was used to plot aggregation plot in R
plotProfile -m ./H2AZ_chromHMM_matrix.gz -out ./H2AZ_chromHMM_new.pdf --numPlotsPerRow 4 --perGroup --outFileNameData ./H2AZ_chromHMM_matrix_enrich.tab

