#!/bin/bash

###chrom_info.txt contains all bam file of feature from control data using for chromHMM training
##10 state defined
java -Xmx4g -jar /path/to/ChromHMM.jar BinarizeBam -paired /path/to/mm10_size.txt /path/of/main/directory/for/bamfile /path/to/chrom_info.txt /path/to/binary/

##emissons_10.txt contains the features probability of each state from chromHMM and control_10_overlap.txt is the genomic meta data of each states
##used these two file for heatmap plotting in R
##combine to get state_statistics.csv
java -Xmx64g -jar /path/to/ChromHMM.jar LearnModel -p 0 -l /path/to/mm10_size.txt -s 1000 /path/to/binary/ /path/to/chrom_10/ 10 mm10

##split different state into different bed file for the aggregation plot of MNase-seq and H2A.Z data
cd /path/to/chrom_10/
for label in $(cut -f4 control_10_segments.bed | sort -u); do
  awk -F'\t' -v L="$label" '$4==L' control_10_segments.bed > "${label}.bed"
done