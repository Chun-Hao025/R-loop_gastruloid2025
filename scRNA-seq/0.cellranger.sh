#!/bin/bash

ml cellranger/6.1.2
### $1 is the folder name for each scRNA-seq samples containing R1, R2, I1, I2.fq.gz
cellranger count --include-introns=True --expect-cells=10000 --id=$1_Results --transcriptome=/path/to/mm10 --fastqs=./$1/ --sample=$1 --localcores=14 --localmem=64