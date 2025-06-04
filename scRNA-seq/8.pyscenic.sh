#!/bin/bash

###bsub -oo scenic.log -q long -W 24:00 -n 16 -R rusage[mem=4G] -R "span[hosts=1]" "sh pyscenic.sh clusterX"
###clusterX would be the name of loom file would like to build the initial GRN by SCENIC and equal to arg1 here.
###download docker image from SCENIC website (https://hub.docker.com/r/aertslab/pyscenic/tags)
arg1=$1
cd /path/to/GRN/
singularity run pyscenic_0.12.1.sif pyscenic grn ${arg1}.loom TF_filter2.txt -o ${arg1}_adj1.csv --num_workers 12
singularity run pyscenic_0.12.1.sif pyscenic ctx ${arg1}_adj1.csv mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl --output ${arg1}_reg1.csv --num_workers 12 --mode dask_multiprocessing --mask_dropouts --expression_mtx_fname ${arg1}.loom
singularity run pyscenic_0.12.1.sif pyscenic aucell ${arg1}.loom ${arg1}_reg1.csv --output ${arg1}_auc1.csv --num_workers 12