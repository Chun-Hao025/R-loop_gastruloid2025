#!/bin/bash
##https://github.com/Roy-lab/scMTNI
ml scmtni/1.0
cd /path/to/GRN/scMTNI/

BASE_DIR="/path/to/GRN/scMTNI/control/subsample"
OUT_DIR="/path/to/GRN/scMTNI/control/log"
indir="/path/to/GRN/scMTNI/"

##For each condition, all control, RHKI control, and RHKI Dox, ID.txt contain file directory of subsampling file
##25 subsamples for one combination
##(-p 0.2, beta1 and beta2: -b -2, -q 2; -b -4, -q 2; -b -4, -q 4) 
for i in `cat ./ID.txt`; do
    Input="${BASE_DIR}/${i}"
    Results="${indir}/scMTNI_control/${i}"

    for file in ${Input}/testdata_config_motifs*.txt; do
        idx=$(basename "${file}" | grep -o '[0-9]\+')
        bsub -J "scMTNI_${i}_${idx}" -oo "${OUT_DIR}/scMTNI_${i}_${idx}.log" -n 2 -R "rusage[mem=1G]" -q short -W 8:00 -R "span[hosts=1]"\
        "scMTNI -f "${file}" -x50 -l "${Input}/TFs_OGs.txt" -n "${Input}/ogids/AllGenes${idx}.txt" -d "${indir}/lineage_tree.txt" -m "${Input}/testdata_ogids.txt" -s "${Input}/celltype_order.txt" -p 0.2 -c yes -b -4 -q 4"
        echo "Submitting job for file: ${file}, idx: ${idx}"
    done
done