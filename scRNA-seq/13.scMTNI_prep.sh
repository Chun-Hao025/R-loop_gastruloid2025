#!/bin/bash

###processing subsample file especially split genes to make it processing parallely
###Gene chunks, eg 500, 1000, depends on cell number of each cell type for each condition.

BASE_DIR="/path/to/GRN/scMTNI/control/subsample"
OUT_DIR="/path/to/GRN/scMTNI/control/log/"
indir="/path/to/GRN/scMTNI/"
net_dir="/path/to/GRN/scMTNI/initial/" ###folder with files of TF-target for each cluster from 10.SCENIC_initial2.R
filelist="${indir}/cluster.txt"
regfile="${indir}/TF_filter2.txt"
for SUB in `cat /path/to/GRN/scMTNI/ID.txt`; do
        Input="${BASE_DIR}/${SUB}"
        mkdir -p "${indir}/scMTNI_control/${SUB}"
        Results="${indir}/scMTNI_control/${SUB}"

        # Define log files for each job
        LOG_FILE="${OUT_DIR}/${SUB}_prep1.log"
        #ERROR_FILE="${OUT_DIR}/${SUB}_prep.err"
        bsub -J "process_${SUB}" -oo "$LOG_FILE" -n 4 -R "rusage[mem=2G]" -q short -W 4:00 -R "span[hosts=1]"\
                "python /path/to/GRN/scMTNI/PreparescMTNIinputfiles.py --filelist $filelist --regfile $regfile --netdir $net_dir --indir $Input --outdir $Results --splitgene 500 --motifs True"

done

BASE_DIR="/path/to/GRN/scMTNI/RH_control/subsample"
OUT_DIR="/path/to/GRN/scMTNI/RH_control/log/"
indir="/path/to/GRN/scMTNI/"
net_dir="/path/to/GRN/scMTNI/prior/"  ### folder with files of TF-target for each cluster from all control with scMTNI pipeline (14-18 on all control data)
filelist="${indir}/cluster.txt"
regfile="${indir}/TF_filter2.txt"
for SUB in `cat /path/to/GRN/scMTNI/ID.txt`; do
        Input="${BASE_DIR}/${SUB}"
        mkdir -p "${indir}/scMTNI_RH_control/${SUB}"
        Results="${indir}/scMTNI_RH_control/${SUB}"

        # Define log files for each job
        LOG_FILE="${OUT_DIR}/${SUB}_prep1.log"
        #ERROR_FILE="${OUT_DIR}/${SUB}_prep.err"
        bsub -J "process_${SUB}" -oo "$LOG_FILE" -n 4 -R "rusage[mem=1G]" -q short -W 4:00 -R "span[hosts=1]"\
                "python /path/to/GRN/scMTNI/PreparescMTNIinputfiles.py --filelist $filelist --regfile $regfile --netdir $net_dir --indir $Input --outdir $Results --splitgene 1000 --motifs True"

done

BASE_DIR="/path/to/GRN/scMTNI/RH_Dox/subsample"
OUT_DIR="/path/to/GRN/scMTNI/RH_Dox/log/"
indir="/path/to/GRN/scMTNI/"
net_dir="/path/to/GRN/scMTNI/prior/" ### folder with files of TF-target for each cluster from all control with scMTNI pipeline (14-18 on all control data)
filelist="${indir}/cluster.txt"
regfile="${indir}/TF_filter2.txt"
for SUB in `cat /path/to/GRN/scMTNI/ID.txt`; do
        Input="${BASE_DIR}/${SUB}"
        mkdir -p "${indir}/scMTNI_RH_Dox/${SUB}"
        Results="${indir}/scMTNI_RH_Dox/${SUB}"

        # Define log files for each job
        LOG_FILE="${OUT_DIR}/${SUB}_prep1.log"
        #ERROR_FILE="${OUT_DIR}/${SUB}_prep.err"
        bsub -J "process_${SUB}" -oo "$LOG_FILE" -n 4 -R "rusage[mem=1G]" -q short -W 4:00 -R "span[hosts=1]"\
                "python /path/to/GRN/scMTNI/PreparescMTNIinputfiles.py --filelist $filelist --regfile $regfile --indir $Input --netdir $net_dir --outdir $Results --splitgene 1000 --motifs True"

done