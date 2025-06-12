#!/bin/bash
#BSUB -n 5  
#BSUB -R rusage[mem=4G] 
#BSUB -W 24:00 
#BSUB -J sub[1-5]
#BSUB -q long   # which queue we want to run in
#BSUB -oo sub_%J_%I.log # log
#BSUB -R "span[hosts=1]" # All hosts on the same chassis"


##all control
indir=/path/to/GRN/scMTNI/
Input=${indir}/control/
filelist=${indir}/cluster.txt

python /path/to/data_subsample.py --filelist $filelist --indir $Input --nseed 10 --fraction 0.5 --start_seed $(((LSB_JOBINDEX-1)*10))


##RHKI control
Input=${indir}/RH_control/
filelist=${indir}/cluster.txt

python /path/to/data_subsample.py --filelist $filelist --indir $Input --nseed 10 --fraction 0.5 --start_seed $(((LSB_JOBINDEX-1)*10))

##RHKI Dox
Input=${indir}/RH_Dox/
filelist=${indir}/cluster.txt

python /path/to/data_subsample.py --filelist $filelist --indir $Input --nseed 10 --fraction 0.5 --start_seed $(((LSB_JOBINDEX-1)*10))
