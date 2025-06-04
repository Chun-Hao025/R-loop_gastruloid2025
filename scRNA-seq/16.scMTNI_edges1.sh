#!/bin/bash

estEdge="/path/to/GRN/scMTNI/Scripts/Postprocessing/estimateedgeconf/estimateEdgeConf"  ### the function in the scMTNI package (https://github.com/Roy-lab/scMTNI)
##The output directory of the scMNTI analysis result for all control, RHKI control or RHKI Dox
indir="/path/to/GRN/scMTNI/scMTNI_control"
#indir="/path/to/GRN/scMTNI/scMTNI_RH_control"
#indir="/path/to/GRN/scMTNI/scMTNI_RH_Dox"

###cluster.txt list all cluster/celltype involved in scMTNI analysis
for i in $(cat ${indir}/cluster.txt); do
	echo ${i}
	mkdir -p ${indir}/new/${i}
	> "${indir}/new/${i}/network_files.txt"
	### To create a list for all result file from all random sample+beta1/beta2 combination for each cluster. Here, 25 random samples result of each beta1/beta2 combination and with total 75 trial
	for rseed in $(cat /path/to/GRN/scMTNI/ID.txt); do
		find ${indir}/${rseed}/ -maxdepth 1 -name "${i}*_var_mb_pw_k50_N1.txt" | grep -E "/${i}[^0-9][^/]*_var_mb_pw_k50_N1.txt$" >> ${indir}/new/${i}/network_files.txt 2>/dev
/null
	done
	for rseed in $(cat /path/to/GRN/scMTNI/ID1.txt); do
		find ${indir}/${rseed}/ -maxdepth 1 -name "${i}*_var_mb_pw_k50_N2.txt" | grep -E "/${i}[^0-9][^/]*_var_mb_pw_k50_N2.txt$" >> ${indir}/new/${i}/network_files.txt 2>/dev
/null
	done
	###combine all TF-target pairs and calculate the fraction of TF-target pair showing in all 75 trials.
	if [[ -s "${indir}/new/${i}/network_files.txt" ]]; then
		${estEdge} "${indir}/new/${i}/network_files.txt" 0 "${indir}/new/${i}/net_" alledges
	fi
done