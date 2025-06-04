This repository contains the scripts used for the following study:

# R-loops promote H2A.Z occupancy and proper differentiation in mouse embryonic stem cells
 
Chun-Hao Chao1 and Thomas G. Fazzio1,2*
1Department of Molecular, Cell, and Cancer Biology, University of Massachusetts Medical School, Worcester, MA 01605, USA.

The raw data can be downloaded at (data will be available after publication)

Note:

Main Root Directory:

Scripts provided in the main root directory cover steps from read alignment, data processing, peak calling, to visualization for CUT&Tag, CUT&RUN, MapR, MNase-seq, and ATAC-seq analyses.

RNA-seq Directory:

Two scripts are available for bulk RNA-seq analysis.

Three BED files correspond to gene classifications based on expression levels, as described in the manuscript.

chromHMM Directory:

chrom_info.txt lists files used for ChromHMM model training.

chrom_new.csv includes the feature probabilities and state annotations derived from the ChromHMM analysis.

chromHMM BED Files:

Ten BED files represent genomic regions corresponding to the ten distinct chromatin states identified through ChromHMM analysis.

scRNA-seq Directory:

Eighteen scripts are provided for the complete analysis pipeline, including initial read alignment via CellRanger, clustering, annotation, trajectory inference, differential gene expression analysis, cell-cell communication (CCC) analysis, and gene regulatory network (GRN) analysis for scRNA-seq data.

GRN Analysis:

Initial GRN construction employs the SCENIC workflow (scripts 7–10 within the scRNA-seq directory).

The motif ranking database required for SCENIC can be directly downloaded from: CisTarget Motif Ranking Database.

scMTNI Analysis:

Reference GRNs and trajectory information derived from control samples across all timepoints (located in the GRN subfolder) were utilized as inputs for scripts 11–18 and other required scripts in the GRN subfolder, constituting the first round of scMTNI analysis (allcontrol_edge.txt.zip, combining edges from all clusters).

A second round of scMTNI analysis, using the combined GRNs from all control conditions as reference (RHKI_edge.txt.zip, combining edges from each cluster and treatment condition), was conducted to evaluate differences between control and treatment groups.

Subsampling Strategy for scMTNI:

For each round of scMTNI analysis, 25 out of 50 subsampled datasets were used per beta1/beta2 parameter combination to reduce computational demands.

Each parameter combination was repeated for a total of 75 trials across three different beta1/beta2 combinations.

Only TF-target pairs consistently observed in ≥80% of the trials were retained as high-confidence edges.

```
├── RNA-seq/
│   ├── [2 scripts for bulk RNA-seq analysis]
│   └── [3 BED files for gene classification]
│
├── chromHMM/
│   ├── chrom_info.txt          # File list for ChromHMM training
│   ├── chrom_new.csv           # Feature probability and state annotation from ChromHMM
│   └── [10 BED files representing chromatin states]
│
├── scRNA-seq/
│   ├── [18 scripts for scRNA-seq analyses: alignment, clustering, annotation, trajectory, DE, CCC, GRN]
│   └── GRN/
│       ├── allcontrol_edge.txt.zip   # Combined edges from all control clusters
│       ├── cluster.txt               # cluster/celltype for the analysis
│       ├── TF_filter2.txt            # TF used for scMTNI and SCENIC analysis
│       ├── lineage_tree.txt          # lineage tree for scMTNI analysis
│       ├── data_subsample.py         # subsample script used in 12.scMTNI_subsample.sh
│       ├── PreparescMTNIinputfiles.py # the script for scMTNI input file prep used in 13.scMTNI_prep.sh
│       └── RHKI_edge.txt.zip         # Combined edges from each cluster and treatment condition
│
├── 1.bowtie_genome.sh          # Genome indexing for alignment
├── 2.bowtie.sh                 # Read alignment
├── 3.chromHMM.sh               # ChromHMM pipeline execution
├── 4.CT_downstream1.sh         # CUT&Tag downstream processing
├── 5.CT_plotting.R             # CUT&Tag result visualization
├── 6.ATAC_downstream.sh        # ATAC-seq downstream processing
├── 7.ATAC_plotting.R           # ATAC-seq result visualization
├── 8.MNase_downstream.sh       # MNase-seq downstream processing
├── 9.MNase-seq_plotting.R      # MNase-seq result visualization
├── LICENSE
└── README.md
```

