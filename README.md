This repository contains the scripts used for the following study:

# R-loops promote H2A.Z occupancy and proper differentiation in mouse embryonic stem cells
 
Chun-Hao Chao1 and Thomas G. Fazzio1,2*

1Department of Molecular, Cell, and Cancer Biology, University of Massachusetts Medical School, Worcester, MA 01605, USA.

The raw data can be downloaded at GEO (data will be available after publication)

For scRNA-seq, we recommand to download the processed seurat object, which contained cell annotation, pseudotime, UMAP, harmony corrected latent space, and all metadata across our data and two published data, from GEO to explore. 

Note:

Main Root Directory:

* Scripts provided in the main root directory cover steps from read alignment, data processing, peak calling, to visualization for CUT&Tag, CUT&RUN, MapR, MNase-seq, and ATAC-seq analyses.

RNA-seq Directory:

* Two scripts are available for bulk RNA-seq analysis.

* Three BED files correspond to gene classifications based on expression levels, as described in the manuscript.

chromHMM Directory:

* chrom_info.txt lists files used for ChromHMM model training.

* chrom_new.csv includes the feature probabilities and state annotations derived from the ChromHMM analysis.

* chromHMM BED Files: Ten BED files represent genomic regions corresponding to the ten distinct chromatin states identified through ChromHMM analysis.

scRNA-seq Directory:

* Eighteen scripts are provided for the complete analysis pipeline, including initial read alignment via CellRanger, clustering, annotation, trajectory inference, differential gene expression analysis, cell-cell communication (CCC) analysis, and gene regulatory network (GRN) analysis for scRNA-seq data.

GRN Analysis:

* Initial GRN construction employs the SCENIC workflow (scripts 7–10 within the scRNA-seq directory).

* The motif ranking database required for SCENIC can be directly downloaded from: CisTarget Motif Ranking Database.

scMTNI Analysis:

* Reference GRNs and trajectory information derived from control samples across all timepoints (located in the GRN subfolder) were utilized as inputs for scripts 11–18 and other required scripts in the GRN subfolder, constituting the first round of scMTNI analysis (allcontrol_edge.txt.zip, combining edges from all clusters).

* A second round of scMTNI analysis, using the combined GRNs from all control conditions as reference (RHKI_edge.txt.zip, combining edges from each cluster and treatment condition), was conducted to evaluate differences between control and treatment groups.

Subsampling Strategy for scMTNI:

* For each round of scMTNI analysis, 25 out of 50 subsampled datasets were used per beta1/beta2 parameter combination to reduce computational demands.

* Each parameter combination was repeated for a total of 75 trials across three different beta1/beta2 combinations.

* Only TF-target pairs consistently observed in ≥80% of the trials were retained as high-confidence edges.

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
```
R version 4.2.0 (2022-04-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 13.6.3

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] tidygraph_1.2.2             data.table_1.14.10          ggraph_2.1.0                CellChat_2.1.1              igraph_1.5.1                Matrix_1.6-4               
 [7] ggbreak_0.1.1               GO.db_3.16.0                biomaRt_2.54.1              cellmatch_0.0.0.9000        harmony_1.0.3               Rcpp_1.0.12                
[13] DoubletFinder_2.0.3         DEGreport_1.34.0            GenomicFeatures_1.50.4      lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1              
[19] purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1                tidyverse_2.0.0             reshape2_1.4.4             
[25] DESeq2_1.38.3               SummarizedExperiment_1.28.0 MatrixGenerics_1.10.0       GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         tximport_1.26.1            
[31] viridis_0.6.2               viridisLite_0.4.2           ComplexHeatmap_2.14.0       circlize_0.4.15             matrixStats_1.1.0           RColorBrewer_1.1-3         
[37] plyr_1.8.9                  ChIPseeker_1.34.1           clusterProfiler_4.6.0       org.Mm.eg.db_3.16.0         AnnotationDbi_1.60.2        IRanges_2.32.0             
[43] S4Vectors_0.36.2            Biobase_2.58.0              BiocGenerics_0.44.0         ggrepel_0.9.5               scCustomize_2.1.2           ggplot2_3.4.4              
[49] patchwork_1.2.0             dplyr_1.1.4                 Seurat_5.0.1                SeuratObject_5.0.1          sp_1.6-0                   

loaded via a namespace (and not attached):
  [1] statnet.common_4.9.0                    svglite_2.1.3                           ica_1.0-3                               Rsamtools_2.14.0                       
  [5] foreach_1.5.2                           lmtest_0.9-40                           crayon_1.5.2                            MASS_7.3-58.2                          
  [9] nlme_3.1-157                            backports_1.4.1                         GOSemSim_2.24.0                         rlang_1.1.3                            
 [13] XVector_0.38.0                          HDO.db_0.99.1                           ROCR_1.0-11                             irlba_2.3.5.1                          
 [17] limma_3.54.2                            filelock_1.0.2                          BiocParallel_1.32.6                     rjson_0.2.21                           
 [21] bit64_4.0.5                             glue_1.7.0                              rngtools_1.5.2                          sctransform_0.4.1                      
 [25] parallel_4.2.0                          vipor_0.4.5                             spatstat.sparse_3.0-0                   dotCall64_1.0-2                        
 [29] DOSE_3.24.2                             spatstat.geom_3.0-6                     tidyselect_1.2.0                        fitdistrplus_1.1-8                     
 [33] XML_3.99-0.14                           zoo_1.8-11                              ggpubr_0.6.0                            GenomicAlignments_1.34.0               
 [37] ggnetwork_0.5.12                        xtable_1.8-4                            RcppHNSW_0.4.1                          magrittr_2.0.3                         
 [41] cli_3.6.2                               zlibbioc_1.44.0                         rstudioapi_0.15.0                       miniUI_0.1.1.1                         
 [45] bslib_0.6.1                             fastmatch_1.1-3                         treeio_1.22.0                           fastDummies_1.7.3                      
 [49] shiny_1.8.0                             xfun_0.41                               clue_0.3-64                             gson_0.0.9                             
 [53] cluster_2.1.4                           caTools_1.18.2                          KEGGREST_1.38.0                         logging_0.10-108                       
 [57] ape_5.6-2                               listenv_0.9.0                           Biostrings_2.66.0                       png_0.1-8                              
 [61] reshape_0.8.9                           future_1.33.1                           withr_3.0.0                             bitops_1.0-7                           
 [65] ggforce_0.4.1                           coda_0.19-4                             pillar_1.9.0                            gplots_3.1.3                           
 [69] GlobalOptions_0.1.2                     cachem_1.0.8                            GetoptLong_1.0.5                        paletteer_1.5.0                        
 [73] vctrs_0.6.5                             ellipsis_0.3.2                          generics_0.1.3                          NMF_0.26                               
 [77] tools_4.2.0                             beeswarm_0.4.0                          munsell_0.5.0                           tweenr_2.0.2                           
 [81] fgsea_1.24.0                            DelayedArray_0.24.0                     fastmap_1.1.1                           compiler_4.2.0                         
 [85] abind_1.4-5                             httpuv_1.6.13                           rtracklayer_1.58.0                      TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
 [89] plotly_4.10.4                           GenomeInfoDbData_1.2.9                  gridExtra_2.3                           edgeR_3.40.2                           
 [93] lattice_0.20-45                         deldir_1.0-6                            utf8_1.2.4                              later_1.3.2                            
 [97] BiocFileCache_2.15.1                    prismatic_1.1.1                         jsonlite_1.8.8                          scales_1.3.0                           
[101] carData_3.0-5                           tidytree_0.4.2                          pbapply_1.7-2                           genefilter_1.80.3                      
[105] lazyeval_0.2.2                          promises_1.2.1                          car_3.1-2                               doParallel_1.0.17                      
[109] goftest_1.2-3                           sna_2.7-2                               spatstat.utils_3.1-3                    reticulate_1.34.0                      
[113] cowplot_1.1.3                           Rtsne_0.16                              downloader_0.4                          uwot_0.1.14                            
[117] survival_3.5-0                          yaml_2.3.8                              plotrix_3.8-2                           systemfonts_1.0.5                      
[121] htmltools_0.5.7                         memoise_2.0.1                           BiocIO_1.8.0                            locfit_1.5-9.7                         
[125] graphlayouts_0.8.4                      digest_0.6.34                           mime_0.12                               rappdirs_0.3.3                         
[129] registry_0.5-1                          spam_2.9-1                              RSQLite_2.3.1                           yulab.utils_0.0.6                      
[133] future.apply_1.11.1                     blob_1.2.4                              splines_4.2.0                           labeling_0.4.3                         
[137] rematch2_2.1.2                          RCurl_1.98-1.12                         broom_1.0.5                             hms_1.1.3                              
[141] colorspace_2.1-0                        ConsensusClusterPlus_1.62.0             BiocManager_1.30.22                     mnormt_2.1.1                           
[145] ggbeeswarm_0.7.1                        shape_1.4.6                             aplot_0.1.9                             ggrastr_1.0.1                          
[149] sass_0.4.8                              RANN_2.6.1                              enrichplot_1.18.3                       fansi_1.0.6                            
[153] tzdb_0.4.0                              parallelly_1.36.0                       R6_2.5.1                                ggridges_0.5.4                         
[157] lifecycle_1.0.4                         ggsignif_0.6.4                          curl_5.2.0                              jquerylib_0.1.4                        
[161] leiden_0.4.3                            snakecase_0.11.0                        qvalue_2.30.0                           RcppAnnoy_0.0.20                       
[165] iterators_1.0.14                        spatstat.explore_3.0-6                  htmlwidgets_1.6.4                       polyclip_1.10-4                        
[169] network_1.18.2                          shadowtext_0.1.2                        timechange_0.3.0                        gridGraphics_0.5-1                     
[173] globals_0.16.2                          spatstat.random_3.1-3                   progressr_0.13.0                        codetools_0.2-19                       
[177] FNN_1.1.4                               gtools_3.9.4                            prettyunits_1.2.0                       psych_2.4.1                            
[181] dbplyr_2.4.0                            gridBase_0.4-7                          RSpectra_0.16-1                         gtable_0.3.4                           
[185] DBI_1.2.1                               ggalluvial_0.12.5                       ggfun_0.0.9                             tensor_1.5                             
[189] httr_1.4.7                              KernSmooth_2.23-20                      presto_1.0.0                            stringi_1.8.3                          
[193] progress_1.2.3                          farver_2.1.1                            annotate_1.76.0                         ggtree_3.6.2                           
[197] xml2_1.3.6                              ggdendro_0.1.23                         boot_1.3-28.1                           BiocNeighbors_1.16.0                   
[201] restfulr_0.0.15                         geneplotter_1.76.0                      ggplotify_0.1.0                         scattermore_1.2                        
[205] bit_4.0.5                               scatterpie_0.1.8                        spatstat.data_3.0-0                     janitor_2.1.0                          
[209] pkgconfig_2.0.3                         ggprism_1.0.4                           rstatix_0.7.2                           knitr_1.45

```
Other package (check the link following these packages to install all dependent packages)
```
anndata 0.10.9
bowtie2 2.5.0
bedtools 2.30.0
cellranger 6.1.2
ChromHMM 1.24 (https://compbio.mit.edu/ChromHMM/)
deepTools 3.5.1 (https://deeptools.readthedocs.io/en/3.0.2/content/installation.html)
hnswlib 0.7.0
homer 4.11 (http://homer.ucsd.edu/homer/introduction/install.html)
Margaret (https://github.com/Zafar-Lab/Margaret)
matplotlib 3.9.2
networkx 3.3
numpy 1.26.4
pandas 2.2.3
picard 2.10.9
pyscenic 0.12.1 (https://github.com/aertslab/pySCENIC)
python 3.10.14
pyvia 0.2.4 (https://pyvia.readthedocs.io/en/latest/Installation.html)
rsem 1.3.3
samtools 1.16.1
scanpy 1.10.2 (https://scanpy.readthedocs.io/en/stable/installation.html)
scikit-learn 1.5.1  
scipy 1.11.4
scmtni (https://github.com/Roy-lab/scMTNI)
star 2.7.10a
trimmomatic 0.39
ucsc-bedgraphtobigwig 472
umap-learn 0.5.6
```
