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
conda enviroment on HPC
```
# packages in environment at /pi/thomas.fazzio-umw/conda/env:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
_r-mutex                  1.0.1               anacondar_1    conda-forge
absl-py                   2.0.0                    pypi_0    pypi
alabaster                 0.7.13                   pypi_0    pypi
alevin-fry                0.8.2                h4ac6f70_0    bioconda
alsa-lib                  1.1.5             h516909a_1002    conda-forge
anndata                   0.8.0                    pypi_0    pypi
annoy                     1.17.3                   pypi_0    pypi
argcomplete               1.12.3                   pypi_0    pypi
argparse                  1.4.0                    pypi_0    pypi
asciitree                 0.3.3                    pypi_0    pypi
astunparse                1.6.3                    pypi_0    pypi
asynctest                 0.13.0                   pypi_0    pypi
babel                     2.14.0                   pypi_0    pypi
backcall                  0.2.0              pyh9f0ad1d_0    conda-forge
backports                 1.0                pyhd8ed1ab_3    conda-forge
backports.functools_lru_cache 1.6.4              pyhd8ed1ab_0    conda-forge
backports.zoneinfo        0.2.1            py37h540881e_5    conda-forge
bamtools                  2.5.2                hdcf5f25_5    bioconda
bedops                    2.4.41               h4ac6f70_1    bioconda
bedtools                  2.30.0               h468198e_3    bioconda
binutils_impl_linux-64    2.39                 he00db2b_1    conda-forge
binutils_linux-64         2.39                h5fc0e48_13    conda-forge
bioconductor-annotate     1.72.0            r41hdfd78af_0    bioconda
bioconductor-annotationdbi 1.56.2            r41hdfd78af_0    bioconda
bioconductor-annotationfilter 1.18.0            r41hdfd78af_0    bioconda
bioconductor-apeglm       1.16.0            r41hc247a5b_2    bioconda
bioconductor-biobase      2.54.0            r41hc0cfd56_2    bioconda
bioconductor-biocfilecache 2.2.0             r41hdfd78af_0    bioconda
bioconductor-biocgenerics 0.40.0            r41hdfd78af_0    bioconda
bioconductor-biocio       1.4.0             r41hdfd78af_0    bioconda
bioconductor-biocparallel 1.28.3            r41hc247a5b_1    bioconda
bioconductor-biomart      2.50.0            r41hdfd78af_0    bioconda
bioconductor-biostrings   2.62.0            r41hc0cfd56_2    bioconda
bioconductor-data-packages 20221112             hdfd78af_0    bioconda
bioconductor-delayedarray 0.20.0            r41hc0cfd56_2    bioconda
bioconductor-deseq2       1.34.0            r41hc247a5b_3    bioconda
bioconductor-ensembldb    2.18.1            r41hdfd78af_0    bioconda
bioconductor-genefilter   1.76.0            r41h38f54d8_2    bioconda
bioconductor-geneplotter  1.72.0            r41hdfd78af_0    bioconda
bioconductor-genomeinfodb 1.30.1            r41hdfd78af_0    bioconda
bioconductor-genomeinfodbdata 1.2.7             r41hdfd78af_2    bioconda
bioconductor-genomicalignments 1.30.0            r41hc0cfd56_2    bioconda
bioconductor-genomicfeatures 1.46.1            r41hdfd78af_0    bioconda
bioconductor-genomicranges 1.46.1            r41hc0cfd56_1    bioconda
bioconductor-iranges      2.28.0            r41hc0cfd56_2    bioconda
bioconductor-keggrest     1.34.0            r41hdfd78af_0    bioconda
bioconductor-matrixgenerics 1.6.0             r41hdfd78af_0    bioconda
bioconductor-protgenerics 1.26.0            r41hdfd78af_0    bioconda
bioconductor-rhtslib      1.26.0            r41hc0cfd56_2    bioconda
bioconductor-rsamtools    2.10.0            r41hc247a5b_2    bioconda
bioconductor-rtracklayer  1.54.0            r41hd029910_0    bioconda
bioconductor-s4vectors    0.32.4            r41hc0cfd56_0    bioconda
bioconductor-summarizedexperiment 1.24.0            r41hdfd78af_0    bioconda
bioconductor-xvector      0.34.0            r41hc0cfd56_2    bioconda
bioconductor-zlibbioc     1.40.0            r41hc0cfd56_2    bioconda
bioframe                  0.7.2                    pypi_0    pypi
biopython                 1.81                     pypi_0    pypi
bismark                   0.24.0               hdfd78af_0    bioconda
bleach                    6.0.0              pyhd8ed1ab_0    conda-forge
blosc                     1.21.4               h0f2a231_0    conda-forge
bokeh                     2.4.3                    pypi_0    pypi
boltons                   24.0.0                   pypi_0    pypi
boost-cpp                 1.74.0               h6cacc03_7    conda-forge
bowtie2                   2.5.0            py37hb24965f_0    bioconda
brotli                    1.0.9                h166bdaf_9    conda-forge
brotli-bin                1.0.9                h166bdaf_9    conda-forge
brunsli                   0.1                  h9c3ff4c_0    conda-forge
bwidget                   1.9.14               ha770c72_1    conda-forge
bzip2                     1.0.8                h7f98852_4    conda-forge
c-ares                    1.18.1               h7f98852_0    conda-forge
c-blosc2                  2.10.2               hb4ffafa_0    conda-forge
ca-certificates           2024.8.30            hbcca054_0    conda-forge
cached-property           1.5.2                hd8ed1ab_1    conda-forge
cached_property           1.5.2              pyha770c72_1    conda-forge
cachetools                5.3.1                    pypi_0    pypi
cairo                     1.16.0            ha12eb4b_1010    conda-forge
cellrank                  1.5.1                    pypi_0    pypi
certifi                   2024.8.30          pyhd8ed1ab_0    conda-forge
cffi                      1.15.1           py37h43b0acd_1    conda-forge
cfitsio                   4.1.0                hd9d235c_0    conda-forge
charls                    2.2.0                h9c3ff4c_0    conda-forge
charset-normalizer        3.2.0                    pypi_0    pypi
click                     8.1.3            py37h89c1867_0    conda-forge
cloudpickle               2.2.1                    pypi_0    pypi
cmake                     3.27.5                   pypi_0    pypi
coloredlogs               15.0.1                   pypi_0    pypi
colormath                 3.0.0                    pypi_0    pypi
comm                      0.2.2              pyhd8ed1ab_0    conda-forge
cooler                    0.9.3                    pypi_0    pypi
cramino                   0.14.5               h5076881_0    bioconda
cryptography              40.0.2                   pypi_0    pypi
curl                      7.86.0               h2283fc2_1    conda-forge
cutadapt                  4.4              py37h8902056_0    bioconda
cycler                    0.11.0             pyhd8ed1ab_0    conda-forge
cykhash                   2.0.1                    pypi_0    pypi
cython                    3.0.2                    pypi_0    pypi
cytoolz                   0.12.3                   pypi_0    pypi
debugpy                   1.7.0                    pypi_0    pypi
decorator                 5.1.1              pyhd8ed1ab_0    conda-forge
deeptools                 3.5.1                      py_0    bioconda
deeptoolsintervals        0.1.9            py37h8902056_4    bioconda
defusedxml                0.7.1              pyhd8ed1ab_0    conda-forge
dill                      0.3.7                    pypi_0    pypi
dnaio                     0.10.0           py37h8902056_1    bioconda
docrep                    0.3.2                    pypi_0    pypi
docutils                  0.19                     pypi_0    pypi
dxpy                      0.363.0                  pypi_0    pypi
editorconfig              0.12.3             pyhd8ed1ab_0    conda-forge
entrypoints               0.4                pyhd8ed1ab_0    conda-forge
et-xmlfile                1.1.0                    pypi_0    pypi
expat                     2.5.0                hcb278e6_1    conda-forge
fastqc                    0.11.9               hdfd78af_1    bioconda
fcsparser                 0.2.4                    pypi_0    pypi
fftw                      3.3.10          mpi_openmpi_h4a81ba8_8    conda-forge
flatbuffers               23.5.26                  pypi_0    pypi
font-ttf-dejavu-sans-mono 2.37                 hab24e00_0    conda-forge
font-ttf-inconsolata      3.000                h77eed37_0    conda-forge
font-ttf-source-code-pro  2.038                h77eed37_0    conda-forge
font-ttf-ubuntu           0.83                 hab24e00_0    conda-forge
fontconfig                2.14.1               hef1e5e3_0  
fonts-conda-ecosystem     1                             0    conda-forge
fonts-conda-forge         1                             0    conda-forge
freetype                  2.12.1               hca18f0e_1    conda-forge
freexl                    1.0.6                h166bdaf_1    conda-forge
fribidi                   1.0.10               h36c2ea0_0    conda-forge
frozendict                2.4.6                    pypi_0    pypi
frozenlist                1.3.3                    pypi_0    pypi
fsspec                    2023.1.0                 pypi_0    pypi
future                    0.18.3                   pypi_0    pypi
gast                      0.4.0                    pypi_0    pypi
gcc                       10.4.0              hb92f740_13    conda-forge
gcc_impl_linux-64         10.4.0              h5231bdf_19    conda-forge
gcc_linux-64              10.4.0              h9215b83_13    conda-forge
gdal                      3.5.0            py37hed81ac4_0    conda-forge
geos                      3.10.2               h9c3ff4c_0    conda-forge
geotiff                   1.7.1                h509b78c_1    conda-forge
gettext                   0.21.1               h27087fc_0    conda-forge
gfortran                  10.4.0              h0c96582_13    conda-forge
gfortran_impl_linux-64    10.4.0              h7d168d2_19    conda-forge
gfortran_linux-64         10.4.0              h69d5af5_13    conda-forge
ghostscript               9.18                          1    bioconda
giflib                    5.2.1                h36c2ea0_2    conda-forge
globus-sdk                3.41.0                   pypi_0    pypi
glpk                      5.0                  h445213a_0    conda-forge
gmp                       6.2.1                h58526e2_0    conda-forge
google-auth               2.35.0                   pypi_0    pypi
google-auth-oauthlib      0.4.6                    pypi_0    pypi
google-pasta              0.2.0                    pypi_0    pypi
graphite2                 1.3.13            h58526e2_1001    conda-forge
grpcio                    1.58.0                   pypi_0    pypi
gsl                       2.7                  he838d99_0    conda-forge
gxx                       10.4.0              hb92f740_13    conda-forge
gxx_impl_linux-64         10.4.0              h5231bdf_19    conda-forge
gxx_linux-64              10.4.0              h6e491c6_13    conda-forge
h5py                      3.7.0           nompi_py37hc97082c_100    conda-forge
harfbuzz                  4.2.0                h40b6f09_0    conda-forge
hdf4                      4.2.15               h9772cbc_5    conda-forge
hdf5                      1.12.1          mpi_openmpi_h41b9b70_4    conda-forge
heapdict                  1.0.1                    pypi_0    pypi
hisat2                    2.2.1                h87f3376_4    bioconda
hnswlib                   0.8.0                    pypi_0    pypi
homer                     4.11            pl5262h9f5acd7_8    bioconda
htslib                    1.10                 h244ad75_0    bioconda
humanfriendly             10.0                     pypi_0    pypi
humanize                  4.6.0                    pypi_0    pypi
hypre                     2.25.0          mpi_openmpi_ha709252_0    conda-forge
icu                       69.1                 h9c3ff4c_0    conda-forge
idna                      3.4                      pypi_0    pypi
igraph                    0.10.8                   pypi_0    pypi
imagecodecs               2021.4.28        py37hd0c323f_0    conda-forge
imagesize                 1.4.1                    pypi_0    pypi
importlib-metadata        4.11.4           py37h89c1867_0    conda-forge
importlib-resources       5.12.0                   pypi_0    pypi
importlib_metadata        4.11.4               hd8ed1ab_0    conda-forge
interlap                  0.2.7                    pypi_0    pypi
ipykernel                 6.16.2             pyh210e3f2_0    conda-forge
ipython                   7.33.0           py37h89c1867_0    conda-forge
ipython_genutils          0.2.0                      py_1    conda-forge
ipywidgets                8.1.2              pyhd8ed1ab_0    conda-forge
isa-l                     2.31.0               h4bc722e_2    conda-forge
java-jdk                  7.0.91                        1    bioconda
jax                       0.3.25                   pypi_0    pypi
jaxlib                    0.3.25                   pypi_0    pypi
jaxopt                    0.8                      pypi_0    pypi
jbig                      2.1               h7f98852_2003    conda-forge
jedi                      0.18.2             pyhd8ed1ab_0    conda-forge
jinja2                    3.0.3                    pypi_0    pypi
joblib                    1.2.0              pyhd8ed1ab_0    conda-forge
jpeg                      9e                   h166bdaf_2    conda-forge
jq                        1.5                           0    bioconda
jsbeautifier              1.14.9             pyhd8ed1ab_0    conda-forge
json-c                    0.16                 hc379101_0    conda-forge
jsonschema                2.6.0                 py37_1002    conda-forge
jupyter                   1.0.0             pyhd8ed1ab_10    conda-forge
jupyter-core              4.12.0                   pypi_0    pypi
jupyter_client            7.4.9              pyhd8ed1ab_0    conda-forge
jupyter_console           6.5.1              pyhd8ed1ab_0    conda-forge
jupyter_core              4.11.1           py37h89c1867_0    conda-forge
jupyterlab_widgets        3.0.10             pyhd8ed1ab_0    conda-forge
jxrlib                    1.1                  hd590300_3    conda-forge
k8                        0.2.5                hd03093a_2    bioconda
kaleido                   0.2.1                    pypi_0    pypi
kealib                    1.4.15               hfe1a663_0    conda-forge
keras                     2.11.0                   pypi_0    pypi
kernel-headers_linux-64   2.6.32              he073ed8_15    conda-forge
keyutils                  1.6.1                h166bdaf_0    conda-forge
kiwisolver                1.4.4            py37h7cecad7_0    conda-forge
krb5                      1.19.3               h08a2579_0    conda-forge
lcms2                     2.12                 hddcbb42_0    conda-forge
ld_impl_linux-64          2.39                 hcc3a1bd_1    conda-forge
leidenalg                 0.10.1                   pypi_0    pypi
lerc                      2.2.1                h9c3ff4c_0    conda-forge
libaec                    1.1.3                h59595ed_0    conda-forge
libblas                   3.9.0           16_linux64_openblas    conda-forge
libbrotlicommon           1.0.9                h166bdaf_9    conda-forge
libbrotlidec              1.0.9                h166bdaf_9    conda-forge
libbrotlienc              1.0.9                h166bdaf_9    conda-forge
libcblas                  3.9.0           16_linux64_openblas    conda-forge
libclang                  16.0.6                   pypi_0    pypi
libcups                   2.3.3                h3e49a29_2    conda-forge
libcurl                   7.86.0               h2283fc2_1    conda-forge
libdap4                   3.20.6               hd7c4107_2    conda-forge
libdeflate                1.7                  h7f98852_5    conda-forge
libedit                   3.1.20191231         he28a2e2_2    conda-forge
libev                     4.33                 h516909a_1    conda-forge
libexpat                  2.5.0                hcb278e6_1    conda-forge
libffi                    3.4.2                h7f98852_5    conda-forge
libgcc                    14.2.0               h77fa898_1    conda-forge
libgcc-devel_linux-64     10.4.0              hd38fd1e_19    conda-forge
libgcc-ng                 14.2.0               h69a702a_1    conda-forge
libgdal                   3.5.0                h3886835_0    conda-forge
libgfortran-ng            13.2.0               h69a702a_6    conda-forge
libgfortran5              13.2.0               h43f5ff8_6    conda-forge
libglib                   2.74.1               h7a41b64_0    conda-forge
libgomp                   14.2.0               h77fa898_1    conda-forge
libhwloc                  2.8.0                h32351e8_1    conda-forge
libiconv                  1.17                 h166bdaf_0    conda-forge
libidn2                   2.3.4                h166bdaf_0    conda-forge
libkml                    1.3.0             h01aab08_1016    conda-forge
liblapack                 3.9.0           16_linux64_openblas    conda-forge
libllvm11                 11.1.0               he0ac6c6_5    conda-forge
libnetcdf                 4.8.1           mpi_openmpi_he7012b2_2    conda-forge
libnghttp2                1.52.0               h61bc06f_0    conda-forge
libnsl                    2.0.0                h7f98852_0    conda-forge
libnuma                   2.0.16               h0b41bf4_1    conda-forge
libopenblas               0.3.21          pthreads_h78a6416_3    conda-forge
libopenssl-static         3.0.7                h0b41bf4_1    conda-forge
libpng                    1.6.43               h2797004_0    conda-forge
libpq                     14.5                 he2d8382_0    conda-forge
librttopo                 1.1.0                hf69c175_9    conda-forge
libsanitizer              10.4.0              h5246dfb_19    conda-forge
libsodium                 1.0.18               h36c2ea0_1    conda-forge
libspatialite             5.0.1               ha867d66_15    conda-forge
libsqlite                 3.40.0               h753d276_0    conda-forge
libssh2                   1.11.0               h0841786_0    conda-forge
libstdcxx                 14.2.0               hc0a3c3a_1    conda-forge
libstdcxx-devel_linux-64  10.4.0              hd38fd1e_19    conda-forge
libstdcxx-ng              13.2.0               h7e041cc_5    conda-forge
libtiff                   4.3.0                hf544144_1    conda-forge
libudunits2               2.2.28               h40f5838_3    conda-forge
libunistring              0.9.10               h7f98852_0    conda-forge
libuuid                   2.38.1               h0b41bf4_0    conda-forge
libuv                     1.44.2               h166bdaf_0    conda-forge
libwebp                   1.2.4                h11a3e52_1  
libwebp-base              1.2.4                h166bdaf_0    conda-forge
libxcb                    1.13              h7f98852_1004    conda-forge
libxcrypt                 4.4.36               hd590300_1    conda-forge
libxml2                   2.9.14               haae042b_4    conda-forge
libxslt                   1.1.33               h0ef7038_3    conda-forge
libzip                    1.10.1               h2629f0a_3    conda-forge
libzlib                   1.2.13               h166bdaf_4    conda-forge
libzopfli                 1.0.3                h9c3ff4c_0    conda-forge
llvmlite                  0.39.1                   pypi_0    pypi
locket                    1.0.0                    pypi_0    pypi
loompy                    2.0.16                     py_0    bioconda
louvain                   0.8.1                    pypi_0    pypi
lz4-c                     1.9.3                h9c3ff4c_1    conda-forge
macs2                     2.2.7.1          py37h8902056_5    bioconda
macs3                     3.0.0a7                  pypi_0    pypi
make                      4.3                  hd18ef5c_1    conda-forge
mappy                     2.26                     pypi_0    pypi
markdown                  3.4.4                    pypi_0    pypi
markdown-it-py            2.2.0                    pypi_0    pypi
markupsafe                2.1.3                    pypi_0    pypi
matplotlib-base           3.4.3            py37h1058ff1_2    conda-forge
matplotlib-inline         0.1.6              pyhd8ed1ab_0    conda-forge
mdurl                     0.1.2                    pypi_0    pypi
mellon                    1.3.1                    pypi_0    pypi
meme                      5.5.2           py37pl5321hea7a43d_0    bioconda
metis                     5.1.0             h59595ed_1007    conda-forge
minimap2                  2.24                 h7132678_1    bioconda
mistune                   0.8.4           pyh1a96a4e_1006    conda-forge
mnnpy                     0.1.9.5                  pypi_0    pypi
mpfr                      4.2.0                hb012696_0    conda-forge
mpi                       1.0                     openmpi    conda-forge
msgpack                   1.0.5                    pypi_0    pypi
multidict                 6.0.5                    pypi_0    pypi
multiprocess              0.70.15                  pypi_0    pypi
multiprocessing-on-dill   3.5.0a4                  pypi_0    pypi
multiqc                   1.19                     pypi_0    pypi
mumps-include             5.2.1               ha770c72_11    conda-forge
mumps-mpi                 5.2.1               hfb3545b_11    conda-forge
mysql-connector-c         6.1.11            h659d440_1008    conda-forge
nanocomp                  1.23.1                   pypi_0    pypi
nanofilt                  2.8.0                    pypi_0    pypi
nanoget                   1.19.1                   pypi_0    pypi
nanolyse                  1.2.1                    pypi_0    pypi
nanomath                  1.3.0                    pypi_0    pypi
nanopack                  1.1.0                    pypi_0    pypi
nanoplot                  1.42.0                   pypi_0    pypi
nanoqc                    0.9.4                    pypi_0    pypi
nanostat                  1.6.0                    pypi_0    pypi
natsort                   8.4.0                    pypi_0    pypi
nbconvert                 5.6.1              pyhd8ed1ab_2    conda-forge
nbformat                  5.7.3              pyhd8ed1ab_0    conda-forge
ncurses                   6.3                  h27087fc_1    conda-forge
nest-asyncio              1.5.8                    pypi_0    pypi
networkx                  2.6.3                    pypi_0    pypi
notebook                  5.7.11           py37h89c1867_0    conda-forge
novoalign                 4.02.02              h82c745c_4    bioconda
nspr                      4.36                 h5888daf_0    conda-forge
nss                       3.89                 he45b914_0    conda-forge
numba                     0.56.4                   pypi_0    pypi
numexpr                   2.8.6                    pypi_0    pypi
numpy                     1.21.6           py37h976b520_0    conda-forge
numpy-groupies            0.9.22                   pypi_0    pypi
numpydoc                  1.5.0                    pypi_0    pypi
nvidia-cublas-cu11        11.10.3.66               pypi_0    pypi
nvidia-cuda-nvrtc-cu11    11.7.99                  pypi_0    pypi
nvidia-cuda-runtime-cu11  11.7.99                  pypi_0    pypi
nvidia-cudnn-cu11         8.5.0.96                 pypi_0    pypi
oauthlib                  3.2.2                    pypi_0    pypi
opencv-python-headless    4.10.0.84                py37_0    fastai
openjdk                   17.0.11              h24d6bf4_0    <unknown>
openjpeg                  2.5.0                h7d73246_0    conda-forge
openmpi                   4.1.5              h414af15_101    conda-forge
openpyxl                  3.1.3                    pypi_0    pypi
openssl                   3.1.7                hb9d3cd8_0    conda-forge
opt-einsum                3.3.0                    pypi_0    pypi
packaging                 23.0               pyhd8ed1ab_0    conda-forge
pairtools                 1.1.0                    pypi_0    pypi
palantir                  1.3.0                    pypi_0    pypi
pandas                    1.3.5                    pypi_0    pypi
pandoc                    2.19.2               h32600fe_1    conda-forge
pandocfilters             1.5.0              pyhd8ed1ab_0    conda-forge
pango                     1.50.7               hbd2fdc8_0    conda-forge
parmetis                  4.0.3             he9a3056_1005    conda-forge
parso                     0.8.3              pyhd8ed1ab_0    conda-forge
partd                     1.4.1                    pypi_0    pypi
patsy                     0.5.3                    pypi_0    pypi
pbzip2                    1.1.13               h1fcc475_2    conda-forge
pcre                      8.45                 h9c3ff4c_0    conda-forge
pcre2                     10.37                hc3806b6_1    conda-forge
perl                      5.32.1          7_hd590300_perl5    conda-forge
perl-app-cpanminus        1.7047          pl5321hd8ed1ab_0    conda-forge
perl-base                 2.23            pl5321hdfd78af_2    bioconda
perl-business-isbn        3.007           pl5321hdfd78af_0    bioconda
perl-business-isbn-data   20210112.006    pl5321hdfd78af_0    bioconda
perl-carp                 1.38            pl5321hdfd78af_4    bioconda
perl-cgi                  4.56            pl5321h031d066_1    bioconda
perl-common-sense         3.75            pl5321hd8ed1ab_0    conda-forge
perl-compress-raw-zlib    2.105           pl5321h87f3376_0    bioconda
perl-config-general       2.65            pl5321hdfd78af_0    bioconda
perl-constant             1.33            pl5321hdfd78af_2    bioconda
perl-data-dumper          2.183           pl5321hec16e2b_1    bioconda
perl-dbi                  1.643           pl5321h166bdaf_0    conda-forge
perl-digest-md5           2.58            pl5321h166bdaf_0    conda-forge
perl-encode               3.19            pl5321hec16e2b_1    bioconda
perl-encode-locale        1.05            pl5321hdfd78af_7    bioconda
perl-exporter             5.72            pl5321hdfd78af_2    bioconda
perl-extutils-makemaker   7.70            pl5321hd8ed1ab_0    conda-forge
perl-file-path            2.18            pl5321hd8ed1ab_0    conda-forge
perl-file-spec            3.48_01         pl5321hdfd78af_2    bioconda
perl-file-temp            0.2304          pl5321hd8ed1ab_0    conda-forge
perl-file-which           1.24            pl5321hd8ed1ab_0    conda-forge
perl-html-parser          3.81            pl5321h4ac6f70_1    bioconda
perl-html-tagset          3.20            pl5321hdfd78af_4    bioconda
perl-html-template        2.97            pl5321hdfd78af_2    bioconda
perl-html-tree            5.07            pl5321hdfd78af_2    bioconda
perl-http-date            6.06            pl5321hdfd78af_0    bioconda
perl-http-message         6.36            pl5321hdfd78af_0    bioconda
perl-inc-latest           0.500           pl5321ha770c72_0    conda-forge
perl-io-html              1.004           pl5321hdfd78af_0    bioconda
perl-json                 4.10            pl5321hdfd78af_1    bioconda
perl-json-xs              4.03            pl5321h4ac6f70_3    bioconda
perl-log-log4perl         1.56            pl5321hd8ed1ab_0    conda-forge
perl-lwp-mediatypes       6.04            pl5321hdfd78af_1    bioconda
perl-math-cdf             0.1             pl5321h031d066_10    bioconda
perl-mime-base64          3.16            pl5321hec16e2b_2    bioconda
perl-module-build         0.4232          pl5321ha770c72_0    conda-forge
perl-parent               0.236           pl5321hdfd78af_2    bioconda
perl-pathtools            3.75            pl5321h166bdaf_0    conda-forge
perl-scalar-list-utils    1.62            pl5321hec16e2b_1    bioconda
perl-sys-info             0.7811          pl5321hdfd78af_1    bioconda
perl-sys-info-base        0.7807          pl5321hdfd78af_1    bioconda
perl-sys-info-driver-linux 0.7905          pl5321hdfd78af_1    bioconda
perl-test-nowarnings      1.06            pl5321hdfd78af_0    bioconda
perl-text-template-simple 0.91            pl5321hdfd78af_1    bioconda
perl-time-local           1.35            pl5321hdfd78af_0    bioconda
perl-timedate             2.33            pl5321hdfd78af_2    bioconda
perl-types-serialiser     1.01            pl5321hdfd78af_0    bioconda
perl-unix-processors      2.046           pl5321h7f98852_1001    conda-forge
perl-uri                  5.12            pl5321hdfd78af_0    bioconda
perl-url-encode           0.03            pl5321h9ee0642_0    bioconda
perl-xml-namespacesupport 1.12            pl5321hdfd78af_1    bioconda
perl-xml-parser           2.44_01         pl5321hc3e0081_1003    conda-forge
perl-xml-sax              1.02            pl5321hdfd78af_1    bioconda
perl-xml-sax-base         1.09            pl5321hdfd78af_1    bioconda
perl-xml-sax-expat        0.51            pl5321hdfd78af_4    bioconda
perl-xml-simple           2.25            pl5321hdfd78af_2    bioconda
perl-yaml                 1.30            pl5321hdfd78af_0    bioconda
petsc                     3.17.3          real_h1c6c8d6_100    conda-forge
petsc4py                  3.17.3          real_h56b374f_101    conda-forge
pexpect                   4.8.0              pyh1a96a4e_2    conda-forge
pickleshare               0.7.5                   py_1003    conda-forge
pigz                      2.8                  h2797004_0    conda-forge
pillow                    9.1.1            py37h44f0d7a_1    conda-forge
pip                       22.3.1             pyhd8ed1ab_0    conda-forge
pixman                    0.40.0               h36c2ea0_0    conda-forge
plotly                    5.13.0             pyhd8ed1ab_0    conda-forge
poppler                   22.04.0              h1434ded_0    conda-forge
poppler-data              0.4.12               hd8ed1ab_0    conda-forge
postgresql                14.5                 ha7cec9f_0    conda-forge
progressbar2              4.2.0                    pypi_0    pypi
proj                      9.0.0                h93bde94_1    conda-forge
prometheus_client         0.16.0             pyhd8ed1ab_0    conda-forge
prompt-toolkit            3.0.36             pyha770c72_0    conda-forge
prompt_toolkit            3.0.36               hd8ed1ab_0    conda-forge
protobuf                  3.19.6                   pypi_0    pypi
psutil                    5.9.6                    pypi_0    pypi
pthread-stubs             0.4               h36c2ea0_1001    conda-forge
ptscotch                  6.0.9                h0a9c416_2    conda-forge
ptyprocess                0.7.0              pyhd3deb0d_0    conda-forge
py2bit                    0.3.0            py37h8902056_6    bioconda
pyaml-env                 1.2.1                    pypi_0    pypi
pyarrow                   12.0.1                   pypi_0    pypi
pyasn1                    0.5.0                    pypi_0    pypi
pyasn1-modules            0.3.0                    pypi_0    pypi
pybigwig                  0.3.18           py37hdc12a6d_2    bioconda
pybktree                  1.1                      pypi_0    pypi
pycparser                 2.21               pyhd8ed1ab_0    conda-forge
pyfaidx                   0.8.1.2                  pypi_0    pypi
pygam                     0.8.0                    pypi_0    pypi
pygments                  2.14.0             pyhd8ed1ab_0    conda-forge
pygpcca                   1.0.4                    pypi_0    pypi
pyjwt                     2.8.0                    pypi_0    pypi
pynndescent               0.5.10                   pypi_0    pypi
pyparsing                 3.0.9              pyhd8ed1ab_0    conda-forge
pysam                     0.16.0           py37ha9a96c6_0    bioconda
python                    3.7.12          hf930737_100_cpython    conda-forge
python-dateutil           2.9.0.post0              pypi_0    pypi
python-deprecated         1.1.0                    pypi_0    pypi
python-fastjsonschema     2.16.2             pyhd8ed1ab_0    conda-forge
python-igraph             0.10.8                   pypi_0    pypi
python-isal               1.1.0            py37h540881e_0    conda-forge
python-tzdata             2023.3             pyhd8ed1ab_0    conda-forge
python-utils              3.5.2                    pypi_0    pypi
python_abi                3.7                     3_cp37m    conda-forge
pytorch-metric-learning   2.4.1                    pypi_0    pypi
pytz                      2024.2             pyhd8ed1ab_0    conda-forge
pytz-deprecation-shim     0.1.0.post0      py37h89c1867_2    conda-forge
pyyaml                    6.0              py37h540881e_4    conda-forge
pyzmq                     24.0.1           py37h0c0c2a8_0    conda-forge
qtconsole-base            5.4.4              pyha770c72_0    conda-forge
qtpy                      2.4.1              pyhd8ed1ab_0    conda-forge
r-abind                   1.4_5           r41hc72bb7e_1004    conda-forge
r-ashr                    2.2_54            r41h7525677_1    conda-forge
r-askpass                 1.1               r41h06615bd_3    conda-forge
r-assertthat              0.2.1             r41hc72bb7e_3    conda-forge
r-backports               1.4.1             r41h06615bd_1    conda-forge
r-base                    4.1.3                hd930d0e_0    conda-forge
r-base64enc               0.1_3           r41h06615bd_1005    conda-forge
r-bbmle                   1.0.25            r41hc72bb7e_1    conda-forge
r-bdsmatrix               1.3_6             r41h06615bd_1    conda-forge
r-beeswarm                0.4.0             r41h06615bd_2    conda-forge
r-bh                      1.81.0_1          r41hc72bb7e_0    conda-forge
r-bit                     4.0.5             r41h06615bd_0    conda-forge
r-bit64                   4.0.5             r41h06615bd_1    conda-forge
r-bitops                  1.0_7             r41h06615bd_1    conda-forge
r-blob                    1.2.4             r41hc72bb7e_0    conda-forge
r-boot                    1.3_28.1          r41hc72bb7e_0    conda-forge
r-broom                   1.0.5             r41hc72bb7e_0    conda-forge
r-bslib                   0.5.0             r41hc72bb7e_0    conda-forge
r-cachem                  1.0.8             r41h57805ef_0    conda-forge
r-cairo                   1.6_0             r41h06615bd_1    conda-forge
r-callr                   3.7.3             r41hc72bb7e_0    conda-forge
r-caret                   6.0_94            r41h133d619_0    conda-forge
r-catools                 1.18.2            r41h7525677_1    conda-forge
r-cellranger              1.1.0           r41hc72bb7e_1005    conda-forge
r-class                   7.3_22            r41h57805ef_0    conda-forge
r-cli                     3.6.1             r41h38f115c_0    conda-forge
r-clipr                   0.8.0             r41hc72bb7e_1    conda-forge
r-clock                   0.7.0             r41ha503ecb_0    conda-forge
r-cluster                 2.1.4             r41h8da6f51_0    conda-forge
r-coda                    0.19_4            r41hc72bb7e_1    conda-forge
r-codetools               0.2_19            r41hc72bb7e_0    conda-forge
r-colorspace              2.1_0             r41h133d619_0    conda-forge
r-commonmark              1.9.0             r41h133d619_0    conda-forge
r-cowplot                 1.1.1             r41hc72bb7e_1    conda-forge
r-cpp11                   0.4.7             r41hc72bb7e_0    conda-forge
r-crayon                  1.5.2             r41hc72bb7e_1    conda-forge
r-crosstalk               1.2.0             r41hc72bb7e_1    conda-forge
r-crul                    1.4.0             r41h785f33e_0    conda-forge
r-curl                    4.3.3             r41h06615bd_1    conda-forge
r-data.table              1.14.8            r41h133d619_0    conda-forge
r-dbi                     1.1.3             r41hc72bb7e_1    conda-forge
r-dbplyr                  2.3.2             r41hc72bb7e_0    conda-forge
r-deldir                  1.0_9             r41h61816a4_0    conda-forge
r-diagram                 1.6.5             r41ha770c72_1    conda-forge
r-digest                  0.6.31            r41h38f115c_0    conda-forge
r-dplyr                   1.1.2             r41ha503ecb_0    conda-forge
r-dqrng                   0.3.0             r41h7525677_1    conda-forge
r-dtplyr                  1.3.1             r41hc72bb7e_0    conda-forge
r-e1071                   1.7_13            r41h38f115c_0    conda-forge
r-ellipsis                0.3.2             r41h06615bd_1    conda-forge
r-emdbook                 1.3.12            r41hc72bb7e_2    conda-forge
r-essentials              1.7.0            r342hf65ed6a_0    R
r-etrunct                 0.1             r41hc72bb7e_1004    conda-forge
r-evaluate                0.21              r41hc72bb7e_0    conda-forge
r-fansi                   1.0.4             r41h133d619_0    conda-forge
r-farver                  2.1.1             r41h7525677_1    conda-forge
r-fastmap                 1.1.1             r41h38f115c_0    conda-forge
r-filelock                1.0.2           r41h133d619_1003    conda-forge
r-fitdistrplus            1.1_11            r41hc72bb7e_0    conda-forge
r-fnn                     1.1.3.2           r41h38f115c_0    conda-forge
r-fontawesome             0.5.1             r41hc72bb7e_0    conda-forge
r-forcats                 1.0.0             r41hc72bb7e_0    conda-forge
r-foreach                 1.5.2             r41hc72bb7e_1    conda-forge
r-foreign                 0.8_84            r41h133d619_0    conda-forge
r-formatr                 1.14              r41hc72bb7e_0    conda-forge
r-fs                      1.6.2             r41ha503ecb_0    conda-forge
r-futile.logger           1.4.3           r41hc72bb7e_1004    conda-forge
r-futile.options          1.0.1           r41hc72bb7e_1003    conda-forge
r-future                  1.32.0            r41hc72bb7e_0    conda-forge
r-future.apply            1.11.0            r41hc72bb7e_0    conda-forge
r-gargle                  1.4.0             r41h785f33e_0    conda-forge
r-generics                0.1.3             r41hc72bb7e_1    conda-forge
r-ggbeeswarm              0.7.2             r41hc72bb7e_0    conda-forge
r-ggplot2                 3.4.2             r41hc72bb7e_0    conda-forge
r-ggrastr                 1.0.2             r41hc72bb7e_0    conda-forge
r-ggrepel                 0.9.3             r41h38f115c_0    conda-forge
r-ggridges                0.5.4             r41hc72bb7e_1    conda-forge
r-gistr                   0.9.0             r41hc72bb7e_1    conda-forge
r-glmnet                  4.1_7             r41haa30946_0    conda-forge
r-globals                 0.16.2            r41hc72bb7e_0    conda-forge
r-glue                    1.6.2             r41h06615bd_1    conda-forge
r-goftest                 1.2_3             r41h06615bd_1    conda-forge
r-googledrive             2.1.0             r41hc72bb7e_0    conda-forge
r-googlesheets4           1.1.0             r41h785f33e_0    conda-forge
r-gower                   1.0.1             r41h133d619_0    conda-forge
r-gplots                  3.1.3             r41hc72bb7e_1    conda-forge
r-gridextra               2.3             r41hc72bb7e_1004    conda-forge
r-gtable                  0.3.3             r41hc72bb7e_0    conda-forge
r-gtools                  3.9.4             r41h06615bd_0    conda-forge
r-hardhat                 1.3.0             r41hc72bb7e_0    conda-forge
r-haven                   2.5.2             r41h38f115c_0    conda-forge
r-here                    1.0.1             r41hc72bb7e_1    conda-forge
r-hexbin                  1.28.3            r41hac0b197_0    conda-forge
r-highr                   0.10              r41hc72bb7e_0    conda-forge
r-hms                     1.1.3             r41hc72bb7e_0    conda-forge
r-htmltools               0.5.5             r41h38f115c_0    conda-forge
r-htmlwidgets             1.6.2             r41hc72bb7e_0    conda-forge
r-httpcode                0.3.0             r41ha770c72_2    conda-forge
r-httpuv                  1.6.11            r41ha503ecb_0    conda-forge
r-httr                    1.4.6             r41hc72bb7e_0    conda-forge
r-ica                     1.0_3             r41hc72bb7e_1    conda-forge
r-ids                     1.0.1             r41hc72bb7e_2    conda-forge
r-igraph                  1.3.5             r41hb34fc8a_0    conda-forge
r-invgamma                1.1               r41hc72bb7e_2    conda-forge
r-ipred                   0.9_14            r41h133d619_0    conda-forge
r-irdisplay               1.1               r41hd8ed1ab_1    conda-forge
r-irkernel                1.3.2             r41h785f33e_0    conda-forge
r-irlba                   2.3.5.1           r41h5f7b363_0    conda-forge
r-isoband                 0.2.7             r41h38f115c_1    conda-forge
r-iterators               1.0.14            r41hc72bb7e_1    conda-forge
r-jquerylib               0.1.4             r41hc72bb7e_1    conda-forge
r-jsonlite                1.8.5             r41h57805ef_0    conda-forge
r-kernsmooth              2.23_21           r41h13b3f57_0    conda-forge
r-knitr                   1.43              r41hc72bb7e_0    conda-forge
r-labeling                0.4.2             r41hc72bb7e_2    conda-forge
r-lambda.r                1.2.4             r41hc72bb7e_2    conda-forge
r-later                   1.3.1             r41ha503ecb_0    conda-forge
r-lattice                 0.21_8            r41h133d619_0    conda-forge
r-lava                    1.7.2.1           r41hc72bb7e_0    conda-forge
r-lazyeval                0.2.2             r41h06615bd_3    conda-forge
r-leiden                  0.4.3             r41hc72bb7e_1    conda-forge
r-lifecycle               1.0.3             r41hc72bb7e_1    conda-forge
r-listenv                 0.9.0             r41hc72bb7e_0    conda-forge
r-lmtest                  0.9_40            r41h8da6f51_1    conda-forge
r-lobstr                  1.1.2             r41h38f115c_2    conda-forge
r-locfit                  1.5_9.7           r41h133d619_0    conda-forge
r-lsei                    1.3_0             r41hc3ea6d6_2    conda-forge
r-lubridate               1.9.2             r41h133d619_1    conda-forge
r-magrittr                2.0.3             r41h06615bd_1    conda-forge
r-maps                    3.4.1             r41h06615bd_1    conda-forge
r-mass                    7.3_58.3          r41h133d619_0    conda-forge
r-matrix                  1.5_4.1           r41h316c678_0    conda-forge
r-matrixstats             1.0.0             r41h57805ef_0    conda-forge
r-memoise                 2.0.1             r41hc72bb7e_1    conda-forge
r-mgcv                    1.8_42            r41he1ae0d6_0    conda-forge
r-mime                    0.12              r41h06615bd_1    conda-forge
r-miniui                  0.1.1.1         r41hc72bb7e_1003    conda-forge
r-mixsqp                  0.3_48            r41h9f5de39_0    conda-forge
r-modelmetrics            1.2.2.2           r41h7525677_2    conda-forge
r-modelr                  0.1.11            r41hc72bb7e_0    conda-forge
r-munsell                 0.5.0           r41hc72bb7e_1005    conda-forge
r-mvtnorm                 1.2_2             r41h61816a4_0    conda-forge
r-nlme                    3.1_162           r41hac0b197_0    conda-forge
r-nnet                    7.3_19            r41h57805ef_0    conda-forge
r-npsurv                  0.5_0             r41hc72bb7e_1    conda-forge
r-numderiv                2016.8_1.1        r41hc72bb7e_4    conda-forge
r-openssl                 2.0.6             r41habfbb5e_0    conda-forge
r-parallelly              1.36.0            r41hc72bb7e_0    conda-forge
r-patchwork               1.1.2             r41hc72bb7e_1    conda-forge
r-pbapply                 1.7_0             r41hc72bb7e_0    conda-forge
r-pbdzmq                  0.3_9             r41hfae1697_0    conda-forge
r-pillar                  1.9.0             r41hc72bb7e_0    conda-forge
r-pkgconfig               2.0.3             r41hc72bb7e_2    conda-forge
r-plogr                   0.2.0           r41hc72bb7e_1004    conda-forge
r-plotly                  4.10.2            r41hc72bb7e_0    conda-forge
r-plyr                    1.8.8             r41h7525677_0    conda-forge
r-png                     0.1_8             r41h10cf519_0    conda-forge
r-polyclip                1.10_4            r41h7525677_0    conda-forge
r-prettyunits             1.1.1             r41hc72bb7e_2    conda-forge
r-proc                    1.18.2            r41ha503ecb_0    conda-forge
r-processx                3.8.1             r41h133d619_0    conda-forge
r-prodlim                 2023.03.31        r41h38f115c_0    conda-forge
r-progress                1.2.2             r41hc72bb7e_3    conda-forge
r-progressr               0.13.0            r41hc72bb7e_0    conda-forge
r-promises                1.2.0.1           r41h7525677_1    conda-forge
r-proxy                   0.4_27            r41h06615bd_1    conda-forge
r-pryr                    0.1.6             r41h38f115c_0    conda-forge
r-ps                      1.7.5             r41h133d619_0    conda-forge
r-purrr                   1.0.1             r41h133d619_0    conda-forge
r-quantmod                0.4.22            r41hc72bb7e_0    conda-forge
r-r6                      2.5.1             r41hc72bb7e_1    conda-forge
r-ragg                    1.2.2             r41hc1f6985_0    conda-forge
r-randomforest            4.7_1.1           r41h8da6f51_1    conda-forge
r-rann                    2.6.1             r41h7525677_3    conda-forge
r-rappdirs                0.3.3             r41h06615bd_1    conda-forge
r-rbokeh                  0.5.2             r41hc72bb7e_2    conda-forge
r-rcolorbrewer            1.1_3             r41h785f33e_1    conda-forge
r-rcpp                    1.0.10            r41h38f115c_0    conda-forge
r-rcppannoy               0.0.20            r41h7525677_0    conda-forge
r-rcpparmadillo           0.12.4.0.0        r41h08d816e_0    conda-forge
r-rcppeigen               0.3.3.9.3         r41h9f5de39_0    conda-forge
r-rcppnumerical           0.5_0             r41h38f115c_0    conda-forge
r-rcppparallel            5.1.6             r41h38f115c_0    conda-forge
r-rcppprogress            0.4.2             r41hc72bb7e_2    conda-forge
r-rcpptoml                0.2.2             r41h38f115c_0    conda-forge
r-rcurl                   1.98_1.12         r41h133d619_0    conda-forge
r-readr                   2.1.4             r41h38f115c_0    conda-forge
r-readxl                  1.4.2             r41h81ef4d7_0    conda-forge
r-recipes                 1.0.6             r41hc72bb7e_0    conda-forge
r-recommended             4.1             r41hd8ed1ab_1005    conda-forge
r-rematch                 1.0.1           r41hc72bb7e_1005    conda-forge
r-rematch2                2.1.2             r41hc72bb7e_2    conda-forge
r-repr                    1.1.6             r41h785f33e_0    conda-forge
r-reprex                  2.0.2             r41hc72bb7e_1    conda-forge
r-reshape2                1.4.4             r41h7525677_2    conda-forge
r-restfulr                0.0.15            r41h73dbb54_0    bioconda
r-reticulate              1.30              r41ha503ecb_1    conda-forge
r-rgeos                   0.5_9             r41hf69c175_1    conda-forge
r-rjson                   0.2.21            r41h7525677_2    conda-forge
r-rlang                   1.1.1             r41ha503ecb_0    conda-forge
r-rmarkdown               2.22              r41hc72bb7e_0    conda-forge
r-rocr                    1.0_11            r41hc72bb7e_2    conda-forge
r-rpart                   4.1.19            r41h06615bd_0    conda-forge
r-rprojroot               2.0.3             r41hc72bb7e_1    conda-forge
r-rspectra                0.16_1            r41h9f5de39_1    conda-forge
r-rsqlite                 2.3.1             r41h38f115c_0    conda-forge
r-rstudioapi              0.14              r41hc72bb7e_1    conda-forge
r-rtsne                   0.16              r41h37cf8d7_1    conda-forge
r-rvest                   1.0.3             r41hc72bb7e_1    conda-forge
r-sass                    0.4.6             r41ha503ecb_0    conda-forge
r-scales                  1.2.1             r41hc72bb7e_1    conda-forge
r-scattermore             1.1               r41ha503ecb_0    conda-forge
r-sctransform             0.3.5             r41h9f5de39_1    conda-forge
r-selectr                 0.4_2             r41hc72bb7e_2    conda-forge
r-seurat                  4.3.0             r41h38f115c_0    conda-forge
r-seuratobject            4.1.3             r41h38f115c_0    conda-forge
r-shape                   1.4.6             r41ha770c72_1    conda-forge
r-shiny                   1.7.4             r41h785f33e_0    conda-forge
r-sitmo                   2.0.2             r41h7525677_1    conda-forge
r-snow                    0.4_4             r41hc72bb7e_1    conda-forge
r-sourcetools             0.1.7_1           r41h38f115c_0    conda-forge
r-sp                      1.6_1             r41h57805ef_0    conda-forge
r-spatial                 7.3_16            r41h133d619_0    conda-forge
r-spatstat.data           3.0_1             r41hc72bb7e_0    conda-forge
r-spatstat.explore        3.2_1             r41h57805ef_0    conda-forge
r-spatstat.geom           3.2_1             r41h57805ef_0    conda-forge
r-spatstat.random         3.1_5             r41ha503ecb_0    conda-forge
r-spatstat.sparse         3.0_1             r41h133d619_0    conda-forge
r-spatstat.utils          3.1_0             r41h2b5f3a1_1    conda-forge
r-squarem                 2021.1            r41hc72bb7e_1    conda-forge
r-stringi                 1.7.6             r41h337692f_1    conda-forge
r-stringr                 1.5.0             r41h785f33e_0    conda-forge
r-survival                3.5_5             r41h133d619_0    conda-forge
r-sys                     3.4.2             r41h57805ef_0    conda-forge
r-systemfonts             1.0.4             r41h0ff29ef_1    conda-forge
r-tensor                  1.5             r41hc72bb7e_1004    conda-forge
r-textshaping             0.3.6             r41ha32badf_1    conda-forge
r-tibble                  3.2.1             r41h133d619_1    conda-forge
r-tidyr                   1.3.0             r41h38f115c_0    conda-forge
r-tidyselect              1.2.0             r41hc72bb7e_0    conda-forge
r-tidyverse               1.3.2             r41hc72bb7e_1    conda-forge
r-timechange              0.2.0             r41h38f115c_0    conda-forge
r-timedate                4022.108          r41hc72bb7e_0    conda-forge
r-tinytex                 0.45              r41hc72bb7e_0    conda-forge
r-triebeard               0.4.1             r41h38f115c_0    conda-forge
r-truncnorm               1.0_9             r41h133d619_0    conda-forge
r-ttr                     0.24.3            r41h06615bd_1    conda-forge
r-tzdb                    0.4.0             r41ha503ecb_0    conda-forge
r-urltools                1.7.3             r41h7525677_3    conda-forge
r-utf8                    1.2.3             r41h133d619_0    conda-forge
r-uuid                    1.1_0             r41h06615bd_1    conda-forge
r-uwot                    0.1.14            r41h7525677_1    conda-forge
r-vctrs                   0.6.2             r41ha503ecb_0    conda-forge
r-vipor                   0.4.5           r41hc72bb7e_1004    conda-forge
r-viridislite             0.4.1             r41hc72bb7e_1    conda-forge
r-vroom                   1.6.3             r41ha503ecb_0    conda-forge
r-withr                   2.5.0             r41hc72bb7e_1    conda-forge
r-xfun                    0.39              r41ha503ecb_0    conda-forge
r-xml                     3.99_0.11         r41h2b86b34_2    conda-forge
r-xml2                    1.3.3             r41h044e5c7_2    conda-forge
r-xtable                  1.8_4             r41hc72bb7e_4    conda-forge
r-xts                     0.13.1            r41h133d619_0    conda-forge
r-yaml                    2.3.7             r41h133d619_0    conda-forge
r-zoo                     1.8_12            r41h133d619_0    conda-forge
rdma-core                 28.9                 h59595ed_1    conda-forge
readline                  8.1.2                h0f457ee_0    conda-forge
regex                     2024.4.16                pypi_0    pypi
requests                  2.28.2                   pypi_0    pypi
requests-oauthlib         1.3.1                    pypi_0    pypi
rhash                     1.4.3                h166bdaf_0    conda-forge
rich                      13.8.0                   pypi_0    pypi
rich-click                1.8.3                    pypi_0    pypi
rpy2                      3.5.1           py37r41hc105733_1    conda-forge
rsa                       4.9                      pypi_0    pypi
samtools                  1.6                 hc3601fc_10    bioconda
scalapack                 2.2.0                h67de57e_1    conda-forge
scanpy                    1.9.3                    pypi_0    pypi
scdml                     0.0.1                    pypi_0    pypi
scikit-learn              1.0.2            py37hf9e9bfc_0    conda-forge
scipy                     1.7.3            py37hf2a6cf1_0    conda-forge
scotch                    6.0.9                hb2e6521_2    conda-forge
scvelo                    0.2.5                    pypi_0    pypi
seaborn                   0.12.2                   pypi_0    pypi
sed                       4.8                  he412f7d_0    conda-forge
send2trash                1.8.0              pyhd8ed1ab_0    conda-forge
session-info              1.0.0                    pypi_0    pypi
setuptools                59.8.0           py37h89c1867_1    conda-forge
simplegeneric             0.8.1                      py_1    conda-forge
simplejson                3.19.3                   pypi_0    pypi
six                       1.16.0             pyh6c4a22f_0    conda-forge
slepc                     3.17.1          real_h11342ed_102    conda-forge
slepc4py                  3.17.1          real_h584801e_103    conda-forge
snappy                    1.1.10               hdb0a2a9_1    conda-forge
snowballstemmer           2.2.0                    pypi_0    pypi
sortedcontainers          2.4.0                    pypi_0    pypi
spectra                   0.0.11                   pypi_0    pypi
sphinx                    5.3.0                    pypi_0    pypi
sphinxcontrib-applehelp   1.0.2                    pypi_0    pypi
sphinxcontrib-devhelp     1.0.2                    pypi_0    pypi
sphinxcontrib-htmlhelp    2.0.0                    pypi_0    pypi
sphinxcontrib-jsmath      1.0.1                    pypi_0    pypi
sphinxcontrib-qthelp      1.0.3                    pypi_0    pypi
sphinxcontrib-serializinghtml 1.1.5                    pypi_0    pypi
sqlite                    3.40.0               h4ff8645_0    conda-forge
sra-tools                 2.8.0                         0    bioconda
sratoolkit                2.5.7                         0    daler
statsmodels               0.13.5                   pypi_0    pypi
stdlib-list               0.9.0                    pypi_0    pypi
suitesparse               5.10.1               h9e50725_1    conda-forge
superlu                   5.2.2                h00795ac_0    conda-forge
superlu_dist              7.2.0                h34f6f4d_0    conda-forge
sysroot_linux-64          2.12                he073ed8_15    conda-forge
tables                    3.7.0                    pypi_0    pypi
tbb                       2021.10.0                pypi_0    pypi
tblib                     2.0.0                    pypi_0    pypi
tenacity                  8.1.0              pyhd8ed1ab_0    conda-forge
tensorboard               2.11.2                   pypi_0    pypi
tensorboard-data-server   0.6.1                    pypi_0    pypi
tensorboard-plugin-wit    1.8.1                    pypi_0    pypi
tensorflow                2.11.0                   pypi_0    pypi
tensorflow-estimator      2.11.0                   pypi_0    pypi
tensorflow-io-gcs-filesystem 0.34.0                   pypi_0    pypi
termcolor                 2.3.0                    pypi_0    pypi
terminado                 0.17.1             pyh41d4057_0    conda-forge
testpath                  0.6.0              pyhd8ed1ab_0    conda-forge
texttable                 1.6.7                    pypi_0    pypi
threadpoolctl             3.1.0              pyh8a188c0_0    conda-forge
tifffile                  2021.7.2           pyhd8ed1ab_0    conda-forge
tiledb                    2.8.3                h3f4058f_1    conda-forge
tk                        0.1.0                    pypi_0    pypi
tktable                   2.10                 hb7b940f_3    conda-forge
toml                      0.10.2             pyhd8ed1ab_0    conda-forge
toolz                     0.12.1                   pypi_0    pypi
torch                     1.13.1                   pypi_0    pypi
torchdiffeq               0.2.3                    pypi_0    pypi
tornado                   6.2              py37h540881e_0    conda-forge
tqdm                      4.66.1                   pypi_0    pypi
traitlets                 5.8.1              pyhd8ed1ab_0    conda-forge
trimmomatic               0.39                 hdfd78af_2    bioconda
typing                    3.10.0.0           pyhd8ed1ab_0    conda-forge
typing-extensions         4.4.0                hd8ed1ab_0    conda-forge
typing_extensions         4.4.0              pyha770c72_0    conda-forge
tzcode                    2024b                hb9d3cd8_0    conda-forge
tzdata                    2023c                h71feb2d_0    conda-forge
tzlocal                   4.2              py37h89c1867_1    conda-forge
ucsc-bedgraphtobigwig     472                  h9b8f530_1    bioconda
ucsc-bigwigmerge          469                  h9b8f530_0    bioconda
ucsc-bigwigtowig          469                  h9b8f530_0    bioconda
ucsc-fasize               469                  h9b8f530_0    bioconda
ucsc-fetchchromsizes      469                  h9b8f530_0    bioconda
ucsc-wigtobigwig          472                  h9b8f530_1    bioconda
ucx                       1.14.1               h64cca9d_5    conda-forge
udunits2                  2.2.28               h40f5838_3    conda-forge
umap-learn                0.5.4                    pypi_0    pypi
umi-tools                 1.1.4                    pypi_0    pypi
unzip                     6.0                  h7f98852_3    conda-forge
uriparser                 0.9.8                hac33072_0    conda-forge
urllib3                   1.26.16                  pypi_0    pypi
varscan                   2.4.6                hdfd78af_0    bioconda
velocity-package          0.2                       dev_0    <develop>
velocyto.py               0.17.17          py37h37892f8_5    bioconda
wcwidth                   0.2.6              pyhd8ed1ab_0    conda-forge
webencodings              0.5.1                      py_1    conda-forge
websocket-client          0.54.0                   pypi_0    pypi
werkzeug                  2.2.3                    pypi_0    pypi
wget                      1.20.3               ha35d2d1_1    conda-forge
wheel                     0.38.4             pyhd8ed1ab_0    conda-forge
widgetsnbextension        4.0.10             pyhd8ed1ab_0    conda-forge
wrapt                     1.15.0                   pypi_0    pypi
xerces-c                  3.2.3                h8ce2273_4    conda-forge
xmltodict                 0.13.0             pyhd8ed1ab_0    conda-forge
xopen                     1.6.0            py37h89c1867_0    conda-forge
xorg-fixesproto           5.0               h7f98852_1002    conda-forge
xorg-inputproto           2.3.2             h7f98852_1002    conda-forge
xorg-kbproto              1.0.7             h7f98852_1002    conda-forge
xorg-libice               1.0.10               h7f98852_0    conda-forge
xorg-libsm                1.2.3             hd9c2040_1000    conda-forge
xorg-libx11               1.7.2                h7f98852_0    conda-forge
xorg-libxau               1.0.9                h7f98852_0    conda-forge
xorg-libxdmcp             1.1.3                h7f98852_0    conda-forge
xorg-libxext              1.3.4                h7f98852_1    conda-forge
xorg-libxfixes            5.0.3             h7f98852_1004    conda-forge
xorg-libxi                1.7.10               h7f98852_0    conda-forge
xorg-libxrender           0.9.10            h7f98852_1003    conda-forge
xorg-libxt                1.2.1                h7f98852_2    conda-forge
xorg-libxtst              1.2.3             h7f98852_1002    conda-forge
xorg-recordproto          1.14.2            h7f98852_1002    conda-forge
xorg-renderproto          0.11.1            h7f98852_1002    conda-forge
xorg-xextproto            7.3.0             h7f98852_1002    conda-forge
xorg-xproto               7.0.31            h7f98852_1007    conda-forge
xz                        5.2.6                h166bdaf_0    conda-forge
yaml                      0.2.5                h7f98852_2    conda-forge
yarl                      1.9.4                    pypi_0    pypi
yq                        3.1.0              pyhd8ed1ab_0    conda-forge
zeromq                    4.3.4                h9c3ff4c_1    conda-forge
zfp                       0.5.5                h9c3ff4c_8    conda-forge
zict                      2.2.0                    pypi_0    pypi
zipp                      3.11.0             pyhd8ed1ab_0    conda-forge
zlib                      1.2.13               h166bdaf_4    conda-forge
zlib-ng                   2.0.7                h0b41bf4_0    conda-forge
zstd                      1.5.2                h3eb15da_6    conda-forge
```
conda environment locally
```
_r-mutex                  1.0.1               anacondar_1    conda-forge
adjusttext                1.3.0                    pypi_0    pypi
aiohappyeyeballs          2.4.4                    pypi_0    pypi
aiohttp                   3.11.11                  pypi_0    pypi
aiosignal                 1.3.2                    pypi_0    pypi
anndata                   0.10.9                   pypi_0    pypi
anyio                     4.4.0                    pypi_0    pypi
appnope                   0.1.4                    pypi_0    pypi
arboreto                  0.1.6                    pypi_0    pypi
argon2-cffi               23.1.0                   pypi_0    pypi
argon2-cffi-bindings      21.2.0                   pypi_0    pypi
array-api-compat          1.8                      pypi_0    pypi
arrow                     1.3.0                    pypi_0    pypi
asttokens                 2.4.1              pyhd8ed1ab_0    conda-forge
async-lru                 2.0.4                    pypi_0    pypi
async-timeout             5.0.1                    pypi_0    pypi
attrs                     24.2.0                   pypi_0    pypi
aws-c-auth                0.7.27               h1e647a1_0    conda-forge
aws-c-cal                 0.7.4                h41e72e7_0    conda-forge
aws-c-common              0.9.27               h99b78c6_0    conda-forge
aws-c-compression         0.2.19               h41e72e7_0    conda-forge
aws-c-event-stream        0.4.3                h79ff00d_0    conda-forge
aws-c-http                0.8.8                h69517e7_1    conda-forge
aws-c-io                  0.14.18              h20e6805_7    conda-forge
aws-c-mqtt                0.10.4              h3e8bf47_18    conda-forge
aws-c-s3                  0.6.4               h5e39592_12    conda-forge
aws-c-sdkutils            0.1.19               h85401af_2    conda-forge
aws-checksums             0.1.18              h85401af_10    conda-forge
aws-crt-cpp               0.28.2               h3482d01_1    conda-forge
aws-sdk-cpp               1.11.379             h8d911dc_8    conda-forge
azure-core-cpp            1.13.0               hd01fc5c_0    conda-forge
azure-identity-cpp        1.8.0                h13ea094_2    conda-forge
azure-storage-blobs-cpp   12.12.0              hfde595f_0    conda-forge
azure-storage-common-cpp  12.7.0               hcf3b6fd_1    conda-forge
azure-storage-files-datalake-cpp 12.11.0              h082e32e_1    conda-forge
babel                     2.16.0                   pypi_0    pypi
bcrypt                    4.2.0                    pypi_0    pypi
beautifulsoup4            4.12.3                   pypi_0    pypi
biopython                 1.84                     pypi_0    pypi
biothings-client          0.4.1                    pypi_0    pypi
bleach                    6.1.0                    pypi_0    pypi
bokeh                     3.4.3                    pypi_0    pypi
boltons                   24.0.0                   pypi_0    pypi
brotli                    1.1.0                hd74edd7_2    conda-forge
brotli-bin                1.1.0                hd74edd7_2    conda-forge
brotli-python             1.1.0           py310hb4ad77e_2    conda-forge
bwidget                   1.10.1               hce30654_0    conda-forge
bzip2                     1.0.8                h99b78c6_7    conda-forge
c-ares                    1.34.4               h5505292_0    conda-forge
ca-certificates           2024.12.14           hf0a4a13_0    conda-forge
cached-property           1.5.2                hd8ed1ab_1    conda-forge
cached_property           1.5.2              pyha770c72_1    conda-forge
cairo                     1.18.2               h6a3b0d2_1    conda-forge
cctools_osx-arm64         1010.6               h3f5b1a0_2    conda-forge
cell2cell                 0.7.4                    pypi_0    pypi
cellrank                  2.0.6              pyhd8ed1ab_1    conda-forge
certifi                   2024.8.30                pypi_0    pypi
cffi                      1.17.1                   pypi_0    pypi
cfgv                      3.4.0                    pypi_0    pypi
charset-normalizer        3.3.2                    pypi_0    pypi
clang                     19.1.6          default_h474c9e2_0    conda-forge
clang-19                  19.1.6          default_hf90f093_0    conda-forge
clang_impl_osx-arm64      19.1.6              h253acbe_23    conda-forge
clang_osx-arm64           19.1.6              h07b0088_23    conda-forge
clangxx                   19.1.6          default_h1ffe849_0    conda-forge
clangxx_impl_osx-arm64    19.1.6              h3e0e5ee_23    conda-forge
clangxx_osx-arm64         19.1.6              h07b0088_23    conda-forge
click                     8.1.7           unix_pyh707e725_0    conda-forge
cloudpickle               3.0.0              pyhd8ed1ab_0    conda-forge
colorama                  0.4.6              pyhd8ed1ab_1    conda-forge
colorcet                  3.1.0                    pypi_0    pypi
comm                      0.2.2              pyhd8ed1ab_0    conda-forge
compiler-rt               19.1.6               hd2aecb6_0    conda-forge
compiler-rt_osx-arm64     19.1.6               h7969c41_0    conda-forge
contourpy                 1.3.0           py310h7306fd8_1    conda-forge
cryptography              43.0.3                   pypi_0    pypi
ctxcore                   0.2.0                    pypi_0    pypi
curl                      8.11.1               h73640d1_0    conda-forge
cycler                    0.12.1             pyhd8ed1ab_1    conda-forge
cython                    0.29.37                  pypi_0    pypi
cytoolz                   0.12.3          py310hd125d64_0    conda-forge
dask                      2024.8.2           pyhd8ed1ab_0    conda-forge
dask-core                 2024.8.2           pyhd8ed1ab_0    conda-forge
dask-expr                 1.1.13             pyhd8ed1ab_0    conda-forge
datashader                0.16.3                   pypi_0    pypi
debugpy                   1.8.5                    pypi_0    pypi
decorator                 5.1.1              pyhd8ed1ab_0    conda-forge
decoupler                 1.8.0                    pypi_0    pypi
defusedxml                0.7.1                    pypi_0    pypi
deprecated                1.2.14                   pypi_0    pypi
dill                      0.3.9                    pypi_0    pypi
distlib                   0.3.9                    pypi_0    pypi
distributed               2024.8.2           pyhd8ed1ab_0    conda-forge
docrep                    0.3.2              pyh44b312d_0    conda-forge
et-xmlfile                2.0.0                    pypi_0    pypi
exceptiongroup            1.2.2              pyhd8ed1ab_0    conda-forge
executing                 2.1.0              pyhd8ed1ab_0    conda-forge
face                      20.1.1                   pypi_0    pypi
fastjsonschema            2.20.0                   pypi_0    pypi
fcsparser                 0.2.8                    pypi_0    pypi
fftw                      3.3.10          mpi_mpich_hd148bab_10    conda-forge
filelock                  3.15.4                   pypi_0    pypi
firefox                   133.0                h286801f_0    conda-forge
font-ttf-dejavu-sans-mono 2.37                 hab24e00_0    conda-forge
font-ttf-inconsolata      3.000                h77eed37_0    conda-forge
font-ttf-source-code-pro  2.038                h77eed37_0    conda-forge
font-ttf-ubuntu           0.83                 h77eed37_3    conda-forge
fontconfig                2.15.0               h1383a14_1    conda-forge
fonts-conda-ecosystem     1                             0    conda-forge
fonts-conda-forge         1                             0    conda-forge
fonttools                 4.53.1                   pypi_0    pypi
fqdn                      1.5.1                    pypi_0    pypi
freetype                  2.12.1               hadb7bae_2    conda-forge
fribidi                   1.0.10               h27ca646_0    conda-forge
frozendict                2.4.6                    pypi_0    pypi
frozenlist                1.5.0                    pypi_0    pypi
fsspec                    2024.9.0           pyhff2d567_0    conda-forge
future                    1.0.0              pyhd8ed1ab_1    conda-forge
gdown                     5.2.0                    pypi_0    pypi
geckodriver               0.35.0               hc069d6b_0    conda-forge
gensim                    4.3.3                    pypi_0    pypi
geos                      3.13.0               hf9b8971_0    conda-forge
gflags                    2.2.2             hc88da5d_1004    conda-forge
gfortran_impl_osx-arm64   13.2.0               h30f4408_1    conda-forge
gfortran_osx-arm64        13.2.0               h57527a5_1    conda-forge
glog                      0.7.1                heb240a5_0    conda-forge
glom                      23.5.0                   pypi_0    pypi
gmp                       6.3.0                h7bae524_2    conda-forge
gprofiler-official        1.0.0                    pypi_0    pypi
graphite2                 1.3.13            hebf3989_1003    conda-forge
graphtools                1.5.3                    pypi_0    pypi
gseapy                    1.0.3                    pypi_0    pypi
gsl                       2.7                  h6e638da_0    conda-forge
h11                       0.14.0                   pypi_0    pypi
h2                        4.1.0              pyhd8ed1ab_0    conda-forge
h5py                      3.11.0                   pypi_0    pypi
harfbuzz                  10.1.0               h9df47df_0    conda-forge
harmonypy                 0.0.10                   pypi_0    pypi
hdf5                      1.14.3          mpi_mpich_h13a04de_8    conda-forge
hnswlib                   0.7.0           py310he53c7d2_2    conda-forge
holoviews                 1.19.1                   pypi_0    pypi
hpack                     4.0.0              pyh9f0ad1d_0    conda-forge
httpcore                  1.0.5                    pypi_0    pypi
httpx                     0.27.2                   pypi_0    pypi
hyperframe                6.0.1              pyhd8ed1ab_0    conda-forge
hypre                     2.32.0          mpi_mpich_h189fe77_1    conda-forge
icu                       75.1                 hfee45f7_0    conda-forge
identify                  2.6.1                    pypi_0    pypi
idna                      3.8                      pypi_0    pypi
igraph                    0.11.6                   pypi_0    pypi
imageio                   2.35.1                   pypi_0    pypi
importlib-metadata        8.4.0              pyha770c72_0    conda-forge
importlib_metadata        8.4.0                hd8ed1ab_0    conda-forge
inflect                   7.4.0                    pypi_0    pypi
interlap                  0.2.7                    pypi_0    pypi
ipykernel                 6.29.5                   pypi_0    pypi
ipython                   8.27.0             pyh707e725_0    conda-forge
ipywidgets                8.1.5              pyhd8ed1ab_0    conda-forge
isl                       0.25                 h9a09cb3_0    conda-forge
isoduration               20.11.0                  pypi_0    pypi
jax                       0.5.0                    pypi_0    pypi
jaxlib                    0.5.0                    pypi_0    pypi
jaxopt                    0.8.3                    pypi_0    pypi
jedi                      0.19.1             pyhd8ed1ab_0    conda-forge
jinja2                    3.1.4              pyhd8ed1ab_0    conda-forge
joblib                    1.4.2              pyhd8ed1ab_1    conda-forge
json5                     0.9.25                   pypi_0    pypi
jsonpointer               3.0.0                    pypi_0    pypi
jsonschema                4.23.0                   pypi_0    pypi
jsonschema-specifications 2023.12.1                pypi_0    pypi
jupyter-client            8.6.2                    pypi_0    pypi
jupyter-core              5.7.2                    pypi_0    pypi
jupyter-events            0.10.0                   pypi_0    pypi
jupyter-lsp               2.2.5                    pypi_0    pypi
jupyter-server            2.14.2                   pypi_0    pypi
jupyter-server-terminals  0.5.3                    pypi_0    pypi
jupyterlab                4.2.5                    pypi_0    pypi
jupyterlab-pygments       0.3.0                    pypi_0    pypi
jupyterlab-server         2.27.3                   pypi_0    pypi
jupyterlab_widgets        3.0.13             pyhd8ed1ab_0    conda-forge
kaleido                   0.2.1                    pypi_0    pypi
kiwisolver                1.4.7           py310h7306fd8_0    conda-forge
kneed                     0.8.5                    pypi_0    pypi
krb5                      1.21.3               h237132a_0    conda-forge
lazy-loader               0.4                      pypi_0    pypi
lcms2                     2.16                 ha0e7c42_0    conda-forge
ld64_osx-arm64            951.9                hb91ea2e_2    conda-forge
legacy-api-wrap           1.4                      pypi_0    pypi
leidenalg                 0.10.2                   pypi_0    pypi
lerc                      4.0.0                h9a09cb3_0    conda-forge
liana                     1.4.0                    pypi_0    pypi
libabseil                 20240116.2      cxx17_h00cdb27_1    conda-forge
libaec                    1.1.3                hebf3989_0    conda-forge
libamd                    3.3.3            ss783_h6dbf161    conda-forge
libarrow                  17.0.0          h20538ec_13_cpu    conda-forge
libarrow-acero            17.0.0          hf9b8971_13_cpu    conda-forge
libarrow-dataset          17.0.0          hf9b8971_13_cpu    conda-forge
libarrow-substrait        17.0.0          hbf8b706_13_cpu    conda-forge
libasprintf               0.22.5               h8414b35_3    conda-forge
libblas                   3.9.0           23_osxarm64_openblas    conda-forge
libbrotlicommon           1.1.0                hd74edd7_2    conda-forge
libbrotlidec              1.1.0                hd74edd7_2    conda-forge
libbrotlienc              1.1.0                hd74edd7_2    conda-forge
libbtf                    2.3.2            ss783_h6c9afe8    conda-forge
libcamd                   3.3.3            ss783_h6c9afe8    conda-forge
libcblas                  3.9.0           23_osxarm64_openblas    conda-forge
libccolamd                3.3.4            ss783_h6c9afe8    conda-forge
libcholmod                5.3.0            ss783_h87d6651    conda-forge
libclang-cpp19.1          19.1.6          default_hf90f093_0    conda-forge
libcolamd                 3.3.4            ss783_h6c9afe8    conda-forge
libcrc32c                 1.1.2                hbdafb3b_0    conda-forge
libcurl                   8.11.1               h73640d1_0    conda-forge
libcxx                    19.1.6               ha82da77_1    conda-forge
libcxx-devel              19.1.6               h6dc3340_1    conda-forge
libdeflate                1.21                 h99b78c6_0    conda-forge
libedit                   3.1.20191231         hc8eb9b7_2    conda-forge
libev                     4.33                 h93a5062_2    conda-forge
libevent                  2.1.12               h2757513_1    conda-forge
libexpat                  2.6.4                h286801f_0    conda-forge
libffi                    3.4.2                h3422bc3_5    conda-forge
libgettextpo              0.22.5               h8414b35_3    conda-forge
libgfortran               5.0.0           13_2_0_hd922786_3    conda-forge
libgfortran-devel_osx-arm64 13.2.0               h5d7a38c_3    conda-forge
libgfortran5              13.2.0               hf226fd6_3    conda-forge
libglib                   2.82.2               h07bd6cf_0    conda-forge
libgoogle-cloud           2.28.0               hfe08963_0    conda-forge
libgoogle-cloud-storage   2.28.0               h1466eeb_0    conda-forge
libgrpc                   1.62.2               h9c18a4f_0    conda-forge
libhwloc                  2.11.2          default_hbce5d74_1001    conda-forge
libiconv                  1.17                 h0d3ecfb_2    conda-forge
libintl                   0.22.5               h8414b35_3    conda-forge
libjpeg-turbo             3.0.0                hb547adb_1    conda-forge
libklu                    2.3.5            ss783_h4a7adf4    conda-forge
liblapack                 3.9.0           23_osxarm64_openblas    conda-forge
libllvm14                 14.0.6               hd1a9a77_4    conda-forge
libllvm19                 19.1.6               hc4b4ae8_0    conda-forge
liblzma                   5.6.3                h39f12f2_1    conda-forge
libnghttp2                1.64.0               h6d7220d_0    conda-forge
libopenblas               0.3.27          openmp_h517c56d_1    conda-forge
libparquet                17.0.0          hf0ba9ef_13_cpu    conda-forge
libpng                    1.6.44               hc14010f_0    conda-forge
libprotobuf               4.25.3               hbfab5d5_0    conda-forge
libptscotch               7.0.6                h1b0e71d_0    conda-forge
libre2-11                 2023.09.01           h7b2c953_2    conda-forge
libscotch                 7.0.6                hf397248_0    conda-forge
libspqr                   4.3.4            ss783_h93d26d6    conda-forge
libsqlite                 3.46.1               hc14010f_0    conda-forge
libssh2                   1.11.1               h9cc3647_0    conda-forge
libsuitesparseconfig      7.8.3            ss783_h714a54a    conda-forge
libthrift                 0.20.0               h64651cc_1    conda-forge
libtiff                   4.7.0                h9c1d414_0    conda-forge
libumfpack                6.3.5            ss783_h852ec90    conda-forge
libutf8proc               2.8.0                h1a8c8d9_0    conda-forge
libwebp-base              1.4.0                h93a5062_0    conda-forge
libxcb                    1.17.0               hdb1d25a_0    conda-forge
libxml2                   2.13.5               h178c5d8_1    conda-forge
libzlib                   1.3.1                hfb2fe0b_1    conda-forge
linkify-it-py             2.0.3                    pypi_0    pypi
llvm-openmp               19.1.6               hdb05f8b_0    conda-forge
llvm-tools                19.1.6               hd2aecb6_0    conda-forge
llvm-tools-19             19.1.6               h87a4c7e_0    conda-forge
llvmlite                  0.43.0                   pypi_0    pypi
locket                    1.0.0              pyhd8ed1ab_0    conda-forge
loompy                    3.0.7                    pypi_0    pypi
louvain                   0.8.2                    pypi_0    pypi
lxml                      5.3.0                    pypi_0    pypi
lz4                       4.3.3           py310hc798581_1    conda-forge
lz4-c                     1.9.4                hb7217d7_0    conda-forge
macs2                     2.2.9.1                  pypi_0    pypi
magic-impute              3.0.0                    pypi_0    pypi
make                      4.4.1                hc9fafa5_2    conda-forge
markdown                  3.7                      pypi_0    pypi
markdown-it-py            3.0.0                    pypi_0    pypi
markupsafe                2.1.5           py310h493c2e1_1    conda-forge
matplotlib                3.9.2                    pypi_0    pypi
matplotlib-base           3.10.0          py310hadbac3a_0    conda-forge
matplotlib-inline         0.1.7              pyhd8ed1ab_0    conda-forge
mdit-py-plugins           0.4.1                    pypi_0    pypi
mdurl                     0.1.2                    pypi_0    pypi
mellon                    1.5.0                    pypi_0    pypi
memento                   0.1.0                    pypi_0    pypi
metis                     5.1.0             h15f6cfe_1007    conda-forge
mistune                   3.0.2                    pypi_0    pypi
mizani                    0.11.4                   pypi_0    pypi
ml-dtypes                 0.5.1                    pypi_0    pypi
more-itertools            10.5.0                   pypi_0    pypi
mpc                       1.3.1                h8f1351a_1    conda-forge
mpfr                      4.2.1                hb693164_3    conda-forge
mpi                       1.0.1                     mpich    conda-forge
mpi4py                    4.0.1           py310h2e9cc40_1    conda-forge
mpich                     4.2.3              h6b0954e_101    conda-forge
mpmath                    1.3.0                    pypi_0    pypi
msgpack-python            1.0.8           py310h7306fd8_1    conda-forge
mudata                    0.3.1                    pypi_0    pypi
multidict                 6.1.0                    pypi_0    pypi
multipledispatch          1.0.0                    pypi_0    pypi
multiprocessing-on-dill   3.5.0a4                  pypi_0    pypi
mumps-include             5.7.3                hce30654_6    conda-forge
mumps-mpi                 5.7.3                h8349c92_6    conda-forge
munkres                   1.1.4              pyh9f0ad1d_0    conda-forge
mygene                    3.2.2                    pypi_0    pypi
natsort                   8.4.0              pyh29332c3_1    conda-forge
nbclient                  0.10.0                   pypi_0    pypi
nbconvert                 7.16.4                   pypi_0    pypi
nbformat                  5.10.4                   pypi_0    pypi
ncurses                   6.5                  h7bae524_1    conda-forge
nest-asyncio              1.6.0                    pypi_0    pypi
networkx                  3.3                      pypi_0    pypi
nodeenv                   1.9.1                    pypi_0    pypi
notebook                  7.2.2                    pypi_0    pypi
notebook-shim             0.2.4                    pypi_0    pypi
nptyping                  2.5.0                    pypi_0    pypi
numba                     0.60.0          py310h0628f0e_0    conda-forge
numba-progress            1.1.0                    pypi_0    pypi
numexpr                   2.10.2                   pypi_0    pypi
numpy                     1.26.4          py310hd45542a_0    conda-forge
numpy_groupies            0.11.2             pyhd8ed1ab_1    conda-forge
omnipath                  1.0.8                    pypi_0    pypi
openjpeg                  2.5.3                h8a3d83b_0    conda-forge
openpyxl                  3.1.5                    pypi_0    pypi
openssl                   3.4.0                h81ee809_1    conda-forge
opt-einsum                3.4.0                    pypi_0    pypi
orc                       2.0.2                h75dedd0_0    conda-forge
outcome                   1.3.0.post0              pypi_0    pypi
overrides                 7.7.0                    pypi_0    pypi
packaging                 24.1               pyhd8ed1ab_0    conda-forge
palantir                  1.3.6                    pypi_0    pypi
pandas                    2.2.3                    pypi_0    pypi
pandocfilters             1.5.1                    pypi_0    pypi
panel                     1.4.5                    pypi_0    pypi
pango                     1.54.0               h73f1e88_4    conda-forge
param                     2.1.1                    pypi_0    pypi
paramiko                  3.5.0                    pypi_0    pypi
parmetis                  4.0.3             ha4b917a_1007    conda-forge
parso                     0.8.4              pyhd8ed1ab_0    conda-forge
partd                     1.4.2              pyhd8ed1ab_0    conda-forge
patsy                     0.5.6                    pypi_0    pypi
pcre2                     10.44                h297a79d_2    conda-forge
pecanpy                   2.0.9                    pypi_0    pypi
petsc                     3.22.2          real_hffa6677_103    conda-forge
petsc4py                  3.22.2          py310h6f8b1d8_0    conda-forge
pexpect                   4.9.0              pyhd8ed1ab_0    conda-forge
phate                     1.0.11                   pypi_0    pypi
phenograph                1.5.7                    pypi_0    pypi
pickleshare               0.7.5                   py_1003    conda-forge
pillow                    11.0.0          py310h530beaf_0    conda-forge
pip                       24.2               pyh8b19718_1    conda-forge
pixman                    0.44.2               h2f9eb0b_0    conda-forge
platformdirs              4.2.2                    pypi_0    pypi
plotly                    5.24.1                   pypi_0    pypi
plotnine                  0.13.6                   pypi_0    pypi
pre-commit                4.0.1                    pypi_0    pypi
progressbar2              4.5.0              pyhd8ed1ab_1    conda-forge
prometheus-client         0.20.0                   pypi_0    pypi
prompt-toolkit            3.0.47             pyha770c72_0    conda-forge
propcache                 0.2.1                    pypi_0    pypi
psutil                    6.0.0           py310ha6dd24b_0    conda-forge
pthread-stubs             0.4               h27ca646_1001    conda-forge
ptyprocess                0.7.0              pyhd3deb0d_0    conda-forge
pure_eval                 0.2.3              pyhd8ed1ab_0    conda-forge
pyarrow                   17.0.0          py310h24597f5_1    conda-forge
pyarrow-core              17.0.0          py310hf3d4daf_1_cpu    conda-forge
pyarrow-hotfix            0.6                pyhd8ed1ab_0    conda-forge
pybind11                  2.13.5                   pypi_0    pypi
pycairo                   1.27.0                   pypi_0    pypi
pycparser                 2.22               pyhd8ed1ab_0    conda-forge
pyct                      0.5.0                    pypi_0    pypi
pycurl                    7.45.3                   pypi_0    pypi
pydeseq2                  0.4.11                   pypi_0    pypi
pygam                     0.9.1              pyhd8ed1ab_1    conda-forge
pygments                  2.18.0             pyhd8ed1ab_0    conda-forge
pygpcca                   1.0.4              pyhd8ed1ab_2    conda-forge
pygsp                     0.5.1                    pypi_0    pypi
pynacl                    1.5.0                    pypi_0    pypi
pynndescent               0.5.13             pyhd8ed1ab_1    conda-forge
pyparsing                 3.1.4                    pypi_0    pypi
pypath-common             0.2.5                    pypi_0    pypi
pypath-omnipath           0.16.17                  pypi_0    pypi
pyreadr                   0.5.2                    pypi_0    pypi
pyscenic                  0.12.1                   pypi_0    pypi
pysftp                    0.2.9                    pypi_0    pypi
pysocks                   1.7.1              pyha2e5f31_6    conda-forge
python                    3.10.14         h2469fbe_0_cpython    conda-forge
python-dateutil           2.9.0.post0              pypi_0    pypi
python-json-logger        2.0.7                    pypi_0    pypi
python-slugify            8.0.4                    pypi_0    pypi
python-tzdata             2024.1             pyhd8ed1ab_0    conda-forge
python-utils              3.8.2                    pypi_0    pypi
python_abi                3.10                    5_cp310    conda-forge
pytz                      2024.1             pyhd8ed1ab_0    conda-forge
pyvia                     0.2.4                    pypi_0    pypi
pyviz-comms               3.0.3                    pypi_0    pypi
pyyaml                    6.0.2           py310h493c2e1_1    conda-forge
pyzmq                     26.2.0                   pypi_0    pypi
qhull                     2020.2               h420ef59_5    conda-forge
r-base                    4.3.3               h527b63a_16    conda-forge
rdata                     0.11.2                   pypi_0    pypi
re2                       2023.09.01           h4cba328_2    conda-forge
readline                  8.2                  h92ec313_1    conda-forge
referencing               0.35.1                   pypi_0    pypi
requests                  2.32.3                   pypi_0    pypi
rfc3339-validator         0.1.4                    pypi_0    pypi
rfc3986-validator         0.1.1                    pypi_0    pypi
rpds-py                   0.20.0                   pypi_0    pypi
rpy2                      3.5.11          py310r43h280b8fa_3    conda-forge
s-gd2                     1.8.1                    pypi_0    pypi
scalapack                 2.2.0                h71a4f75_4    conda-forge
scanpro                   0.3.2                    pypi_0    pypi
scanpy                    1.10.2                   pypi_0    pypi
scikit-image              0.24.0                   pypi_0    pypi
scikit-learn              1.5.1                    pypi_0    pypi
scipy                     1.11.4                   pypi_0    pypi
scprep                    1.2.3                    pypi_0    pypi
scvelo                    0.3.3              pyhd8ed1ab_0    conda-forge
seaborn                   0.13.2               hd8ed1ab_3    conda-forge
seaborn-base              0.13.2             pyhd8ed1ab_3    conda-forge
selenium                  4.27.1                   pypi_0    pypi
send2trash                1.8.3                    pypi_0    pypi
session-info              1.0.0              pyhd8ed1ab_0    conda-forge
setuptools                73.0.1             pyhd8ed1ab_0    conda-forge
sh                        2.2.1                    pypi_0    pypi
shapely                   2.0.6           py310h6b3522b_2    conda-forge
sigtool                   0.1.3                h44b9a77_0    conda-forge
simplegeneric             0.8.1              pyhd8ed1ab_2    conda-forge
six                       1.16.0             pyh6c4a22f_0    conda-forge
slepc                     3.22.2          real_hf0b33f3_300    conda-forge
slepc4py                  3.22.2          py310h61b750a_0    conda-forge
smart-open                7.0.4                    pypi_0    pypi
snappy                    1.2.1                hd02b534_0    conda-forge
sniffio                   1.3.1                    pypi_0    pypi
sortedcontainers          2.4.0              pyhd8ed1ab_0    conda-forge
soupsieve                 2.6                      pypi_0    pypi
sqlparse                  0.5.1                    pypi_0    pypi
stack-data                0.6.3                    pypi_0    pypi
stack_data                0.6.2              pyhd8ed1ab_0    conda-forge
statannotations           0.4.4                    pypi_0    pypi
statsmodels               0.14.2                   pypi_0    pypi
stdlib-list               0.10.0                   pypi_0    pypi
stream                    1.1                      pypi_0    pypi
superlu                   5.2.2                hc615359_0    conda-forge
superlu_dist              9.1.0                h89afcdd_0    conda-forge
sympy                     1.13.1                   pypi_0    pypi
tabulate                  0.9.0                    pypi_0    pypi
tapi                      1300.6.5             h03f4b80_0    conda-forge
tasklogger                1.2.0                    pypi_0    pypi
tbb                       2022.0.0             h0cbf7ec_0    conda-forge
tblib                     3.0.0              pyhd8ed1ab_0    conda-forge
tenacity                  9.0.0                    pypi_0    pypi
tensorly                  0.8.1                    pypi_0    pypi
termcolor                 2.4.0                    pypi_0    pypi
terminado                 0.18.1                   pypi_0    pypi
text-unidecode            1.3                      pypi_0    pypi
texttable                 1.7.0                    pypi_0    pypi
threadpoolctl             3.5.0              pyhc1e730c_0    conda-forge
tifffile                  2024.8.30                pypi_0    pypi
timeloop                  1.0.2                    pypi_0    pypi
tinycss2                  1.3.0                    pypi_0    pypi
tk                        8.6.13               h5083fa2_1    conda-forge
tktable                   2.10                 h1e387b8_6    conda-forge
toml                      0.10.2                   pypi_0    pypi
tomli                     2.0.1                    pypi_0    pypi
toolz                     0.12.1             pyhd8ed1ab_0    conda-forge
torch                     2.5.0                    pypi_0    pypi
torchaudio                2.5.0                    pypi_0    pypi
torchvision               0.20.0                   pypi_0    pypi
tornado                   6.4.1           py310h493c2e1_1    conda-forge
tqdm                      4.66.5                   pypi_0    pypi
traitlets                 5.14.3             pyhd8ed1ab_0    conda-forge
trio                      0.27.0                   pypi_0    pypi
trio-websocket            0.11.1                   pypi_0    pypi
typeguard                 4.3.0                    pypi_0    pypi
types-python-dateutil     2.9.0.20240821           pypi_0    pypi
typing_extensions         4.12.2             pyha770c72_0    conda-forge
tzdata                    2024a                h8827d51_1    conda-forge
tzlocal                   5.2             py310hbe9552e_1    conda-forge
uc-micro-py               1.0.3                    pypi_0    pypi
umap-learn                0.5.6                    pypi_0    pypi
unicodedata2              15.1.0          py310hf9df320_1    conda-forge
uri-template              1.3.0                    pypi_0    pypi
urllib3                   2.2.2              pyhd8ed1ab_1    conda-forge
virtualenv                20.26.6                  pypi_0    pypi
wcwidth                   0.2.13             pyhd8ed1ab_0    conda-forge
webcolors                 24.8.0                   pypi_0    pypi
webencodings              0.5.1                    pypi_0    pypi
websocket-client          1.8.0                    pypi_0    pypi
wget                      3.2                      pypi_0    pypi
wheel                     0.44.0             pyhd8ed1ab_0    conda-forge
widgetsnbextension        4.0.13             pyhd8ed1ab_0    conda-forge
wrapt                     1.16.0                   pypi_0    pypi
wsproto                   1.2.0                    pypi_0    pypi
xarray                    2024.7.0                 pypi_0    pypi
xlrd                      2.0.1                    pypi_0    pypi
xorg-libxau               1.0.11               hb547adb_0    conda-forge
xorg-libxdmcp             1.1.3                h27ca646_0    conda-forge
xyzservices               2024.9.0           pyhd8ed1ab_0    conda-forge
xz                        5.2.6                h57fd34a_0    conda-forge
yaml                      0.2.5                h3422bc3_2    conda-forge
yarl                      1.18.3                   pypi_0    pypi
zict                      3.0.0              pyhd8ed1ab_0    conda-forge
zipp                      3.20.1             pyhd8ed1ab_0    conda-forge
zlib                      1.3.1                hfb2fe0b_1    conda-forge
zstandard                 0.23.0          py310h2665a74_1    conda-forge
zstd                      1.5.6                hb46c0d2_0    conda-forge
```
Other package

star 2.7.10a cellranger 6.1.2 rsem 1.3.3 picard 2.10.9 ChromHMM 1.24 pyscenic 0.12.1
