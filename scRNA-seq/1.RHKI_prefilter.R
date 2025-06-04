library(Seurat)
library(DoubletFinder)
library(harmony)
library(dplyr)
library(patchwork)
library(scCustomize)
library(ggplot2)
library(cellmatch)
library(biomaRt)
library(GO.db)

###gene excluded
ribo.genes <- c(grep(pattern = "^Rpl", x = rownames(RH@assays$RNA@data), value = TRUE),grep(pattern = "^Rps", x = rownames(RH@assays$RNA@data), value = TRUE))
mito.genes = c(grep(pattern = "^mt", x = rownames(RH@assays$RNA@data), value = TRUE), grep(pattern = "^Mrpl", x = rownames(RH@assays$RNA@data), value = TRUE), grep(pattern = "^Mrps", x = rownames(RH@assays$RNA@data), value = TRUE))
Gm.genes = c(grep(pattern = "Gm", x = rownames(RH@assays$RNA@data), value = TRUE))
GOBPOFFSPRING[["GO:0003735"]]
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
View(ensembl@filters)
gene.data <- getBM(attributes=c('ensembl_gene_id','external_gene_name', 'go_id'), filters = 'go', values = 'GO:0003735', mart = ensembl)
exclude.gene = c(mito.genes, Gm.genes, unique(gene.data$external_gene_name))
exclude.gene = unique(exclude.gene)

###cell cycle genes
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv")
cell_cycle_genes <- read.csv(text = cc_file)
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"), ignore.case = TRUE)
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
edb <- ah[[id]]
annotations <- genes(edb, return.type = "data.frame")
annotations <- annotations %>% dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
s_genes <- cell_cycle_markers %>% dplyr::filter(phase == "S") %>% pull("gene_name")
g2m_genes <- cell_cycle_markers %>% dplyr::filter(phase == "G2/M") %>% pull("gene_name")

###doublet prediction
dirs = list.files("/path/to/cellranger/outputfolder/", pattern = "Results")
samples = gsub("_Results", "", dirs)
for (i in 1:length(samples)) {
  tmp = Seurat::Read10X(paste0("/path/to/cellranger/outputfolder/", dirs[i], "/outs/filtered_feature_bc_matrix")) ### contains *_barcodes.tsv.gz, *_feature.tsv.gz and *_matrix.mtx.gz files
  tmp = CreateSeuratObject(tmp) %>% NormalizeData() %>% FindVariableFeatures()%>%
    ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
  
  sweep.res.obj = paramSweep_v3(tmp, PCs = 1:10, sct = F)
  sweep.stats_obj <- summarizeSweep(sweep.res.obj, GT = FALSE)
  bcmvn <- find.pK(sweep.stats_obj)
  barplot(bcmvn$BCmetric, names.arg =  bcmvn$pK, las=2)
  x = as.numeric(levels(bcmvn$pK)[which.max(bcmvn$BCmetric)])
  nExp = ncol(tmp)*0.08 # estimated doublet rate
  obj_db = doubletFinder_v3(tmp, PCs = 1:10, pN = 0.25, pK = x, nExp = nExp, sct =F)
  metadf = obj_db@meta.data
  df = metadf[,grep("pANN|DF",colnames(metadf))]
  write.csv(df, paste0("/path/to/Doublet/",samples[i],".csv"),row.names = T)
  
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^mt-")
  tmp[["percent.Rb"]] <- PercentageFeatureSet(tmp, pattern = "^Rpl|Rps")
  tmp <- CellCycleScoring(tmp,g2m.features = g2m_genes,s.features = s_genes)
  tmp$CC.Difference = tmp$S.Score - tmp$G2M.Score
  saveRDS(tmp, paste0("/path/to/prefilter/",samples[i],".rds"))
}

### prefilter individual samples
### For detail sample to batch information, read scRNA-seq_info.csv
### remove small cluster which are also with a higher mitochondria counts or high ribosome counts
### percent.MT <10 and nFeature_RNA > 1800 as cutoff
### Since we are going to merge our data with two other public data and the sequencing depth is different among different dataset, we do give a upper bound of nFeature_RNA to make everything in the same cutoff range.

### Rnaseh1_D0_Dox_rep2####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D0_Dox_rep2.rds")
tmpdb = read.csv("/path/to//Doublet/Rnaseh1_D0_Dox_rep2.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.03_835.84
tmp$doublet = tmpdb$DF.classifications_0.25_0.03_835.84
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents = c(0,2,10), invert = T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]##doublet!!
tmp_metadf$orig.ident = "Rnaseh1_D0_Dox_rep2"
tmp_metadf$replicate = "rep8"
tmp_metadf$treatment="Dox"
tmp_metadf$cell="RHKI"
tmp_metadf$resource = "exp"
tmp_metadf$time = "0h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D0_Dox_rep2:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D0_Dox_rep2:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D0_Dox_rep2.rds")

### Rnaseh1_D0_Dox_rep4 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D0_Dox_rep4.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D0_Dox_rep4.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.18_924.56
tmp$doublet = tmpdb$DF.classifications_0.25_0.18_924.56
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D0_Dox_rep4"
tmp_metadf$replicate = "rep10"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="Dox"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "0h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D0_Dox_rep4:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D0_Dox_rep4:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D0_Dox_rep4.rds")

### Rnaseh1_D0_N_rep2 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D0_N_rep2.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D0_N_rep2.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.3_949.6
tmp$doublet = tmpdb$DF.classifications_0.25_0.3_949.6
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(0,1,9,11), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]##doublet!
tmp_metadf$orig.ident = "Rnaseh1_D0_N_rep2"
tmp_metadf$replicate = "rep8"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="control"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "0h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D0_N_rep2:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D0_N_rep2:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D0_N_rep2.rds")

### Rnaseh1_D0_N_rep4 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D0_N_rep4.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D0_N_rep4.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.29_877.68
tmp$doublet = tmpdb$DF.classifications_0.25_0.29_877.68
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D0_N_rep4"
tmp_metadf$replicate = "rep10"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="control"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "0h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D0_N_rep4:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D0_N_rep4:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D0_N_rep4.rds")

### Rnaseh1_D3_Dox_rep1 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D3_Dox_rep1.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D3_Dox_rep1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.17_860.88
tmp$doublet = tmpdb$DF.classifications_0.25_0.17_860.88
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(5,11), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]##doublet!
tmp_metadf$orig.ident = "Rnaseh1_D3_Dox_rep1"
tmp_metadf$replicate = "rep7"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="Dox"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D3_Dox_rep1:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D3_Dox_rep1:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D3_Dox_rep1.rds")



### Rnaseh1_D3_Dox_rep4 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D3_Dox_rep4.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D3_Dox_rep4.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.26_640.56
tmp$doublet = tmpdb$DF.classifications_0.25_0.26_640.56
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]##doublet!
tmp_metadf$orig.ident = "Rnaseh1_D3_Dox_rep4"
tmp_metadf$replicate = "rep10"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="Dox"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D3_Dox_rep4:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D3_Dox_rep4:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D3_Dox_rep4.rds")

### Rnaseh1_D3_Dox_rep5 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D3_Dox_rep5.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D3_Dox_rep5.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.18_680.96
tmp$doublet = tmpdb$DF.classifications_0.25_0.18_680.96
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D3_Dox_rep5"
tmp_metadf$replicate = "rep11"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="Dox"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D3_Dox_rep5:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D3_Dox_rep5:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D3_Dox_rep5.rds")

### Rnaseh1_D3_Dox_rep6 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D3_Dox_rep6.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D3_Dox_rep6.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.12_595.12
tmp$doublet = tmpdb$DF.classifications_0.25_0.12_595.12
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(10), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]##doublet!
tmp_metadf$orig.ident = "Rnaseh1_D3_Dox_rep6"
tmp_metadf$replicate = "rep12"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="Dox"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D3_Dox_rep6:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D3_Dox_rep6:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D3_Dox_rep6.rds")

### Rnaseh1_D3_N_rep1 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D3_N_rep1.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D3_N_rep1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.3_514.72
tmp$doublet = tmpdb$DF.classifications_0.25_0.3_514.72
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(7,9), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D3_N_rep1"
tmp_metadf$replicate = "rep7"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="control"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D3_N_rep1:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D3_N_rep1:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D3_N_rep1.rds")



### Rnaseh1_D3_N_rep4 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D3_N_rep4.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D3_N_rep4.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.15_677.92
tmp$doublet = tmpdb$DF.classifications_0.25_0.15_677.92
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D3_N_rep4"
tmp_metadf$replicate = "rep10"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="control"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D3_N_rep4:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D3_N_rep4:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D3_N_rep4.rds")

### Rnaseh1_D3_N_rep5 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D3_N_rep5.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D3_N_rep5.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.25_610.16
tmp$doublet = tmpdb$DF.classifications_0.25_0.25_610.16
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(13), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D3_N_rep5"
tmp_metadf$replicate = "rep11"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="control"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D3_N_rep5:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D3_N_rep5:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D3_N_rep5.rds")

### Rnaseh1_D3_N_rep6 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D3_N_rep6.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D3_N_rep6.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.09_360.24
tmp$doublet = tmpdb$DF.classifications_0.25_0.09_360.24
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(11), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D3_N_rep6"
tmp_metadf$replicate = "rep12"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="control"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D3_N_rep6:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D3_N_rep6:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D3_N_rep6.rds")

### Rnaseh1_D5_Dox_rep1 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D5_Dox_rep1.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D5_Dox_rep1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.29_982.72
tmp$doublet = tmpdb$DF.classifications_0.25_0.29_982.72
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(6,8,13), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D5_Dox_rep1"
tmp_metadf$replicate = "rep7"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="Dox"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D5_Dox_rep1:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D5_Dox_rep1:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D5_Dox_rep1.rds")

### Rnaseh1_D5_Dox_rep4 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D5_Dox_rep4.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D5_Dox_rep4.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.16_593.6
tmp$doublet = tmpdb$DF.classifications_0.25_0.16_593.6
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(10), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]##doublet!
tmp_metadf$orig.ident = "Rnaseh1_D5_Dox_rep4"
tmp_metadf$replicate = "rep10"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="Dox"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D5_Dox_rep4:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D5_Dox_rep4:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D5_Dox_rep4.rds")



### Rnaseh1_D5_Dox_rep6 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D5_Dox_rep6.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D5_Dox_rep6.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.18_761.12
tmp$doublet = tmpdb$DF.classifications_0.25_0.18_761.12
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D5_Dox_rep6"
tmp_metadf$replicate = "rep12"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="Dox"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D5_Dox_rep6:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D5_Dox_rep6:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D5_Dox_rep6.rds")

### Rnaseh1_D5_Dox_rep7 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D5_Dox_rep7.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D5_Dox_rep7.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.3_644.16
tmp$doublet = tmpdb$DF.classifications_0.25_0.3_644.16
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D5_Dox_rep7"
tmp_metadf$replicate = "rep13"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="Dox"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D5_Dox_rep7:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D5_Dox_rep7:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D5_Dox_rep7.rds")

### Rnaseh1_D5_N_rep1 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D5_N_rep1.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D5_N_rep1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.01_858.16
tmp$doublet = tmpdb$DF.classifications_0.25_0.01_858.16
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(6,13,14), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D5_N_rep1"
tmp_metadf$replicate = "rep7"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="control"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D5_N_rep1:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D5_N_rep1:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D5_N_rep1.rds")

### Rnaseh1_D5_N_rep3 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D5_N_rep3.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D5_N_rep3.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.06_685.12
tmp$doublet = tmpdb$DF.classifications_0.25_0.06_685.12
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D5_N_rep3"
tmp_metadf$replicate = "rep9"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="control"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D5_N_rep3:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D5_N_rep3:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D5_N_rep3.rds")

### Rnaseh1_D5_N_rep4 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D5_N_rep4.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D5_N_rep4.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.28_674.32
tmp$doublet = tmpdb$DF.classifications_0.25_0.28_674.32
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(13), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D5_N_rep4"
tmp_metadf$replicate = "rep10"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="control"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D5_N_rep4:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D5_N_rep4:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D5_N_rep4.rds")


### Rnaseh1_D5_N_rep6 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D5_N_rep6.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D5_N_rep6.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.07_712.4
tmp$doublet = tmpdb$DF.classifications_0.25_0.07_712.4
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D5_N_rep6"
tmp_metadf$replicate = "rep12"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="control"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D5_N_rep6:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D5_N_rep6:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D5_N_rep6.rds")


### Rnaseh1_D5_N_rep7 ####
tmp = readRDS("/path/to/prefilter/Rnaseh1_D5_N_rep7.rds")
tmpdb = read.csv("/path/to/Doublet/Rnaseh1_D5_N_rep7.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.25_753.6
tmp$doublet = tmpdb$DF.classifications_0.25_0.25_753.6
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "doublet")
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "Rnaseh1_D5_N_rep7"
tmp_metadf$replicate = "rep13"
tmp_metadf$resource = "exp"
tmp_metadf$treatment="control"
tmp_metadf$cell="RHKI"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "Rnaseh1_D5_N_rep7:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "Rnaseh1_D5_N_rep7:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/post-filter/Rnaseh1_D5_N_rep7.rds")

###merge all our dataset
for (i in 1:length(samples)) {
  tmp = readRDS(paste0("/path/to/post-filter/",samples[i],".rds"))
  if (i==1) {
    obj = tmp} else {
      obj = merge(obj, y = tmp, merge.data=T)
    }
}
saveRDS(tmp, "/path/to/post-filter/RHKI_merge.rds")
