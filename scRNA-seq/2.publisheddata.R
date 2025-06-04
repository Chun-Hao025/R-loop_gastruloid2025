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
dirs = list.files("/path/to/cellranger/publisheddata/output/folder/", pattern = "Results")
samples = gsub("_Results", "", dirs)
for (i in 1:length(samples)) {
  tmp = Seurat::Read10X(paste0("/path/to/cellranger/publisheddata/output/folder/", dirs[i], "/outs/filtered_feature_bc_matrix")) ### contains *_barcodes.tsv.gz, *_feature.tsv.gz and *_matrix.mtx.gz files
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
  write.csv(df, paste0("/path/to/publisheddata/Doublet/",samples[i],".csv"),row.names = T)
  
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^mt-")
  tmp[["percent.Rb"]] <- PercentageFeatureSet(tmp, pattern = "^Rpl|Rps")
  tmp <- CellCycleScoring(tmp,g2m.features = g2m_genes,s.features = s_genes)
  tmp$CC.Difference = tmp$S.Score - tmp$G2M.Score
  saveRDS(tmp, paste0("/path/to/publisheddata/prefilter/",samples[i],".rds"))
}

### prefilter individual samples
### For detail sample to batch information, read scRNA-seq_info.csv
### remove small cluster which are also with a higher mitochondria counts or high ribosome counts
### percent.MT <10 and nFeature_RNA > 1800 as cutoff
### Since we are going to merge our data with two other published data and the sequencing depth is different among different dataset, we do give a upper bound of nFeature_RNA to make everything in the same cutoff range.
### 0 hr data from Suppinger et. al. (denoted as ref here) is not the real 0 hr data (express late stage marker and overlap with cells from late timepoints) so I exclude it

### ref_24h ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_24h.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_24h.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.09_414.08
tmp$doublet = tmpdb$DF.classifications_0.25_0.09_414.08
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents = c(9), invert = T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp1=subset(tmp1, idents=c(7), invert=T)
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_24h"
tmp_metadf$replicate = "rep1"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "24h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_24h_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_24h_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_24h.rds")

### ref_36h ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_36h.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_36h.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.12_694.88
tmp$doublet = tmpdb$DF.classifications_0.25_0.12_694.88
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
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,group.by = "doublet")
VlnPlot(tmp1, features = c("percent.mt"))
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_36h"
tmp_metadf$replicate = "rep1"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "36h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_36h_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_36h_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_36h.rds")

### ref_48h_1 ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_48h_1.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_48h_1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.04_488.72
tmp$doublet = tmpdb$DF.classifications_0.25_0.04_488.72
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
VlnPlot(tmp1, features = c("percent.mt"))
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_48h_1"
tmp_metadf$replicate = "rep1"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "48h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_48h_1_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_48h_1_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_48h_1.rds")

### ref_48h_2 ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_48h_2.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_48h_2.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.3_459.68
tmp$doublet = tmpdb$DF.classifications_0.25_0.3_459.68
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
VlnPlot(tmp1, features = c("percent.mt"))
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_48h_2"
tmp_metadf$replicate = "rep2"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "48h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_48h_2_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_48h_2_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_48h_2.rds")

### ref_52h ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_52h.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_52h.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.26_830.96
tmp$doublet = tmpdb$DF.classifications_0.25_0.26_830.96
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
VlnPlot(tmp1, features = c("percent.mt"))
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_52h"
tmp_metadf$replicate = "rep1"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "52h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_52h_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_52h_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_52h.rds")

### ref_56h ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_56h.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_56h.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.06_487.04
tmp$doublet = tmpdb$DF.classifications_0.25_0.06_487.04
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
VlnPlot(tmp1, features = c("percent.mt"))
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp1=subset(tmp1, idents=c(9), invert=T)
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_56h"
tmp_metadf$replicate = "rep1"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "56h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_56h_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_56h_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_56h.rds")

### ref_60h ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_60h.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_60h.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.3_434.32
tmp$doublet = tmpdb$DF.classifications_0.25_0.3_434.32
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(8), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_60h"
tmp_metadf$replicate = "rep1"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "60h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_60h_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_60h_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_60h.rds")

### ref_72h ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_72h.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_72h.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.12_428.72
tmp$doublet = tmpdb$DF.classifications_0.25_0.12_428.72
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
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_72h"
tmp_metadf$replicate = "rep1"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_72h_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_72h_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_72h.rds")

### ref_84h_1 ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_84h_1.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_84h_1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.14_362.48
tmp$doublet = tmpdb$DF.classifications_0.25_0.14_362.48
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
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_84h_1"
tmp_metadf$replicate = "rep1"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "84h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_84h_1_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_84h_1_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_84h_1.rds")

### ref_84h_2 ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_84h_2.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_84h_2.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.23_417.04
tmp$doublet = tmpdb$DF.classifications_0.25_0.23_417.04
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
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_84h_2"
tmp_metadf$replicate = "rep2"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "84h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_84h_2_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_84h_2_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_84h_2.rds")

### ref_96h ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_96h.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_96h.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.3_444.8
tmp$doublet = tmpdb$DF.classifications_0.25_0.3_444.8
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
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_96h"
tmp_metadf$replicate = "rep1"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "96h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_96h_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_96h_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_96h.rds")

### ref_108h_1 ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_108h_1.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_108h_1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.24_419.84
tmp$doublet = tmpdb$DF.classifications_0.25_0.24_419.84
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
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_108h_1"
tmp_metadf$replicate = "rep1"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "108h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_108h_1_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_108h_1_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_108h_1.rds")

### ref_108h_2 ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_108h_2.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_108h_2.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.13_398.88
tmp$doublet = tmpdb$DF.classifications_0.25_0.13_398.88
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
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_108h_2"
tmp_metadf$replicate = "rep2"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "108h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_108h_2_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_108h_2_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_108h_2.rds")

### ref_120h ####
tmp = readRDS("/path/to/publisheddata/prefilter/ref_120h.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/ref_120h.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.26_307.52
tmp$doublet = tmpdb$DF.classifications_0.25_0.26_307.52
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
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "ref_120h"
tmp_metadf$replicate = "rep1"
tmp_metadf$resource = "cell_resource"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "ref_120h_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "ref_120h_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/publisheddata/post-filter/ref_120h.rds")

### biox_72h_1 ####
tmp = readRDS("/path/to/publisheddata/prefilter/biox_72h_1.rds")
tmpdb = read.csv("/path/to/publisheddata/doublet/biox_72h_1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.07_456.8
tmp$doublet = tmpdb$DF.classifications_0.25_0.07_456.8
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 500)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(7,8,11), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_72h_1"
tmp_metadf$replicate = "rep3"
tmp_metadf$resource = "biox"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_72h_1_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_72h_1_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_72h_1.rds")

### biox_72h_2 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_72h_2.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_72h_2.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.21_683.2
tmp$doublet = tmpdb$DF.classifications_0.25_0.21_683.2
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
VlnPlot(tmp1, features = c("percent.mt"))
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1, group.by="doublet")
DimPlot(tmp1,label = T)+NoLegend()
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_72h_2"
tmp_metadf$replicate = "rep4"
tmp_metadf$resource = "biox"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_72h_2_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_72h_2_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_72h_2.rds")

### biox_72h_3 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_72h_3.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_72h_3.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.06_788.96
tmp$doublet = tmpdb$DF.classifications_0.25_0.06_788.96
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
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_72h_3"
tmp_metadf$replicate = "rep5"
tmp_metadf$resource = "biox"
tmp_metadf$time = "72h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_72h_3_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_72h_3_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_72h_3.rds")

### biox_90h_1 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_90h_1.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_90h_1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.04_629.2
tmp$doublet = tmpdb$DF.classifications_0.25_0.04_629.2
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
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_90h_1"
tmp_metadf$replicate = "rep3"
tmp_metadf$resource = "biox"
tmp_metadf$time = "90h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_90h_1_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_90h_1_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_90h_1.rds")

### biox_90h_2 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_90h_2.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_90h_2.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.16_552.96
tmp$doublet = tmpdb$DF.classifications_0.25_0.16_552.96
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
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_90h_2"
tmp_metadf$replicate = "rep4"
tmp_metadf$resource = "biox"
tmp_metadf$time = "90h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_90h_2_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_90h_2_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_90h_2.rds")

### biox_90h_3 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_90h_3.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_90h_3.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.25_621.12
tmp$doublet = tmpdb$DF.classifications_0.25_0.25_621.12
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
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_90h_3"
tmp_metadf$replicate = "rep5"
tmp_metadf$resource = "biox"
tmp_metadf$time = "90h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_90h_3_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_90h_3_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_90h_3.rds")

### biox_90h_4 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_90h_4.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_90h_4.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.3_582.64
tmp$doublet = tmpdb$DF.classifications_0.25_0.3_582.64
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_90h_4"
tmp_metadf$replicate = "rep6"
tmp_metadf$resource = "biox"
tmp_metadf$time = "90h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_90h_4_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_90h_4_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_90h_4.rds")

### biox_96h ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_96h.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_96h.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.14_379.52
tmp$doublet = tmpdb$DF.classifications_0.25_0.14_379.52
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(4), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_96h"
tmp_metadf$replicate = "rep3"
tmp_metadf$resource = "biox"
tmp_metadf$time = "96h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_96h_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_96h_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_96h.rds")


### biox_102h_1 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_102h_1.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_102h_1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.27_385.44
tmp$doublet = tmpdb$DF.classifications_0.25_0.27_385.44
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_102h_1"
tmp_metadf$replicate = "rep3"
tmp_metadf$resource = "biox"
tmp_metadf$time = "102h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_102h_1_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_102h_1_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_102h_1.rds")


### biox_102h_2 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_102h_2.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_102h_2.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.07_366.16
tmp$doublet = tmpdb$DF.classifications_0.25_0.07_366.16
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_102h_2"
tmp_metadf$replicate = "rep4"
tmp_metadf$resource = "biox"
tmp_metadf$time = "102h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_102h_2_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_102h_2_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_102h_2.rds")


### biox_102h_3 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_102h_3.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_102h_3.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.1_450.72
tmp$doublet = tmpdb$DF.classifications_0.25_0.1_450.72
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_102h_3"
tmp_metadf$replicate = "rep5"
tmp_metadf$resource = "biox"
tmp_metadf$time = "102h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_102h_3_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_102h_3_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_102h_3.rds")


### biox_102h_4 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_102h_4.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_102h_4.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.19_435.44
tmp$doublet = tmpdb$DF.classifications_0.25_0.19_435.44
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(13), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
DimPlot(tmp1,label = T)+NoLegend()
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,group.by = "doublet")
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
VlnPlot(tmp1, features = c("percent.mt"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_102h_4"
tmp_metadf$replicate = "rep6"
tmp_metadf$resource = "biox"
tmp_metadf$time = "102h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_102h_4_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_102h_4_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_102h_4.rds")
### biox_114h_1 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_114h_1.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_114h_1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.18_607.68
tmp$doublet = tmpdb$DF.classifications_0.25_0.18_607.68
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,group.by = "doublet")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_114h_1"
tmp_metadf$replicate = "rep3"
tmp_metadf$resource = "biox"
tmp_metadf$time = "114h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_114h_1_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_114h_1_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_114h_1.rds")

### biox_114h_2 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_114h_2.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_114h_2.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.07_633.84
tmp$doublet = tmpdb$DF.classifications_0.25_0.07_633.84
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,group.by = "doublet")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_114h_2"
tmp_metadf$replicate = "rep4"
tmp_metadf$resource = "biox"
tmp_metadf$time = "114h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_114h_2_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_114h_2_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_114h_2.rds")


### biox_114h_3 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_114h_3.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_114h_3.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.19_621.28
tmp$doublet = tmpdb$DF.classifications_0.25_0.19_621.28
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,group.by = "doublet")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_114h_3"
tmp_metadf$replicate = "rep5"
tmp_metadf$resource = "biox"
tmp_metadf$time = "114h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_114h_3_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_114h_3_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_114h_3.rds")


### biox_120h_1 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_120h_1.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_120h_1.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.29_541.44
tmp$doublet = tmpdb$DF.classifications_0.25_0.29_541.44
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, idents=c(8,9,12), invert=T)
tmp1 = subset(tmp1, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,group.by = "doublet")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_120h_1"
tmp_metadf$replicate = "rep3"
tmp_metadf$resource = "biox"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_120h_1_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_120h_1_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_120h_1.rds")

### biox_120h_2 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_120h_2.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_120h_2.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.26_629.52
tmp$doublet = tmpdb$DF.classifications_0.25_0.26_629.52
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,group.by = "doublet")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_120h_2"
tmp_metadf$replicate = "rep4"
tmp_metadf$resource = "biox"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_120h_2_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_120h_2_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_120h_2.rds")

### biox_120h_3 ####
tmp = readRDS("/path/to/cellranger/publisheddata/prefilter/biox_120h_3.rds")
tmpdb = read.csv("/path/to/cellranger/publisheddata/doublet/biox_120h_3.csv")
all(colnames(tmp)==tmpdb$X)
tmp$db_score = tmpdb$pANN_0.25_0.3_465.92
tmp$doublet = tmpdb$DF.classifications_0.25_0.3_465.92
DimPlot(tmp, label = T) + NoLegend()
VlnPlot(tmp, features = c("percent.mt"))+geom_hline(yintercept = 7.5)
DimPlot(tmp,group.by = "Phase")
DimPlot(tmp,group.by = "doublet")
VlnPlot(tmp, features = c("nCount_RNA"))+geom_hline(yintercept = 1000)+ geom_hline(yintercept = 30000)
VlnPlot(tmp, features = c("nFeature_RNA"))+geom_hline(yintercept = 1800)+ geom_hline(yintercept = 7000)
FeatureScatter(tmp , "nCount_RNA", "nFeature_RNA")
summary(tmp@meta.data)
Idents(tmp) = tmp$seurat_clusters
tmp1 = subset(tmp, subset = percent.mt < 10 & nFeature_RNA > 1800)
tmp1 = FindVariableFeatures(tmp1) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA() %>% RunUMAP(dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()
tmp1[["percent.Rb"]] <- PercentageFeatureSet(tmp1, pattern = "^Rpl|Rps")
VlnPlot(tmp1, features = c("percent.mt"))
DimPlot(tmp1,group.by = "Phase")
DimPlot(tmp1,group.by = "doublet")
DimPlot(tmp1,label = T)+NoLegend()
VlnPlot(tmp1, features = c("nCount_RNA", "nFeature_RNA"))
tmp_metadf = tmp1@meta.data
tmp_metadf = tmp_metadf[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score","Phase","db_score", "doublet", "percent.Rb")]
tmp_metadf$orig.ident = "biox_120h_3"
tmp_metadf$replicate = "rep5"
tmp_metadf$resource = "biox"
tmp_metadf$time = "120h"
tmp2 = tmp1@assays$RNA
colnames(tmp2@counts) = gsub("^", "biox_120h_3_Results:", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("^", "biox_120h_3_Results:", colnames(tmp2@data))
colnames(tmp2@counts) = gsub("-1", "x", colnames(tmp2@counts))
colnames(tmp2@data) = gsub("-1", "x", colnames(tmp2@data))
rownames(tmp_metadf) = colnames(tmp2)
tmp2 = CreateSeuratObject(tmp2)
tmp2@meta.data = tmp_metadf
saveRDS(tmp2, "/path/to/cellranger/publisheddata/post-filter/biox_120h_3.rds")

##merge all published data together
for (i in 1:length(samples)) {
  tmp = readRDS(paste0("/path/to/cellranger/publisheddata/post-filter/",samples[i],".rds"))
  if (i==1) {
    obj = tmp} else {
      obj = merge(obj, y = tmp, merge.data=T)
    }
}
saveRDS(tmp, "/path/to/cellranger/publisheddata/post-filter/ref_merge.rds")
library(GO.db)
library(biomaRt)
mito.genes = c(grep(pattern = "^mt", x = rownames(obj@assays$RNA@data), value = TRUE), grep(pattern = "^Mrpl", x = rownames(obj@assays$RNA@data), value = TRUE), grep(pattern = "^Mrps", x = rownames(obj@assays$RNA@data), value = TRUE))
Gm.genes = c(grep(pattern = "^Gm", x = rownames(obj@assays$RNA@data), value = TRUE))
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
View(ensembl@filters)
gene.data <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'go', values = 'GO:0003735', mart = ensembl)
exclude.gene = c(mito.genes, Gm.genes, unique(gene.data$external_gene_name))
exclude.gene = unique(exclude.gene)
mat1 = obj@assays$RNA@counts
tmp_adf=obj@meta.data
mat1 = mat1[!(rownames(mat1)%in%exclude.gene),]
ref=CreateSeuratObject(mat1)
ref@meta.data=tmp_adf
ref=subset(ref, subset= doublet %in% c("Singlet"))
ref = NormalizeData(ref) %>% FindVariableFeatures(nfeatures = 2000)
var=VariableFeatures(ref)
var1=c()
for (i in 1:length(samples)){
  tmp=readRDS(paste0("/Users/ChaoC1/Dropbox (UMass Medical School)/scRNA-seq-home/pub_cell_resource/intron/post-filter/", samples[i],".rds"))
  tmp = NormalizeData(tmp) %>% FindVariableFeatures(nfeatures = 500)
  tmp2=VariableFeatures(tmp)
  if (i==1) {
    var1 = tmp2} else {
      var1 = unique(c(var1, tmp2))
    }
}

####
set.seed(1234)
ref=subset(ref, subset = cell %in% c("E14Tg2a","CGR8","RHKI"))
ref=NormalizeData(ref) %>% FindVariableFeatures(nfeatures = 2000)
ref = ScaleData(ref, vars.to.regress =c("CC.Difference")) %>% RunPCA()
ElbowPlot(ref, ndims = 50)
ref=RunHarmony(ref, "replicate", reduction.save = "harmony_rep", max_iter=30)
ref = RunUMAP(ref, dims = 1:30, reduction = "harmony_rep") %>% FindNeighbors(dims = 1:30, reduction = "harmony_rep") %>% FindClusters()

harm=ref@reductions$harmony_rep@cell.embeddings
umap=ref@reductions$umap@cell.embeddings
meta_adf=ref@meta.data
write.csv(meta_adf, "/Users/ChaoC1/Dropbox (UMass Medical School)/scRNA-seq-home/20231216scRNA-seq_Rnaseh/via/via_meta_new2.csv")
write.csv(harm, "/Users/ChaoC1/Dropbox (UMass Medical School)/scRNA-seq-home/20231216scRNA-seq_Rnaseh/via/via_harmony_new2.csv")
write.csv(umap, "/Users/ChaoC1/Dropbox (UMass Medical School)/scRNA-seq-home/20231216scRNA-seq_Rnaseh/via/via_umap_new2.csv")
meta_adf=read.csv("/Users/ChaoC1/Dropbox (UMass Medical School)/scRNA-seq-home/20231216scRNA-seq_Rnaseh/via/via_meta_new2.csv", row.names = 1)

ref$PARC_noZX_1=meta_adf2$PARC_noZX1 #KNN=40
ref$PARC_noZX_2=meta_adf2$PARC_noZX2 #KNN=30
ref$PARC_noZX_3=meta_adf2$PARC_noZX2 #KNN=30
ref$PARC_noZX_1=paste0("PARC_", ref$PARC_noZX_1)
ref$PARC_noZX_2=paste0("PARC_", ref$PARC_noZX_2)
ref$PARC_noZX_3=paste0("PARC_", ref$PARC_noZX_3)

library(cellmatch)
ref_bulk1 = cellmatch::AggregateByCluster(expm1(Seurat::GetAssayData(ref, "RNA")), ref$PARC_noZX_1, method = "average")
ref_bulk2 = cellmatch::AggregateByCluster(expm1(Seurat::GetAssayData(ref, "RNA")), ref$PARC_noZX_2, method = "average")
ref_bulk3 = cellmatch::AggregateByCluster(expm1(Seurat::GetAssayData(ref, "RNA")), ref$ann_new_noZX_KNN40, method = "average")
ref_bulk1=ref_bulk1*100
ref_bulk2=ref_bulk2*100
ref_bulk3=ref_bulk3*100
embryo_bulk = readRDS("/Users/ChaoC1/Dropbox (UMass Medical School)/scRNA-seq-home/embryo/embryo_all_aggre.rds")

setwd("/Users/ChaoC1/Dropbox (UMass Medical School)/scRNA-seq-home/20231216scRNA-seq_Rnaseh/")
temp_dir1 = getwd() %>% file.path("PARC_noZX_1_2500")
temp_dir2 = getwd() %>% file.path("PARC_noZX_1_2000")
temp_dir3 = getwd() %>% file.path("PARC_noZX_2_2500")
temp_dir4 = getwd() %>% file.path("PARC_noZX_2_2000")
temp_dir5 = getwd() %>% file.path("ann_noZX_KNN40")
dir.create(temp_dir1)
dir.create(temp_dir2)
dir.create(temp_dir3)
dir.create(temp_dir4)
dir.create(temp_dir5)
matching_output1 = RunCellMatch(query = ref_bulk1,reference = embryo_bulk,results_path = temp_dir1,K = 2500,num_init = 50)
matching_output2 = RunCellMatch(query = ref_bulk1,reference = embryo_bulk,results_path = temp_dir2,K = 2000,num_init = 50)
matching_output3 = RunCellMatch(query = ref_bulk2,reference = embryo_bulk,results_path = temp_dir3,K = 2500,num_init = 50)
matching_output4 = RunCellMatch(query = ref_bulk2,reference = embryo_bulk,results_path = temp_dir4,K = 2000,num_init = 50)
matching_output5 = RunCellMatch(query = ref_bulk3,reference = embryo_bulk,results_path = temp_dir5,K = 2500,num_init = 50)