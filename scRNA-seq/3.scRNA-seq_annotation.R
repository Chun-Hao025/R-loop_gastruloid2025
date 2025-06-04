library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)
library(scCustomize)
library(ggplot2)
library(cellmatch)
library(biomaRt)
library(GO.db)
library(ggbreak)
library(Matrix)
library(reshape2)

set.seed(1234)
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

###combine ourdata with publised dataset
ref=readRDS("/path/to/cellranger/publisheddata/post-filter/ref_merge.rds")
RH=readRDS("/path/to/post-filter/RHKI_merge.rds")

tmp= merge(ref, y = RH, merge.data=T)
rm(RH, ref)

mat1 = tmp@assays$RNA@counts
tmp_adf=tmp@meta.data
mat1 = mat1[!(rownames(mat1)%in%exclude.gene),]
ref=CreateSeuratObject(mat1)
ref@meta.data=tmp_adf
ref=subset(ref, subset= doublet %in% c("Singlet"))
ref=NormalizeData(ref) %>% FindVariableFeatures(nfeatures = 2000)

###regression out CC.Difference as the trade-off between the heterogeniety of proliferation potency in different timepoints of gastruloid differentiation and skewing by cell proliferation 
ref = ScaleData(ref, vars.to.regress =c("CC.Difference")) %>% RunPCA()
ElbowPlot(ref, ndims = 50)

###batch effect correction through harmony
ref=RunHarmony(ref, "replicate", reduction.save = "harmony_rep", max_iter=30)
ref = RunUMAP(ref, dims = 1:30, reduction = "harmony_rep") %>% FindNeighbors(dims = 1:30, reduction = "harmony_rep") %>% FindClusters()

###export harmony latent space for cell clustering through pyVIA in python (check 4.trajectory_inference.py)
harm=ref@reductions$harmony_rep@cell.embeddings
umap=ref@reductions$umap@cell.embeddings
meta_adf=ref@meta.data
mat=ref@assays$RNA@counts
write.csv(meta_adf, "/path/to/trajectory/via_meta_new.csv")
write.csv(harm, "/path/to/trajectory/via_harmony_new.csv")
write.csv(umap, "/path/to/trajectory/via_umap_new.csv")
writeMM(mat, file="/path/to/trajectory/ref_RH_count.mtx")

###PARC_KNN40 is the cluster from pyVIA 
meta_adf=read.csv("/path/to/trajectory/via_meta_new1.csv", row.names = 1)
ref@meta.data=meta_adf
ref$PARC_KNN40=paste0("PARC_", ref$PARC_KNN40)

###markers of clusters
Idents(ref)=ref$PARC_KNN40
marker=FindAllMarkers(ref, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
marker%>%group_by(cluster)%>%slice_max(n=2, order_by=avg_log2FC)

ref_bulk = cellmatch::AggregateByCluster(expm1(Seurat::GetAssayData(ref, "RNA")), ref$PARC_KNN40, method = "average")

###embryo data and annotation from https://tome.gs.washington.edu/
###combination of multiple dataset (GEO: GSE100597, GEO: GSE109071, Atlas-E-MTAB-6967, GEO: GSE186069, GEO: GSE186068)
embryo=readRDS("/path/to/embryo/embryo_all.rds")
DefaultAssay(embryo)="RNA"
##exclude extra embryonic and cell type related to some terminal differentiation status such as blood cell ...
embryo = subset(embryo, subset= cell_type %in% c("Allantois", "Amniochorionic mesoderm", "Amniochorionic mesoderm A", "Amniochorionic mesoderm B", "Anterior floor plate", "Blood progenitors", "Embryonic visceral endoderm", "Extraembryonic ectoderm",
                                                 "Extraembryonic mesoderm", "Extraembryonic visceral endoderm", "First heart field", "Hypoblast","Parietal endoderm", "Posterior floor plate", "Primitive erythroid cells","Second heart field", "Visceral endoderm"), invert=T)
embryo = NormalizeData(embryo)
embryo_bulk = cellmatch::AggregateByCluster(expm1(Seurat::GetAssayData(embryo, "data")), embryo$cell_type, method = "average")

setwd("/path/to/")
temp_dir1 = getwd() %>% file.path("PARC_KNN40_2000")
dir.create(temp_dir1)
matching_output1 = RunCellMatch(query = ref_bulk,reference = embryo_bulk, results_path = temp_dir1,K = 2000,num_init = 50)

###project to the processed data (with annotation) from GSE212050
biox=readRDS("path/to/publisheddata/GSE212050_seurat_final")

Anchor = FindTransferAnchors(reference = biox, query = ref, dim=1:30, reference.reduction = "pca")

ref = MapQuery(anchorset = Anchor, reference = biox, query = ref, reference.reduction = "pca", refdata=biox$cell_type, reduction.model="umap")
ref$biox_ann = ref$predicted.id
ref$biox_ann_score=ref$predicted.id.score

##based on the similarity between cluster in the merge dataset and embryo dataset, assign the annotation.
##PARC_4 is a awkard cluster which is with different cell number distribution across all timepoints comparing to PARC_11 and PARC_23 but these three are assigned to caudal neuroectoderm by embryo data
##PARC_4 is with the similar score to Neuromesodermal progenitor by embryo data.
##PARC_4 is majorly annotated to late spinal cord by the GSE212050
##PARC_4 shares feature among Caudal neuroectoderm, NMPs and spinal cord but the top marker genes such as Crapb2, Hoxc8, Hes3 and Irx3 are more similiar to NMPs and spinal cord.
##Therefore, PARC_4 are here annotated to NMP2 instead of the annotation from cellmatch result.
ann=rep(NA, nrow(meta_adf))##KNN40
for (i in 1:nrow(meta_adf)){
  if(meta_adf$PARC_noZX1[i] %in% c("0")){
    ann[i] = "Naive_Pluripotent_Cells"
  }else if(meta_adf$PARC_KNN40[i] %in% c("2")){
    ann[i] = "Primitive_Streak"
  }else if(meta_adf$PARC_KNN40[i] %in% c("7","10","16")){
    ann[i] = "Epiblasts"
  }else if(meta_adf$PARC_KNN40[i] %in% c("18")){
    ann[i] = "Placodal_area"
  }else if(meta_adf$PARC_KNN40[i] %in% c("11", "23")){
    ann[i] = "Caudal_Neuroectoderm1"
  }else if(meta_adf$PARC_KNN40[i] %in% c("4")){
    ann[i] = "Neuromesodermal_Progenitors2"
  }else if(meta_adf$PARC_KNN40[i] %in% c("1","13","20")){
    ann[i] = "Neuromesodermal_Progenitors1"
  }else if(meta_adf$PARC_KNN40[i] %in% c("5","15")){
    ann[i] = "Neuromesodermal_Progenitors2"
  }else if(meta_adf$PARC_KNN40[i] %in% c("24")){
    ann[i] = "Neuromesodermal_Progenitors3"
  }else if(meta_adf$PARC_KNN40[i] %in% c("3")){
    ann[i] = "Spinal_Cord"
  }else if(meta_adf$PARC_KNN40[i] %in% c("26")){
    ann[i] = "Neuron-like_Cells"
  }else if(meta_adf$PARC_KNN40[i] %in% c("22")){
    ann[i] = "Definitive_Endoderm_Guts"
  }else if(meta_adf$PARC_KNN40[i] %in% c("19")){
    ann[i] = "Nascent_Mesoderm"
  }else if(meta_adf$PARC_KNN40[i] %in% c("12","14","17")){
    ann[i] = "Paraxial_Mesoderm_B"
  }else if(meta_adf$PARC_KNN40[i] %in% c("6","8")){
    ann[i] = "Paraxial_Mesoderm_A"
  }else if(meta_adf$PARC_KNN40[i] %in% c("21", "9")){
    ann[i] = "Splanchnic_Mesoderm"
  }else if(meta_adf$PARC_KNN40[i] %in% c("25")){
    ann[i] = "Endothelium"
  }else {
    ann[i] = "unknown"
  }
}
ref$ann_KNN40=ann

##rerun the similarity analysis for Sup Fig 5C
temp_dir2 = getwd() %>% file.path("ann")
dir.create(temp_dir2)
matching_output1 = RunCellMatch(query = ref_bulk1,reference = embryo_bulk, results_path = temp_dir2,K = 2000,num_init = 50)

###
ref$ann_F=factor(ref$ann_KNN40, levels = c("Naive_Pluripotent_Cells", "Epiblasts", "Primitive_Streak", "Caudal_Neuroectoderm1", "Neuromesodermal_Progenitors1", "Neuromesodermal_Progenitors2", "Neuromesodermal_Progenitors3",
                                                                         "Spinal_Cord", "Neuron-like_Cells","Nascent_Mesoderm","Paraxial_Mesoderm_A","Paraxial_Mesoderm_B","Splanchnic_Mesoderm","Endothelium","Definitive_Endoderm_Guts","Placodal_area"))

ref_bulk1 = cellmatch::AggregateByCluster(expm1(Seurat::GetAssayData(ref, "RNA")), ref$ann_new_noZX_KNN40, method = "average")

###lineage assign
cat=rep(NA, nrow(meta_adf))
for (i in 1:nrow(meta_adf)){
  if (meta_adf$ann__KNN40[i] %in% c("Nascent_Mesoderm", "Paraxial_Mesoderm_A", "Paraxial_Mesoderm_B", "Splanchnic_Mesoderm", "Endothelium")){
    cat[i] = "Mesoderm"
  }else if (meta_adf$ann_KNN40[i] %in% c("Caudal_Neuroectoderm1", "Caudal_Neuroectoderm2", "Neuromesodermal_Progenitors1", "Neuromesodermal_Progenitors2", "Neuromesodermal_Progenitors3", "Spinal_Cord", "Neuron-like_Cells")){
    cat[i] = "Ectoderm"
  }else if (meta_adf$ann_KNN40[i] %in% c("Definitive_Endoderm_Guts")){
    cat[i] = "Endoderm"
  }else{
    cat[i] = meta_adf$ann_KNN40[i]
  }
}
ref$ann_cat=cat
###fix color code
color=rep(NA, nrow(meta_adf))
for (i in 1:nrow(meta_adf)){
  if(meta_adf$ann_KNN40[i] %in% c("Naive_Pluripotent_Cells")){
    color[i] = "#5A5156FF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Epiblast")){
    color[i] = "#E4E1E3FF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Primitive_Streak")){
    color[i] = "#F6222EFF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Caudal_Neuroectoderm1")){
    color[i] = "#FE00FAFF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Neuromesodermal_Progenitors1")){
    color[i] = "#16FF32FF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Neuromesodermal_Progenitors2")){
    color[i] = "#3283FEFF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Neuromesodermal_Progenitors3")){
    color[i] = "#FEAF16FF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Spinal_Cord")){
    color[i] = "#B00068FF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Neuron-like_Cells")){
    color[i] = "#1CFFCEFF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Nascent_Mesoderm")){
    color[i] = "#90AD1CFF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Paraxial_Mesoderm_A")){
    color[i] = "#2ED9FFFF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Paraxial_Mesoderm_B")){
    color[i] = "#DEA0FDFF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Splanchnic_Mesoderm")){
    color[i] = "#AA0DFEFF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Endothelium")){
    color[i] = "#F8A19FFF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Definitive_Endoderm_Guts")){
    color[i] = "#325A9BFF"
  }else if(meta_adf$ann_KNN40[i] %in% c("Placodal_area")){
    color[i] = "#C4451CFF"
  }else {
    color[i] = "unknown"
  }
}
ref$ann_color=color
###you can download the final seurat object from GEO page to explore it (ref_RH_final.rds)
###generate DimPlot by time, data resource and cluster
ref$time_F=factor(ref$time, levels = c("0h", "24h", "36h", "48h", "52h", "56h", "60h", "72h", "84h", "90h", "96h", "102h", "108h", "114h", "120h"))
DimPlot_scCustom(ref, label=T, group.by="ann_F", raster = F,repel = T,label.size = 10, aspect_ratio = 1)+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.title.x = element_text(size=20, color='black'),axis.title.y = element_text(size=20, color='black'),axis.text = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),title =element_text(size=0))
DimPlot_scCustom(ref, label=F, group.by="time_F", raster = F,repel = T,label.size = 10, aspect_ratio = 1)+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.title.x = element_text(size=20, color='black'),axis.title.y = element_text(size=20, color='black'),axis.text = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),title =element_text(size=0))
DimPlot_scCustom(ref, label=F, group.by="resource", raster = F,repel = T,label.size = 10, aspect_ratio = 1)+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.title.x = element_text(size=20, color='black'),axis.title.y = element_text(size=20, color='black'),axis.text = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),title =element_text(size=0))

###cell proportion barplot
meta_adf=subset(ref@meta.data, subset = cell %in% c("RHKI"))
meta_adf$ID=paste0(meta_adf$time, "_", meta_adf$treatment)
time1=table(meta_adf$ID, meta_adf$ann_new_noZX_cat)
time2=table(meta_adf$ID, meta_adf$ann_new_noZX_F)
time1=data.frame(time1)
time2=data.frame(time2)
time1=dcast(time1, Var1 ~ Var2, value.var = 'Freq')
time2=dcast(time2, Var1 ~ Var2, value.var = 'Freq')
rownames(time1)=time1$Var1
rownames(time2)=time2$Var1
time1=time1[,-1]
time2=time2[,-1]
time_perc1 <- apply(time1, 2,function(x) x/rowSums(time1))
time_perc2 <- apply(time2, 2,function(x) x/rowSums(time2))
time_perc1=as.data.frame(time_perc1)
time_perc2=as.data.frame(time_perc2)
time_perc1$ID=rownames(time_perc1)
time_perc2$ID=rownames(time_perc2)
time_perc1$time=c("0hr", "0hr", "120hr", "120hr", "72hr", "72hr")
time_perc1$time=factor(time_perc1$time, levels=c("0hr", "72hr", "120hr"))
time_perc1$treatment=c("control", "Dox","control", "Dox","control", "Dox")
time_perc2$time=c("0hr", "0hr", "120hr", "120hr", "72hr", "72hr")
time_perc2$time=factor(time_perc2$time, levels=c("0hr", "72hr", "120hr"))
time_perc2$treatment=c("control", "Dox","control", "Dox","control", "Dox")

time_T1 <- time_perc1 %>%
  melt(id.var=c("ID", "time", "treatment")) %>%
  arrange(ID, variable)
time_T2 <- time_perc2 %>%
  melt(id.var=c("ID", "time", "treatment")) %>%
  arrange(ID, variable)
time_T1$variable=factor(time_T1$variable, levels = c("Naive_Pluripotent_Cells", "Epiblasts", "Primitive_Streak", "Ectoderm", "Mesoderm", "Endoderm", "Placodal_area"))

##lineage check following time
ggplot(time_T1) +
  aes(x = treatment, fill=variable, y=value)+scale_fill_manual(values=c(q1$colour[1:3], "darkorange", "deepskyblue1", "deeppink", q1$colour[16]))+
  scale_y_continuous(labels = scales::percent, expand = c(0.005, 0.005))+scale_x_discrete(expand = c(0.25, 0.25))+
  geom_bar(position = "fill", stat="identity")+facet_wrap(~time)+coord_fixed(ratio = 2)+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.title.x = element_text(size=20, color='black'),axis.title.y = element_text(size=20, color='black'),axis.text = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),title =element_text(size=0))

ggplot(time_T2) +
  aes(x = "", fill=treatment, y=value)+
  scale_y_continuous(labels = scales::percent, expand = c(0.005, 0.005))+
  geom_bar(data=subset(sleep.long3, subset = time %in% c("120hr")),position = "fill", stat="identity", width=1)+coord_polar("y", start=0)+
  facet_wrap(~celltype)+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=0, color='black'), axis.ticks = element_line(size = 0), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"))

ggplot(data=subset(time_T2, subset = time %in% c("72hr","120hr") &variable %in% c("Neuromesodermal_Progenitors1"))) +
  aes(x = treatment, y=value, fill=treatment)+scale_y_continuous(labels = scales::percent, expand = c(0.00001, 0.00001))+scale_fill_manual(values=c('black', 'red'))+scale_y_break(breaks = c(0.01,0.15), scales = 0.5)+
  geom_bar(stat="identity", width=0.8)+facet_grid(~time)+ggtitle("Neuromesodermal_Progenitors1")+NoLegend()+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white"), panel.grid = element_line(size=0), plot.title = element_text(size=20, color='black',hjust = 0.5),axis.line = element_line(size=0), axis.title.x = element_text(size=0, color='black'),axis.title.y = element_text(size=0, color='black'),axis.text = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),title =element_text(size=0))

ggplot(data=subset(time_T2, subset = time %in% c("72hr","120hr") &variable %in% c("Neuromesodermal_Progenitors2"))) +
  aes(x = treatment, y=value, fill=treatment)+scale_y_continuous(labels = scales::percent, expand = c(0.00001, 0.00005))+scale_fill_manual(values=c('black', 'red'))+scale_y_break(breaks = c(0.02,0.1), scales = 0.5)+
  geom_bar(stat="identity", width=0.8)+facet_grid(~time)+ggtitle("Neuromesodermal_Progenitors2")+NoLegend()+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white"), panel.grid = element_line(size=0), plot.title = element_text(size=20, color='black',hjust = 0.5),axis.line = element_line(size=0), axis.title.x = element_text(size=0, color='black'),axis.title.y = element_text(size=0, color='black'),axis.text = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),title =element_text(size=0))

ggplot(data=subset(time_T2, subset = time %in% c("72hr", "120hr") &variable %in% c("Spinal_Cord"))) +
  aes(x = treatment, y=value, fill=treatment)+scale_y_continuous(labels = scales::percent, expand = c(0.00001, 0.00005))+scale_fill_manual(values=c('black', 'red'))+scale_y_break(breaks = c(0.02,0.1), scales = 0.5)+
  geom_bar(stat="identity", width=0.8)+facet_grid(~time)+ggtitle("Spinal_Cord")+NoLegend()+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white"), panel.grid = element_line(size=0), plot.title = element_text(size=20, color='black',hjust = 0.5),axis.line = element_line(size=0), axis.title.x = element_text(size=0, color='black'),axis.title.y = element_text(size=0, color='black'),axis.text = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),title =element_text(size=0))

###find marker for annotated celltype
Idents(ref)=ref$ann_KNN40
marker1=FindAllMarkers(ref, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
marker1%>%group_by(cluster)%>%slice_max(n=2, order_by=avg_log2FC)

##marker gene violin plot for Sup. Fig. 5D
genes=c("Zfp42","Nanog", "Dppa5a", "Pou3f1", "Pim2","T", "Fst","Epcam","Cdx2", "Nkx1-2", "Ptn", "Crabp2", "Hoxaas3", "Hoxc8","Irx3", "Hes5","Tubb3", "Mixl1", "Tbx6","Rspo3" ,"Cdh11","Aldh1a2","Pcdh19","Meox1","Tcf15","Pax3","Tgfb2", "Twist1", "Kdr","Etv2", "Trh", "Krt8")
Stacked_VlnPlot(ref, genes, pt.size=0, group.by = "ann_F")+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white"), panel.grid = element_line(size=0),axis.line = element_line(size=0), axis.title.x = element_text(size=0, color='black'),axis.title.y = element_text(size=0, color='black'),axis.text.x = element_text(size=0, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),title =element_text(size=0))

saveRDS(ref, "/path/to/ref_RH.rds")
