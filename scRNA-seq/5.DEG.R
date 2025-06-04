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

ref=readRDS("/path/to/ref_RH.rds")

###differential expressed gene between control and Dox
###Using MAST method to identify differentially expressed genes
###Most of cell type include subsample process to exclude the bias due to cell number difference between condition but some cell types with relative few cells would skip the downsample step. And this is shown in the following code for each cell type.
###We also take care of replicate variable to have a very rigid DEG calling
###Combine all these files to get Table S2
RH=subset(ref, subset= cell %in% c("RHKI"))
rm(ref)
#Naive Pluripotent Cells
tmp=subset(RH, subset= ann_KNN40 %in% c("Naive_Pluripotent_Cells"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', max.cells.per.ident = 10000, random.seed = 1234, test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Naive_Pluripotent_Cells.csv")
#Epiblasts
tmp=subset(RH, subset= ann_KNN40 %in% c("Epiblasts"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', max.cells.per.ident = 1000, random.seed = 1234, test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Epiblasts.csv")
#Primitive_Streak
tmp=subset(RH, subset= ann_KNN40 %in% c("Primitive_Streak"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Primitive_Streak.csv")
#Caudal_Neuroectoderm1
tmp=subset(RH, subset= ann_KNN40 %in% c("Caudal_Neuroectoderm1"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', max.cells.per.ident = 3000, random.seed = 1234, test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Caudal_Neuroectoderm1.csv")
#Neuromesodermal_Progenitors1
tmp=subset(RH, subset= ann_KNN40 %in% c("Neuromesodermal_Progenitors1"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', max.cells.per.ident = 7500, random.seed = 1234, test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Neuromesodermal_Progenitors1.csv")
#Neuromesodermal_Progenitors2
tmp=subset(RH, subset= ann_KNN40 %in% c("Neuromesodermal_Progenitors2"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', max.cells.per.ident = 5500, random.seed = 1234, test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Neuromesodermal_Progenitors2.csv")
#Neuromesodermal_Progenitors3
tmp=subset(RH, subset= ann_KNN40 %in% c("Neuromesodermal_Progenitors3"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Neuromesodermal_Progenitors3.csv")
#Spinal_Cord
tmp=subset(RH, subset= ann_KNN40 %in% c("Spinal_Cord"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', max.cells.per.ident = 2500, random.seed = 1234, test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Spinal_Cord.csv")
#Neuron-like_Cells
tmp=subset(RH, subset= ann_KNN40 %in% c("Neuron-like_Cells"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Neuron-like_Cells.csv")
#Nascent_Mesoderm
tmp=subset(RH, subset= ann_KNN40 %in% c("Nascent_Mesoderm"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', max.cells.per.ident = 300, random.seed = 1234, test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Nascent_Mesoderm.csv")
#Paraxial_Mesoderm_A
tmp=subset(RH, subset= ann_KNN40 %in% c("Paraxial_Mesoderm_A"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', max.cells.per.ident = 6000, random.seed = 1234, test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Paraxial_Mesoderm_A.csv")
#Paraxial_Mesoderm_B
tmp=subset(RH, subset= ann_KNN40 %in% c("Paraxial_Mesoderm_B"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', max.cells.per.ident = 1000, random.seed = 1234, test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Paraxial_Mesoderm_B.csv")
#Splanchnic_Mesoderm
tmp=subset(RH, subset= ann_KNN40 %in% c("Splanchnic_Mesoderm"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', max.cells.per.ident = 2000, random.seed = 1234, test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Splanchnic_Mesoderm.csv")
#Endothelium
tmp=subset(RH, subset= ann_KNN40 %in% c("Endothelium"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Endothelium.csv")
#Definitive_Endoderm_Guts
tmp=subset(RH, subset= ann_KNN40 %in% c("Definitive_Endoderm_Guts"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Definitive_Endoderm_Guts.csv")
#Placodal_area
tmp=subset(RH, subset= ann_KNN40 %in% c("Placodal_area"))
DEG=FindMarkers(tmp, ident.1 = "control", ident.2 = 'Dox', min.pct = 0.25, logfc.threshold = 0.25, group.by = 'treatment', max.cells.per.ident = 1000, random.seed = 1234, test.use = "MAST", latent.vars="replicate")
DEG=subset(DEG, subset= p_val_adj <= 0.05)
write.csv(DEG, "/path/to/DEG/Placodal_area.csv")

###For visualization
###For Fig 5E, Sup Fig 5E, the |avg_log2FC| >= 0.4 is used as criteria to peak genes from DEG list above
### features=c(genes would like to show), ident=c(cell type would like to show)
Stacked_VlnPlot(RH, features=c(), group.by = "ann_F", split.by = "treatment", pt.size = 0, x_lab_rotate = 45, colors_use = c("black", "red"), idents = c())+NoLegend()

###For specific set of genes acrossing multiple cell type in Fig 5F, 5G, Sup Fig. 5F, 5G(use *** to denote statistic significance, p_val_adj)
##critical genes for early decision points
Stacked_VlnPlot(RH, c("Nanog","Pou3f1", "Sox2","T","Nkx1-2", "Cdx2", "Hes3", "Hes5", "Pax6"), group.by = "ann_F", split.by = "treatment", pt.size = 0, x_lab_rotate = 45, colors_use = c("black", "red"), idents = c('Epiblasts', 'Primitive_Streak','Caudal_Neuroectoderm1', 'Neuromesodermal_Progenitors1', 'Neuromesodermal_Progenitors2', 'Spinal_Cord', 'Neuron-like_Cells'))+NoLegend()
##RA
Stacked_VlnPlot(RH,c("Crabp2", "Rara", "Rarb","Rxrb","Rxrg","Aldh1a2", "Cyp26a1"), group.by = "ann_F", split.by = "treatment", pt.size = 0, colors_use = c("black", "red"),idents = c('Neuromesodermal_Progenitors1', 'Neuromesodermal_Progenitors2', 'Spinal_Cord', 'Paraxial_Mesoderm_B', 'Paraxial_Mesoderm_A', 'Splanchnic_Mesoderm'))+NoLegend()
##other cell signaling gene
Stacked_VlnPlot(RH, c("Fgfr1", "Fgf8", "Tgfbr3", "Dlk1","Notch1", "Notch2", "Notch3"), group.by = "ann_new_noZX_F", split.by = "treatment", pt.size = 0, x_lab_rotate = 45, colors_use = c("black", "red"), idents = c(c('Primitive_Streak', 'Nascent_Mesoderm', 'Paraxial_Mesoderm_B', 'Paraxial_Mesoderm_A', 'Splanchnic_Mesoderm')))+NoLegend()
Stacked_VlnPlot(RH, c("Tcf7l1", "Lef1"), group.by = "ann_F", split.by = "treatment", pt.size = 0, x_lab_rotate = 45, colors_use = c("black", "red"))+NoLegend()
