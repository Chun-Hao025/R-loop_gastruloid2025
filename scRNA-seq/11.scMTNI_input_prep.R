library(Seurat)
library(dplyr)
library(patchwork)

ref=readRDS("/path/to/ref_RH.rds")
cluster=c("cluster0", "cluster1", "cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10",
          "cluster11","cluster12","cluster13","cluster14","cluster15")
celltype=c('Caudal_Neuroectoderm1','Definitive_Endoderm_Guts','Endothelium','Epiblasts','Naive_Pluripotent_Cells','Nascent_Mesoderm',
           'Neuromesodermal_Progenitors1','Neuromesodermal_Progenitors2','Neuromesodermal_Progenitors3','Neuron-like_Cells','Paraxial_Mesoderm_A','Paraxial_Mesoderm_B',
           'Placodal_area','Primitive_Streak','Spinal_Cord','Splanchnic_Mesoderm')
ann=data.frame(cluster, celltype)
tmp=subset(ref, subset=treatment %in% c("control"))
rm(ref)
gc()

for (i in 1:nrow(ann)){
  tmp1=subset(tmp, subset= treatment %in% c("control") & ann_new_noZX_KNN40 %in% celltype[i])
  tmp2=subset(tmp, subset= cell %in% c("RHKI") & treatment %in% c("control") & ann_new_noZX_KNN40 %in% celltype[i])
  tmp3=subset(tmp, subset= cell %in% c("RHKI") & treatment %in% c("Dox") & ann_new_noZX_KNN40 %in% celltype[i])
  data1=tmp1@assays$RNA@data
  data1=as.data.frame(data1)
  data2=tmp2@assays$RNA@data
  data2=as.data.frame(data2)
  data3=tmp3@assays$RNA@data
  data3=as.data.frame(data3)
  write.table(data1, paste0("./GRN/scMTNI/control/", ann$cluster[i], "_count.tsv"), sep = "\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
  rm(tmp1, data)
  gc()
  write.table(data2, paste0("./GRN/scMTNI/RH_control/", ann$cluster[i], "_count.tsv"), sep = "\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
  rm(tmp1, data)
  gc()
  write.table(data3, paste0("./GRN/scMTNI/RH_Dox/", ann$cluster[i], "_count.tsv"), sep = "\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
  rm(tmp1, data)
  gc()
}

