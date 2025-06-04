library(data.table)
library(dplyr)
library(patchwork)


total=read.csv("/path/to/GRN/total_rss.csv", row.names = 1)
colnames(total) <- gsub("...", "", colnames(total), fixed = TRUE)
total$cell=rownames(total)
total2=melt(total, id.var=c("cell"))
total2_sorted <- total2 %>%
  arrange(value) %>%    # sorts from low to high
  mutate(rank = row_number())
ggplot(total2, aes(x=value))+
  stat_ecdf(geom = "step") +
  facet_wrap(~ cell)

ggplot(total2_sorted, aes(x = rank, y = value)) +
  geom_point() +
  facet_wrap(~ cell, scales = "free") +
  labs(x = "Rank (low to high)", y = "Value") +
  theme_minimal()

cell_T=unique(total2$cell)
celltype=c('Caudal_Neuroectoderm1', 'Definitive_Endoderm_Guts', 'Endothelium', 'Epiblasts','Naive_Pluripotent_Cells', 'Nascent_Mesoderm', 'Neuromesodermal_Progenitors1',
           'Neuromesodermal_Progenitors2', 'Neuromesodermal_Progenitors3', 'Neuron-like_Cells', 'Paraxial_Mesoderm_A', 'Paraxial_Mesoderm_B', 'Placodal_area', 'Primitive_Streak', 'Spinal_Cord', 'Splanchnic_Mesoderm')
cluster=c("cluster0", "cluster1", "cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10",
         "cluster11","cluster12","cluster13","cluster14","cluster15")
##Define cell type specific TF by rss score above 60 quantile for each cell type and add these TF-target to the cell type specific TF-target pair.

TFs=list()
for (i in 1:length(cell_T)){
  tmp=total2[total2$cell %in% cell_T[i],]
  thre=quantile(tmp$value,probs = c(0.6))
  tmp=tmp[tmp$value >= thre,]
  TFs[[i]]=as.character(tmp$variable)
}
names(TFs)=cluster
for (i in 1:length(cluster)){
  tmp=read.csv(paste0("/path/to/GRN/", cluster[i], "_reg_TF_target.csv"), row.names = 1)
  total1=read.csv("/path/to/GRN/total_reg_TF_target.csv", row.names = 1)
  tmp <- tmp %>%
    mutate(weight_percentile = (rank(weight, ties.method = "average") / n()))
  total1=subset(total1, subset = TF %in% TFs[[i]])
  total1 <- total1 %>%
    mutate(weight_percentile = (rank(weight, ties.method = "average") / n()))
  tmp1=rbind(tmp, total1)
  tmp1 <- tmp1 %>%
    group_by(TF, target) %>%
    slice_max(order_by = weight_percentile, n = 1, with_ties = FALSE) %>%
    ungroup()
  tmp1=tmp1[,c(1,2,4)]
  tmp1$TF=paste0(tmp1$TF, "_", cluster[i])
  tmp1$target=paste0(tmp1$target, "_", cluster[i])
  write.table(tmp1, paste0("/path/to/GRN", cluster[i], "_network.txt"), sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)
}

