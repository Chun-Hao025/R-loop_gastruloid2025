library(tximport)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ComplexHeatmap)
library(reshape2)
library(dplyr)
library(tidyverse)
library(GenomicFeatures)
library(DEGreport)
library(viridis)
library(patchwork)

###geneInfo.tab was from STAR output folder and contained ensemble ID to gene symbol information
Gene=read.table("/path/to/geneInfo.tab", skip = 1)
colnames(Gene)=c("ensembl_gene_id_version", 'external_gene_name', 'type')
dup=duplicated(Gene$V2)

for (i in 1:length(dup)) {
  if (dup[i]){
    Gene$V2[i]=paste0(Gene$V2[i], ".1")
  }
}
dup=duplicated(Gene$V2)

for (i in 1:length(dup)) {
  if (dup[i]){
    Gene$V2[i]=gsub(".1", "", Gene$V2[i])
    Gene$V2[i]=paste0(Gene$V2[i], ".2")
  }
}
##remove version index in gene_id
Gene$ensembl_gene_id=gsub("\\..*", "", Gene$ensembl_gene_id_version)

#### Day0 control gene expression level classification ####
txdb <- makeTxDbFromGFF(file = "/path/to/gencode.vM25.annotation.gtf", dataSource = "gencode", organism = "Mus musculus",)
g = genes(txdb)
g_chr=as.data.frame(g@seqnames)
g_chr=cbind(g_chr, as.data.frame(g@ranges))
colnames(g_chr)[1]="chr"
g_chr=cbind(g_chr, as.data.frame(g@strand))
colnames(g_chr)[6]="strand"
K = g_chr[match(Gene$ensembl_gene_id_version, g_chr[,5]),]
K = cbind(K, Gene)
K = K[,c(-5)]
#read all control 0hr RSEM output
control1=read.delim("/path/to/RSEM/output/control_0hr_rep1.RSEM.genes.results")
control2=read.delim("/path/to/RSEM/output/control_0hr_rep2.RSEM.genes.results")
control3=read.delim("/path/to/RSEM/output/control_0hr_rep3.RSEM.genes.results")
control1=control1[match(K$ensembl_gene_id_version,control1$gene_id),]
control2=control2[match(K$ensembl_gene_id_version,control2$gene_id),]
control3=control3[match(K$ensembl_gene_id_version,control3$gene_id),]
K$N_FPKM_rep1=control1$FPKM
K$N_FPKM_rep2=control2$FPKM
K$N_FPKM_rep3=control3$FPKM
K$mean_FPKM=rowMeans(K[,c("N_FPKM_rep1", "N_FPKM_rep2", "N_FPKM_rep3")])
#exclude non-expressed gene first since it's a big number and would dilute out the high expression gene class.
K1=subset(K, subset= mean_FPKM > 0)
#Get 75 percentile and 50 percentile value for classifying
quantile(K1, c(0.5,0.75))

for (i in 1:nrow(K)){
  if (K$mean_FPKM[i] <= 1.2){
    K$class[i]="class3"
  }else if (K$mean_FPKM[i] > 1.2 & K$mean_FPKM[i] <= 13.5){
    K$class[i]="class2"
  }else{
    K$class[i]="class1"
  }
}
K=K[order(K$mean_FPKM, decreasing = TRUE),]
#heatmap for average log2(FPKM+1) for each class
exp_level=c("class1", "class2", "class3")
exp_matrix=data.frame(c(0,0,0),row.names = exp_level)
colnames(exp_matrix)=c("MeanFPKM")
for (i in 1:length(exp_level)){
  tmp=K[K$class %in% exp_level[i],]
  exp_matrix$MeanFPKM[i]=mean(tmp$mean_FPKM)
  exp_matrix$MeanLog[i]=log2(exp_matrix$MeanFPKM[i]+1)
}
exp_matrix <- as.matrix(exp_matrix)
exp_matrix <- apply(exp_matrix, 2, as.numeric)
Heatmap(exp_matrix[,2], 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        col = viridis(n = 10, alpha = 0.5, begin = 0, end = 1, option = "D"), 
        border = TRUE, 
        row_labels = c('High', 'Medium', 'Low'), 
        row_names_side = "left", 
        border_gp = gpar(lwd = 2, col = "black"))
class1=subset(K, subset=class %in% c("class1"))
class2=subset(K, subset=class %in% c("class2"))
class3=subset(K, subset=class %in% c("class3"))
write.table(class1, "/path/to/class1.bed", row.names=F, col.names=F, quote=F, sep="\t")
write.table(class2, "/path/to/class2.bed", row.names=F, col.names=F, quote=F, sep="\t")
write.table(class3, "/path/to/class3.bed", row.names=F, col.names=F, quote=F, sep="\t")

###RSEM output contains raw gene count
dirs = list.files("/path/to/RSEM/output/", pattern = "RSEM.genes.results")
files = gsub("^", "/path/to/RSEM/output/", dirs)
samples = gsub(".RSEM.genes.results", "", dirs)
names(files) = samples
txi.rsem=tximport(files, type = "rsem", txIn=FALSE, txOut=FALSE)
Gene=Gene[match(rownames(txi.rsem$counts),Gene$V1),]
rownames(txi.rsem$counts)=Gene$V2
rownames(txi.rsem$abundance)=Gene$V2
rownames(txi.rsem$length)=Gene$V2
colnames(txi.rsem$counts)=samples
colnames(txi.rsem$abundance)=samples
colnames(txi.rsem$length)=samples
zero_length_and_unexpressed = (apply(txi.rsem$counts, 1, max) == 0) &
  (apply(txi.rsem$length, 1, min) == 0)

txi.rsem$length = txi.rsem$length[!zero_length_and_unexpressed,]
txi.rsem$abundance = txi.rsem$abundance[!zero_length_and_unexpressed,]
txi.rsem$counts = txi.rsem$counts[!zero_length_and_unexpressed,]
counts=txi.rsem$counts
counts = round(counts)

#colData
dss_time = c(rep("Day0", 6), rep("Day3", 6), rep("Day5", 6))
dss_rep = c(rep(c("rep1", "rep2", "rep3"),6))
dss_cell= c(rep("RHKI", 18))
dss_treatment=c(rep(c(rep("Dox", 3), rep("control", 3)),3))
dss_time=factor(dss_time, levels = c("Day0", "Day3", "Day5"))
dss_treatment=factor(dss_treatment, levels = c("control", "Dox"))
dss_rep=factor(dss_rep, levels = c("rep1", "rep2", "rep3"))
colData=DataFrame(dss_cell, dss_time, dss_treatment, dss_rep, row.names = colnames(counts))

#PCA analysis
dds = DESeqDataSetFromMatrix(counts, colData = colData, design = ~dss_treatment+dss_time+dss_treatment:dss_time )
dds <- DESeq(dds)
vsd = vst(dds, blind=FALSE)
p1=plotPCA(vsd, intgroup=c("dss_treatment", "dss_time"))
pca=p1$data
ggplot(pca, aes(x=PC1, y=PC2, color=dss_treatment, label=name))+geom_point(size=3)+scale_color_manual(values=c('black', 'red', 'blue', 'grey', 'orange', 'purple'))+
  geom_label_repel(aes(label = name), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', force = 30)+
  xlab(paste0("PC1(77%)")) +ylab(paste0("PC2(22%)")) +coord_fixed(ratio = 2)+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"))

##split data by timepoints
#Day0 (0hr)
counts_Day0=counts[,c(1:6)]
colData_Day0=colData[c(1:6),]
dds_Day0 = DESeqDataSetFromMatrix(counts_Day0, colData = colData_Day0, design = ~dss_treatment)
dds_Day0 = DESeq(dds_Day0)

Day0=as.data.frame(results(dds_Day0))
Day0=Day0[,c(1,2,6)]
colnames(Day0)=c("RHKI_baseMean","RHKI_LFC", "RHKI_padj")
Day0[,3]=ifelse(is.na(Day0[,3]), 1, Day0[,3])
Day0[,2]=ifelse(is.na(Day0[,2]), 0, Day0[,2])
tmp=rep(NA, nrow(Day0))
tmp1=rep(NA, nrow(Day0))
adj=Day0[,3]
LFC=Day0[,2]

#significant by adjusted p-value <= 0.05
for(i in 1:length(adj)){
  if(adj[i] <= 0.05){
    tmp[i] = "sig"
  }else{
    tmp[i] = "none"
  }
}
#Up and downregulation with 0.58 cutoff (1.5 fold change)
for(i in 1:length(LFC)){
  if(adj[i] <= 0.05 & LFC[i] >= 0.58){
    tmp1[i] = "UP"
  }else if (adj[i] <= 0.05 & LFC[i] <= -0.58){
    tmp1[i] = "Down"
  }else{
    tmp1[i] = "none"
  }
}
Day0=cbind(Day0, tmp, tmp1)
colnames(Day0)[4:5]=c("sig", "DEG")

Day0$delabel=rep(NA, nrow(Day0))
Day0$symbol=rownames(Day0)
Day0_highlight=c("Rnaseh1", "Cdx1", "Hhip", "S100a1")
for (i in 1:nrow(Day0)){
  if (rownames(Day0)[i] %in% Day0_highlight) {
    Day0$delabel[i] = Day0$symbol[i]
  }
}

ggplot(Day0, aes(x=RHKI_LFC, y=-log10(RHKI_padj)))+xlim(-7.5,7.5)+ylim(0,100)+
  geom_point(data=subset(Day0, subset= sig %in% c("none")), size=4, alpha=0.3, col='grey')+
  geom_point(data=subset(Day0, subset= sig %in% c("sig")), col="wheat", size=4)+
  geom_point(data=subset(Day0, subset= DEG %in% c("UP", "Down")), aes(col=DEG), size=4)+scale_color_manual(values=c('darkorange','deepskyblue'))+
  ###geom_point(data=subset(Day0_RHZ, subset= DEG_RHKI %in% c("UP", "Down")), aes(col=DEG_RHKI))+###
  geom_vline(xintercept = 0, color="black")+
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed') +geom_vline(xintercept = c(-0.58,0.58), linetype = 'dashed') +scale_y_break(c(15, 70), scales = 0.2)+
  geom_text_repel(aes(label=delabel), max.overlaps = 50, color="black",box.padding = 1, size=10)+
  xlab("log2(Dox/control)")+ylab("-log10(adjusted p-value)")+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"), axis.text = element_text(size = 12,colour = "black"))
genes_Day0=c("Pou5f1", "Sox2", "Rnaseh1", "Hhip", "Cdx1", "S100a1")
normal_count_Day0=t(assay(rlog(dds_Day0, blind = FALSE)))
normal_count_Day0=as.data.frame(normal_count_Day0)
normal_count_Day0$cell=colData_Day0$dss_cell
normal_count_Day0$treatment=colData_Day0$dss_treatment
normal_count_Day0$ID=rownames(normal_count_Day0)

plots = list()
for (i in 1:length(genes_Day0)){
  gene=genes_Day0[i]
  plots[[i]]=ggplot(normal_count_Day0, aes(x=treatment, y=.data[[gene]], group=treatment, fill=treatment)) + geom_violin(draw_quantiles = c(0.5), aes(group = interaction(cell, treatment), fill = treatment))+scale_fill_manual(values=c('grey', 'red', 'bisque2', 'deepskyblue', 'orange', 'purple','darkred'))+
    geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size = 4)+scale_color_manual(values=c('black', 'red', 'blue', 'grey', 'orange', 'purple', 'darkred'))+
    theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"), axis.text = element_text(size = 12,colour = "black"))
}
wrap_plots(plots, nrow=3)
#Day3 (72hr)
counts_Day3=counts[,c(7:12)]
colData_Day3=colData[c(7:12),]
dds_Day3 = DESeqDataSetFromMatrix(counts_Day3, colData = colData_Day3, design = ~dss_treatment) ###
dds_Day3 = DESeq(dds_Day3)

Day3=as.data.frame(results(dds_Day3))
Day3=Day3[,c(1,2,6)]
colnames(Day3)=c("RHKI_baseMean","RHKI_LFC", "RHKI_padj")
Day3[,3]=ifelse(is.na(Day3[,3]), 1, Day3[,3])
Day3[,2]=ifelse(is.na(Day3[,2]), 0, Day3[,2])
tmp=rep(NA, nrow(Day3))
tmp1=rep(NA, nrow(Day3))
adj=Day3[,3]
LFC=Day3[,2]

#significant by adjusted p-value <=0.05
for(i in 1:length(adj)){
  if(adj[i] <= 0.05){
    tmp[i] = "sig"
  }else{
    tmp[i] = "none"
  }
}
#up and downregulation with LFC cutoff 0.58 (1.5 fold change)
for(i in 1:length(LFC)){
  if(adj[i] <= 0.05 & LFC[i] >= 0.58){
    tmp1[i] = "UP"
  }else if (adj[i] <= 0.05 & LFC[i] <= -0.58){
    tmp1[i] = "Down"
  }else{
    tmp1[i] = "none"
  }
}
Day3=cbind(Day3, tmp, tmp1)
colnames(Day3)[4:5]=c("sig", "DEG")

Day3$delabel=rep(NA, nrow(Day3))
Day3$symbol=rownames(Day3)
Day3_highlight=c("Rnaseh1")
for (i in 1:nrow(Day3)){
  if (rownames(Day3)[i] %in% Day3_highlight) {
    Day3$delabel[i] = Day3$symbol[i]
  }
}

ggplot(Day3, aes(x=RHKI_LFC, y=-log10(RHKI_padj)))+xlim(-7.5,7.5)+ylim(0,140)+
  geom_point(data=subset(Day3, subset= sig %in% c("none")), size=4, alpha=0.3, col='grey')+
  geom_point(data=subset(Day3, subset= sig %in% c("sig")), col="wheat", size=4)+
  geom_point(data=subset(Day3, subset= DEG %in% c("UP", "Down")), aes(col=DEG), size=4)+scale_color_manual(values=c('darkorange','deepskyblue'))+
  geom_vline(xintercept = 0, color="black")+
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed') +geom_vline(xintercept = c(-0.58,0.58), linetype = 'dashed') +scale_y_break(c(15, 70), scales = 0.2)+
  geom_text_repel(aes(label=delabel), max.overlaps = 50, color="black",box.padding = 1, size=10)+
  xlab("log2(Dox/control)")+ylab("-log10(adjusted p-value)")+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"), axis.text = element_text(size = 12,colour = "black"))

normal_count_Day3=t(assay(rlog(dds_Day3, blind = FALSE)))
normal_count_Day3=as.data.frame(normal_count_Day3)
normal_count_Day3$cell=colData_Day3$dss_cell
normal_count_Day3$treatment=colData_Day3$dss_treatment
normal_count_Day3$ID=rownames(normal_count_Day3)

#Day5 (120 hrs)
counts_Day5=counts[,c(13:18)]
colData_Day5=colData[c(13:18),]
dds_Day5 = DESeqDataSetFromMatrix(counts_Day5, colData = colData_Day5, design = ~dss_treatment) ###
dds_Day5 = DESeq(dds_Day5)

Day5=as.data.frame(results(dds_Day5))
Day5=Day5[,c(1,2,6)]
colnames(Day5)=c("RHKI_baseMean","RHKI_LFC", "RHKI_padj")
Day5[,3]=ifelse(is.na(Day5[,3]), 1, Day5[,3])
Day5[,2]=ifelse(is.na(Day5[,2]), 0, Day5[,2])
tmp=rep(NA, nrow(Day5))
tmp1=rep(NA, nrow(Day5))
adj=Day5[,3]
LFC=Day5[,2]

#signficant defined by adjusted p-value <=0.05
for(i in 1:length(adj)){
  if(adj[i] <= 0.05){
    tmp[i] = "sig"
  }else{
    tmp[i] = "none"
  }
}
#up and downregulated gene is with LFC cutoff 0.58 (1.5 fold change)
for(i in 1:length(LFC)){
  if(adj[i] <= 0.05 & LFC[i] >= 0.58){
    tmp1[i] = "UP"
  }else if (adj[i] <= 0.05 & LFC[i] <= -0.58){
    tmp1[i] = "Down"
  }else{
    tmp1[i] = "none"
  }
}
Day5=cbind(Day5, tmp, tmp1)
colnames(Day5)[4:5]=c("sig", "DEG")

Day5$delabel=rep(NA, nrow(Day5))
Day5$symbol=rownames(Day5)
Day5_highlight=c("Cdh1", "Fgfr4", "Bmp7", "Tgfbr3","Nkx3-1", "Cdx1", "Cdx2", "Hoxa1", "Nkx6-1", "Lif", "Foxa2", "Meox2","T", "Oligo2", "Elavl3", "Map2", "Crabp1","Rnaseh1", "Sox8","Sox5", "Neurog3", "Pcdh17", "Elavl4", "Aldh1a2")
for (i in 1:nrow(Day5)){
  if (rownames(Day5)[i] %in% Day5_highlight) {
    Day5$delabel[i] = Day5$symbol[i]
  }
}

ggplot(Day5, aes(x=RHKI_LFC, y=-log10(RHKI_padj)))+xlim(-7.5,7.5)+##ylim(0,140)+
  geom_point(data=subset(Day5, subset= sig %in% c("none")), size=4, alpha=0.3, col='grey')+
  geom_point(data=subset(Day5, subset= sig %in% c("sig")), col="wheat", size=4)+
  geom_point(data=subset(Day5, subset= DEG %in% c("UP", "Down")), aes(col=DEG), size=4, alpha=0.8)+scale_color_manual(values=c('darkorange','deepskyblue'))+
  geom_vline(xintercept = 0, color="black")+
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed') +geom_vline(xintercept = c(-0.58,0.58), linetype = 'dashed') +##scale_y_break(c(30, 50), scales = 0.2)+
  geom_text_repel(aes(label=delabel), max.overlaps = 50, color="black",box.padding = 1, size=10)+
  xlab("log2(Dox/control)")+ylab("-log10(adjusted p-value)")+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"), axis.text = element_text(size = 12,colour = "black"))
genes_Day5=c("Rnaseh1", "Neurog3", "Sox8","Ascl1","Pcdh17", "Crabp1", "Aldh1a2", "Foxa2", "Cdx2", "T", "Nkx3-1", "Hoxa1", "Tgfbr3", "Fgfr4", "Wnt5a")
normal_count_Day5=t(assay(rlog(dds_Day5, blind = FALSE)))
normal_count_Day5=as.data.frame(normal_count_Day5)
normal_count_Day5$cell=colData_Day5$dss_cell
normal_count_Day5$treatment=colData_Day5$dss_treatment
normal_count_Day5$ID=rownames(normal_count_Day5)

# heatmap of differntially expressed genes of interested
UP=rownames(Day5[Day5$DEG %in% c("UP"),])
genes=grep("^Hox", UP, value=TRUE)
genes_Day5=c("Rnaseh1", "Neurog3", "Sox8","Ascl1","Pcdh17", "Hes3", "Crabp1", "Foxa2", "Cdx2", "T", "Tgfbr3", "Fgfr4", "Wnt5a", genes)
test=normal_count_Day5[rownames(Day5) %in% genes_Day5,]
pheatmap(test,cluster_rows = T, scale = "row", cluster_cols=F,color = viridis(n=8, alpha = 1, begin=0, end=1, option = "B"), row_km=2)

#GO analysis on differential expressed gene
Day5UP1_GO = enrichGO(rownames(Day5[Day5$DEG %in% c("UP"),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
Day5DN1_GO = enrichGO(rownames(Day5[Day5$DEG %in% c("Down"),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

p1=dotplot(Day5UP1_GO, showCategory=10)+ggtitle("Upregulation")+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), legend.direction = 'horizontal', legend.position = c(0.8, 0.1), axis.line = element_line(size=0), legend.background = element_rect(fill = "transparent"), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15, angle = 45), legend.key.size = unit(0.4,"cm"))
p2=dotplot(Day5DN1_GO, showCategory=10)+ggtitle("Downregulation")+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), legend.direction = 'horizontal', legend.position = c(0.8, 0.1), axis.line = element_line(size=0),legend.background = element_rect(fill = "transparent"), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15, angle = 45), legend.key.size = unit(0.4,"cm"))

plot=list(p1,p2)
wrap_plots(plot, nrow=1)

##gene expression trend following timepoints of differentially expressed genes 
ma = assay(rlog(dds))
UP=rownames(Day5[Day5$DEG %in% c("UP"),])
UP=setdiff(UP, c("Rnaseh1", "Hprt"))
DN=rownames(Day5[Day5$DEG %in% c("Down"),])
ma2=ma[rownames(ma) %in% UP,]
ma3=ma[rownames(ma) %in% DN,]

clusters <- degPatterns(ma2, metadata = colData, time = "dss_time", col = "dss_treatment")
clusters1 <- degPatterns(ma3, metadata = colData, time = "dss_time", col = "dss_treatment")

h=clusters$plot$data
h1=clusters1$plot$data
pheatmap(ma3,cluster_rows = T, scale = "row", cluster_cols=F,color = viridis(n=8, alpha = 1, begin=0, end=1, option = "B"), row_km=5)
p1=ggplot(h, aes(x=dss_time, y=value, group=merge, fill=dss_treatment)) + ylim(-2,2) +geom_violin(draw_quantiles = c(0.5))+xlab("")+ylab("z-score expression")+ggtitle("Upregulation(168-59)")+
  geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size = 1)+geom_smooth(method = "gam",se = TRUE,aes(x=dss_time, y=value, group=dss_treatment, color = dss_treatment), formula = y ~ s(x, bs = "cs", k = 3), linetype='dashed', alpha=0.3)+
  scale_color_manual(values=c('black','red'))+scale_fill_manual(values=c('grey','red'))+facet_wrap(~cluster, scales = "fixed", nrow = 1)+
  theme(strip.placement = "outside", panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"), legend.position = "none", axis.text = element_text(size = 12,colour = "black"))
p2=ggplot(h1, aes(x=dss_time, y=value, group=merge, fill=dss_treatment))+ ylim(-2,2) + geom_violin(draw_quantiles = c(0.5))+
  geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size = 1)+geom_smooth(method = "gam",se = TRUE,aes(x=dss_time, y=value, group=dss_treatment, color = dss_treatment), formula = y ~ s(x, bs = "cs", k = 3), linetype='dashed', alpha=0.3)+
  scale_color_manual(values=c('black','red'))+scale_fill_manual(values=c('grey','red'))+facet_wrap(~cluster, scales = "fixed", nrow = 1)+xlab("")+ylab("z-score expression")+ggtitle("downregulation(118-72)")+
  theme(strip.placement = "outside", panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"),axis.text = element_text(size = 12,colour = "black"))
plots=list(p1,p2)
wrap_plots(plots, nrow=2)


write.csv(Day0, '/Users/ChaoC1/Documents/figure/Day0_LFC.csv')
write.csv(Day3, '/Users/ChaoC1/Documents/figure/Day3_LFC.csv')
write.csv(Day5, '/Users/ChaoC1/Documents/figure/Day5_LFC.csv')

write.csv(normal_count_Day0, '/Users/ChaoC1/Documents/figure/Day0_rlog.csv')
write.csv(normal_count_Day3, '/Users/ChaoC1/Documents/figure/Day3_rlog.csv')
write.csv(normal_count_Day5, '/Users/ChaoC1/Documents/figure/Day5_rlog.csv')

saveRDS(dds_Day0, '/Users/ChaoC1/Documents/figure/bulk_Day0.rds')
saveRDS(dds_Day3, '/Users/ChaoC1/Documents/figure/bulk_Day3.rds')
saveRDS(dds_Day5, '/Users/ChaoC1/Documents/figure/bulk_Day5.rds')
saveRDS(dds, '/Users/ChaoC1/Documents/figure/bulk_all.rds')
