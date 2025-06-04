library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(ChIPseeker)
library(plyr)
library(graphics)
library(patchwork)
library(RColorBrewer)
library(matrixStats)

#meta of differential peak of H2A.Z
meta = read.csv("/path/to/meta_analysis.csv", header=T)
meta_sub = subset(meta, type %in% c("genome","upregulation", "downregulation"))
ggplot(meta, aes(fill=Annotation, y=value, x=type)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette="Set3")+
  scale_x_discrete(limits=c("genome","upregulation", "downregulation"))+scale_y_continuous(labels = scales::percent)+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"))

#GO of differntial peaks of H2A.Z
down_ann = read.table("/path/to/H2AZ/down_annotation.txt", header = T)
UP_ann = read.table("/path/to/H2AZ/up_annotation.txt", header = T)
egoU = enrichGO(gene = UP_ann$Gene.Name, keyType = "SYMBOL", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
egoD = enrichGO(gene = down_ann$Gene.Name, keyType = "SYMBOL", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
p1 = dotplot(egoU, showCategory=12)+ggtitle("upregulation")+
  theme(plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "snow3", size=2))
p2 = dotplot(egoD, showCategory=12)+ggtitle("downregulation")+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"))
plots = list(p1, p2)
names(plots) = c("upregulation", "downregulation")
wrap_plots(plots)

##H2AZ chromHMM plotting
H2AZ=read.delim("/path/to/H2AZ_chromHMM_matrix_enrich.tab", sep="\t", header=FALSE, skip = 2)
H2AZ=as.data.frame(H2AZ)
H2AZ=H2AZ[order(H2AZ$V2),]
colnames(H2AZ)=c("ID", "HMM", 1:200)
HMM=unique(H2AZ$HMM)
for (i in 1:length(HMM)){
  tmp=H2AZ[H2AZ$HMM %in% HMM[i],]
  tmp2=t(tmp)
  tmp2=as.data.frame(tmp2)
  colnames(tmp2)=tmp$ID
  tmp2=tmp2[c(-1,-2),]
  tmp2=sapply(tmp2, function(x) as.numeric(as.character(x)))
  K=c(rowMeans(tmp2[,c(1,2,3,4)]),rowMeans(tmp2[,c(5,6,7,8)]))
  L=c(rowSds(tmp2[,c(1,2,3,4)])/sqrt(4),rowSds(tmp2[,c(5,6,7,8)])/sqrt(4))
  tmp3=data.frame(K,L)
  colnames(tmp3)=c("avg", "err")
  tmp3$treatment=c(rep("control", 200), rep("Dox", 200))
  tmp3$position=rep(c(1:200), 2)
  tmp3$region=rep(HMM[i], length(tmp3))
  if(i==1){
    obj=tmp3
  }else{
    obj=rbind(obj, tmp3)
  }
}
ggplot(obj, aes(x=position, y=avg, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= avg - err, ymax = avg + err), alpha=0.5, linewidth=0.5, width=0.1, linetype="solid")+
  geom_line()+geom_vline(xintercept = c(50, 150), linetype="dotted")+xlim(0,200)+
  scale_colour_manual(values=c("black", "red"))+facet_wrap(~region)+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"))

file=c("H2AZ_TSS.txt", "H2AZ_DHS.txt", "H2AZ_DHS_noTSS.txt")
type=c("TSS","DHS", "DHS(noTSS)")
for (i in 1:length(file)){
  tmp=read.delim(paste0("/path/to/", file[i]), sep="\t", header=TRUE, row.names = 1)
  tmp=tmp[,c(1,4,7,10,13,16,19,22)]
  tmp1=sapply(tmp, function(x) as.numeric(as.character(x)))
  K=c(rowMeans(tmp1[,c(1,2,3,4)]),rowMeans(tmp1[,c(5,6,7,8)]))
  L=c(rowSds(tmp1[,c(1,2,3,4)], useNames = FALSE)/sqrt(4),rowSds(tmp1[,c(5,6,7,8)], useNames = FALSE)/sqrt(4))
  tmp2=data.frame(K,L)
  colnames(tmp2)=c("avg", "err")
  tmp2$treatment=c(rep("control", 201), rep("Dox", 201))
  tmp2$position=rep(c(rownames(tmp)),2)
  tmp2$region=rep(type[i], nrow(tmp2))
  if(i==1){
    obj=tmp2
  }else{
    obj=rbind(obj, tmp2)
  }
}
obj$region=factor(obj$region, levels=c("TSS","DHS", "DHS(noTSS)"))
obj$position=as.numeric(obj$position)
ggplot(obj, aes(x=position, y=avg, color=treatment, group=treatment)) + 
  geom_errorbar(aes(ymin= avg - err, ymax = avg + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+facet_wrap(~region, nrow=1)+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"))

###aggregation plot from HOMER annotatePeak.pl
files=c("CT_K4me3_TSS.csv", "CT_K27ac_peakcenter.csv", "CT_K27me3_peakcenter.csv", "CT_K9me3_peakcenter.csv", "CT_H2AZ_peakcenter.csv", "CT_H4ac_TSS.csv",
        "CT_K36me3_meta.csv", "CT_RNAP2S2P_meta.csv", "CR_CTCF_TSS.csv", "CT_H33_TSS.csv", "CT_ZX_H2AZ_peakcenter.csv")
dirs=paste0("/path/to/", files)

K4=read.csv(dirs[1], header=T)
colnames(K4)=c("position", "control_rep1", "control_rep2", "Dox_rep1", "Dox_rep2")
K=c(rowMeans(K4[,c(2,3)]), rowMeans(K4[,c(4,5)]))
L=c(rowSds(as.matrix(K4[,c(2,3)]),useNames = FALSE)/sqrt(2), rowSds(as.matrix(K4[,c(4,5)]),useNames = FALSE)/sqrt(2))
K4_L=data.frame(rep(K4$position,2),K,L,c(rep("control", nrow(K4)), rep("Dox", nrow(K4))))
colnames(K4_L)=c("position", "mean", "err", "treatment")
p1=ggplot(K4_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("TSS")+ggtitle("H3K4me3")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),legend.position = "none")
  
K27ac=read.csv(dirs[2], header=T)
colnames(K27ac)=c("position", "control_rep1", "control_rep2", "Dox_rep1", "Dox_rep2")
K=c(rowMeans(K27ac[,c(2,3)]), rowMeans(K27ac[,c(4,5)]))
L=c(rowSds(as.matrix(K27ac[,c(2,3)]),useNames = FALSE)/sqrt(2), rowSds(as.matrix(K27ac[,c(4,5)]),useNames = FALSE)/sqrt(2))
K27ac_L=data.frame(rep(K27ac$position,2),K,L,c(rep("control", nrow(K27ac)), rep("Dox", nrow(K27ac))))
colnames(K27ac_L)=c("position", "mean", "err", "treatment")
p2=ggplot(K27ac_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("peakcenter")+ggtitle("H3K27ac")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),legend.position = "none")

K27me=read.csv("/path/to/CT_K27me3_peakcenter1.csv", header=T)
colnames(K27me)[1]="position"
colnames(K27me)=c("position", "control_rep1", "control_rep2", "control_rep3", "control_rep4", "Dox_rep1", "Dox_rep2", "Dox_rep3", "Dox_rep4")
K=c(rowMeans(K27me[,c(2,3,4,5)]), rowMeans(K27me[,c(6,7,8,9)]))
L=c(rowSds(as.matrix(K27me[,c(2,3,4,5)]),useNames = FALSE)/sqrt(4), rowSds(as.matrix(K27me[,c(6,7,8,9)]),useNames = FALSE)/sqrt(4))
K27me_L=data.frame(rep(K27me$position,2),K,L,c(rep("control", nrow(K27me)), rep("Dox", nrow(K27me))))
colnames(K27me_L)=c("position", "mean", "err", "treatment")
p3=ggplot(K27me_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("peakcenter")+ggtitle("H3K27me3")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),legend.position = "none")
###
K27me_ST=K27me %>% melt(id.var=c("position"))
K27me_ST$variable=as.character(K27me_ST$variable)
for (i in 1:nrow(K27me_ST)){
  tmp=strsplit(K27me_ST$variable[i], split = "_")[[1]]
  K27me_ST$treatment[i]=tmp[1]
  K27me_ST$replicate[i]=tmp[2]
}
K27me_ST$treatment=factor(K27me_ST$treatment, levels = c("control", "Dox"))
model1 <- lmer(value ~ treatment + (1 | position/replicate), data = K27me_ST)
summary(model1)
emm1 <- emmeans(model1, ~ treatment)
pairs(emm1)
summary(emm1)
plot(emm1)
###
K9me=read.csv(dirs[4], header=T)
colnames(K9me)=c("position", "control_rep1", "control_rep2", "Dox_rep1", "Dox_rep2")
K=c(rowMeans(K9me[,c(2,3)]), rowMeans(K9me[,c(4,5)]))
L=c(rowSds(as.matrix(K9me[,c(2,3)]),useNames = FALSE)/sqrt(2), rowSds(as.matrix(K9me[,c(4,5)]),useNames = FALSE)/sqrt(2))
K9me_L=data.frame(rep(K9me$position,2),K,L,c(rep("control", nrow(K9me)), rep("Dox", nrow(K9me))))
colnames(K9me_L)=c("position", "mean", "err", "treatment")
p4=ggplot(K9me_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("peakcenter")+ggtitle("H3K9me3")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),legend.position = "none")

H2A=read.csv(dirs[5], header=T)
colnames(H2A)=c("position", "control_rep1", "control_rep2", "control_rep3", "control_rep4", "Dox_rep1", "Dox_rep2", "Dox_rep3", "Dox_rep4")
K=c(rowMeans(H2A[,c(2,3,4,5)]), rowMeans(H2A[,c(6,7,8,9)]))
L=c(rowSds(as.matrix(H2A[,c(2,3,4,5)]),useNames = FALSE)/sqrt(4), rowSds(as.matrix(H2A[,c(6,7,8,9)]),useNames = FALSE)/sqrt(4))
H2A_L=data.frame(rep(H2A$position,2),K,L,c(rep("control", nrow(H2A)), rep("Dox", nrow(H2A))))
colnames(H2A_L)=c("position", "mean", "err", "treatment")
p5=ggplot(H2A_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("peakcenter")+ggtitle("H2A.Z")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),legend.position = "none")

H4ac=read.csv(dirs[6], header=T)
colnames(H4ac)=c("position", "control_rep1", "control_rep2", "Dox_rep1", "Dox_rep2")
K=c(rowMeans(H4ac[,c(2,3)]), rowMeans(H4ac[,c(4,5)]))
L=c(rowSds(as.matrix(H4ac[,c(2,3)]),useNames = FALSE)/sqrt(2), rowSds(as.matrix(H4ac[,c(4,5)]),useNames = FALSE)/sqrt(2))
H4ac_L=data.frame(rep(H4ac$position,2),K,L,c(rep("control", nrow(H4ac)), rep("Dox", nrow(H4ac))))
colnames(H4ac_L)=c("position", "mean", "err", "treatment")
p6=ggplot(H4ac_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("TSS")+ggtitle("H4ac")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),legend.position = "none")

K36me=read.csv(dirs[7], header=T)
colnames(K36me)=c("position", "control_rep1", "control_rep2", "control_rep3", "control_rep4", "Dox_rep1", "Dox_rep2", "Dox_rep3", "Dox_rep4")
K=c(rowMeans(K36me[,c(2,3,4,5)]), rowMeans(K36me[,c(6,7,8,9)]))
L=c(rowSds(as.matrix(K36me[,c(2,3,4,5)]),useNames = FALSE)/sqrt(4), rowSds(as.matrix(K36me[,c(6,7,8,9)]),useNames = FALSE)/sqrt(4))
K36me_L=data.frame(rep(K36me$position,2),K,L,c(rep("control", nrow(K36me)), rep("Dox", nrow(K36me))))
colnames(K36me_L)=c("position", "mean", "err", "treatment")
p7=ggplot(K36me_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("TSS")+ggtitle("H3K36me3")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),legend.position = "none")

RNAP=read.csv(dirs[8], header=T)
colnames(RNAP)=c("position", "control_rep1", "control_rep2", "Dox_rep1", "Dox_rep2")
K=c(rowMeans(RNAP[,c(2,3)]), rowMeans(RNAP[,c(4,5)]))
L=c(rowSds(as.matrix(RNAP[,c(2,3)]),useNames = FALSE)/sqrt(2), rowSds(as.matrix(RNAP[,c(4,5)]),useNames = FALSE)/sqrt(2))
RNAP_L=data.frame(rep(RNAP$position,2),K,L,c(rep("control", nrow(RNAP)), rep("Dox", nrow(RNAP))))
colnames(RNAP_L)=c("position", "mean", "err", "treatment")
p8=ggplot(RNAP_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("TSS")+ggtitle("RNAP2-S2P")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),legend.position = "none")

CTCF=read.csv(dirs[9], header=T)
colnames(CTCF)=c("position", "control_rep1", "control_rep2", "Dox_rep1", "Dox_rep2")
K=c(rowMeans(CTCF[,c(2,3)]), rowMeans(CTCF[,c(4,5)]))
L=c(rowSds(as.matrix(CTCF[,c(2,3)]),useNames = FALSE)/sqrt(2), rowSds(as.matrix(CTCF[,c(4,5)]),useNames = FALSE)/sqrt(2))
CTCF_L=data.frame(rep(CTCF$position,2),K,L,c(rep("control", nrow(CTCF)), rep("Dox", nrow(CTCF))))
colnames(CTCF_L)=c("position", "mean", "err", "treatment")
p9=ggplot(CTCF_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("TSS")+ggtitle("CTCF")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),legend.position = "none")

H33=read.csv(dirs[10], header=T)
colnames(H33)=c("position", "control_rep1", "control_rep2", "Dox_rep1", "Dox_rep2")
K=c(rowMeans(H33[,c(2,3)]), rowMeans(H33[,c(4,5)]))
L=c(rowSds(as.matrix(H33[,c(2,3)]),useNames = FALSE)/sqrt(2), rowSds(as.matrix(H33[,c(4,5)]),useNames = FALSE)/sqrt(2))
H33_L=data.frame(rep(H33$position,2),K,L,c(rep("control", nrow(H33)), rep("Dox", nrow(H33))))
colnames(H33_L)=c("position", "mean", "err", "treatment")
p10=ggplot(H33_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("TSS")+ggtitle("H3.3")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),legend.position = "none")

ZX=read.csv(dirs[11], header=T)
colnames(ZX)=c("position", "control_rep1", "control_rep2","control_rep4", "control_rep5", "control_rep6", "Dox_rep1", "Dox_rep2", "Dox_rep3", "Dox_rep4", "Dox_rep4")
K=c(rowMeans(ZX[,c(2,3,4,5,6)]), rowMeans(ZX[,c(7,8,9,10,11)]))
L=c(rowSds(as.matrix(ZX[,c(2,3,4,5,6)]),useNames = FALSE)/sqrt(5), rowSds(as.matrix(ZX[,c(7,8,9,10,11)]),useNames = FALSE)/sqrt(5))
ZX_L=data.frame(rep(ZX$position,2),K,L,c(rep("control", nrow(ZX)), rep("Dox", nrow(ZX))))
colnames(ZX_L)=c("position", "mean", "err", "treatment")
p11=ggplot(ZX_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("peakcenter")+ggtitle("ZX1 H2A.Z")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"))

K4P=read.csv("/path/to/CT_K4me3_peakcenter.csv", header=T)
colnames(K4P)=c("position", "control_rep1", "control_rep2", "Dox_rep1", "Dox_rep2")
K=c(rowMeans(K4P[,c(2,3)]), rowMeans(K4P[,c(4,5)]))
L=c(rowSds(as.matrix(K4P[,c(2,3)]),useNames = FALSE)/sqrt(2), rowSds(as.matrix(K4P[,c(4,5)]),useNames = FALSE)/sqrt(2))
K4P_L=data.frame(rep(K4P$position,2),K,L,c(rep("control", nrow(K4P)), rep("Dox", nrow(K4P))))
colnames(K4P_L)=c("position", "mean", "err", "treatment")
p12=ggplot(K4P_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("peakcenter")+ggtitle("H3K4me3")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"),legend.position = "none")

plot=list(p12, p2, p3, p4,p5,p6,p7,p8,p9,p10,p11)
wrap_plots(plot, nrow=3)

RL=read.csv("/path/to/MapR_TSS.csv", header=T)
colnames(RL)=c("position", "control_rep1", "control_rep2", "control_rep3", "control_rep4", "RA_rep1", "RA_rep2", "RA_rep3", "RA_rep4", "pA_rep1", "pA_rep2", "pA_rep3", "pA_rep4")
K=c(rowMeans(RL[,c(2,3,4,5)]), rowMeans(RL[,c(6,7,8,9)]), rowMeans(RL[,c(10,11,12,13)]))
L=c(rowSds(as.matrix(RL[,c(2,3,4,5)]),useNames = FALSE)/sqrt(4), rowSds(as.matrix(RL[,c(6,7,8,9)]),useNames = FALSE)/sqrt(4), rowSds(as.matrix(RL[,c(10,11,12,13)]),useNames = FALSE)/sqrt(4))
RL_L=data.frame(rep(RL$position,3),K,L,c(rep("MapR", nrow(RL)), rep("MapR+RNaseA", nrow(RL)), rep("pAMN", nrow(RL))))
colnames(RL_L)=c("position", "mean", "err", "treatment")
ggplot(RL_L, aes(x=position, y=mean, group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin= mean - err, ymax = mean + err), alpha=0.5, linewidth=0.5, width=2, linetype="solid")+
  geom_line(linewidth=1, alpha=1)+geom_vline(xintercept = c(0), linetype="dotted")+ylab("Mean of CPM")+xlab("TSS")+ggtitle("MapR")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", "darkorange"))+
  theme(panel.background = element_rect(color = "black", size=1, fill = "white",), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.text = element_text(size=20, color='black'), axis.text.x = element_text(angle = 0, hjust = 0.5),axis.title.x = element_text(size=0, color='black'), axis.title.y = element_text(size=20, color='black'), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=20), text = element_text(size=20),legend.key.size = unit(0.8,"cm"))


