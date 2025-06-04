library(circlize)
library(matrixStats)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(viridis)


###ChromHMM state
state=read.csv("/path/to/state_statistics.csv", header=T, row.names = 1)
col_fun = colorRamp2(c(0,2,10,60), c("white","yellow", "red","violetred"))

#reorder and name state based on the epigenetic features
#state=c("E2", "E9", "E3", "E10", "E5", "E8", "E6", "E1", "E7")
label=c("Active TSS", "Bivalent TSS", "Enhancer", "H2A.Z high", "Active Transcribed", "H3K27me3 heterochromatin", "H3K9me3 heterochromatin", "Insulator", "Intergenic region", "Intergenic region")

ht1=Heatmap(state[1:10,1:12], cluster_rows = F, cluster_columns = F, col=viridis(n=10, alpha = 0.5, begin=0, end=1, option = "D"), border = T, row_labels = label, row_names_side = "left", row_split = c(rep("A", 5), rep("B",2), rep("C",3)),border_gp = gpar(lwd = 2,col="black"), )
ht2 = Heatmap(state[1:10,14:16], cluster_rows = F, cluster_columns = F, col=col_fun,border = T,border_gp = gpar(lwd = 2,col="black"),cell_fun = function(j, i, x, y, width, height, fill) {
  if(state[1:10,14:16][i, j] > 1.5)
    grid.text(sprintf("%.1f", state[1:10,14:16][i, j]), x, y, gp = gpar(fontsize = 10))
})

ht1+ht2

##plot MNase-seq result from deeptools
## MNase-seq split by expression level
sub_exp=read.delim("/path/to/MN_sub_expressionlevel_enrich_matrix.tab", sep="\t", header=FALSE, skip = 2)
sub_exp=as.data.frame(sub_exp)
sub_exp=sub_exp[order(sub_exp$V2),]
colnames(sub_exp)=c("ID", "peak", 1:400)
peak=unique(sub_exp$peak)
for (i in 1:length(peak)){
  tmp=sub_exp[sub_exp$peak %in% peak[i],]
  tmp2=t(tmp)
  tmp2=as.data.frame(tmp2)
  colnames(tmp2)=tmp$ID
  tmp2=tmp2[c(-1,-2),]
  tmp2=sapply(tmp2, function(x) as.numeric(as.character(x)))
  K=c(rowMeans(tmp2[,c(1,2)]),rowMeans(tmp2[,c(3,4)]))
  L=c(rowSds(tmp2[,c(1,2)])/sqrt(2),rowSds(tmp2[,c(3,4)])/sqrt(2))
  tmp3=data.frame(K,L)
  colnames(tmp3)=c("avg", "err")
  tmp3$treatment=c(rep("control", 400), rep("Dox", 400))
  tmp3$position=rep(c(1:400), 2)
  tmp3$region=rep(peak[i], length(tmp3))
  if(i==1){
    obj=tmp3
  }else{
    obj=rbind(obj, tmp3)
  }
}
obj$region=factor(obj$region, levels = c("high", "medium", "low"))
ggplot(obj, aes(x=position, y=avg, color=treatment)) + 
  geom_errorbar(aes(ymin= avg - err, ymax = avg + err), alpha=0.5, linewidth=0.5, width=0.1, linetype="solid")+
  geom_line()+geom_vline(xintercept = c(200), linetype="dotted")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1"))+facet_wrap(~region)+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"))

mono_exp=read.delim("/path/to/MN_mono_expressionlevel_enrich_matrix.tab", sep="\t", header=FALSE, skip = 2)
mono_exp=as.data.frame(mono_exp)
mono_exp=mono_exp[order(mono_exp$V2),]
colnames(mono_exp)=c("ID", "peak", 1:400)
peak=unique(mono_exp$peak)
for (i in 1:length(peak)){
  tmp=mono_exp[mono_exp$peak %in% peak[i],]
  tmp2=t(tmp)
  tmp2=as.data.frame(tmp2)
  colnames(tmp2)=tmp$ID
  tmp2=tmp2[c(-1,-2),]
  tmp2=sapply(tmp2, function(x) as.numeric(as.character(x)))
  K=c(rowMeans(tmp2[,c(1,2)]),rowMeans(tmp2[,c(3,4)]))
  L=c(rowSds(tmp2[,c(1,2)])/sqrt(2),rowSds(tmp2[,c(3,4)])/sqrt(2))
  tmp3=data.frame(K,L)
  colnames(tmp3)=c("avg", "err")
  tmp3$treatment=c(rep("control", 400), rep("Dox", 400))
  tmp3$position=rep(c(1:400), 2)
  tmp3$region=rep(peak[i], length(tmp3))
  if(i==1){
    obj=tmp3
  }else{
    obj=rbind(obj, tmp3)
  }
}
obj$region=factor(obj$region, levels = c("high", "medium", "low"))
ggplot(obj, aes(x=position, y=avg, color=treatment)) + 
  geom_errorbar(aes(ymin= avg - err, ymax = avg + err), alpha=0.5, linewidth=0.5, width=0.1, linetype="solid")+
  geom_line()+geom_vline(xintercept = c(200), linetype="dotted")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1"))+facet_wrap(~region)+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"))

##MNase-seq split by chromHMM state
sub=read.delim("/path/to/MN_sub_chromHMM_matrix_enrich.tab", sep="\t", header=FALSE, skip = 2)
sub=as.data.frame(sub)
sub=sub[order(sub$V2),]
colnames(sub)=c("ID", "HMM", 1:200)
HMM=unique(sub$HMM)
for (i in 1:length(HMM)){
  tmp=sub2[sub2$HMM %in% HMM[i],]
  tmp2=t(tmp)
  tmp2=as.data.frame(tmp2)
  colnames(tmp2)=tmp$ID
  tmp2=tmp2[c(-1,-2),]
  tmp2=sapply(tmp2, function(x) as.numeric(as.character(x)))
  K=c(rowMeans(tmp2[,c(1,2)]),rowMeans(tmp2[,c(3,4)]))
  L=c(rowSds(tmp2[,c(1,2)])/sqrt(2),rowSds(tmp2[,c(3,4)])/sqrt(2))
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

mono=read.delim("/path/to/MN_mono_chromHMM_matrix_enrich.tab", sep="\t", header=FALSE, skip = 2)
mono=as.data.frame(mono)
mono=mono[order(mono$V2),]
colnames(mono)=c("ID", "HMM", 1:200)
HMM=unique(mono$HMM)
for (i in 1:length(HMM)){
  tmp=mono2[mono2$HMM %in% HMM[i],]
  tmp2=t(tmp)
  tmp2=as.data.frame(tmp2)
  colnames(tmp2)=tmp$ID
  tmp2=tmp2[c(-1,-2),]
  tmp2=sapply(tmp2, function(x) as.numeric(as.character(x)))
  K=c(rowMeans(tmp2[,c(1,2)]),rowMeans(tmp2[,c(3,4)]))
  L=c(rowSds(tmp2[,c(1,2)])/sqrt(2),rowSds(tmp2[,c(3,4)])/sqrt(2))
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
  geom_line()+geom_vline(xintercept = c(50, 150), linetype="dashed")+xlim(0,200)+
  scale_colour_manual(values=c("black", "red"))+facet_wrap(~region)+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"))
