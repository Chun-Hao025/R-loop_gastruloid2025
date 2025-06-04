library(ggplot2)
library(ggrepel)
library(patchwork)
library(matrixStats)

ATAC_exp=read.delim("/path/to/ATAC_expressionlevel_enrich_matrix.tab", sep="\t", header=FALSE, skip = 2)
ATAC_exp=as.data.frame(ATAC_exp)
ATAC_exp=ATAC_exp[order(ATAC_exp$V2),]
colnames(ATAC_exp)=c("ID", "peak", 1:400)
peak=unique(ATAC_exp$peak)
for (i in 1:length(peak)){
  tmp=ATAC_exp[ATAC_exp$peak %in% peak[i],]
  tmp2=t(tmp)
  tmp2=as.data.frame(tmp2)
  colnames(tmp2)=tmp$ID
  tmp2=tmp2[c(-1,-2),]
  tmp2=sapply(tmp2, function(x) as.numeric(as.character(x)))
  K=c(rowMeans(tmp2[,c(1,2,3,4)]),rowMeans(tmp2[,c(5,6,7,8)]))
  L=c(rowSds(tmp2[,c(1,2,3,4)], useNames = F)/sqrt(4),rowSds(tmp2[,c(5,6,7,8)], useNames = F)/sqrt(4))
  tmp3=data.frame(K,L)
  colnames(tmp3)=c("avg", "err")
  tmp3$treatment=c(rep("control", 400), rep("Dox", 400))
  tmp3$position=rep(c(1:400), 2)
  tmp3$region=rep(peak[i], length(tmp3))
  if(i==1){
    obj5=tmp3
  }else{
    obj5=rbind(obj5, tmp3)
  }
}
obj5$region=factor(obj5$region, levels = c("high", "medium", "low"))
ggplot(obj5, aes(x=position, y=avg, color=treatment, group=treatment)) + 
  geom_errorbar(aes(ymin= avg - err, ymax = avg + err), alpha=0.5, linewidth=0.5, width=0.1, linetype="solid")+
  geom_line()+geom_vline(xintercept = c(200), linetype="dotted")+
  scale_colour_manual(values=c("black", "red", "deepskyblue1", 'orange'))+facet_wrap(~region)+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15), legend.key.size = unit(0.8,"cm"))
