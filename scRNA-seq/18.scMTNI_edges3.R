library(data.table)
library(ComplexHeatmap)
library(ggplot2)
library(ggraph)
library(viridis)
library(patchwork)
library(org.Mm.eg.db)
library(clusterProfiler)
library(reshape2)
library(dplyr)
library(tidygraph)

### This script was used for filtering the TF-targets pair that identified >= 80% of the trial and defined it as the high confident edges.
cluster=c("cluster0", "cluster1", "cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8", "cluster9","cluster10",
          "cluster11","cluster12","cluster13","cluster14","cluster15")
cell=c('Caudal_Neuroectoderm1', 'Definitive_Endoderm_Guts', 'Endothelium', 'Epiblasts', 
              'Naive_Pluripotent_Cells', 'Nascent_Mesoderm', 'Neuromesodermal_Progenitors1', 
              'Neuromesodermal_Progenitors2', 'Neuromesodermal_Progenitors3', 'Neuron-like_Cells',
              'Paraxial_Mesoderm_A', 'Paraxial_Mesoderm_B', 'Placodal_area', 'Primitive_Streak', 
              'Spinal_Cord', 'Splanchnic_Mesoderm')
combine=NULL
for (i in 1:length(cluster)) {
    tmp1=fread(paste0("/path/to/GRN/scMTNI/scMTNI_control", cluster[i], "/net_alledge_origcoef.txt"))  ### change to the correct file corresponded to all control, RHKI control or RHKI Dox from scRNA-seq_scMTNI_edges2.py
    colnames(tmp1)=c("TF", "target", "weight", "orig_coef")
    tmp1$TF=sapply(strsplit(tmp1$TF, "_"), function(x) x[c(1)])
    tmp1$target=sapply(strsplit(tmp1$target, "_"), function(x) x[c(1)])
    tmp1=tmp1[tmp1$weight >=0.8,]
    tmp1$ann=cell[i]
    tmp1$ID=cluster[i]
    if(i ==1){
        combine=tmp1
    }else{
        combine=rbind(combine, tmp1)
    }
}
write.table(combine, file=paste0("/path/to/GRN/scMTNI/scMTNI_control/combine_edge_orig.txt"), sep="\t", row.names=F, col.names=T, quote=F)

###all control: allcontrol_edge.txt
###RHKI (control and Dox): RHKI_edge.txt
samples=c("scMTNI_control", "scMTNI_RH_control", "scMTNI_RH_Dox")

for (i in 1:length(samples)){
    tmp=fread(paste0("/path/to/GRN/scMTNI/", samples[i], "/combine_edge_orig.txt"))
    if (i ==1){
        obj=tmp
    }else{
        obj=rbind(obj, tmp)
    }
}
write.table(obj, file=paste0("/path/to/GRN/scMTNI/combine_edge_orig.txt"), sep="\t", row.names=F, col.names=T, quote=F)

###network plot
###
###filter the TF and cell type desired to plot
tmp=subset(obj, subset = TF %in% c("Rarb") & ID %in% c("cluster14_control", "cluster14_Dox"))
tmp=tmp[,c("TF", "target", "orig_coef", "ann")]
colnames(tmp)=c("source", "target", "orig_coef", "ann")

##define unique edges
edges_status = tmp %>%
  group_by(source, target) %>%
  mutate(n_anns = n_distinct(ann)) %>%
  ungroup() %>%
  mutate(status = if_else(n_anns > 1, "shared", ann))

g_all = edges_status %>%
  dplyr::select(source, target, orig_coef, ann, status) %>%
  as_tbl_graph(directed=TRUE)
g_all = g_all %>%
  activate(nodes) %>%
  mutate(
    degree    = centrality_degree(mode="all"),
    gene_type = case_when(
      name %in% unique(edges_status$source) ~ "TF",
      TRUE                                 ~ "Gene"
    )
  )

###depends on how many edges, choose different layout for better visualization
lay = create_layout(g_all,layout  = "sparse_stress",pivots  = 1000,weights = NA)
#lay = create_layout(g_all,layout  = "nicely")

### genes would like to highlight on the plot
label=c("Rarb", "Pax6", "Tenm3", "Lrig1", "Nrcam", "Fgfr2", "Dlg2", "Crabp2")
lay$highlight <- lay$name %in% label

p <- ggraph(lay) +
  # edges: width ~ abs(Gain), color ~ sign(Cor), alpha fixed
  #geom_edge_link(aes(edge_width = abs(Gain),color      = sign(Cor),linetype   = status),
  geom_edge_link(aes(edge_width = abs(orig_coef),color      = sign(orig_coef),linetype   = status),
                 arrow = arrow(type="closed", length = unit(1.5, "mm")),end_cap = circle(5, "pt"),show.legend = TRUE, alpha=0.6) +
  #geom_edge_fan(aes(edge_width = abs(orig_coef),color      = sign(orig_coef),linetype   = status),
  #              arrow    = arrow(type="closed", length = unit(1.5, "mm")),end_cap  = circle(5, "pt"),show.legend = TRUE, strength = 6) +
  # nodes: shape for TF vs Gene; size by degree
  geom_node_point(aes(
    shape = gene_type,
    size  = degree
  ),
  fill = "white",
  stroke = 0.5
  ) +
  # labels for TFs only
  #geom_node_text(aes(label = ifelse(gene_type=="TF", name, NA)),repel = TRUE,size = 10) +
  geom_node_text(data    = filter(lay, name %in% label),aes(x = x, y = y, label = name),repel= TRUE,size= 5,nudge_y = 0.05)+
  #geom_mark_hull(data = filter(lay2, !is.na(early_K10)),aes(x = x, y = y, group = early_K10, fill = early_K10),concavity = 5, alpha = 0.2,label.buffer = unit(2, "mm"), label.fontsize = 4) +
  geom_node_point(data = subset(lay, highlight),aes(x = x, y = y, size = degree),shape = 16,color = "red",alpha = 0.8) +
  # diverging color scale for Cor (activation vs repression)
  scale_edge_color_gradient2(low= "darkorange",mid= "white",high= "dodgerblue",midpoint= 0,name = "Cor") +
  # line‐type for shared vs unique
  scale_edge_linetype_manual(
    values = c(shared="solid", Caudal_Neuroectoderm1_control = "dashed", Caudal_Neuroectoderm1_Dox = "dashed", Primitive_Streak_control="dashed", Primitive_Streak_Dox="dashed", Epiblasts_control="dashed", Epiblasts_Dox="dashed", Splanchnic_Mesoderm_control="dashed", Splanchnic_Mesoderm_Dox="dashed", 
               Neuromesodermal_Progenitors1_control="dashed", Neuromesodermal_Progenitors1_Dox="dashed", Neuromesodermal_Progenitors2_control="dashed", Neuromesodermal_Progenitors2_Dox="dashed",Spinal_Cord_control="dashed", Spinal_Cord_Dox="dashed", Paraxial_Mesoderm_A_control="dashed", Paraxial_Mesoderm_A_Dox="dashed",
               Paraxial_Mesoderm_B_control="dashed", Paraxial_Mesoderm_B_Dox="dashed", Nascent_Mesoderm_control="dashed", Nascent_Mesoderm_Dox="dashed", Caudal_Neuroectoderm1 = "dashed", Primitive_Streak="dashed", Epiblasts="dashed"),
    name   = "Edge status"
  ) +
  # control the edge width range
  scale_edge_width_continuous(range = c(0.05, 1.5),name  = "|Gain|") +
  # node shapes and sizes
  scale_shape_manual(values = c(TF=21, Gene=16)) +
  guides(
    edge_alpha = "none",
    size       = guide_legend("Node degree"),
    shape      = guide_legend("Node type")
  ) +
  # facet edges by the original condition
  facet_edges(~ ann) +
  theme_graph(base_family="sans") +
  ##ggtitle("TF–Target Networks: Condition A vs Condition B (shared edges solid)")

print(p)

###Upset plot
###split by either positive correlation or negative correlation edges
genelist=list(control_SP_P=subset(tmp, subset= ID %in% c("cluster14_control") & orig_coef >0)$target, Dox_SP_P=subset(tmp, subset=ID %in% c("cluster14_Dox") & orig_coef >0)$target,
               control_SP_N=subset(tmp, subset= ID %in% c("cluster14_control") & orig_coef <0)$target, Dox_SP_N=subset(tmp, subset=ID %in% c("cluster14_Dox") & orig_coef <0)$target)

m=make_comb_mat(genelist)
UpSet(m,set_order = c("control_SP_P", "Dox_SP_P", "control_SP_N", "Dox_SP_N"),
    top_annotation = upset_top_annotation(m, add_numbers = TRUE,annotation_name_rot = 0,annotation_name_side = "right"),
    right_annotation = upset_right_annotation(m,add_numbers = TRUE )
    )
###The category of overlap targets for GO analysis
ego <- enrichGO(gene = extract_comb(m, comb_name = "0100"), keyType = "SYMBOL", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
dotplot(ego, showCategory=5)+
  theme(panel.background = element_rect(colour = "black", size=1, fill = "white"), panel.grid = element_line(size=0), axis.line = element_line(size=0), axis.ticks = element_line(size = 1), axis.title = element_text(size = 15), legend.text = element_text(size=15, angle = 45), legend.key.size = unit(0.8,"cm"), legend.position = "bottom", legend.direction = "horizontal")