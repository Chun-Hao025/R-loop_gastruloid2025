library(CellChat)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

##cell cell commnunication analysis
##exclude Placodal_area
RH=readRDS("/path/to/RHKI_new.rds")
RH=subset(RH, subset= ann_new_noZX_KNN40 %in% c('Placodal_area'), invert=T) ## exclude Placodal area and resistant type
RH$ann_new_noZX_F1=factor(RH$ann_new_noZX_KNN40, levels=c('Naive_Pluripotent_Cells', 'Epiblasts', 'Primitive_Streak', 'Caudal_Neuroectoderm1','Neuromesodermal_Progenitors1','Neuromesodermal_Progenitors2','Neuromesodermal_Progenitors3',
                                                'Spinal_Cord', 'Neuron-like_Cells','Nascent_Mesoderm','Paraxial_Mesoderm_A','Paraxial_Mesoderm_B',
                                                'Splanchnic_Mesoderm','Endothelium','Definitive_Endoderm_Guts'))
level=c('Naive Pluripotent Cells', 'Epiblasts', 'Primitive Streak', 'Caudal Neuroectoderm 1','Caudal Neuroectoderm 2','Neuromesodermal Progenitors 1','Neuromesodermal Progenitors 2',
        'Spinal Cord', 'Neuron-like Cells','Nascent Mesoderm','Paraxial Mesoderm A1','Paraxial Mesoderm A2','Paraxial Mesoderm A3','Paraxial Mesoderm B',
        'Splanchnic Mesoderm','Endothelium','Definitive Endoderm/Guts')
RH$treatment=factor(RH$treatment, levels=c("control", "Dox"))
RH=SplitObject(RH, split.by = "treatment")
RH=RH[c(2,1)]

##correct the miss annotation in CellChatDB 
for (i in 1:nrow(CellChatDB$interaction)){
  if (CellChatDB$interaction$ligand[i] == "RA-ALDH1A2"){
    CellChatDB$interaction$interaction_name_2[i] = gsub(pattern = "^RA-ALDH1A1", replacement = "RA-ALDH1A2", x = CellChatDB$interaction$interaction_name_2[i])
  }else if (CellChatDB$interaction$ligand[i] == "RA-ALDH1A3") {
    CellChatDB$interaction$interaction_name_2[i] = gsub(pattern = "^RA-ALDH1A1", replacement = "RA-ALDH1A3", x = CellChatDB$interaction$interaction_name_2[i])
  }else{
    CellChatDB$interaction$interaction_name_2[i] = CellChatDB$interaction$interaction_name_2[i]
  }
}
cellchat_RH=list()
for (i in 1:length(RH)){
  data.input = RH[[i]]@assays$RNA@data# normalized data matrix
  meta = RH[[i]]@meta.data # a dataframe with rownames containing cell mata data
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "ann_new_noZX_F1")
  cellchat <- setIdent(cellchat, ident.use = "ann_new_noZX_F1") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents))
  
  ###CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat, features=NULL) # This step is necessary even if using the whole database
  future::plan("multisession", workers = 4)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  #droplevel(cellchat@idents), remove the unused idents
  future::plan("multisession", workers = 4) # do parallel
  options(future.globals.maxSize = 3 * 1024^3)
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use=FALSE, population.size=TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat_RH[[i]]=cellchat
}
names(cellchat_RH)=c("control", "Dox")
cellchat <- mergeCellChat(cellchat_RH, add.names = names(cellchat_RH))

##RA signaling
# 1) extract the communication‐probability matrix for RA
commProb <- subsetCommunication(
  cellchat_RH[[1]],
  signaling    = "RA",
  slot.name    = "netP"       # this returns the aggregated probability matrix
)
commProb1 <- subsetCommunication(
  cellchat_RH[[2]],
  signaling    = "RA",
  slot.name    = "netP"       # this returns the aggregated probability matrix
)
# commProb is a square matrix with rownames=source cell types
# and colnames=target cell types

# 2) subset to your cells of interest
commSub <- subset(commProb, subset = source %in% c('Paraxial_Mesoderm_A','Paraxial_Mesoderm_B','Splanchnic_Mesoderm') & target %in% c('Neuromesodermal_Progenitors1','Neuromesodermal_Progenitors2', 'Spinal_Cord', 'Paraxial_Mesoderm_A','Paraxial_Mesoderm_B','Splanchnic_Mesoderm'))
commSub1 <- subset(commProb1, subset = source %in% c('Paraxial_Mesoderm_A','Paraxial_Mesoderm_B','Splanchnic_Mesoderm') & target %in% c('Neuromesodermal_Progenitors1','Neuromesodermal_Progenitors2', 'Spinal_Cord', 'Paraxial_Mesoderm_A','Paraxial_Mesoderm_B','Splanchnic_Mesoderm'))

# 2) pivot to wide
mat_df <- dcast(
  commSub,
  source ~ target,
  value.var = "prob",
  fun.aggregate = sum,     # there should be only one entry per pair
  fill = 0
)
mat_df1 <- dcast(
  commSub1,
  source ~ target,
  value.var = "prob",
  fun.aggregate = sum,     # there should be only one entry per pair
  fill = 0
)

# 3) set rownames and convert to numeric matrix
rownames(mat_df) <- mat_df$source
mat_df$source   <- NULL
mat_df[3,]=c(0,0,0,0,0,0)
rownames(mat_df) = rownames(mat_df1)
rownames(mat_df1) <- mat_df1$source
mat_df1$source   <- NULL

mat <- as.matrix(mat_df[,c(1,2,5)])
mat1 <- as.matrix(mat_df[,c(3,4,6)])
mat2 <- as.matrix(mat_df1[,c(1,2,5)])
mat3 <- as.matrix(mat_df1[,c(3,4,6)])
# 4) draw with ComplexHeatmap
ht1=Heatmap(mat, name = "RA prob",col= colorRamp2(c(0, max(mat)), c("white","firebrick")),
  cluster_rows= FALSE,cluster_columns = FALSE,column_title= "ectoderm", row_title = "control",row_title_side = 'left',column_title_side = 'top',
  show_row_names= TRUE,show_column_names = TRUE,row_names_side = "left", border = T)
ht2=Heatmap(mat1,name = "RA prob",col = colorRamp2(c(0, max(mat)), c("white","firebrick")),
  cluster_rows= FALSE,cluster_columns = FALSE,column_title = "mesoderm",column_title_side = 'top',
  show_row_names = TRUE, show_column_names = TRUE,row_names_side = "left", border = T)
ht3=Heatmap(mat2,name = "RA prob",col = colorRamp2(c(0, max(mat)), c("white","firebrick")),
            cluster_rows= FALSE,cluster_columns = FALSE,row_title = "control",row_title_side = 'left',
            show_row_names = TRUE, show_column_names = TRUE,row_names_side = "left", border = T)
ht4=Heatmap(mat3,name = "RA prob",col = colorRamp2(c(0, max(mat)), c("white","firebrick")),
            cluster_rows= FALSE,cluster_columns = FALSE,
            show_row_names = TRUE, show_column_names = TRUE,row_names_side = "left", border = T)
ht1 + ht2 + ht3 + ht4  # this builds a HeatmapList

#specific ligand-receptor pair in RA signaling
pairLR=data_frame(interaction_name=c("RetinoicAcid-RA-ALDH1A2_RARB_CRABP2", "RetinoicAcid-RA-ALDH1A2_RXRA_CRABP2", "RetinoicAcid-RA-ALDH1A2_RXRB_CRABP2"))
p1=netVisual_bubble(cellchat, sources.use = c(11:13), targets.use = c(5,6,8,11:13),  comparison = c(1, 2), angle.x = 45,
                 remove.isolate = F, dot.size.min = 5, grid.on = F, pairLR.use = pairLR, sort.by.target = TRUE)+scale_color_viridis_c(option="A")


###other signaling

commProb2 <- subsetCommunication(
  cellchat_RH[[1]],
  signaling    = c("WNT", "ncWNT", "TGFb", "BMP", "NOTCH", "FGF"),
  slot.name    = "netP"       # this returns the aggregated probability matrix
)
commProb3 <- subsetCommunication(
  cellchat_RH[[2]],
  signaling    = c("WNT", "ncWNT", "TGFb", "BMP", "NOTCH", "FGF"),
  slot.name    = "netP"       # this returns the aggregated probability matrix
)
pathway=c("WNT", "ncWNT", "TGFb", "NOTCH", "FGF")
plot=list()
my_row_order = c('Neuromesodermal_Progenitors1','Neuromesodermal_Progenitors2', 'Spinal_Cord','Paraxial_Mesoderm_A','Paraxial_Mesoderm_B','Splanchnic_Mesoderm')
for (i in 1:length(pathway)){
  commProb2 <- subsetCommunication(cellchat_RH[[1]], signaling = pathway[i],slot.name = "netP")
  commSub2 <- subset(commProb2, subset = source %in% c('Neuromesodermal_Progenitors1','Neuromesodermal_Progenitors2', 'Spinal_Cord','Paraxial_Mesoderm_A','Paraxial_Mesoderm_B','Splanchnic_Mesoderm') & target %in% c('Neuromesodermal_Progenitors1','Neuromesodermal_Progenitors2', 'Spinal_Cord', 'Paraxial_Mesoderm_A','Paraxial_Mesoderm_B','Splanchnic_Mesoderm'))
  mat_df2 <- dcast(commSub2, source ~ target,value.var = "prob",fun.aggregate = sum, fill = 0)
  rownames(mat_df2) <- mat_df2$source
  mat_df2$source   <- NULL
  missing_rows <- setdiff(my_row_order, rownames(mat_df2))
  if (length(missing_rows)>0){
    zero_df <- as.data.frame(
      matrix(0,nrow = length(missing_rows),ncol = ncol(mat_df2),dimnames = list(missing_rows, colnames(mat_df2)))
    )
    mat_df2 <- rbind(mat_df2, zero_df)
  }
  mat_df2 <- mat_df2[my_row_order, , drop = FALSE]
  mat_df2 <- mat_df2[, my_row_order, drop = FALSE]
  mat4 <- as.matrix(mat_df2)
  plot[[i]]=Heatmap(mat4, name = pathway[i], col = colorRamp2(c(0, max(mat)), c("white","firebrick")),
                    row_order = my_row_order,column_order = my_row_order,
          cluster_rows= FALSE,cluster_columns = FALSE,row_title = "control",row_title_side = 'left',
          show_row_names = TRUE, show_column_names = TRUE,row_names_side = "left", border = T)
}
for (i in 1:length(pathway)){
  j=i+5
  commProb2 <- subsetCommunication(cellchat_RH[[2]], signaling = pathway[i],slot.name = "netP")
  commSub2 <- subset(commProb2, subset = source %in% c('Neuromesodermal_Progenitors1','Neuromesodermal_Progenitors2', 'Spinal_Cord','Paraxial_Mesoderm_A','Paraxial_Mesoderm_B','Splanchnic_Mesoderm') & target %in% c('Neuromesodermal_Progenitors1','Neuromesodermal_Progenitors2', 'Spinal_Cord', 'Paraxial_Mesoderm_A','Paraxial_Mesoderm_B','Splanchnic_Mesoderm'))
  mat_df2 <- dcast(commSub2, source ~ target,value.var = "prob",fun.aggregate = sum, fill = 0)
  rownames(mat_df2) <- mat_df2$source
  mat_df2$source   <- NULL
  missing_rows <- setdiff(my_row_order, rownames(mat_df2))
  if (length(missing_rows)>0){
    zero_df <- as.data.frame(
      matrix(0,nrow = length(missing_rows),ncol = ncol(mat_df2),dimnames = list(missing_rows, colnames(mat_df2)))
    )
    mat_df2 <- rbind(mat_df2, zero_df)
  }
  mat_df2 <- mat_df2[my_row_order, , drop = FALSE]
  mat_df2 <- mat_df2[, my_row_order, drop = FALSE]
  mat4 <- as.matrix(mat_df2)
  plot[[j]]=Heatmap(mat4, name = pathway[i], col = colorRamp2(c(0, max(mat)), c("white","firebrick")),
                    row_order = my_row_order,column_order = my_row_order,
                    cluster_rows= FALSE,cluster_columns = FALSE,row_title = "control",row_title_side = 'left',
                    show_row_names = TRUE, show_column_names = TRUE,row_names_side = "left", border = T)
}
ht_list <- do.call(c, plot)
ht_all <- Reduce(`+`, plot)
draw(ht_all,
     ncol   = 5,                  # number of columns
     ht_gap = unit(0.5, "cm"))    # gap between heatmaps
draw(plot[[1]], plot[[2]],plot[[3]],plot[[4]],plot[[5]])
