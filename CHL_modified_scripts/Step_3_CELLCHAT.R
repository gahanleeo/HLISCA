##############################################################################
# Script information                                                      
# Title: Step_3_continue
# Author: Erping Long
# Date: 2022-02-10
# Description: None
##############################################################################
# Date: 02182025
# ChiaHan's run

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)
library(GenomeInfoDb)
library(tidyverse)
library(grid)
library(readr)
library(harmony)
library(GenomicRanges)
library(SeuratData)


# load RDS 
setwd('/data/Choi_lung/ChiaHan/')

# lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds')
 
# Idents(lung_4) <- "orig.ident"
# 
# lung_4 = RenameIdents(lung_4, 'M1' = 'MN1', 'M3' = 'MS1', 'M4' = 'MS2', 'M5' = 'MS3', 'M6' = 'MS4',
#                           'M7' = 'MN2', 'M8' = 'MN3', 'M9' = 'MN4', 'N17_20' = 'FN4', 'NCI_115' = 'FN3', 
#                           'NCI_118' = 'FS2', 'NCI_13' = 'FN1', 'NCI_28' = 'FN2', 'NCI_55' = 'FS1', 
#                           'S1' = 'FS3', 'S2' = 'FS4')
# 
# lung_4$SampleID <- Idents(lung_4)
# # add sex and smoking status
# 
# Idents(lung_4) <- "orig.ident"
# 
# lung_4 = RenameIdents(lung_4, 'M1' = 'non_smoker', 'M3' = 'smoker', 'M4' = 'smoker', 'M5' = 'smoker', 'M6' = 'smoker',
#                           'M7' = 'non_smoker', 'M8' = 'non_smoker', 'M9' = 'non_smoker', 'N17_20' = 'non_smoker', 'NCI_115' = 'non_smoker', 
#                           'NCI_118' = 'smoker', 'NCI_13' = 'non_smoker', 'NCI_28' = 'non_smoker', 'NCI_55' = 'smoker', 
#                           'S1' = 'smoker', 'S2' = 'smoker')
# 
# lung_4$Smoking <- Idents(lung_4)
# 
# Idents(lung_4) <- "orig.ident"
# 
# lung_4 = RenameIdents(lung_4, 'M1' = 'Male', 'M3' = 'Male', 'M4' = 'Male', 'M5' = 'Male', 'M6' = 'Male',
#                           'M7' = 'Male', 'M8' = 'Male', 'M9' = 'Male', 'N17_20' = 'Female', 'NCI_115' = 'Female', 
#                           'NCI_118' = 'Female', 'NCI_13' = 'Female', 'NCI_28' = 'Female', 'NCI_55' = 'Female', 
#                           'S1' = 'Female', 'S2' = 'Female')
# 
# lung_4$Gender <- Idents(lung_4)
# 
# saveRDS(lung_4,'CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds')


##############################################################################
# Description: CellCHAT
##############################################################################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(celldex)
library(cowplot)
library(tidyverse)
library(CellChat)


# 
# split_list <- SplitObject(lung_4, split.by = "SampleID")
# 
# list_ns <- split_list[c("MN1", "MN2", "MN3", "MN4", "FN1", "FN2", "FN3", "FN4")]
# list_s <- split_list[c("MS1", "MS2", "MS3", "MS4", "FS1", "FS2", "FS3", "FS4")]
# 
# lung_ns <- merge(list_ns$MN1, 
#                  y = c(list_ns$MN2, list_ns$MN3, list_ns$MN4, list_ns$FN1, list_ns$FN2, list_ns$FN3, list_ns$FN4))
# 
# lung_s <- merge(list_s$MS1, 
#                 y = c(list_s$MS2, list_s$MS3, list_s$MS4, list_s$FS1, list_s$FS2, list_s$FS3, list_s$FS4))
# 
# saveRDS(lung_ns,"/data/leec20/Jiyeon_lab_single_cell_quest/Scripts/Tmp_RDS/lung_ns.rds")
# saveRDS(lung_s,"/data/leec20/Jiyeon_lab_single_cell_quest/Scripts/Tmp_RDS/lung_S.rds")
# 
# print('saved RDS')
lung_ns <- readRDS('/data/leec20/Jiyeon_lab_single_cell_quest/Scripts/Tmp_RDS/lung_ns.rds') 
lung_s <- readRDS('/data/leec20/Jiyeon_lab_single_cell_quest/Scripts/Tmp_RDS/lung_S.rds') 


#create cellchat object

DefaultAssay(lung_ns) <- "RNA"
DefaultAssay(lung_s) <- "RNA"


cellchat_ns <- createCellChat(lung_ns, meta = lung_ns@meta.data, group.by = "CellType")
cellchat_s <- createCellChat(lung_s, meta = lung_s@meta.data, group.by = "CellType")

cellchat_ns@DB <- CellChatDB.human
cellchat_s@DB <- CellChatDB.human

#calculate the communication probability

cellchat_ns <- subsetData(cellchat_ns)
cellchat_ns <- identifyOverExpressedGenes(cellchat_ns)
cellchat_ns <- identifyOverExpressedInteractions(cellchat_ns)
cellchat_ns <- computeCommunProb(cellchat_ns)
cellchat_ns <- computeCommunProbPathway(cellchat_ns)
cellchat_ns <- aggregateNet(cellchat_ns)

cellchat_s <- subsetData(cellchat_s)
cellchat_s <- identifyOverExpressedGenes(cellchat_s)
cellchat_s <- identifyOverExpressedInteractions(cellchat_s)
cellchat_s <- computeCommunProb(cellchat_s)
cellchat_s <- computeCommunProbPathway(cellchat_s)
cellchat_s <- aggregateNet(cellchat_s)

cellchat_ns <- netAnalysis_computeCentrality(object = cellchat_ns, slot.name = "netP")
cellchat_s <- netAnalysis_computeCentrality(object = cellchat_s, slot.name = "netP")

object.list <- list(ns = cellchat_ns, s = cellchat_s)

cellchat_com <- mergeCellChat(object.list, add.names = names(object.list))

#save the probability

ns.netp <- subsetCommunication(cellchat_ns, slot.name = 'netP')
write.csv(ns.netp, "Cellchat_result/pathway_nonsmoke_CHL.csv")

s.netp <- subsetCommunication(cellchat_s, slot.name = 'netP')
write.csv(s.netp, "Cellchat_result/pathway_smoke_CHL.csv")

ns.net <- subsetCommunication(cellchat_ns, slot.name = 'net')
write.csv(ns.net, "Cellchat_result/L_R_nonsmoke_CHL.csv")

s.net <- subsetCommunication(cellchat_s, slot.name = 'net')
write.csv(s.net, "Cellchat_result/L_R_smoke_CHL.csv")

### DONE ####
