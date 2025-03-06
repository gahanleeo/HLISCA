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

lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds')
# 
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


###############################################################
## DEG
###############################################################

library(DESeq2)
library(SingleCellExperiment)

library(Matrix.utils)

library(magrittr)
library(purrr)
library(scater)
library(cowplot)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(ashr)
library(MAST)


#generate pseudobulk data
##create single cell experiment object

counts <- lung_4@assays$RNA@counts
metadata <- lung_4@meta.data

metadata$cluster_id <- factor(lung_4$CellType)

sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

sce$sample_id <- as.factor(sce$SampleID)

groups <- colData(sce)[, c("cluster_id", "sample_id")]


##generate sample level metadata

kids <- purrr::set_names(levels(sce$cluster_id))
nk <- length(kids)
sids <- purrr::set_names(levels(sce$sample_id))
ns <- length(sids)

n_cells <- as.numeric(table(sce$sample_id))

m <- match(sids, sce$sample_id)
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% dplyr::select(-"cluster_id")




##aggregate the counts per sample_id and cluster_id to generate psudobulk data 

groups <- colData(sce)[, c("cluster_id", "sample_id")]


pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 


# THIS ONE IS PROBLEMIC CAN NOT DIFFERENT NK or NK_T, if not changing the sampleId

splitf <- gsub('_MN[1-9]+|_MS[1-9]+|_FN[1-9]+|_FS[1-9]+','',rownames(pb))

#splitf <- gsub('_M[1-9]+|_NCI.+|_N17_20|_S[1-2]','',rownames(pb))

# splitf <- sapply(stringr::str_split(rownames(pb),
#                                     pattern = "_",
#                                     n = 2), `[`, 1)
# 

# pb <- split.data.frame(pb,
#  factor(splitf)) %>%
#   lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "[MNS].*")))

# 
# pb <- split.data.frame(pb, factor(splitf)) %>% 
#   lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "[MS]\\d+|NCI_\\d+|N17_20")))
# 

pb <- split.data.frame(pb, factor(splitf)) %>% 
  lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "MN[1-9]+|MS[1-9]+|FN[1-9]+|FS[1-9]+")))



class(pb)
str(pb)

options(width = 100)

table(sce$cluster_id, sce$sample_id)

##differential expression analysis

get_sample_ids <- function(x){
  pb[[x]] %>%colnames()}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()


gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

ei$group_id <- factor(c("N", "S", "S", "S", "S", "N", "N", "N", "N", "N", "S", "N", "N", "S", "S", "S"))

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")]) 

metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id) 


# list of cell type 
clusters <- levels(as.factor(metadata$cluster_id))

###################################### Begin loop ###########################################

for( i in c(1:23)){
  # i <- 2
  # for this example, clusters[1] is Artery
  cluster_metadata <- metadata[which(metadata$cluster_id == clusters[i]), ]
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  counts <- pb[[clusters[i]]]
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  # check
  all(rownames(cluster_metadata) == colnames(cluster_counts))
  dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata, 
                                design = ~ group_id)
  
  dds <- DESeq(dds)
  plotDispEsts(dds)
  
  contrast <- c("group_id", levels(cluster_metadata$group_id)[2], levels(cluster_metadata$group_id)[1])

  res <- results(dds,
                 contrast = contrast,
                 alpha = 0.05)
  
  
  # maybe not for now, or try different "normal" or "apeglm"
  # res <- lfcShrink(dds, 
  #                  contrast =  contrast,
  #                  res=res, type = 'ashr')
  
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  write.csv(res_tbl,
            paste0('DEG_results/DESeq2_default_result/',clusters[i], "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  
  rm(cluster_metadata)
  rm(counts)
  rm(cluster_counts)
  
  print(paste0('done',' ',clusters[i]))
}


######################### end of loop #########################


################################New way to aggreate #######################################################
##############################################################################################################################

####Find differential expression genes

Pesudo_lung_4 <- AggregateExpression(lung_4, assays = "RNA", return.seurat = T, group.by = c("Smoking", "SampleID", "CellType"))

tail(Cells(Pesudo_lung_4))

Pesudo_lung_4$celltype.stim <- paste(Pesudo_lung_4$CellType, Pesudo_lung_4$Smoking, sep = "_")

Idents(Pesudo_lung_4) <- "celltype.stim"


celllist <- Pesudo_lung_4$CellType %>% unique() %>% as.character()


for (i in 1:23){
  # i <- 3
  bulk.de <- FindMarkers(object = Pesudo_lung_4, 
                              ident.1 = paste0(celllist[i],"_smoker"), 
                              ident.2 = paste0(celllist[i],"_non-smoker"),
                              test.use = "DESeq2")
  
  filename <- paste0("DEG_results/Pseudobulking/",celllist[i],"_smoking.DEG.DESEQ2.csv")
  write.csv(bulk.de,filename)
  
  rm(bulk.de)
}

########################################################################################
########CUT OFF, testing rsults files###################################################
########################################################################################

l1 <- list.files('DEG_results/DESeq2_default_result',full.names = T)
lp <- list.files('DEG_results/Pseudobulking_result',full.names = T)


set1 <- data.frame()
setp <- data.frame()

for (i in l1){
  #i <- l1[1]
  dd <- read.csv(i)
  dd$Cell <- paste0(gsub('DEG_results/DESeq2_default_result/|_S_vs_N_all_genes.csv','',i))
  set1 <- rbind(set1,dd)
  
}

set1 <- set1 %>% filter(padj <= 0.05)



for (i in lp){
  #i <- lp[1]
  ddd <- read.csv(i)
  ddd$Cell <- paste0(gsub('DEG_results/Pseudobulking_result/|_smoking.DEG.DESEQ2.csv','',i))
  setp <- rbind(setp,ddd)
  
}

setp <- setp %>% filter(p_val_adj <= 0.05)

