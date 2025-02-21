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
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")


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

  res <- lfcShrink(dds, 
                   contrast =  contrast,
                   res=res, type = 'ashr')
  
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  write.csv(res_tbl,
            paste0('DEG_results/',clusters[i], "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  
  rm(cluster_metadata)
  rm(counts)
  rm(cluster_counts)
  
  print(paste0('done',' ',clusters[i]))
}


######################### end of loop #########################
# 
# 
# # for this example, clusters[1] is Artery
# cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
# rownames(cluster_metadata) <- cluster_metadata$sample_id
# 
# counts <- pb[[clusters[1]]]
# 
# cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
# 
# # check
# all(rownames(cluster_metadata) == colnames(cluster_counts))
# 
# # # make factor
# # cluster_metadata$sample_id <- as.factor(cluster_metadata$sample_id)
# 
# # This will generate error... 
# # Since the group ID is highly correlated with smoking ID 
# 
# # dds <- DESeqDataSetFromMatrix(cluster_counts, 
# #                               colData = cluster_metadata, 
# #                               design = ~ group_id + sample_id)
# 
# dds <- DESeqDataSetFromMatrix(cluster_counts, 
#                               colData = cluster_metadata, 
#                               design = ~ group_id)
# 
# dds <- DESeq(dds)
# plotDispEsts(dds)
# 
# contrast <- c("group_id", levels(cluster_metadata$group_id)[2], levels(cluster_metadata$group_id)[1])
# 
# res <- results(dds, 
#                contrast = contrast,
#                alpha = 0.05)
# 
# res <- lfcShrink(dds, 
#                  contrast =  contrast,
#                  res=res, type = 'ashr')
# 
# res_tbl <- res %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>%
#   as_tibble()
# 
# write.csv(res_tbl,
#           paste0('DEG_results/',clusters[1], "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_all_genes.csv"),
#           quote = FALSE, 
#           row.names = FALSE)


#################################################################################
#################################################################################
# DAR
#################################################################################
##############################################################################
# Script information                                                      
# Title: DAR
# Author: Bolun Li
# Date: 2024-06-01
# Description: None
##############################################################################

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)
library(tidyverse)
library(grid)
library(GenomicRanges)
library(SeuratData)
library(ChIPseeker)

####ATAC analysis


# lung_4$CellType <- factor(lung_4$CellType, levels = c('AT1', 'AT2', 'AT1_AT2','AT2_pro',
#                                                   'Club','Ciliated','Goblet','Basal',
#                                                   'Artery','Vein','Capillary', "Lymphatic",
#                                                   'SMC','MyoFib','Fibroblast','Mesothelial',
#                                                   'Macrophage','Monocyte','Dendritic','NK','NK_T', 'B', 'T'))

DefaultAssay(lung_4) <- "ATAC"

#dividing each cell type into smoker and non_smoker clusters

Idents(lung_4) <- lung_4$Smoking

lung_4$Celltype_ss <- paste(lung_4$CellType, lung_4$Smoking, sep = "_")

table(lung_4$Celltype_ss)

###calling peaks using divided 46 clusters

# call peaks using MACS2

peaks <- CallPeaks(lung_4, group.by = "Celltype_ss",
                   macs2.path = "/data/leec20/conda/envs/MACS3/bin/macs3")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# save the peak as rds file
saveRDS(peaks, file = "peaks.rds")

# quantify counts in each peak
peaks <- readRDS(file = "peaks.rds")
macs2_counts <- FeatureMatrix(
  fragments = Fragments(lung_4),
  features = peaks,
  cells = colnames(lung_4)
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"


# create a new assay using the MACS2 peak set and add it to the Seurat object
lung_4[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  sep = c(":", "-"),
  fragments = Fragments(lung_4),
  annotation = annotations
)

DefaultAssay(lung_4) <- "peaks"

lung_4 <- FindTopFeatures(lung_4, min.cutoff = 5)
lung_4 <- RunTFIDF(lung_4)
lung_4 <- RunSVD(lung_4)

saveRDS(lung_4, file = "peak_called_by_cell_types_smoking.rds")
lung_4 <- readRDS("peak_called_by_cell_types_smoking.rds")


# run FindMarkers by parallel
library(future)
availableCores()
plan("multisession", workers = 16)

# Signac DAR
celltypes <- levels(lung_4$CellType)
DefaultAssay(lung_4)
table(Idents(lung_4))
Idents(lung_4) <- "Celltype_ss"
DAR.list <- list()
for (i in celltypes){
  cluster1 <- paste0(i,"_smoker")
  cluster2 <- paste0(i,"_non_smoker")
  marker_i <- FindMarkers(lung_4, ident.1 = cluster1, ident.2 = cluster2, 
                          test.use = 'LR',latent.vars = 'nCount_peaks')
  DAR.list[[i]] <- marker_i
  filename <- paste0(i,".smoking.DAR.LR.csv")
  write.csv(marker_i,filename)
}

# select significant DAR with adjusted p value < 0.05
DAR.list.sig <- lapply(DAR.list, function(x) {
  x$peak <- rownames(x)
  x <- subset(x, p_val_adj < 0.05 )
})


# add annotation info for DAR
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(org.Hs.eg.db)
peaks <- lung_4@assays$peaks@ranges
peakAnno <- annotatePeak(peaks,
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

write.csv(peakAnno,"project/AT2_co_CS/peaks_annotation_by_ChIPseeker.csv")
peaks_annotation_by_ChIPseeker <- read.csv("project/AT2_co_CS/peaks_annotation_by_ChIPseeker.csv")
peaks_annotation_by_ChIPseeker <- peaks_annotation_by_ChIPseeker[,-1]
rownames(peaks_annotation_by_ChIPseeker) <- paste(peaks_annotation_by_ChIPseeker$seqnames, 
                                                  peaks_annotation_by_ChIPseeker$start,
                                                  peaks_annotation_by_ChIPseeker$end, sep = "-")

DAR.list.new <-list()
celltype <- names(DAR.list.sig)
for (i in 1:23) {
  if (nrow(DAR.list.sig[[i]]) > 0) {
    markers <- DAR.list.sig[[i]]
    markers$cluster <- celltype[i]
    markers_annot <- cbind(markers, peaks_annotation_by_ChIPseeker[markers$peak,])
    DAR.list.new[[i]] <- markers_annot
  }else{
    DAR.list.new[[i]] <- NULL
  }
}

DAR.new <- Reduce(rbind, DAR.list.new)
write.table(DAR.new, file = "project/AT2_co_CS/DAR_LR_smoking_w_CT.txt",
            sep = "\t",quote = FALSE)



##############################################################################
##############Cell chat ##########
##############################################################################







