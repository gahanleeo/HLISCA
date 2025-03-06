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
library(ChIPseeker)


# load RDS 
setwd('/data/Choi_lung/ChiaHan/')

lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds')

#################################################################################
#################################################################################
# DAR
#################################################################################


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
saveRDS(peaks, file = "DAR_results/peaks_smoking.rds")
# quantify counts in each peak
peaks <- readRDS(file = "DAR_results/peaks_smoking.rds")

print('peak saved')

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

# saveRDS(lung_4, file = "DAR_results/peak_called_by_cell_types_smoking_CHL.rds")
# lung_4 <- readRDS("peak_called_by_cell_types_smoking.rds")


print('saved RDS')

# # run FindMarkers by parallel
# library(future)
# availableCores()
# plan("multisession", workers = 16)

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
  filename <- paste0('DAR_results/',i,".smoking.DAR.LR.csv")
  write.csv(marker_i,filename)
}

print('Done here')


########## Start 02212025 ##########

lung_4 <- readRDS('DAR_results/peak_called_by_cell_types_smoking_CHL.rds')

DefaultAssay(lung_4) <- "peaks"
Idents(lung_4) <- "Celltype_ss"


celltypes <- levels(lung_4$CellType)

DAR.list <- list()

for (i in celltypes){
  #i <- celltypes[1]
  ff <- read.csv(paste0('DAR_results/',i,".smoking.DAR.LR.csv"))
  DAR.list[[i]] <- ff
  rm(ff)
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

# start ChIPseeker function
peakAnno <- annotatePeak(peaks,
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

write.csv(peakAnno,"DAR_results/smoking_peaks_annotation_by_ChIPseeker_CHL.csv")
peaks_annotation_by_ChIPseeker <- read.csv("DAR_results/smoking_peaks_annotation_by_ChIPseeker_CHL.csv")
peaks_annotation_by_ChIPseeker <- peaks_annotation_by_ChIPseeker[,-1]
rownames(peaks_annotation_by_ChIPseeker) <- paste(peaks_annotation_by_ChIPseeker$seqnames,
                                                  peaks_annotation_by_ChIPseeker$start,
                                                  peaks_annotation_by_ChIPseeker$end, sep = "-")

peaks_annotation_by_ChIPseeker$X <- rownames(peaks_annotation_by_ChIPseeker)

DAR.list.new <-list()

celltype <- names(DAR.list.sig)

for (i in 1:23) {
  # i <- 3
  if (nrow(DAR.list.sig[[i]]) > 0) {
    markers <- DAR.list.sig[[i]]
    markers$cluster <- celltype[i]
    markers_annot <- left_join(markers, peaks_annotation_by_ChIPseeker,by = "X")
    DAR.list.new[[i]] <- markers_annot
  }else{
    DAR.list.new[[i]] <- NULL
  }
}

DAR.new <- Reduce(rbind, DAR.list.new)

write.table(DAR.new, file = "DAR_results/SData7_DAR_LR_smoking_w_CT.tsv",
            sep = "\t",quote = FALSE,row.names = F,col.names = T)

# t1 <- read.delim('DAR_results/SData7_DAR_LR_smoking_w_CT.tsv',header = T)


#########################################################################################################
#### CUT OFFF testing files
#########################################################################################################

darnew <- read.delim('DAR_results/SData7_DAR_LR_smoking_w_CT.tsv',header = T)
darnew <- darnew[,c(1,9:12)]


Original <- read.csv('/data/leec20/Jiyeon_lab_single_cell_quest/Original_DAR_list.csv',header = F)

# 174 same position 
intersect(darnew$X,Original$V1)


# Split into three columns

Original <- Original %>%
  separate(V1, into = c("chr", "start", "end"), sep = "-", convert = TRUE)


Original$width <- Original$end - Original$start

# something is wired: like chr2-208960769-20896097614 is way beyond end of chr2

Original <- Original[Original$width < 1000000, ]



# since it's kind of few base different, maybe use overlap

library(GenomicRanges)

# Compute widths for both datasets
darnew$width <- darnew$end - darnew$start

# Convert to GRanges
darnew_gr <- GRanges(seqnames = darnew$seqnames, 
                     ranges = IRanges(start = darnew$start, end = darnew$end))

Original_gr <- GRanges(seqnames = Original$chr, 
                      ranges = IRanges(start = Original$start, end = Original$end))

########################################
# findOverlaps(query, subject,...)
########################################

# Find close matches (≤500bp)
proximate_hits <- findOverlaps(darnew_gr, Original_gr, maxgap = 500, type = "any")

# Extract matched rows
matched_darnew <- darnew[queryHits(proximate_hits), ]
matched_Original <- Original[subjectHits(proximate_hits), ]

mtlist <- data.frame(proximate_hits)

########################################
########################################

# the uniq region, some of the region is repeated 
length(darnew$X %>% unique())
# The number that match around the Original data  
mtlist$queryHits %>% unique() %>% length()
# around 602 new region 
2448 - 1846


