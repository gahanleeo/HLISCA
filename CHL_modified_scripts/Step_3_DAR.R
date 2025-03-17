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
# lung_4 <- readRDS("DAR_results/peak_called_by_cell_types_smoking_CHL.rds")


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

dar.new <- read.delim('DAR_results/SData7_DAR_LR_smoking_w_CT.tsv',header = T)
dar.new <- dar.new[,c(1,9:12)]

Ori.list <- list.files('/data/Choi_lung/lbl/SHARE-seq/R/Seurat_DAR_LR',full.names = T)

# Original <- read.csv('/data/leec20/Jiyeon_lab_single_cell_quest/Original_DAR_list.csv',header = F)

Ori_f <- data.frame()

for (j in Ori.list){
  # j <- Ori.list[23]
  df <- read.csv(j)
  df <- df %>% filter(p_val_adj <= 0.05)
  if (nrow(df) == 0){
   print('ohoh')
  }else{
  df$Cell <- gsub('/data/Choi_lung/lbl/SHARE-seq/R/Seurat_DAR_LR/|.smoking.DAR.LR.adjseqdepth.csv','',j)
  Ori_f <- rbind(Ori_f,df)
  }
}



# 174 same position 
intersect(dar.new$X,Ori_f$X)


# Split into three columns

Ori_f <- Ori_f %>%
  separate(X, into = c("chr", "start", "end"), sep = "-", convert = TRUE)

# since it's kind of few base different, maybe use overlap

library(GenomicRanges)

# Convert to GRanges
dar.new_gr <- GRanges(seqnames = dar.new$seqnames, 
                     ranges = IRanges(start = dar.new$start, end = dar.new$end))

Ori_f_gr <- GRanges(seqnames = Ori_f$chr, 
                      ranges = IRanges(start = Ori_f$start, end = Ori_f$end))

########################################
# findOverlaps(query, subject,...)
########################################

# Find close matches (≤500bp)
proximate_hits <- findOverlaps(dar.new_gr, Ori_f_gr, maxgap = 500, type = "any")

mtlist <- data.frame(proximate_hits)
nrow(dar.new) - unique(mtlist$queryHits) %>% length()


# Do summary table

Ori_f <- data.frame()
for (j in Ori.list){
  # j <- Ori.list[23]
  df <- read.csv(j)
  df <- df %>% filter(p_val_adj <= 0.05)
  if (nrow(df) == 0){
    print('ohoh')
  }else{
    df$Cell <- gsub('/data/Choi_lung/lbl/SHARE-seq/R/Seurat_DAR_LR/|.smoking.DAR.LR.adjseqdepth.csv','',j)
    Ori_f <- rbind(Ori_f,df)
  }
}


# get the unique match X 
matched_X_dar_new <- unique(dar.new$X[mtlist$queryHits])
matched_X_Ori_f <- unique(Ori_f$X[mtlist$subjectHits])

summary_s1 <- data.frame(
  DataSet = "DAR",
  Chia_DAR_total_nb = nrow(dar.new),
  Chia_DAR_unique_total_nb = length(unique(dar.new$X)),
  EP_total_nb = nrow(Ori_f),
  EP_unique_total_nb = length(unique(Ori_f$X)),
  Chia_DAR_new_discovery = length(unique(dar.new$X)) - length(matched_X_dar_new)
  # EP_DAR_diff = length(unique(Ori_f$X)) - length(matched_X_Ori_f)
  
)

#############################################
#Compare the peaks results
#############################################

peaks_by_cell_Chl <- read.csv('Peak_and_motifBreaker_results/peaks_by_cell_types_CHL_chranno.csv')
peaks_by_cell_Chl$combine_chr <- paste0(peaks_by_cell_Chl$seqnames,"-",peaks_by_cell_Chl$start,"-",peaks_by_cell_Chl$end)
  
# vs EP
peaks_by_cell_EP <- read.csv("/vf/users/Choi_lung/SHARE-seq/Seurat_SHARE/peaks_analysis/peaks_by_cell_types_combined_22.csv")


peaks_by_cell_EP$chr <- paste0('chr',peaks_by_cell_EP$chr)
peaks_by_cell_EP$combine_chr <- paste0(peaks_by_cell_EP$chr,"-",peaks_by_cell_EP$start,"-",peaks_by_cell_EP$end)

library(GenomicRanges)

# Convert to GRanges

peaks_by_cell_Chl_gr <- GRanges(seqnames = peaks_by_cell_Chl$seqnames, 
                      ranges = IRanges(start = peaks_by_cell_Chl$start, end = peaks_by_cell_Chl$end))

peaks_by_cell_EP_gr <- GRanges(seqnames = peaks_by_cell_EP$chr, 
                    ranges = IRanges(start = peaks_by_cell_EP$start, end = peaks_by_cell_EP$end))


# Find close matches (≤500bp)
proximate_hits <- findOverlaps(peaks_by_cell_Chl_gr, peaks_by_cell_EP_gr, maxgap = 500, type = "any")

# Extract matched rows
# matched_dar.new <- peaks_by_cell_Chl_gr[queryHits(proximate_hits), ]
# matched_Ori_f <- peaks_by_cell_EP_gr[subjectHits(proximate_hits), ]

mtlist <- data.frame(proximate_hits)
# The key point is: the new finding on my analysis vs EP's 
nrow(peaks_by_cell_Chl) - unique(mtlist$queryHits) %>% length()


matched_X_peak <- unique(peaks_by_cell_Chl$combine_chr[mtlist$queryHits])


summary_s2 <- data.frame(
  DataSet = "Peak_by_celltype",
  Chia_Peak_total_nb = nrow(peaks_by_cell_Chl),
  Chia_Peak_unique_total_nb = length(unique(peaks_by_cell_Chl$combine_chr)),
  
  EP_total_nb = nrow(peaks_by_cell_EP),
  EP_unique_total_nb = length(unique(peaks_by_cell_EP$combine_chr)),
  Chia_DAR_new_discovery = length(unique(peaks_by_cell_Chl$combine_chr)) - length(matched_X_peak)

)

# 11408 new peaks on cell type in my result
# make table

##############################
##############################

peak_Smoking_CHL <- readRDS('DAR_results/peaks_smoking.rds')
# vs EP
peak_Smoking_EP <- readRDS('/data/Choi_lung/lbl/SHARE-seq/peak_called_by_cell_types_smoking.rds')
peak_Smoking_EP_pk <- peak_Smoking_EP[["peaks"]]@ranges

proximate_hits <- findOverlaps(peak_Smoking_CHL, peak_Smoking_EP_pk, maxgap = 500, type = "any")

mtlist <- data.frame(proximate_hits)

matched_X_peak <- unique(peak_Smoking_CHL[mtlist$queryHits])

summary_s3 <- data.frame(
  DataSet = "Peak_by_cell_With_SMOKING",
  Chia_Peak_total_nb = length(peak_Smoking_CHL),
  Chia_Peak_unique_total_nb = length(unique(peak_Smoking_CHL)),
  
  EP_total_nb = length(peak_Smoking_EP_pk),
  EP_unique_total_nb = length(unique(peak_Smoking_EP_pk)),
  
  Chia_DAR_new_discovery = length(unique(peak_Smoking_CHL)) - length(matched_X_peak)
  
)



317837 - unique(mtlist$queryHits) %>% length()
# 12889 differenct in cell type by smoking in my result

summary_tot <- data.frame()

names(summary_s1) <- c("DataSet","Chia_total","Chia_total_unique","EP_total","EP_total_unique","Chia_new_discovery")
names(summary_s2) <- c("DataSet","Chia_total","Chia_total_unique","EP_total","EP_total_unique","Chia_new_discovery")
names(summary_s3) <- c("DataSet","Chia_total","Chia_total_unique","EP_total","EP_total_unique","Chia_new_discovery")


summary_tot <- rbind(summary_s1,summary_s2,summary_s3)



