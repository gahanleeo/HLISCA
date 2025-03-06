
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

# ###############################################################
# # TF frootprint
# ###############################################################

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

setwd('/data/Choi_lung/ChiaHan/')

lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds')

DefaultAssay(lung_4) <- "peaks"


pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

lung_4 <- AddMotifs(lung_4, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)


# The TF list is from motifbreaker
# run per SNPs
ids <- colnames(Motifs(lung_4))
motif_in_lung4 <- ConvertMotifID(object = lung_4, id = ids)

motif_B <- read.table('Peak_and_motifBreaker_results/motifBreakR_results_03032025.tsv',header = T)
motif_list <- motif_B$geneSymbol %>% sort() %>% unique()

# find common TF list 

targe_vec <- intersect(motif_list , motif_in_lung4)


for (i in seq(1, length(targe_vec), by = 25)) {
  end <- min(i + 24, length(targe_vec))  # Ensure the last group is correctly handled
  
  print(paste("Processed:", i, "to", end))
  
  lung_4 <- Footprint(
    object = lung_4,
    motif.name = targe_vec[i:end],
    genome = BSgenome.Hsapiens.UCSC.hg38
  )

  print(paste("Processed:", i, "to", end))
  saveRDS(lung_4, paste0("/data/leec20/Jiyeon_lab_single_cell_quest/TF_foot_RDS/","lung_",end,".rds"))
  print('done save')
  
}

# Idents(lung_4) <- "CellType"

# 
# 
# p2 <- PlotFootprint(lung_4, features = c("BHLHE40"))
# p2 + patchwork::plot_layout(ncol = 1)
# #
# ggsave("footprint_CHL.pdf",
#        p2+ patchwork::plot_layout(ncol = 1), width = 10, height = 30)

lung_4 <- readRDS('/data/leec20/Jiyeon_lab_single_cell_quest/TF_foot_RDS/lung_50.rds')


