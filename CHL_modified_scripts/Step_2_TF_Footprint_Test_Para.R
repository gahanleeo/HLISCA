
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

#lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds')

lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_Linkage1MB_TF_03142025.rds')
# 
# DefaultAssay(lung_4) <- "peaks"
# 
# # The TF list is from motifbreaker
# # run per SNPs
# ids <- colnames(Motifs(lung_4))
# motif_in_lung4 <- ConvertMotifID(object = lung_4, id = ids)
# 
# motif_B <- read.table('Peak_and_motifBreaker_results/motifBreakR_results_03032025.tsv',header = T)
# motif_list <- motif_B$geneSymbol %>% sort() %>% unique()
# 
# # find common TF list 
# 
# targe_vec <- intersect(motif_list , motif_in_lung4)
# 
# 
# # called 250g mem 
# 
# # set to 200GB
# options(future.globals.maxSize = 350 * 1024^3)
# 
# library(future)
# 
# plan("multicore", workers = 16)
# plan()
# 
# 
# print('start')
# 
# for (i in seq(1, length(targe_vec), by = 25)) {
#   end <- min(i + 24, length(targe_vec))  # Ensure the last group is correctly handled
#   
#   print(paste("Processed:", i, "to", end))
#   
#   lung_4 <- Footprint(
#     object = lung_4,
#     motif.name = targe_vec[i:end],
#     genome = BSgenome.Hsapiens.UCSC.hg38
#   )
#   
#   print(paste("Done Processed:", i, "to", end))
#   
# }
# 
# 
# 
# print('All loop fined')
# 
# saveRDS(lung_4, paste0("/data/leec20/Jiyeon_lab_single_cell_quest/TF_foot_RDS/","lung_","All_loop",".rds"))
# 
# print('saved')


###################### ######################
############### Plot the TFFootprint ########
###################### ######################

Idents(lung_4) <- "CellType"

ids <- colnames(Motifs(lung_4))
motif_in_lung4 <- ConvertMotifID(object = lung_4, id = ids)

motif_B <- read.table('Peak_and_motifBreaker_results/motifBreakR_results_03032025.tsv',header = T)
motif_list <- motif_B$geneSymbol %>% sort() %>% unique()

# find common TF list

targe_vec <- intersect(motif_list , motif_in_lung4)

library(patchwork)
library(ggplot2)


# Create a sequence to iterate in chunks of 4

motifs_per_plot <- 4

for (i in seq(1, length(targe_vec), by = 4)) {
  # Select up to 4 motifs at a time
  # i <- 1
  
  selected_motifs <- targe_vec[i:min(i + motifs_per_plot - 1, length(targe_vec))]
  
  # Generate footprint plots
  # label.top = 3 maybe can be use ?
  p <- PlotFootprint(lung_4, features = selected_motifs) 
  # Arrange them in a 2x2 grid
  p_combined <- p + patchwork::plot_layout(ncol = 1, nrow = 4)
  # Save each set of 4 plots
  ggsave(paste0("TF_Footprint_plot/TF_footprint_plot_", i, ".pdf"), p_combined, width = 10, height = 16)
  print('save one figure!')
}




# Generate footprint plots
p2 <- PlotFootprint(lung_4, features = c(""))
# Arrange the plots in a 2x2 grid
p2_combined <- p2 + patchwork::plot_layout(ncol = 1, nrow = 2)
# Save the plot
ggsave("footprint_CHL.pdf", p2_combined, width = 10, height = 10)



