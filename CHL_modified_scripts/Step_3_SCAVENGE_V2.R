
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

library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(viridis)

setwd('/data/Choi_lung/ChiaHan/')

lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds')

#storing the matrix in a RangedSummarizedExperiment object

# SE_gvar <- SummarizedExperiment(assays = list(counts = peakbycellmat),
#                                 rowRanges = rowRanges(proj_PeakMatrix),
#                                 colData = DataFrame(names = colnames(peakbycellmat)))
# 



SE_gvar <- SummarizedExperiment(assays = list(counts = lung_4@assays[["peaks"]]@counts),
                                rowRanges = lung_4@assays[["peaks"]]@ranges,
                                colData = DataFrame(names = colnames(lung_4@assays[["peaks"]])))


assayNames(SE_gvar) <- "counts"

SE_gvar <- addGCBias(SE_gvar, genome = BSgenome.Hsapiens.UCSC.hg38)
SE_gvar_bg <- getBackgroundPeaks(SE_gvar, niterations=200)
#save(SE_gvar, SE_gvar_bg, file="NCI_M1_SE_gvar.rda")

#######################################################################################
# add embed
#######################################################################################
# Extract UMAP embeddings and convert to a data frame
umap_coords <- Embeddings(lung_4, reduction = "wnn.umap") 
umap_df <- data.frame(
  x = umap_coords[,1],  # UMAP1
  y = umap_coords[,2]   # UMAP2
)

# Extract CellType metadata
cell_types <- lung_4@meta.data$CellType

# Directly assign to colData(SE_gvar) since cell names match
colData(SE_gvar)$UMAP1 <- umap_df$x
colData(SE_gvar)$UMAP2 <- umap_df$y

names(cell_types) <- rownames(lung_4@meta.data)  # Assign cell barcodes to the vector
# Assign CellType to SE_gvar, ensuring correct barcode matching
colData(SE_gvar)$CellType <- cell_types[colnames(SE_gvar)]


#######################################################################################
#######################################################################################

#### SCAVENGE #####

library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)


trait_file <- "/data/leec20/Jiyeon_lab_single_cell_quest/Overall_PP_SHARE_Flat.bed"

###gchromVAR analysis

SE_lung5k <- addGCBias(SE_gvar, genome = BSgenome.Hsapiens.UCSC.hg38)
SE_lung5k_bg <- getBackgroundPeaks(SE_lung5k, niterations=200)

trait_import <- importBedScore(rowRanges(SE_lung5k), trait_file, colidx=5)

print('start run SE')

SE_lung5k_DEV <- computeWeightedDeviations(SE_lung5k, trait_import, background_peaks = SE_lung5k_bg)

print('done SE')

###Reformat results
z_score_mat <- data.frame(colData(SE_lung5k), z_score=t(assays(SE_lung5k_DEV)[["z"]]) %>% c())
head(z_score_mat)

### Generate the seed cell index (using the top 5% if too many cells are eligible)
seed_idx <- seedindex(z_score_mat$z_score, 0.05)

###calculate scale factor
scale_factor <- cal_scalefactor(z_score=z_score_mat$z_score, 0.01)

###Construct m-knn graph
peak_by_cell_mat <- assay(SE_lung5k)
tfidf_mat <- tfidf(bmat=peak_by_cell_mat, mat_binary=TRUE, TF=TRUE, log_TF=TRUE)

###Calculate lsi-mat
lsi_mat <- do_lsi(tfidf_mat, dims=30)

###Calculate m-knn graph
mutualknn30 <- getmutualknn(lsi_mat, 30)

###Network propagation
np_score <- randomWalk_sparse(intM=mutualknn30, rownames(mutualknn30)[seed_idx], gamma=0.05)

###Trait relevant score (TRS) with scaled and normalized
omit_idx <- np_score==0
mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]
TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
TRS <- TRS * scale_factor
mono_mat <- data.frame(z_score_mat[!omit_idx, ], seed_idx[!omit_idx], np_score, TRS)

write.csv(mono_mat,"Arrow_SCAVENGE_results/mono_mat_results_CHL.csv")

print('done CSV')



# mono_mat <- read.csv('Arrow_SCAVENGE_results/results_CHL.csv',row.names = 1)

#############################################
### Those umap plot not sure how to use...###
#############################################
# 
# ###UMAP plots of cell type annotation and cell-to-cell graph
# p <- ggplot(data=mono_mat, aes(UMAP1, UMAP2, color=CellType)) + geom_point(size=1, na.rm = TRUE) +
#   pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
# p
# 
# ###Comparsion before and after SCAVENGE analysis
# 
# ###scatter plot
# p1 <- ggplot(data=mono_mat, aes(UMAP1, UMAP2, color=z_score)) + geom_point(size=1, na.rm = TRUE, alpha = 0.6) +
#   scale_color_viridis_c() + scale_alpha()+
#   pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
# p1
# 
# ###bar plot
# pp1 <- ggplot(data=mono_mat,  aes(x=CellType, y=z_score))  +
#   geom_boxplot(aes(fill=CellType, color=CellType), outlier.shape=NA) +
#   guides(fill=FALSE) + pretty_plot(fontsize = 10) +
#   stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")+
#   ylim(c(0,5))
# 
# 
# pp1
# 
# pp2 <- ggplot(data=mono_mat,  aes(x=CellType, y=TRS))  +
#   geom_boxplot(aes(fill=CellType, color=CellType), outlier.shape=NA) + 
#   guides(fill=FALSE) + 
#   pretty_plot(fontsize = 10) +
#   stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")
# 
# pp2


############
############

print('start permu')
### Trait relevant cell determination from permutation test, reuqire a lot memory ~225G
mono_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30, seed_idx=mono_mat$seed_idx, topseed_npscore=mono_mat$np_score, permutation_times=1000, true_cell_significance=0.05, rda_output=F, mycores=8, rw_gamma=0.05)
mono_mat2 <- data.frame(mono_mat, mono_permu)
print('done mono_mat2')

write.csv(mono_mat2,"Arrow_SCAVENGE_results/mono_mat2_results_CHL.csv")


### Look at the distribution of statistically significant phenotypically enriched and depleted cells
##Enriched cells
# 
# mono_mat2 |>
#   dplyr::group_by(CellType) |> 
#   dplyr::summarise(enriched_cell=sum(true_cell_top_idx)) |> 
#   ggplot(aes(x=CellType, y=enriched_cell, fill=CellType)) + 
#   geom_bar(stat="identity") + 
#   theme_classic()
# 
# 
# 
# ##Depleted cells
# 
# mono_mat2$rev_true_cell_top_idx <- !mono_mat2$true_cell_top_idx
# 
# mono_mat2 |>
#   dplyr::group_by(CellType) |> 
#   dplyr::summarise(depleted_cell=sum(rev_true_cell_top_idx)) |> 
#   ggplot(aes(x=CellType, y=depleted_cell, fill=CellType)) + 
#   geom_bar(stat="identity") + 
#   theme_classic()


