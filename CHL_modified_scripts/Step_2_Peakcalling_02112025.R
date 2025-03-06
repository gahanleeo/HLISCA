##############################################################################
# Script information                                                      
# Title: Peak calling
# Author: Erping Long
# Date: 2022-02-10
# Description: None
##############################################################################
# Date: 02102025
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


#RDS file: /data/Choi_lung/lbl/SHARE-seq/correction/pbmc_integrate_all_removal_combined_res1_annot.rds

####ATAC analysis
setwd('/data/Choi_lung/ChiaHan/')

# lung_4 <- readRDS(file = "/data/Choi_lung/lbl/SHARE-seq/correction/pbmc_integrate_all_removal_combined_res1_annot.rds")

# check plot
# DimPlot(lung_4,reduction = 'wnn.umap',raster = F)

DefaultAssay(lung_4) <- "ATAC"

# call peaks using MACS2

peaks <- CallPeaks(lung_4, group.by = "CellType",macs2.path = '/data/leec20/conda/envs/MACS3/bin/macs3')

print('fin peak call')

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

write.csv(peaks,"./peaks_by_cell_types_CHL_chranno.csv")
print('done save csv')

# SAVE now so that it will be benficiant downstream
saveRDS(peaks,'./peaks_CHL_chranno.rds')
print('done save peak rds')

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(lung_4),
  features = peaks,
  cells = colnames(lung_4)
)
# 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

##################### Anno chrmomsome problem #####################

# https://github.com/stuart-lab/signac/issues/870

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


# saveRDS(lung_4, "./CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno.rds")

# print('done save full rds')

# if finished, run Bolun's script "Seurat_ATAC.R" start from 91 to 117. 
# the line 93 file is in my folder in Choi_lung
# the old_rst is the old result of CCV, need to compare to new result and run motifbreaker for new snps
# after peak calling, do Peak-linkage (https://github.com/pumclyy/HLISCA/blob/main/11.cCREs_analysis/Linkage.R) and ChIPseeker



##############################################################################
##############################################################################
####### THIS IS LINE 91 to 117 from Seurat_ATAC.R ##########
##############################################################################
##############################################################################

# # load files 
lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno.rds')
peaks <- readRDS('peaks_CHL_chranno.rds')

Idents(lung_4) <- "CellType"
# 
lung_GWAS <- read.table("Lung_GWAS_51loci.txt",header = TRUE, sep = "\t")

dim(lung_GWAS)
head(peaks@ranges)

# loaded the peaks result of mine
peak_info <- read.csv("Peak_and_motifBreaker_results/peaks_by_cell_types_CHL_chranno.csv", header = TRUE, row.names = 1)

peak_info$peak <- paste(peak_info$seqnames, peak_info$start, peak_info$end, sep = "-")
lung_GWAS$POS <- as.numeric(str_split_fixed(lung_GWAS$Variant_ID, pattern = ":", n = 2)[,2])
lung_GWAS$POS[1482] = 118255679
lung_GWAS$POS[1562] = 118255679
lung_GWAS$colocalized_cCRE <- NA
for(i in 1:nrow(lung_GWAS)){
  cCRE <- peak_info$peak[which(peak_info$seqnames == lung_GWAS$CHR[i] &
                                 as.numeric(peak_info$start) <= as.numeric(lung_GWAS$POS[i]) &
                                 as.numeric(peak_info$end) >= as.numeric(lung_GWAS$POS[i]))]
  if (length(cCRE) == 1) {
    lung_GWAS$colocalized_cCRE[i] <- cCRE
  }
}

GWAS_subset <- lung_GWAS[which(!is.na(lung_GWAS$colocalized_cCRE)),]

old_rst <- read.csv("/data/Choi_lung/SHARE-seq/Seurat_SHARE/peaks_analysis/test_SNP/merge.csv", header = FALSE)
old_ccv <- str_split_fixed(old_rst$V1, pattern = "\\.", n = 2)[,1]
length(intersect(old_ccv, GWAS_subset$RSNUM))
setdiff(GWAS_subset$RSNUM , old_ccv)

tmp <- lung_GWAS[which(lung_GWAS$RSNUM %in% setdiff(GWAS_subset$RSNUM , old_ccv)),]

# ##############################################################################
# ##############################################################################
# ##############################################################################
# ##############################################################################
# 
# ##############################################################################
# # Since we have new SNPs, the new SNPs need to do motifbreaker
# ##############################################################################
library(motifbreakR)
library(readr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(stringr)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(MotifDb)
library(ggplot2)

#snp <- read.csv("snp_list.csv")
snp <- tmp$RSNUM

snps.fromlist <- snps.from.rsid(rsid = snp,dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38,search.genome = BSgenome.Hsapiens.UCSC.hg38)

data("hocomoco")

results <- motifbreakR(snpList = snps.fromlist, filterp = TRUE,
                       pwmList = hocomoco,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())

results <- calculatePvalue(results)

# 1. Convert GRanges to data.frame without using row names
df <- as.data.frame(results, row.names = NULL)
# 2. Preserve duplicate GRanges names in a new column
df$variantName <- names(results)
# 3. Flatten any list columns by collapsing them into strings
#    (Example: the 'motifPos' column)
df$motifPos <- vapply(df$motifPos, function(x) paste(x, collapse = ","), character(1))

write.table(df,"motifBreakR_results_CHL_chranno.tsv",row.names = F,col.names = T,quote = F,sep = '\t')

 
# ##############################################################################
# # ChIPSeeker
# ##############################################################################
 
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#read the results of peak calling

data <- read.csv("peaks_by_cell_types_CHL_chranno.csv")

#change the result to GRanges object

peaks <- with(data, GRanges(seqnames = seqnames, IRanges(start, end), strand))

#annotate peaks

peakAnno <- annotatePeak(peaks, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

write.csv(peakAnno,"peaks_annotation_by_ChIPseeker_CHL_chranno.csv")

####################################################################################
# Linkage
####################################################################################

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
library(BSgenome.Hsapiens.UCSC.hg38)


####Linking peaks to genes

# can do parallelization, but not support in Seuratv5
#library(future)
# in Rstudio
#plan("multisession", workers = 10)
# in R
#plan("multicore", workers = 10)
# set 90g memory
#options(future.globals.maxSize = 200 * 1024^3)
#plan()

DefaultAssay(lung_4) <- "peaks"

# first compute the GC content for each peak
lung_4 <- RegionStats(lung_4, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
lung_4 <- LinkPeaks(
  object = lung_4,
  peak.assay = "peaks",
  expression.assay = "SCT",
  distance = 1e+06,
  min.cells = 10,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  verbose = TRUE
)

write.csv(Links(lung_4[["peaks"]]),"./links_1Mb_CHL.csv")

saveRDS(lung_4, file = "lung_4_linkage_1Mb_CHL.rds")

print('save 1MB')

# 2Mb

lung_4 <- RegionStats(lung_4, genome = BSgenome.Hsapiens.UCSC.hg38)

lung_4 <- LinkPeaks(
  object = lung_4,
  peak.assay = "peaks",
  expression.assay = "SCT",
  distance = 2e+06
)

saveRDS(lung_4, "lung_4_linkage_2Mb_CHL.rds")

write.csv(Links(lung_4[["peaks"]]), "./links_2Mb_CHL.csv")

print('save 2MB')

# 5Mb

lung_4 <- RegionStats(lung_4, genome = BSgenome.Hsapiens.UCSC.hg38)

lung_4 <- LinkPeaks(
  object = lung_4,
  peak.assay = "peaks",
  expression.assay = "SCT",
  distance = 5e+06
)

saveRDS(lung_4, "lung_4_linkage_5Mb_CHL.rds")

write.csv(Links(lung_4[["peaks"]]), "./links_5Mb_CHL.csv")

print('save 5MB')

###############################################################
# TF frootprint
###############################################################

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

# lung_peaks <- readRDS(file = "lung.rds")
DefaultAssay(lung_4) <- "peaks"

# The TF list is from motifbreaker 
# run per SNPs 
motif_B <- read.table('motifBreakR_results_CHL_chranno.tsv',header = T)

# 
# snps <- motif_B$SNP_id %>% unique() %>% sort()
# # motif list per SNPs

# for (i in snps){
#   #i <- snps[1]
#   print(i)
#   mf <- motif_B %>% filter(SNP_id == i)
#   mflist <- mf$geneSymbol %>% unique() %>% sort()
#   print(mflist)
# }


#### I added according to https://stuartlab.org/signac/articles/footprint ####
# extract position frequency matrices for the motifs

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information

lung_4 <- AddMotifs(lung_4, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)


#### 

# No NR0B1,ZNF219
rs104_target <- c("CTCF","MAZ","ZNF148")
rs124_target <- c("SOX18")
rs347_target <- c("TCF7")
rs680_target <- c("KLF15")
rs727_target <- c("GLIS3","SRF")

#### TF Footprinting 
lung_4 <- Footprint(
  object = lung_4,
  motif.name = c("CTCF","MAZ","ZNF148","SOX18","TCF7","KLF15","GLIS3","SRF"),
  genome = BSgenome.Hsapiens.UCSC.hg38)


saveRDS(lung_4, "CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED.rds")

print('save!')


p2 <- PlotFootprint(lung_4, features = c("CTCF","MAZ","ZNF148","SOX18","TCF7","KLF15","GLIS3","SRF"))
p2 + patchwork::plot_layout(ncol = 1)
# 
ggsave("footprint_CHL.pdf", 
        p2+ patchwork::plot_layout(ncol = 1), width = 10, height = 30)


####################################################################################
      #####################      run Cicero    ##########################
####################################################################################

library(SeuratWrappers)
library(cicero)
library(monocle3)


DefaultAssay(lung_4) <- "peaks"

# convert to CellDataSet format and make the cicero object
# wnn.umap has been changed to WNN.UMAP

lung_4.cds <- as.cell_data_set(lung_4)
lung_4.cicero <- make_cicero_cds(lung_4.cds, reduced_coordinates = reducedDims(lung_4.cds)$WNN.UMAP) 

#use a customized hg38 file

hg38_genome <- read.csv("../SHARE-seq/Seurat_SHARE/hg38_sequence_length.csv")

## Usually run with sample_num = 100 ##
conns <- run_cicero(lung_4.cicero, hg38_genome, sample_num = 100)

write.csv(conns,"./conns.csv")

ccans <- generate_ccans(conns)

write.csv(ccans,"./ccans.csv")

links <- ConnectionsToLinks(conns = conns, ccans = ccans)

write.csv(links,"./links.csv")


####################################################################################
# CUTOFF TESTING AREA
####################################################################################

lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED.rds')
peaks <- readRDS('peaks_CHL_chranno.rds')

link1 <- read.csv('Linkage_result/links_1Mb_CHL.csv')
link2 <- read.csv('Linkage_result/links_2Mb_CHL.csv')
link5 <- read.csv('Linkage_result/links_5Mb_CHL.csv')

CONNS <- read.csv('Cicero_result/conns_CHL.csv')
CCANs <- read.csv('Cicero_result/ccans_CHL.csv')
C_links <- read.csv('Cicero_result/links_CHL.csv')

# download and unzip
temp <- tempfile()
download.file("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)

# rename some columns to match requirements
gene_anno$chromosome <- gene_anno$seqid
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name


plot_connections(CONNS, "chr2", 9773451, 9848598,
                 gene_model = gene_anno, 
                 coaccess_cutoff = .25, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )

### Add Linkage links result to master RDS file

#lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED.rds')

ling_1MB_links <- readRDS('Linkage_result/lung_4_linkage_1Mb_CHL.rds')

Links(lung_4) <- Links(ling_1MB_links[["peaks"]])

# saveRDS(lung_4,'CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds')
lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds')


######

SC1 <- read.csv('Arrow_SCAVENGE_results/mono_mat_results_CHL.csv') 
SC2 <- read.csv('Arrow_SCAVENGE_results/mono_mat_with_permu_results_CHL.csv')

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




