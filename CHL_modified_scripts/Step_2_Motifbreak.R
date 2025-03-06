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


lung_4 <- readRDS('CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds')

DefaultAssay(lung_4) <- "peaks"

peaks <- readRDS('Peak_and_motifBreaker_results/peaks_CHL_chranno.rds')

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
#library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

#snp <- read.csv("snp_list.csv")
#snp <- tmp$RSNUM
# Now do all the SNPs

snp <- GWAS_subset$RSNUM
# rename chr14:35342071 to rs543000516
snp[269] <- "rs543000516"

snps.fromlist <- snps.from.rsid(rsid = snp,dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38,search.genome = BSgenome.Hsapiens.UCSC.hg38)
# snps.fromlist <- snps.from.rsid(rsid = snp,dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38,search.genome = BSgenome.Hsapiens.UCSC.hg38)

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

write.table(df,"motifBreakR_results_03032025.tsv",row.names = F,col.names = T,quote = F,sep = '\t')

print('done save')

