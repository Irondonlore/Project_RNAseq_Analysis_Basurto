################################################################################
#                   DESEQ2 NORMALIZATION FOR GSEA ANALYSIS                    #
################################################################################

#| GSEA does not normalize the data. If the gene expression data is derived from 
#| RNAseq it must be normalized for between-sample comparisons using an external 
#| normalization procedure (e.g. those in DESeq2 or Voom).

#| By using DESeq2 normalization we take into account: gene count comparisons between
#|samples and for DE analysis; NOT for within sample comparisons

#| See: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

################################################################################


################################################################################
################################ LIBRARIES #####################################
################################################################################
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
################################################################################



################################################################################
################################## DATA ########################################
################################################################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_ENRICHMENT/GSEA Analysis PCa vs BPH/Data/"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DESEQ2_ANALYSIS/Results/Tables/sample_info_extracted.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/3_STEP_RAW_COUNTS/Results/FullCounts_AC-45_RNAseq-FFPE.txt"
setwd(dir.proj)

#| Table Full counts
counts_data <- read.table(counts.file, sep = "\t", stringsAsFactors = TRUE)

#| Sample informatio data
sample_info <- read.csv(info.file, sep ="\t")
names(sample_info)
################################################################################


################################################################################
##########################   DESeq2 NORMALIZATION   ############################
################################################################################

#| Defining the design formula
design <- ~ 1 + DV200_Zscore + Diagnostico

#| Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = sample_info, design = design)

#| Estimating the size factors
dds <- estimateSizeFactors(dds)

#| Computing normalization
normalized_counts <- counts(dds, normalized=TRUE)

#| Adding the "description" column
description <- data.frame(rep("NA", times= length(rownames(normalized_counts))))
colnames(description) <- "description"
normalized_counts <- cbind(description, normalized_counts)
#| Saving the counts data
write.table(normalized_counts, file="normalized_counts_dds_for_GSEA_analysis.txt", sep="\t", col.names=NA)
normalized_counts

length(rownames(normalized_counts))
dim(normalized_counts)

f <- as.data.frame(sample_info$Diagnostico)

write.table(t(f), "sample_info_GSEA.cls", col.names=F, row.names = F, sep ="\t")
